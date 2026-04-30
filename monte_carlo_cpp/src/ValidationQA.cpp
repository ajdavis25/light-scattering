#include "ValidationQA.hpp"

#include "MonteCarloDriver.hpp"
#include "PhaseFunctions.hpp"
#include "Polarization.hpp"
#include "Vec3.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace {
constexpr double PI = 3.14159265358979323846;

struct ReferencePoint
{
    std::size_t original_index = 0;
    double tangent_altitude_km = 0.0;
    bool has_tangent_altitude = false;
    double zenith_deg = 0.0;
    double azimuth_deg = 0.0;
    double relative_azimuth_deg = 0.0;
    bool uses_relative_azimuth = false;
    double solar_zenith_deg = 0.0;
    double solar_azimuth_deg = 0.0;
    bool has_solar_geometry = false;
    double intensity = 0.0;
    bool has_intensity = false;
    double q = 0.0;
    double u = 0.0;
    double dop = 0.0;
    bool has_dop = false;
    double aop_deg = 0.0;
    bool has_aop = false;
    std::string neutral_point_label;
    bool has_neutral_point_label = false;
};

struct MeasurementModelCalibration
{
    double intensity_gain = 1.0;
    double dolp_scale = 1.0;
    double aop_offset_deg = 0.0;
};

std::string trim(const std::string &value)
{
    const std::string whitespace = " \t\r\n";
    const std::size_t begin = value.find_first_not_of(whitespace);
    if (begin == std::string::npos) {
        return "";
    }
    const std::size_t end = value.find_last_not_of(whitespace);
    return value.substr(begin, end - begin + 1);
}

std::vector<std::string> splitCsvLine(const std::string &line)
{
    std::vector<std::string> tokens;
    std::stringstream stream(line);
    std::string token;
    while (std::getline(stream, token, ',')) {
        tokens.push_back(trim(token));
    }
    return tokens;
}

double normalizeAzimuthDeg(double value)
{
    while (value < 0.0) {
        value += 360.0;
    }
    while (value >= 360.0) {
        value -= 360.0;
    }
    return value;
}

double wrapHalfTurnDeg(double value)
{
    while (value <= -90.0) {
        value += 180.0;
    }
    while (value > 90.0) {
        value -= 180.0;
    }
    return value;
}

double angularDifferenceAopDeg(double lhs, double rhs)
{
    return std::abs(wrapHalfTurnDeg(lhs - rhs));
}

std::vector<std::string> splitConfigList(const std::string &value)
{
    std::vector<std::string> items;
    std::stringstream stream(value);
    std::string token;
    while (std::getline(stream, token, ';')) {
        const std::string stripped = trim(token);
        if (!stripped.empty()) {
            items.push_back(stripped);
        }
    }
    return items;
}

std::string sanitizeMetricToken(const std::string &value)
{
    std::string token;
    token.reserve(value.size());
    bool lastUnderscore = false;
    for (unsigned char ch : value) {
        if (std::isalnum(ch)) {
            token.push_back(static_cast<char>(std::tolower(ch)));
            lastUnderscore = false;
        } else if (!lastUnderscore) {
            token.push_back('_');
            lastUnderscore = true;
        }
    }

    while (!token.empty() && token.front() == '_') {
        token.erase(token.begin());
    }
    while (!token.empty() && token.back() == '_') {
        token.pop_back();
    }

    if (token.empty()) {
        return "case";
    }
    return token;
}

std::string caseMetricPrefix(const std::string &kind, const SimulationConfig &config)
{
    return kind + "_" + sanitizeMetricToken(config.output.case_id) + "_";
}

std::vector<ReferencePoint> loadReferenceCsv(const std::string &filename)
{
    std::ifstream input(filename);
    if (!input) {
        return {};
    }

    std::string headerLine;
    if (!std::getline(input, headerLine)) {
        return {};
    }

    const std::vector<std::string> header = splitCsvLine(headerLine);
    std::map<std::string, std::size_t> columnIndex;
    for (std::size_t index = 0; index < header.size(); ++index) {
        columnIndex[header[index]] = index;
    }

    const auto getValue = [&](const std::vector<std::string> &values, const std::string &key, double defaultValue) {
        const auto iterator = columnIndex.find(key);
        if (iterator == columnIndex.end() || iterator->second >= values.size() || values[iterator->second].empty()) {
            return defaultValue;
        }
        return std::stod(values[iterator->second]);
    };

    std::vector<ReferencePoint> points;
    std::string line;
    while (std::getline(input, line)) {
        const std::string stripped = trim(line);
        if (stripped.empty()) {
            continue;
        }
        const std::vector<std::string> values = splitCsvLine(stripped);
        ReferencePoint point {};
        point.original_index = points.size();
        if (columnIndex.count("original_index")) {
            point.original_index = static_cast<std::size_t>(
                getValue(values, "original_index", static_cast<double>(point.original_index))
            );
        } else if (columnIndex.count("index")) {
            point.original_index = static_cast<std::size_t>(
                getValue(values, "index", static_cast<double>(point.original_index))
            );
        }
        if (columnIndex.count("tangent_altitude_km")) {
            point.tangent_altitude_km = getValue(values, "tangent_altitude_km", 0.0);
            point.has_tangent_altitude = true;
        }
        if (columnIndex.count("zenith_deg")) {
            point.zenith_deg = getValue(values, "zenith_deg", 0.0);
        } else if (columnIndex.count("altitude_deg")) {
            point.zenith_deg = 90.0 - getValue(values, "altitude_deg", 0.0);
        } else {
            continue;
        }

        if (columnIndex.count("relative_azimuth_deg")) {
            point.relative_azimuth_deg = getValue(values, "relative_azimuth_deg", 0.0);
            point.uses_relative_azimuth = true;
        } else if (columnIndex.count("azimuth_deg")) {
            point.azimuth_deg = getValue(values, "azimuth_deg", 0.0);
        }

        const auto solarZenithIterator = columnIndex.find("solar_zenith_deg");
        const auto solarAzimuthIterator = columnIndex.find("solar_azimuth_deg");
        if (solarZenithIterator != columnIndex.end() &&
            solarAzimuthIterator != columnIndex.end() &&
            solarZenithIterator->second < values.size() &&
            solarAzimuthIterator->second < values.size() &&
            !values[solarZenithIterator->second].empty() &&
            !values[solarAzimuthIterator->second].empty()) {
            point.solar_zenith_deg = std::stod(values[solarZenithIterator->second]);
            point.solar_azimuth_deg = normalizeAzimuthDeg(std::stod(values[solarAzimuthIterator->second]));
            point.has_solar_geometry = true;
        }

        const auto intensityIterator = columnIndex.find("intensity");
        if (intensityIterator != columnIndex.end() &&
            intensityIterator->second < values.size() &&
            !values[intensityIterator->second].empty()) {
            point.intensity = std::stod(values[intensityIterator->second]);
            point.has_intensity = true;
        }
        point.q = getValue(values, "q", 0.0);
        point.u = getValue(values, "u", 0.0);
        if (columnIndex.count("dop")) {
            point.dop = getValue(values, "dop", 0.0);
            point.has_dop = true;
        } else if (point.intensity > 0.0 && (std::abs(point.q) > 0.0 || std::abs(point.u) > 0.0)) {
            point.dop = std::sqrt(point.q * point.q + point.u * point.u) / point.intensity;
            point.has_dop = true;
        }
        if (columnIndex.count("aop_deg")) {
            point.aop_deg = getValue(values, "aop_deg", 0.0);
            point.has_aop = true;
        }
        const auto neutralLabelIterator = columnIndex.find("neutral_point_label");
        if (neutralLabelIterator != columnIndex.end() &&
            neutralLabelIterator->second < values.size() &&
            !values[neutralLabelIterator->second].empty()) {
            point.neutral_point_label = values[neutralLabelIterator->second];
            point.has_neutral_point_label = true;
        }

        points.push_back(point);
    }
    return points;
}

std::map<std::size_t, MeasurementModelCalibration> loadMeasurementModelCalibrationCsv(const std::string &filename)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("Unable to open measurement model calibration CSV: " + filename);
    }

    std::string headerLine;
    while (std::getline(input, headerLine)) {
        if (!trim(headerLine).empty() && trim(headerLine).front() != '#') {
            break;
        }
    }
    if (trim(headerLine).empty()) {
        throw std::runtime_error("Measurement model calibration CSV is missing a header: " + filename);
    }

    const std::vector<std::string> header = splitCsvLine(headerLine);
    std::map<std::string, std::size_t> columnIndex;
    for (std::size_t index = 0; index < header.size(); ++index) {
        columnIndex[header[index]] = index;
    }
    for (const std::string &required : {"index", "intensity_gain", "dolp_scale", "aop_offset_deg"}) {
        if (!columnIndex.count(required)) {
            throw std::runtime_error(
                "Measurement model calibration CSV is missing required column '" + required + "': " + filename
            );
        }
    }

    const auto getValue = [&](const std::vector<std::string> &values, const std::string &key) {
        const auto iterator = columnIndex.find(key);
        if (iterator == columnIndex.end() || iterator->second >= values.size() || values[iterator->second].empty()) {
            throw std::runtime_error(
                "Measurement model calibration CSV has an empty required value for '" + key + "': " + filename
            );
        }
        return std::stod(values[iterator->second]);
    };

    std::map<std::size_t, MeasurementModelCalibration> calibration;
    std::string line;
    while (std::getline(input, line)) {
        const std::string stripped = trim(line);
        if (stripped.empty() || stripped.front() == '#') {
            continue;
        }
        const std::vector<std::string> values = splitCsvLine(stripped);
        const std::size_t index = static_cast<std::size_t>(getValue(values, "index"));
        calibration[index] = {
            std::max(0.0, getValue(values, "intensity_gain")),
            std::max(0.0, getValue(values, "dolp_scale")),
            getValue(values, "aop_offset_deg"),
        };
    }
    return calibration;
}

StokesVector applyMeasurementModelCalibration(
    const StokesVector &raw,
    const MeasurementModelCalibration &calibration
)
{
    const double angle = 2.0 * calibration.aop_offset_deg * PI / 180.0;
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    const double rotatedQ = raw.Q * cosine - raw.U * sine;
    const double rotatedU = raw.Q * sine + raw.U * cosine;
    const double polarizationGain = calibration.intensity_gain * calibration.dolp_scale;
    return {
        raw.I * calibration.intensity_gain,
        rotatedQ * polarizationGain,
        rotatedU * polarizationGain,
        raw.V * polarizationGain,
    };
}

void applyMeasurementModelCalibration(
    SkyBinResult &bin,
    const MeasurementModelCalibration &calibration
)
{
    bin.mean = applyMeasurementModelCalibration(bin.mean, calibration);
    bin.first_order = applyMeasurementModelCalibration(bin.first_order, calibration);
    bin.second_order = applyMeasurementModelCalibration(bin.second_order, calibration);
    bin.higher_order = applyMeasurementModelCalibration(bin.higher_order, calibration);
    bin.second_rr = applyMeasurementModelCalibration(bin.second_rr, calibration);
    bin.second_ar = applyMeasurementModelCalibration(bin.second_ar, calibration);
    bin.second_ra = applyMeasurementModelCalibration(bin.second_ra, calibration);
    bin.second_aa = applyMeasurementModelCalibration(bin.second_aa, calibration);
}

void pushMetric(
    ValidationReport &report,
    const std::string &name,
    double value,
    double threshold,
    bool pass
)
{
    report.metrics.push_back({name, value, threshold, pass});
    report.overall_pass = report.overall_pass && pass;
}

double percentile(std::vector<double> values, double fraction)
{
    if (values.empty()) {
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    const std::size_t index = static_cast<std::size_t>(std::clamp(
        fraction * static_cast<double>(values.size() - 1),
        0.0,
        static_cast<double>(values.size() - 1)
    ));
    return values[index];
}

std::vector<SkyDirection> absoluteDirections(
    const std::vector<ReferencePoint> &reference,
    double sunAzimuthDeg
)
{
    std::vector<SkyDirection> directions;
    directions.reserve(reference.size());
    for (const ReferencePoint &point : reference) {
        directions.push_back({
            point.zenith_deg,
            point.uses_relative_azimuth
                ? normalizeAzimuthDeg(sunAzimuthDeg + point.relative_azimuth_deg)
                : normalizeAzimuthDeg(point.azimuth_deg),
        });
    }
    return directions;
}

SkyDirection absoluteDirectionForPoint(const ReferencePoint &point, double sunAzimuthDeg)
{
    return {
        point.zenith_deg,
        point.uses_relative_azimuth
            ? normalizeAzimuthDeg(sunAzimuthDeg + point.relative_azimuth_deg)
            : normalizeAzimuthDeg(point.azimuth_deg),
    };
}

std::vector<SkyBinResult> sampleReferenceDirections(
    const SimulationConfig &baseConfig,
    const std::vector<ReferencePoint> &reference
)
{
    const bool usesPerPointSolarGeometry = std::any_of(
        reference.begin(),
        reference.end(),
        [](const ReferencePoint &point) { return point.has_solar_geometry; }
    );

    if (!usesPerPointSolarGeometry) {
        const std::vector<SkyDirection> directions = absoluteDirections(reference, baseConfig.solar.azimuth_deg);
        return sampleSkyDirections(baseConfig, directions);
    }

    std::map<std::pair<long long, long long>, std::vector<std::size_t>> groupedIndices;
    for (std::size_t index = 0; index < reference.size(); ++index) {
        const double solarZenithDeg = reference[index].has_solar_geometry
            ? reference[index].solar_zenith_deg
            : baseConfig.solar.zenith_deg;
        const double solarAzimuthDeg = reference[index].has_solar_geometry
            ? normalizeAzimuthDeg(reference[index].solar_azimuth_deg)
            : normalizeAzimuthDeg(baseConfig.solar.azimuth_deg);
        groupedIndices[{
            static_cast<long long>(std::llround(solarZenithDeg * 1000.0)),
            static_cast<long long>(std::llround(solarAzimuthDeg * 1000.0)),
        }].push_back(index);
    }

    std::vector<SkyBinResult> sampled(reference.size());
    for (const auto &[solarKey, indices] : groupedIndices) {
        SimulationConfig groupedConfig = baseConfig;
        groupedConfig.solar.use_explicit_angles = true;
        groupedConfig.solar.zenith_deg = static_cast<double>(solarKey.first) / 1000.0;
        groupedConfig.solar.azimuth_deg = static_cast<double>(solarKey.second) / 1000.0;

        std::vector<ReferencePoint> groupedReference;
        groupedReference.reserve(indices.size());
        for (std::size_t index : indices) {
            groupedReference.push_back(reference[index]);
        }

        const std::vector<SkyDirection> groupedDirections = absoluteDirections(
            groupedReference,
            groupedConfig.solar.azimuth_deg
        );
        const std::vector<SkyBinResult> groupedResults = sampleSkyDirections(groupedConfig, groupedDirections);
        for (std::size_t groupedIndex = 0; groupedIndex < indices.size(); ++groupedIndex) {
            sampled[indices[groupedIndex]] = groupedResults[groupedIndex];
        }
    }

    return sampled;
}

void writeBenchmarkComparisonCsv(
    const SimulationConfig &baseConfig,
    const SimulationConfig &benchmarkConfig,
    const std::vector<ReferencePoint> &reference,
    const std::vector<SkyBinResult> &model
)
{
    const std::filesystem::path outputDir =
        std::filesystem::path(baseConfig.output.output_dir) / "validation";
    std::filesystem::create_directories(outputDir);
    const std::string filename = benchmarkConfig.output.case_id.empty()
        ? "benchmark_comparison.csv"
        : benchmarkConfig.output.case_id + "_comparison.csv";
    std::ofstream output(outputDir / filename);
    output
        << "index,tangent_altitude_km,zenith_deg,azimuth_deg,solar_zenith_deg,solar_azimuth_deg,"
        << "reference_intensity,model_intensity,normalized_reference,normalized_model,relative_intensity_error,"
        << "reference_dolp,model_dolp,absolute_dolp_error\n";

    double referencePeak = 0.0;
    double modelPeak = 0.0;
    for (std::size_t index = 0; index < reference.size() && index < model.size(); ++index) {
        referencePeak = std::max(referencePeak, reference[index].intensity);
        modelPeak = std::max(modelPeak, model[index].mean.I);
    }

    for (std::size_t index = 0; index < reference.size() && index < model.size(); ++index) {
        const double referenceDolp = reference[index].has_dop ? reference[index].dop : 0.0;
        const double modelDolp = degreeOfLinearPolarization(model[index].mean);
        const double normalizedReference = reference[index].intensity / std::max(1.0e-12, referencePeak);
        const double normalizedModel = model[index].mean.I / std::max(1.0e-12, modelPeak);
        const double relativeIntensityError = std::abs(normalizedModel - normalizedReference) /
            std::max(1.0e-12, normalizedReference);
        const double absoluteDolpError = reference[index].has_dop
            ? std::abs(modelDolp - referenceDolp)
            : 0.0;

        output
            << index << ","
            << (reference[index].has_tangent_altitude ? reference[index].tangent_altitude_km : 0.0) << ","
            << reference[index].zenith_deg << ","
            << reference[index].azimuth_deg << ","
            << reference[index].solar_zenith_deg << ","
            << reference[index].solar_azimuth_deg << ","
            << reference[index].intensity << ","
            << model[index].mean.I << ","
            << normalizedReference << ","
            << normalizedModel << ","
            << relativeIntensityError << ","
            << referenceDolp << ","
            << modelDolp << ","
            << absoluteDolpError << "\n";
    }
}

void writeMeasurementComparisonCsv(
    const SimulationConfig &baseConfig,
    const SimulationConfig &measurementConfig,
    const std::vector<ReferencePoint> &reference,
    const std::vector<SkyDirection> &referenceDirections,
    const std::vector<SkyBinResult> &model
)
{
    const std::filesystem::path outputDir =
        std::filesystem::path(baseConfig.output.output_dir) / "validation";
    std::filesystem::create_directories(outputDir);
    const std::string filename = measurementConfig.output.case_id.empty()
        ? "measurement_comparison.csv"
        : measurementConfig.output.case_id + "_comparison.csv";
    std::ofstream output(outputDir / filename);
    output
        << "index,zenith_deg,relative_azimuth_deg,absolute_azimuth_deg,"
        << "reference_intensity,model_intensity,normalized_reference,normalized_model,"
        << "reference_dop,model_dop,dop_abs_error,reference_aop_deg,model_aop_deg,aop_abs_error_deg,"
        << "signed_dolp_bias,first_frac,second_frac,higher_frac,"
        << "second_rr_frac,second_ar_frac,second_ra_frac,second_aa_frac\n";

    double referencePeak = 0.0;
    double modelPeak = 0.0;
    for (std::size_t index = 0; index < reference.size() && index < model.size(); ++index) {
        if (reference[index].has_intensity) {
            referencePeak = std::max(referencePeak, reference[index].intensity);
        }
        modelPeak = std::max(modelPeak, model[index].mean.I);
    }

    for (std::size_t index = 0; index < reference.size() && index < model.size(); ++index) {
        const double sunAzimuthDeg = reference[index].has_solar_geometry
            ? reference[index].solar_azimuth_deg
            : measurementConfig.solar.azimuth_deg;
        const double relativeAzimuthDeg = reference[index].uses_relative_azimuth
            ? reference[index].relative_azimuth_deg
            : normalizeAzimuthDeg(referenceDirections[index].azimuth_deg - sunAzimuthDeg);
        const double normalizedReference = reference[index].has_intensity
            ? reference[index].intensity / std::max(1.0e-12, referencePeak)
            : 0.0;
        const double normalizedModel = model[index].mean.I / std::max(1.0e-12, modelPeak);
        const double modelDolp = degreeOfLinearPolarization(model[index].mean);
        const double modelAopDeg = angleOfLinearPolarizationRad(model[index].mean) * 180.0 / PI;
        const double dopAbsError = reference[index].has_dop
            ? std::abs(modelDolp - reference[index].dop)
            : 0.0;
        const double aopAbsErrorDeg = reference[index].has_aop
            ? angularDifferenceAopDeg(modelAopDeg, reference[index].aop_deg)
            : 0.0;
        const double signedDolpBias = reference[index].has_dop
            ? modelDolp - reference[index].dop
            : 0.0;
        const double totalModelI = std::max(1.0e-12, model[index].mean.I);
        const double secondTotalI = std::max(1.0e-12, model[index].second_order.I);

        output
            << index << ","
            << reference[index].zenith_deg << ","
            << relativeAzimuthDeg << ","
            << referenceDirections[index].azimuth_deg << ","
            << (reference[index].has_intensity ? reference[index].intensity : 0.0) << ","
            << model[index].mean.I << ","
            << normalizedReference << ","
            << normalizedModel << ","
            << (reference[index].has_dop ? reference[index].dop : 0.0) << ","
            << modelDolp << ","
            << dopAbsError << ","
            << (reference[index].has_aop ? reference[index].aop_deg : 0.0) << ","
            << modelAopDeg << ","
            << aopAbsErrorDeg << ","
            << signedDolpBias << ","
            << model[index].first_order.I / totalModelI << ","
            << model[index].second_order.I / totalModelI << ","
            << model[index].higher_order.I / totalModelI << ","
            << model[index].second_rr.I / secondTotalI << ","
            << model[index].second_ar.I / secondTotalI << ","
            << model[index].second_ra.I / secondTotalI << ","
            << model[index].second_aa.I / secondTotalI << "\n";
    }
}

Vec3 unitVectorFromZenithAzimuth(double zenithDeg, double azimuthDeg)
{
    const double zenithRad = zenithDeg * PI / 180.0;
    const double azimuthRad = azimuthDeg * PI / 180.0;
    return {
        std::sin(zenithRad) * std::cos(azimuthRad),
        std::sin(zenithRad) * std::sin(azimuthRad),
        std::cos(zenithRad),
    };
}

double angularSeparationDeg(
    double zenithA,
    double azimuthA,
    double zenithB,
    double azimuthB
)
{
    const Vec3 a = unitVectorFromZenithAzimuth(zenithA, azimuthA);
    const Vec3 b = unitVectorFromZenithAzimuth(zenithB, azimuthB);
    return std::acos(std::clamp(dot(a, b), -1.0, 1.0)) * 180.0 / PI;
}

double relativeAzimuthFromAbsolute(double absoluteAzimuthDeg, double sunAzimuthDeg)
{
    return normalizeAzimuthDeg(absoluteAzimuthDeg - sunAzimuthDeg);
}

void evaluateConvergence(
    ValidationReport &report,
    const SimulationConfig &baseConfig
)
{
    SimulationConfig coarseConfig = baseConfig;
    coarseConfig.monte_carlo.photons_per_bin = std::max(64, baseConfig.monte_carlo.photons_per_bin / 2);
    coarseConfig.output.case_id = baseConfig.output.case_id + "_validation_coarse";

    SimulationConfig fineConfig = baseConfig;
    fineConfig.output.case_id = baseConfig.output.case_id + "_validation_fine";

    const SkyResult coarse = runMonteCarloSimulation(coarseConfig);
    const SkyResult fine = runMonteCarloSimulation(fineConfig);

    const double intensityChange = std::abs(fine.peak_intensity - coarse.peak_intensity) / std::max(1.0e-12, fine.peak_intensity);
    const double dopChange = std::abs(fine.peak_dolp - coarse.peak_dolp);
    const double fluxChange = std::abs(fine.hemispheric_flux_estimate - coarse.hemispheric_flux_estimate)
        / std::max(1.0e-12, fine.hemispheric_flux_estimate);

    pushMetric(report, "convergence_peak_intensity_rel", intensityChange, 0.05, intensityChange <= 0.05);
    pushMetric(report, "convergence_peak_dolp_abs", dopChange, 0.02, dopChange <= 0.02);
    pushMetric(report, "convergence_flux_rel", fluxChange, 0.02, fluxChange <= 0.02);
}

bool evaluateBenchmarkConfig(
    ValidationReport &report,
    const SimulationConfig &baseConfig,
    const std::string &configPath
)
{
    const SimulationConfig benchmarkConfig = loadSimulationConfig(configPath);
    const std::string prefix = caseMetricPrefix("benchmark", benchmarkConfig);
    if (benchmarkConfig.output.benchmark_reference_csv.empty()) {
        pushMetric(report, prefix + "reference_present", 1.0, 0.0, false);
        report.notes.push_back(
            "Benchmark validation failed for case " + benchmarkConfig.output.case_id +
            " because benchmark_reference_csv was not set."
        );
        return false;
    }

    const std::vector<ReferencePoint> reference = loadReferenceCsv(benchmarkConfig.output.benchmark_reference_csv);
    if (reference.empty()) {
        pushMetric(report, prefix + "reference_loaded", 1.0, 0.0, false);
        report.notes.push_back(
            "Benchmark validation failed for case " + benchmarkConfig.output.case_id +
            " because the benchmark reference CSV was empty or unreadable."
        );
        return false;
    }

    const std::vector<SkyBinResult> model = sampleReferenceDirections(benchmarkConfig, reference);
    writeBenchmarkComparisonCsv(baseConfig, benchmarkConfig, reference, model);

    double referencePeak = 0.0;
    double modelPeak = 0.0;
    for (std::size_t index = 0; index < reference.size(); ++index) {
        referencePeak = std::max(referencePeak, reference[index].intensity);
        modelPeak = std::max(modelPeak, model[index].mean.I);
    }

    std::vector<double> intensityErrors;
    std::vector<double> dopErrors;
    for (std::size_t index = 0; index < reference.size(); ++index) {
        const double normalizedReference = reference[index].intensity / std::max(1.0e-12, referencePeak);
        if (normalizedReference < baseConfig.validation.benchmark_mask_fraction_of_peak) {
            continue;
        }
        const double normalizedModel = model[index].mean.I / std::max(1.0e-12, modelPeak);
        intensityErrors.push_back(
            std::abs(normalizedModel - normalizedReference) / std::max(1.0e-12, normalizedReference)
        );

        if (reference[index].has_dop) {
            dopErrors.push_back(std::abs(degreeOfLinearPolarization(model[index].mean) - reference[index].dop));
        }
    }

    if (intensityErrors.empty()) {
        pushMetric(report, prefix + "points_in_mask", 1.0, 0.0, false);
        report.notes.push_back(
            "Benchmark validation failed for case " + benchmarkConfig.output.case_id +
            " because no benchmark points passed the benchmark intensity mask."
        );
        return false;
    }

    const double medianIntensity = percentile(intensityErrors, 0.5);
    const double p95Intensity = percentile(intensityErrors, 0.95);
    pushMetric(
        report,
        prefix + "median_intensity_rel",
        medianIntensity,
        baseConfig.validation.median_intensity_error_limit,
        medianIntensity <= baseConfig.validation.median_intensity_error_limit
    );
    pushMetric(
        report,
        prefix + "p95_intensity_rel",
        p95Intensity,
        baseConfig.validation.p95_intensity_error_limit,
        p95Intensity <= baseConfig.validation.p95_intensity_error_limit
    );

    if (!dopErrors.empty()) {
        const double medianDolp = percentile(dopErrors, 0.5);
        const double p95Dolp = percentile(dopErrors, 0.95);
        pushMetric(
            report,
            prefix + "median_dolp_abs",
            medianDolp,
            baseConfig.validation.median_dolp_abs_error_limit,
            medianDolp <= baseConfig.validation.median_dolp_abs_error_limit
        );
        pushMetric(
            report,
            prefix + "p95_dolp_abs",
            p95Dolp,
            baseConfig.validation.p95_dolp_abs_error_limit,
            p95Dolp <= baseConfig.validation.p95_dolp_abs_error_limit
        );
        return true;
    }

    report.notes.push_back(
        "Benchmark reference for case " + benchmarkConfig.output.case_id +
        " does not include polarization, so DoLP benchmark metrics were not evaluated."
    );
    return false;
}

void evaluateBenchmarkCase(
    ValidationReport &report,
    const SimulationConfig &baseConfig
)
{
    if (baseConfig.output.benchmark_case_config.empty()) {
        pushMetric(report, "benchmark_case_present", 1.0, 0.0, false);
        report.notes.push_back("Benchmark validation failed because benchmark_case_config was not set.");
        return;
    }

    const std::vector<std::string> benchmarkCases = splitConfigList(baseConfig.output.benchmark_case_config);
    if (benchmarkCases.empty()) {
        pushMetric(report, "benchmark_case_present", 1.0, 0.0, false);
        report.notes.push_back("Benchmark validation failed because benchmark_case_config did not contain any usable paths.");
        return;
    }

    bool polarizationReferencePresent = false;
    for (const std::string &configPath : benchmarkCases) {
        polarizationReferencePresent = evaluateBenchmarkConfig(report, baseConfig, configPath) || polarizationReferencePresent;
    }

    pushMetric(
        report,
        "benchmark_polarization_reference_present",
        polarizationReferencePresent ? 0.0 : 1.0,
        0.0,
        polarizationReferencePresent
    );
}

bool evaluateMeasurementConfig(
    ValidationReport &report,
    const SimulationConfig &baseConfig,
    const std::string &configPath
)
{
    const SimulationConfig measurementConfig = loadSimulationConfig(configPath);
    const std::string prefix = caseMetricPrefix("measurement", measurementConfig);
    if (measurementConfig.output.measurement_reference_csv.empty()) {
        pushMetric(report, prefix + "reference_present", 1.0, 0.0, false);
        report.notes.push_back(
            "Measurement validation failed for case " + measurementConfig.output.case_id +
            " because measurement_reference_csv was not set."
        );
        return false;
    }

    const std::vector<ReferencePoint> reference = loadReferenceCsv(measurementConfig.output.measurement_reference_csv);
    if (reference.empty()) {
        pushMetric(report, prefix + "reference_loaded", 1.0, 0.0, false);
        report.notes.push_back(
            "Measurement validation failed for case " + measurementConfig.output.case_id +
            " because the measurement reference CSV was empty or unreadable."
        );
        return false;
    }

    std::vector<SkyBinResult> model = sampleReferenceDirections(measurementConfig, reference);
    if (!measurementConfig.output.measurement_model_calibration_csv.empty()) {
        const std::map<std::size_t, MeasurementModelCalibration> calibration =
            loadMeasurementModelCalibrationCsv(measurementConfig.output.measurement_model_calibration_csv);
        for (std::size_t index = 0; index < reference.size(); ++index) {
            const auto iterator = calibration.find(reference[index].original_index);
            if (iterator == calibration.end()) {
                throw std::runtime_error(
                    "Measurement model calibration is missing reference index " +
                    std::to_string(reference[index].original_index) +
                    " for case " + measurementConfig.output.case_id
                );
            }
            applyMeasurementModelCalibration(model[index], iterator->second);
        }
        pushMetric(
            report,
            prefix + "model_calibration_applied",
            1.0,
            1.0,
            true
        );
        report.notes.push_back(
            "Measurement model calibration applied for case " + measurementConfig.output.case_id +
            " from " + measurementConfig.output.measurement_model_calibration_csv +
            ". This is a calibrated model-quality closure, not an independent first-principles pass."
        );
    }
    std::vector<SkyDirection> referenceDirections;
    referenceDirections.reserve(reference.size());
    for (const ReferencePoint &point : reference) {
        const double sunAzimuthDeg = point.has_solar_geometry
            ? point.solar_azimuth_deg
            : measurementConfig.solar.azimuth_deg;
        referenceDirections.push_back(absoluteDirectionForPoint(point, sunAzimuthDeg));
    }
    writeMeasurementComparisonCsv(baseConfig, measurementConfig, reference, referenceDirections, model);

    const bool hasIntensityReference = std::any_of(
        reference.begin(),
        reference.end(),
        [](const ReferencePoint &point) { return point.has_intensity; }
    );

    double referencePeak = 0.0;
    double modelPeak = 0.0;
    if (hasIntensityReference) {
        for (std::size_t index = 0; index < reference.size(); ++index) {
            if (!reference[index].has_intensity) {
                continue;
            }
            referencePeak = std::max(referencePeak, reference[index].intensity);
            modelPeak = std::max(modelPeak, model[index].mean.I);
        }
    }

    double sumSquared = 0.0;
    double count = 0.0;
    std::vector<double> dopErrors;
    std::vector<double> aopErrorsDeg;
    std::vector<double> solarVerticalSignedDolpBias;
    std::size_t brightestReferenceIndex = 0;
    std::size_t brightestModelIndex = 0;
    std::size_t brightestRegionModelIndex = 0;
    double brightestReferenceValue = -1.0;
    double brightestModelValue = -1.0;
    double brightestRegionModelValue = -1.0;
    bool hasBrightestRegionModel = false;
    const double brightestRegionReferenceFraction = std::clamp(
        baseConfig.validation.brightest_region_reference_fraction_of_peak,
        0.0,
        1.0
    );

    for (std::size_t index = 0; index < reference.size(); ++index) {
        const double modelDolp = degreeOfLinearPolarization(model[index].mean);
        const double modelAopDeg = angleOfLinearPolarizationRad(model[index].mean) * 180.0 / PI;
        const double sunAzimuthDeg = reference[index].has_solar_geometry
            ? reference[index].solar_azimuth_deg
            : measurementConfig.solar.azimuth_deg;
        const double relativeAzimuthDeg = reference[index].uses_relative_azimuth
            ? reference[index].relative_azimuth_deg
            : relativeAzimuthFromAbsolute(referenceDirections[index].azimuth_deg, sunAzimuthDeg);

        if (hasIntensityReference && reference[index].has_intensity) {
            const double normalizedReference = reference[index].intensity / std::max(1.0e-12, referencePeak);
            const double normalizedModel = model[index].mean.I / std::max(1.0e-12, modelPeak);

            if (normalizedReference > brightestReferenceValue) {
                brightestReferenceValue = normalizedReference;
                brightestReferenceIndex = index;
            }
            if (normalizedModel > brightestModelValue) {
                brightestModelValue = normalizedModel;
                brightestModelIndex = index;
            }
            if (brightestRegionReferenceFraction > 0.0 &&
                normalizedReference >= brightestRegionReferenceFraction &&
                normalizedModel > brightestRegionModelValue) {
                brightestRegionModelValue = normalizedModel;
                brightestRegionModelIndex = index;
                hasBrightestRegionModel = true;
            }

            if (normalizedReference >= baseConfig.validation.measurement_mask_fraction_of_peak) {
                const double diff = normalizedModel - normalizedReference;
                sumSquared += diff * diff;
                count += 1.0;
                if (reference[index].has_dop) {
                    const double signedDolpBias = modelDolp - reference[index].dop;
                    dopErrors.push_back(std::abs(signedDolpBias));
                    if (reference[index].zenith_deg >= 30.0 &&
                        reference[index].zenith_deg <= 70.0 &&
                        relativeAzimuthDeg >= 240.0 &&
                        relativeAzimuthDeg <= 300.0) {
                        solarVerticalSignedDolpBias.push_back(signedDolpBias);
                    }
                }
                if (reference[index].has_aop && reference[index].has_dop && reference[index].dop >= 0.15) {
                    aopErrorsDeg.push_back(angularDifferenceAopDeg(modelAopDeg, reference[index].aop_deg));
                }
            }
        } else if (!hasIntensityReference && reference[index].has_dop) {
            const double signedDolpBias = modelDolp - reference[index].dop;
            dopErrors.push_back(std::abs(signedDolpBias));
            if (reference[index].zenith_deg >= 30.0 &&
                reference[index].zenith_deg <= 70.0 &&
                relativeAzimuthDeg >= 240.0 &&
                relativeAzimuthDeg <= 300.0) {
                solarVerticalSignedDolpBias.push_back(signedDolpBias);
            }
            if (reference[index].has_aop && reference[index].has_dop && reference[index].dop >= 0.15) {
                aopErrorsDeg.push_back(angularDifferenceAopDeg(modelAopDeg, reference[index].aop_deg));
            }
        }
    }

    if (hasIntensityReference && count <= 0.0) {
        pushMetric(report, prefix + "points_in_mask", 1.0, 0.0, false);
        report.notes.push_back(
            "Measurement validation failed for case " + measurementConfig.output.case_id +
            " because no measurement points passed the intensity mask."
        );
        return false;
    }

    if (hasIntensityReference) {
        const double rmse = std::sqrt(sumSquared / count);
        pushMetric(
            report,
            prefix + "normalized_rmse",
            rmse,
            baseConfig.validation.normalized_rmse_limit,
            rmse <= baseConfig.validation.normalized_rmse_limit
        );

        const std::size_t validationBrightestModelIndex =
            hasBrightestRegionModel ? brightestRegionModelIndex : brightestModelIndex;
        const double brightestLocationError = angularSeparationDeg(
            referenceDirections[brightestReferenceIndex].zenith_deg,
            referenceDirections[brightestReferenceIndex].azimuth_deg,
            model[validationBrightestModelIndex].zenith_deg,
            model[validationBrightestModelIndex].azimuth_deg
        );
        pushMetric(
            report,
            prefix + "brightest_location_deg",
            brightestLocationError,
            baseConfig.validation.brightest_region_deg_limit,
            brightestLocationError <= baseConfig.validation.brightest_region_deg_limit
        );
    } else {
        report.notes.push_back(
            "Measurement reference for case " + measurementConfig.output.case_id +
            " does not include intensity, so only polarization metrics were evaluated."
        );
    }

    if (!dopErrors.empty()) {
        const double medianDolp = percentile(dopErrors, 0.5);
        const double p95Dolp = percentile(dopErrors, 0.95);
        pushMetric(
            report,
            prefix + "median_dolp_abs",
            medianDolp,
            baseConfig.validation.median_dolp_abs_error_limit,
            medianDolp <= baseConfig.validation.median_dolp_abs_error_limit
        );
        pushMetric(
            report,
            prefix + "p95_dolp_abs",
            p95Dolp,
            baseConfig.validation.p95_dolp_abs_error_limit,
            p95Dolp <= baseConfig.validation.p95_dolp_abs_error_limit
        );
    } else {
        report.notes.push_back(
            "Measurement reference for case " + measurementConfig.output.case_id +
            " does not include polarization, so DoLP measurement metrics were not evaluated."
        );
    }

    if (!aopErrorsDeg.empty()) {
        const double medianAop = percentile(aopErrorsDeg, 0.5);
        const double p95Aop = percentile(aopErrorsDeg, 0.95);
        pushMetric(
            report,
            prefix + "median_aop_deg",
            medianAop,
            baseConfig.validation.median_aop_error_deg_limit,
            medianAop <= baseConfig.validation.median_aop_error_deg_limit
        );
        pushMetric(
            report,
            prefix + "p95_aop_deg",
            p95Aop,
            baseConfig.validation.p95_aop_error_deg_limit,
            p95Aop <= baseConfig.validation.p95_aop_error_deg_limit
        );
    } else {
        report.notes.push_back(
            "Measurement reference for case " + measurementConfig.output.case_id +
            " does not include AoP values above the DoLP threshold, so AoP metrics were not evaluated."
        );
    }

    if (!solarVerticalSignedDolpBias.empty()) {
        const double meanBias = std::accumulate(
            solarVerticalSignedDolpBias.begin(),
            solarVerticalSignedDolpBias.end(),
            0.0
        ) / static_cast<double>(solarVerticalSignedDolpBias.size());
        pushMetric(
            report,
            prefix + "solar_vertical_signed_dolp_bias",
            std::abs(meanBias),
            baseConfig.validation.solar_vertical_signed_dolp_bias_limit,
            std::abs(meanBias) <= baseConfig.validation.solar_vertical_signed_dolp_bias_limit
        );
    } else {
        report.notes.push_back(
            "Measurement reference for case " + measurementConfig.output.case_id +
            " did not provide enough solar-vertical points for the signed DoLP bias metric."
        );
    }

    if (std::any_of(reference.begin(), reference.end(), [](const ReferencePoint &point) {
        return point.has_neutral_point_label;
    })) {
        std::vector<double> neutralPointErrorsDeg;
        for (std::size_t index = 0; index < reference.size(); ++index) {
            if (!reference[index].has_neutral_point_label) {
                continue;
            }
            double bestError = std::numeric_limits<double>::infinity();
            for (std::size_t candidateIndex = 0; candidateIndex < model.size(); ++candidateIndex) {
                const double candidateError = angularSeparationDeg(
                    referenceDirections[index].zenith_deg,
                    referenceDirections[index].azimuth_deg,
                    model[candidateIndex].zenith_deg,
                    model[candidateIndex].azimuth_deg
                );
                if (degreeOfLinearPolarization(model[candidateIndex].mean) <= 0.1) {
                    bestError = std::min(bestError, candidateError);
                }
            }
            if (std::isfinite(bestError)) {
                neutralPointErrorsDeg.push_back(bestError);
            }
        }
        if (!neutralPointErrorsDeg.empty()) {
            const double p95NeutralPoint = percentile(neutralPointErrorsDeg, 0.95);
            pushMetric(
                report,
                prefix + "neutral_point_location_deg",
                p95NeutralPoint,
                baseConfig.validation.neutral_point_location_deg_limit,
                p95NeutralPoint <= baseConfig.validation.neutral_point_location_deg_limit
            );
        }
    }

    return !dopErrors.empty();
}

void evaluateMeasurementCase(
    ValidationReport &report,
    const SimulationConfig &baseConfig
)
{
    if (baseConfig.output.measurement_case_config.empty()) {
        pushMetric(report, "measurement_case_present", 1.0, 0.0, false);
        report.notes.push_back("Measurement validation failed because measurement_case_config was not set.");
        return;
    }

    const std::vector<std::string> measurementCases = splitConfigList(baseConfig.output.measurement_case_config);
    if (measurementCases.empty()) {
        pushMetric(report, "measurement_case_present", 1.0, 0.0, false);
        report.notes.push_back("Measurement validation failed because measurement_case_config did not contain any usable paths.");
        return;
    }

    bool polarizationReferencePresent = false;
    for (const std::string &configPath : measurementCases) {
        polarizationReferencePresent = evaluateMeasurementConfig(report, baseConfig, configPath) || polarizationReferencePresent;
    }

    if (!baseConfig.output.paper_primary_measurement_case_config.empty()) {
        const bool frozen = baseConfig.output.paper_primary_measurement_frozen;
        pushMetric(
            report,
            "paper_primary_measurement_frozen",
            frozen ? 0.0 : 1.0,
            0.0,
            frozen
        );
        const bool present = std::find(
            measurementCases.begin(),
            measurementCases.end(),
            baseConfig.output.paper_primary_measurement_case_config
        ) != measurementCases.end();
        pushMetric(
            report,
            "paper_primary_measurement_case_present",
            present ? 0.0 : 1.0,
            0.0,
            present
        );
        if (!frozen) {
            report.notes.push_back(
                "Paper validation is blocked because the configured primary twilight full-sky measurement case is still marked as interim and not frozen."
            );
        }
    }

    pushMetric(
        report,
        "measurement_polarization_reference_present",
        polarizationReferencePresent ? 0.0 : 1.0,
        0.0,
        polarizationReferencePresent
    );
}
}

ValidationReport runValidationSuite(const std::string &config_path)
{
    ValidationReport report;

    const int integrationSamples = 4000;
    double rayleighIntegral = 0.0;
    for (int index = 0; index < integrationSamples; ++index) {
        const double mu = -1.0 + (2.0 * index + 1.0) / integrationSamples;
        rayleighIntegral += 2.0 * PI * rayleighPhase(mu) * (2.0 / integrationSamples);
    }
    pushMetric(report, "rayleigh_phase_normalization", std::abs(rayleighIntegral - 1.0), 1.0e-3, std::abs(rayleighIntegral - 1.0) <= 1.0e-3);

    const StokesVector ninetyDegree = applyRayleighMueller(initUnpolarized(1.0), 0.0);
    const double dop90 = degreeOfLinearPolarization(ninetyDegree);
    pushMetric(report, "rayleigh_90deg_dolp_error", std::abs(dop90 - 1.0), 1.0e-6, std::abs(dop90 - 1.0) <= 1.0e-6);
    pushMetric(report, "rayleigh_physical_stokes", isPhysicallyValid(ninetyDegree) ? 0.0 : 1.0, 0.0, isPhysicallyValid(ninetyDegree));

    const SimulationConfig baseConfig = loadSimulationConfig(config_path);
    evaluateConvergence(report, baseConfig);
    evaluateBenchmarkCase(report, baseConfig);
    evaluateMeasurementCase(report, baseConfig);

    return report;
}

void writeValidationReport(const ValidationReport &report, const std::string &output_dir)
{
    std::filesystem::create_directories(output_dir);
    std::ofstream output(std::filesystem::path(output_dir) / "validation_report.json");
    output << "{\n";
    output << "  \"overall_pass\": " << (report.overall_pass ? "true" : "false") << ",\n";
    output << "  \"metrics\": [\n";
    for (std::size_t index = 0; index < report.metrics.size(); ++index) {
        const ValidationMetric &metric = report.metrics[index];
        output << "    {\"name\": \"" << metric.name << "\", \"value\": " << metric.value
               << ", \"threshold\": " << metric.threshold
               << ", \"pass\": " << (metric.pass ? "true" : "false") << "}";
        output << (index + 1 < report.metrics.size() ? ",\n" : "\n");
    }
    output << "  ],\n";
    output << "  \"notes\": [\n";
    for (std::size_t index = 0; index < report.notes.size(); ++index) {
        output << "    \"" << report.notes[index] << "\"";
        output << (index + 1 < report.notes.size() ? ",\n" : "\n");
    }
    output << "  ]\n";
    output << "}\n";
}
