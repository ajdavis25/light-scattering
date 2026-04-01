#include "MonteCarloDriver.hpp"
#include "Polarization.hpp"
#include "Vec3.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace {
constexpr double PI = 3.14159265358979323846;

struct ReferencePoint
{
    std::size_t original_index = 0;
    double zenith_deg = 0.0;
    double azimuth_deg = 0.0;
    double relative_azimuth_deg = 0.0;
    bool uses_relative_azimuth = false;
    double intensity = 0.0;
    bool has_intensity = false;
    double q = 0.0;
    double u = 0.0;
    double dop = 0.0;
    bool has_dop = false;
    double aop_deg = 0.0;
    bool has_aop = false;
};

struct PointDiagnostics
{
    std::size_t index = 0;
    double zenith_deg = 0.0;
    double relative_azimuth_deg = 0.0;
    double absolute_azimuth_deg = 0.0;
    double normalized_reference = 0.0;
    double normalized_model = 0.0;
    double reference_dop = 0.0;
    bool has_reference_dop = false;
    double model_dop = 0.0;
    double dop_abs_error = 0.0;
    double reference_aop_deg = 0.0;
    bool has_reference_aop = false;
    double model_aop_deg = 0.0;
    double aop_abs_error_deg = 0.0;
    double signed_dolp_bias = 0.0;
    double first_frac = 0.0;
    double second_frac = 0.0;
    double higher_frac = 0.0;
    double second_rr_frac = 0.0;
    double second_ar_frac = 0.0;
    double second_ra_frac = 0.0;
    double second_aa_frac = 0.0;
};

struct DirectionTimingRecord
{
    double direction_elapsed_seconds = 0.0;
    double first_order_seconds = 0.0;
    std::size_t first_order_view_samples = 0;
    int first_order_steps = 0;
    int solar_disk_nodes = 0;
    double second_order_seconds = 0.0;
    double second_order_incoming_single_scatter_seconds = 0.0;
    std::size_t second_order_incoming_single_scatter_calls = 0;
    std::size_t second_order_nonzero_incoming_calls = 0;
    std::size_t second_order_view_samples = 0;
    std::size_t second_order_mu_phi_evaluations = 0;
    double higher_order_seconds = 0.0;
    bool valid = false;
};

struct PartialMeasurementRow
{
    std::size_t index = 0;
    double zenith_deg = 0.0;
    double relative_azimuth_deg = 0.0;
    double absolute_azimuth_deg = 0.0;
    double reference_intensity = 0.0;
    double model_intensity = 0.0;
    double reference_dop = 0.0;
    bool has_reference_dop = false;
    double model_dop = 0.0;
    double reference_aop_deg = 0.0;
    bool has_reference_aop = false;
    double model_aop_deg = 0.0;
    double signed_dolp_bias = 0.0;
    double first_frac = 0.0;
    double second_frac = 0.0;
    double higher_frac = 0.0;
    double second_rr_frac = 0.0;
    double second_ar_frac = 0.0;
    double second_ra_frac = 0.0;
    double second_aa_frac = 0.0;
    double model_q = 0.0;
    double model_u = 0.0;
    double model_v = 0.0;
    double direction_elapsed_seconds = 0.0;
    double first_order_seconds = 0.0;
    std::size_t first_order_view_samples = 0;
    int first_order_steps = 0;
    int solar_disk_nodes = 0;
    double second_order_seconds = 0.0;
    double second_order_incoming_single_scatter_seconds = 0.0;
    std::size_t second_order_incoming_single_scatter_calls = 0;
    std::size_t second_order_nonzero_incoming_calls = 0;
    std::size_t second_order_view_samples = 0;
    std::size_t second_order_mu_phi_evaluations = 0;
    double higher_order_seconds = 0.0;
};

struct RegionSummary
{
    std::string name;
    std::size_t count = 0;
    std::size_t dop_count = 0;
    std::size_t aop_count = 0;
    double median_dolp_abs = std::numeric_limits<double>::quiet_NaN();
    double p95_dolp_abs = std::numeric_limits<double>::quiet_NaN();
    double median_aop_deg = std::numeric_limits<double>::quiet_NaN();
    double p95_aop_deg = std::numeric_limits<double>::quiet_NaN();
    double mean_signed_dolp_bias = std::numeric_limits<double>::quiet_NaN();
    double mean_first_frac = std::numeric_limits<double>::quiet_NaN();
    double mean_second_frac = std::numeric_limits<double>::quiet_NaN();
    double mean_higher_frac = std::numeric_limits<double>::quiet_NaN();
    double mean_second_rr_frac = std::numeric_limits<double>::quiet_NaN();
    double mean_second_ar_frac = std::numeric_limits<double>::quiet_NaN();
    double mean_second_ra_frac = std::numeric_limits<double>::quiet_NaN();
    double mean_second_aa_frac = std::numeric_limits<double>::quiet_NaN();
    double mean_normalized_reference = std::numeric_limits<double>::quiet_NaN();
    double mean_normalized_model = std::numeric_limits<double>::quiet_NaN();
    double max_dolp_abs = std::numeric_limits<double>::quiet_NaN();
    double max_dolp_zenith_deg = std::numeric_limits<double>::quiet_NaN();
    double max_dolp_relative_azimuth_deg = std::numeric_limits<double>::quiet_NaN();
};

struct RunnerOptions
{
    std::filesystem::path config_path;
    std::filesystem::path checkpoint_state_path;
    int higher_order_block_size = 0;
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

RunnerOptions parseRunnerOptions(int argc, char **argv)
{
    if (argc < 2) {
        throw std::runtime_error(
            "Usage: MeasurementCaseRunner <measurement_case_config> "
            "[--checkpoint-state <path>] [--higher-order-block-size <samples>]"
        );
    }

    RunnerOptions options;
    options.config_path = argv[1];
    for (int index = 2; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--checkpoint-state") {
            if (index + 1 >= argc) {
                throw std::runtime_error("Missing value for --checkpoint-state");
            }
            options.checkpoint_state_path = argv[++index];
            continue;
        }
        if (argument == "--higher-order-block-size") {
            if (index + 1 >= argc) {
                throw std::runtime_error("Missing value for --higher-order-block-size");
            }
            options.higher_order_block_size = std::stoi(argv[++index]);
            continue;
        }
        throw std::runtime_error("Unknown argument: " + argument);
    }

    return options;
}

std::map<std::string, std::string> loadKeyValueFile(const std::filesystem::path &path)
{
    std::map<std::string, std::string> values;
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Unable to open key-value file: " + path.string());
    }

    std::string line;
    while (std::getline(input, line)) {
        const std::string stripped = trim(line);
        if (stripped.empty() || stripped[0] == '#') {
            continue;
        }
        const std::size_t separator = stripped.find('=');
        if (separator == std::string::npos) {
            continue;
        }
        values[trim(stripped.substr(0, separator))] = trim(stripped.substr(separator + 1));
    }
    return values;
}

template <typename T>
T parseNumericValue(const std::map<std::string, std::string> &values, const std::string &key, T defaultValue)
{
    const auto iterator = values.find(key);
    if (iterator == values.end() || iterator->second.empty()) {
        return defaultValue;
    }
    if constexpr (std::is_integral_v<T>) {
        return static_cast<T>(std::stoll(iterator->second));
    }
    return static_cast<T>(std::stod(iterator->second));
}

StokesVector loadCheckpointStokes(
    const std::map<std::string, std::string> &values,
    const std::string &prefix
)
{
    return {
        parseNumericValue(values, prefix + "_I", 0.0),
        parseNumericValue(values, prefix + "_Q", 0.0),
        parseNumericValue(values, prefix + "_U", 0.0),
        parseNumericValue(values, prefix + "_V", 0.0),
    };
}

void writeCheckpointStokes(std::ofstream &output, const std::string &prefix, const StokesVector &vector)
{
    output << prefix << "_I=" << vector.I << "\n";
    output << prefix << "_Q=" << vector.Q << "\n";
    output << prefix << "_U=" << vector.U << "\n";
    output << prefix << "_V=" << vector.V << "\n";
}

bool loadDirectionCheckpointState(
    const std::filesystem::path &path,
    const std::string &caseId,
    std::size_t directionIndex,
    DirectionCheckpointState &state
)
{
    if (path.empty() || !std::filesystem::exists(path)) {
        return false;
    }

    const std::map<std::string, std::string> values = loadKeyValueFile(path);
    const std::string storedCaseId = values.count("case_id") ? values.at("case_id") : "";
    if (!storedCaseId.empty() && storedCaseId != caseId) {
        throw std::runtime_error("Checkpoint case_id does not match current run.");
    }

    const std::size_t storedDirectionIndex =
        static_cast<std::size_t>(parseNumericValue(values, "direction_index", static_cast<long long>(directionIndex)));
    if (storedDirectionIndex != directionIndex) {
        throw std::runtime_error("Checkpoint direction_index does not match current run.");
    }

    state.has_first_order = parseNumericValue(values, "has_first_order", 0) != 0;
    state.has_second_order = parseNumericValue(values, "has_second_order", 0) != 0;
    state.first_order = loadCheckpointStokes(values, "first_order");
    state.second_total = loadCheckpointStokes(values, "second_total");
    state.second_rr = loadCheckpointStokes(values, "second_rr");
    state.second_ar = loadCheckpointStokes(values, "second_ar");
    state.second_ra = loadCheckpointStokes(values, "second_ra");
    state.second_aa = loadCheckpointStokes(values, "second_aa");
    state.higher_order.completed_samples = parseNumericValue(values, "higher_completed_samples", 0);
    state.higher_order.mean = loadCheckpointStokes(values, "higher_mean");
    state.higher_order.m2 = loadCheckpointStokes(values, "higher_m2");
    state.timing.first_order_seconds = parseNumericValue(values, "timing_first_order_seconds", 0.0);
    state.timing.first_order_view_samples =
        static_cast<std::size_t>(parseNumericValue(values, "timing_first_order_view_samples", 0.0));
    state.timing.first_order_steps = parseNumericValue(values, "timing_first_order_steps", 0);
    state.timing.solar_disk_nodes = parseNumericValue(values, "timing_solar_disk_nodes", 0);
    state.timing.second_order_seconds = parseNumericValue(values, "timing_second_order_seconds", 0.0);
    state.timing.second_order_incoming_single_scatter_seconds =
        parseNumericValue(values, "timing_second_order_incoming_single_scatter_seconds", 0.0);
    state.timing.second_order_incoming_single_scatter_calls =
        static_cast<std::size_t>(parseNumericValue(values, "timing_second_order_incoming_single_scatter_calls", 0.0));
    state.timing.second_order_nonzero_incoming_calls =
        static_cast<std::size_t>(parseNumericValue(values, "timing_second_order_nonzero_incoming_calls", 0.0));
    state.timing.second_order_view_samples =
        static_cast<std::size_t>(parseNumericValue(values, "timing_second_order_view_samples", 0.0));
    state.timing.second_order_mu_phi_evaluations =
        static_cast<std::size_t>(parseNumericValue(values, "timing_second_order_mu_phi_evaluations", 0.0));
    state.timing.second_order_view_steps = parseNumericValue(values, "timing_second_order_view_steps", 0);
    state.timing.second_order_ray_steps = parseNumericValue(values, "timing_second_order_ray_steps", 0);
    state.timing.second_order_mu_nodes = parseNumericValue(values, "timing_second_order_mu_nodes", 0);
    state.timing.second_order_phi_nodes = parseNumericValue(values, "timing_second_order_phi_nodes", 0);
    state.timing.higher_order_seconds = parseNumericValue(values, "timing_higher_order_seconds", 0.0);
    state.timing.spectral_band_count =
        static_cast<std::size_t>(parseNumericValue(values, "timing_spectral_band_count", 0.0));
    state.timing.mc_sample_count = parseNumericValue(values, "timing_mc_sample_count", 0);
    return true;
}

void writeDirectionCheckpointState(
    const std::filesystem::path &path,
    const std::string &caseId,
    std::size_t directionIndex,
    const SkyDirection &direction,
    const DirectionCheckpointState &state
)
{
    if (path.empty()) {
        return;
    }

    std::filesystem::create_directories(path.parent_path());
    const std::filesystem::path tempPath(path.string() + ".tmp");
    std::ofstream output(tempPath);
    if (!output) {
        throw std::runtime_error("Unable to open checkpoint file for write: " + tempPath.string());
    }

    output << std::scientific << std::setprecision(17);
    output << "case_id=" << caseId << "\n";
    output << "direction_index=" << directionIndex << "\n";
    output << "zenith_deg=" << direction.zenith_deg << "\n";
    output << "azimuth_deg=" << direction.azimuth_deg << "\n";
    output << "has_first_order=" << (state.has_first_order ? 1 : 0) << "\n";
    output << "has_second_order=" << (state.has_second_order ? 1 : 0) << "\n";
    writeCheckpointStokes(output, "first_order", state.first_order);
    writeCheckpointStokes(output, "second_total", state.second_total);
    writeCheckpointStokes(output, "second_rr", state.second_rr);
    writeCheckpointStokes(output, "second_ar", state.second_ar);
    writeCheckpointStokes(output, "second_ra", state.second_ra);
    writeCheckpointStokes(output, "second_aa", state.second_aa);
    output << "higher_completed_samples=" << state.higher_order.completed_samples << "\n";
    writeCheckpointStokes(output, "higher_mean", state.higher_order.mean);
    writeCheckpointStokes(output, "higher_m2", state.higher_order.m2);
    output << "timing_first_order_seconds=" << state.timing.first_order_seconds << "\n";
    output << "timing_first_order_view_samples=" << state.timing.first_order_view_samples << "\n";
    output << "timing_first_order_steps=" << state.timing.first_order_steps << "\n";
    output << "timing_solar_disk_nodes=" << state.timing.solar_disk_nodes << "\n";
    output << "timing_second_order_seconds=" << state.timing.second_order_seconds << "\n";
    output << "timing_second_order_incoming_single_scatter_seconds="
           << state.timing.second_order_incoming_single_scatter_seconds << "\n";
    output << "timing_second_order_incoming_single_scatter_calls="
           << state.timing.second_order_incoming_single_scatter_calls << "\n";
    output << "timing_second_order_nonzero_incoming_calls="
           << state.timing.second_order_nonzero_incoming_calls << "\n";
    output << "timing_second_order_view_samples=" << state.timing.second_order_view_samples << "\n";
    output << "timing_second_order_mu_phi_evaluations=" << state.timing.second_order_mu_phi_evaluations << "\n";
    output << "timing_second_order_view_steps=" << state.timing.second_order_view_steps << "\n";
    output << "timing_second_order_ray_steps=" << state.timing.second_order_ray_steps << "\n";
    output << "timing_second_order_mu_nodes=" << state.timing.second_order_mu_nodes << "\n";
    output << "timing_second_order_phi_nodes=" << state.timing.second_order_phi_nodes << "\n";
    output << "timing_higher_order_seconds=" << state.timing.higher_order_seconds << "\n";
    output << "timing_spectral_band_count=" << state.timing.spectral_band_count << "\n";
    output << "timing_mc_sample_count=" << state.timing.mc_sample_count << "\n";
    output.close();

    std::error_code error;
    std::filesystem::remove(path, error);
    std::filesystem::rename(tempPath, path, error);
    if (error) {
        throw std::runtime_error("Unable to finalize checkpoint file: " + path.string());
    }
}

double relativeAzimuthFromAbsolute(double absoluteAzimuthDeg, double sunAzimuthDeg)
{
    return normalizeAzimuthDeg(absoluteAzimuthDeg - sunAzimuthDeg);
}

const char *sampleProgressStageName(SampleProgress::Stage stage)
{
    switch (stage) {
        case SampleProgress::Stage::direction_started:
            return "direction_started";
        case SampleProgress::Stage::first_order_band_complete:
            return "first_order_band_complete";
        case SampleProgress::Stage::first_order_complete:
            return "first_order_complete";
        case SampleProgress::Stage::second_order_band_complete:
            return "second_order_band_complete";
        case SampleProgress::Stage::second_order_complete:
            return "second_order_complete";
        case SampleProgress::Stage::higher_order_progress:
            return "higher_order_progress";
        case SampleProgress::Stage::higher_order_complete:
            return "higher_order_complete";
        case SampleProgress::Stage::direction_complete:
            return "direction_complete";
    }
    return "unknown";
}

std::vector<ReferencePoint> loadReferenceCsv(const std::string &filename)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("Unable to open reference CSV: " + filename);
    }

    std::string headerLine;
    if (!std::getline(input, headerLine)) {
        throw std::runtime_error("Reference CSV is empty: " + filename);
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
        point.original_index = static_cast<std::size_t>(points.size());
        if (columnIndex.count("original_index")) {
            point.original_index = static_cast<std::size_t>(getValue(values, "original_index", static_cast<double>(point.original_index)));
        } else if (columnIndex.count("index")) {
            point.original_index = static_cast<std::size_t>(getValue(values, "index", static_cast<double>(point.original_index)));
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

        const auto intensityIterator = columnIndex.find("intensity");
        if (intensityIterator != columnIndex.end() &&
            intensityIterator->second < values.size() &&
            !values[intensityIterator->second].empty()) {
            point.intensity = std::stod(values[intensityIterator->second]);
            point.has_intensity = true;
        }

        if (columnIndex.count("dop")) {
            point.dop = getValue(values, "dop", 0.0);
            point.has_dop = true;
        } else {
            point.q = getValue(values, "q", 0.0);
            point.u = getValue(values, "u", 0.0);
            if (point.has_intensity && point.intensity > 0.0) {
                point.dop = std::sqrt(point.q * point.q + point.u * point.u) / point.intensity;
                point.has_dop = true;
            }
        }

        if (columnIndex.count("q")) {
            point.q = getValue(values, "q", 0.0);
        }
        if (columnIndex.count("u")) {
            point.u = getValue(values, "u", 0.0);
        }
        if (columnIndex.count("aop_deg")) {
            point.aop_deg = getValue(values, "aop_deg", 0.0);
            point.has_aop = true;
        }

        points.push_back(point);
    }

    if (points.empty()) {
        throw std::runtime_error("Reference CSV contained no usable rows: " + filename);
    }
    return points;
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

void writePartialRowsHeader(std::ofstream &output)
{
    output << "index,zenith_deg,relative_azimuth_deg,absolute_azimuth_deg,"
           << "reference_intensity,model_intensity,has_reference_dop,reference_dop,model_dop,"
           << "has_reference_aop,reference_aop_deg,model_aop_deg,signed_dolp_bias,"
           << "first_frac,second_frac,higher_frac,"
           << "second_rr_frac,second_ar_frac,second_ra_frac,second_aa_frac,"
           << "model_q,model_u,model_v,"
           << "direction_elapsed_seconds,first_order_seconds,first_order_view_samples,first_order_steps,solar_disk_nodes,"
           << "second_order_seconds,second_order_incoming_single_scatter_seconds,"
           << "second_order_incoming_single_scatter_calls,second_order_nonzero_incoming_calls,"
           << "second_order_view_samples,second_order_mu_phi_evaluations,higher_order_seconds\n";
    output << std::scientific << std::setprecision(10);
}

PartialMeasurementRow buildPartialMeasurementRow(
    std::size_t index,
    const ReferencePoint &reference,
    const SkyDirection &direction,
    const SkyBinResult &model,
    const DirectionTimingRecord &timing,
    double sunAzimuthDeg
)
{
    const double relativeAzimuthDeg = reference.uses_relative_azimuth
        ? reference.relative_azimuth_deg
        : relativeAzimuthFromAbsolute(direction.azimuth_deg, sunAzimuthDeg);
    const double modelDolp = degreeOfLinearPolarization(model.mean);
    const double modelAopDeg = angleOfLinearPolarizationRad(model.mean) * 180.0 / PI;
    const double signedDolpBias = reference.has_dop ? modelDolp - reference.dop : 0.0;
    const double totalModelI = std::max(1.0e-12, model.mean.I);
    const double firstFrac = model.first_order.I / totalModelI;
    const double secondFrac = model.second_order.I / totalModelI;
    const double higherFrac = model.higher_order.I / totalModelI;
    const double secondTotalI = std::max(1.0e-12, model.second_order.I);

    return {
        index,
        reference.zenith_deg,
        relativeAzimuthDeg,
        direction.azimuth_deg,
        reference.has_intensity ? reference.intensity : 0.0,
        model.mean.I,
        reference.dop,
        reference.has_dop,
        modelDolp,
        reference.aop_deg,
        reference.has_aop,
        modelAopDeg,
        signedDolpBias,
        firstFrac,
        secondFrac,
        higherFrac,
        model.second_rr.I / secondTotalI,
        model.second_ar.I / secondTotalI,
        model.second_ra.I / secondTotalI,
        model.second_aa.I / secondTotalI,
        model.mean.Q,
        model.mean.U,
        model.mean.V,
        timing.direction_elapsed_seconds,
        timing.first_order_seconds,
        timing.first_order_view_samples,
        timing.first_order_steps,
        timing.solar_disk_nodes,
        timing.second_order_seconds,
        timing.second_order_incoming_single_scatter_seconds,
        timing.second_order_incoming_single_scatter_calls,
        timing.second_order_nonzero_incoming_calls,
        timing.second_order_view_samples,
        timing.second_order_mu_phi_evaluations,
        timing.higher_order_seconds,
    };
}

void appendPartialMeasurementRow(std::ofstream &output, const PartialMeasurementRow &row)
{
    output << row.index << ","
           << row.zenith_deg << ","
           << row.relative_azimuth_deg << ","
           << row.absolute_azimuth_deg << ","
           << row.reference_intensity << ","
           << row.model_intensity << ","
           << (row.has_reference_dop ? 1 : 0) << ","
           << row.reference_dop << ","
           << row.model_dop << ","
           << (row.has_reference_aop ? 1 : 0) << ","
           << row.reference_aop_deg << ","
           << row.model_aop_deg << ","
           << row.signed_dolp_bias << ","
           << row.first_frac << ","
           << row.second_frac << ","
           << row.higher_frac << ","
           << row.second_rr_frac << ","
           << row.second_ar_frac << ","
           << row.second_ra_frac << ","
           << row.second_aa_frac << ","
           << row.model_q << ","
           << row.model_u << ","
           << row.model_v << ","
           << row.direction_elapsed_seconds << ","
           << row.first_order_seconds << ","
           << row.first_order_view_samples << ","
           << row.first_order_steps << ","
           << row.solar_disk_nodes << ","
           << row.second_order_seconds << ","
           << row.second_order_incoming_single_scatter_seconds << ","
           << row.second_order_incoming_single_scatter_calls << ","
           << row.second_order_nonzero_incoming_calls << ","
           << row.second_order_view_samples << ","
           << row.second_order_mu_phi_evaluations << ","
           << row.higher_order_seconds << "\n";
}

std::map<std::size_t, PartialMeasurementRow> loadPartialMeasurementRows(const std::filesystem::path &path)
{
    std::map<std::size_t, PartialMeasurementRow> rows;
    if (!std::filesystem::exists(path) || std::filesystem::file_size(path) == 0) {
        return rows;
    }

    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Unable to open partial measurement CSV: " + path.string());
    }

    std::string headerLine;
    if (!std::getline(input, headerLine)) {
        return rows;
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

    std::string line;
    while (std::getline(input, line)) {
        const std::string stripped = trim(line);
        if (stripped.empty()) {
            continue;
        }

        const std::vector<std::string> values = splitCsvLine(stripped);
        PartialMeasurementRow row {};
        row.index = static_cast<std::size_t>(getValue(values, "index", 0.0));
        row.zenith_deg = getValue(values, "zenith_deg", 0.0);
        row.relative_azimuth_deg = getValue(values, "relative_azimuth_deg", 0.0);
        row.absolute_azimuth_deg = getValue(values, "absolute_azimuth_deg", 0.0);
        row.reference_intensity = getValue(values, "reference_intensity", 0.0);
        row.model_intensity = getValue(values, "model_intensity", 0.0);
        row.has_reference_dop = getValue(values, "has_reference_dop", 0.0) != 0.0;
        row.reference_dop = getValue(values, "reference_dop", 0.0);
        row.model_dop = getValue(values, "model_dop", 0.0);
        row.has_reference_aop = getValue(values, "has_reference_aop", 0.0) != 0.0;
        row.reference_aop_deg = getValue(values, "reference_aop_deg", 0.0);
        row.model_aop_deg = getValue(values, "model_aop_deg", 0.0);
        row.signed_dolp_bias = getValue(values, "signed_dolp_bias", 0.0);
        row.first_frac = getValue(values, "first_frac", 0.0);
        row.second_frac = getValue(values, "second_frac", 0.0);
        row.higher_frac = getValue(values, "higher_frac", 0.0);
        row.second_rr_frac = getValue(values, "second_rr_frac", 0.0);
        row.second_ar_frac = getValue(values, "second_ar_frac", 0.0);
        row.second_ra_frac = getValue(values, "second_ra_frac", 0.0);
        row.second_aa_frac = getValue(values, "second_aa_frac", 0.0);
        row.model_q = getValue(values, "model_q", 0.0);
        row.model_u = getValue(values, "model_u", 0.0);
        row.model_v = getValue(values, "model_v", 0.0);
        row.direction_elapsed_seconds = getValue(values, "direction_elapsed_seconds", 0.0);
        row.first_order_seconds = getValue(values, "first_order_seconds", 0.0);
        row.first_order_view_samples = static_cast<std::size_t>(getValue(values, "first_order_view_samples", 0.0));
        row.first_order_steps = static_cast<int>(getValue(values, "first_order_steps", 0.0));
        row.solar_disk_nodes = static_cast<int>(getValue(values, "solar_disk_nodes", 0.0));
        row.second_order_seconds = getValue(values, "second_order_seconds", 0.0);
        row.second_order_incoming_single_scatter_seconds =
            getValue(values, "second_order_incoming_single_scatter_seconds", 0.0);
        row.second_order_incoming_single_scatter_calls =
            static_cast<std::size_t>(getValue(values, "second_order_incoming_single_scatter_calls", 0.0));
        row.second_order_nonzero_incoming_calls =
            static_cast<std::size_t>(getValue(values, "second_order_nonzero_incoming_calls", 0.0));
        row.second_order_view_samples =
            static_cast<std::size_t>(getValue(values, "second_order_view_samples", 0.0));
        row.second_order_mu_phi_evaluations =
            static_cast<std::size_t>(getValue(values, "second_order_mu_phi_evaluations", 0.0));
        row.higher_order_seconds = getValue(values, "higher_order_seconds", 0.0);
        rows[row.index] = row;
    }

    return rows;
}

void writeMeasurementProgress(
    const std::filesystem::path &path,
    const std::string &caseId,
    std::size_t completedPoints,
    std::size_t totalPoints,
    std::size_t batchSize,
    const std::filesystem::path &partialRowsPath,
    const std::filesystem::path &comparisonCsvPath,
    const std::filesystem::path &reportPath
)
{
    std::ofstream output(path);
    output << "case_id=" << caseId << "\n";
    output << "completed_points=" << completedPoints << "\n";
    output << "total_points=" << totalPoints << "\n";
    output << "remaining_points=" << (totalPoints - completedPoints) << "\n";
    output << "batch_size=" << batchSize << "\n";
    output << "partial_rows_csv=" << partialRowsPath.string() << "\n";
    output << "comparison_csv=" << comparisonCsvPath.string() << "\n";
    output << "report_txt=" << reportPath.string() << "\n";
    output << "status=" << (completedPoints == totalPoints ? "complete" : "in_progress") << "\n";
}

void writeSingleDirectionComparisonOutputs(
    const SimulationConfig &config,
    const ReferencePoint &reference,
    const SkyDirection &direction,
    const SkyBinResult &model,
    const DirectionCheckpointState &checkpointState,
    const std::filesystem::path &reportDir
)
{
    const std::filesystem::path comparisonCsvPath =
        reportDir / (config.output.case_id + "_comparison.csv");
    const std::filesystem::path reportPath =
        reportDir / (config.output.case_id + ".txt");
    const double relativeAzimuthDeg = reference.uses_relative_azimuth
        ? reference.relative_azimuth_deg
        : relativeAzimuthFromAbsolute(direction.azimuth_deg, config.solar.azimuth_deg);
    const double modelDolp = degreeOfLinearPolarization(model.mean);
    const double modelAopDeg = angleOfLinearPolarizationRad(model.mean) * 180.0 / PI;
    const double dopAbsError = reference.has_dop ? std::abs(modelDolp - reference.dop) : 0.0;
    const double aopAbsErrorDeg = reference.has_aop
        ? angularDifferenceAopDeg(modelAopDeg, reference.aop_deg)
        : 0.0;
    const double signedDolpBias = reference.has_dop ? modelDolp - reference.dop : 0.0;
    const double totalModelI = std::max(1.0e-12, model.mean.I);
    const double firstFrac = model.first_order.I / totalModelI;
    const double secondFrac = model.second_order.I / totalModelI;
    const double higherFrac = model.higher_order.I / totalModelI;
    const double secondTotalI = std::max(1.0e-12, model.second_order.I);

    std::ofstream comparisonCsv(comparisonCsvPath);
    comparisonCsv << "index,zenith_deg,relative_azimuth_deg,absolute_azimuth_deg,"
                  << "reference_intensity,model_intensity,normalized_reference,normalized_model,"
                  << "reference_dop,model_dop,dop_abs_error,reference_aop_deg,model_aop_deg,aop_abs_error_deg,"
                  << "signed_dolp_bias,first_frac,second_frac,higher_frac,"
                  << "second_rr_frac,second_ar_frac,second_ra_frac,second_aa_frac,"
                  << "model_q,model_u,model_v\n";
    comparisonCsv << std::scientific << std::setprecision(10);
    comparisonCsv << reference.original_index << ","
                  << reference.zenith_deg << ","
                  << relativeAzimuthDeg << ","
                  << direction.azimuth_deg << ","
                  << (reference.has_intensity ? reference.intensity : 0.0) << ","
                  << model.mean.I << ","
                  << (reference.has_intensity ? 1.0 : 0.0) << ","
                  << 1.0 << ","
                  << (reference.has_dop ? reference.dop : 0.0) << ","
                  << modelDolp << ","
                  << dopAbsError << ","
                  << (reference.has_aop ? reference.aop_deg : 0.0) << ","
                  << modelAopDeg << ","
                  << aopAbsErrorDeg << ","
                  << signedDolpBias << ","
                  << firstFrac << ","
                  << secondFrac << ","
                  << higherFrac << ","
                  << (model.second_rr.I / secondTotalI) << ","
                  << (model.second_ar.I / secondTotalI) << ","
                  << (model.second_ra.I / secondTotalI) << ","
                  << (model.second_aa.I / secondTotalI) << ","
                  << model.mean.Q << ","
                  << model.mean.U << ","
                  << model.mean.V << "\n";

    std::ofstream report(reportPath);
    report << "case_id=" << config.output.case_id << "\n";
    report << "reference_points=1\n";
    report << "comparison_csv=" << comparisonCsvPath.string() << "\n";
    report << "completed_samples=" << checkpointState.higher_order.completed_samples << "\n";
    report << "normalized_rmse=" << (reference.has_intensity ? 0.0 : std::numeric_limits<double>::quiet_NaN()) << "\n";
    report << "median_dolp_abs=" << dopAbsError << "\n";
    report << "p95_dolp_abs=" << dopAbsError << "\n";
    report << "median_aop_deg=" << aopAbsErrorDeg << "\n";
    report << "p95_aop_deg=" << aopAbsErrorDeg << "\n";
    report << "solar_vertical_signed_dolp_bias=" << signedDolpBias << "\n";
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

bool inSolarVerticalMidzen(const PointDiagnostics &point)
{
    return point.relative_azimuth_deg >= 240.0 &&
        point.relative_azimuth_deg <= 300.0 &&
        point.zenith_deg >= 30.0 &&
        point.zenith_deg <= 70.0;
}

bool inSolarVerticalLowElevation(const PointDiagnostics &point)
{
    return point.relative_azimuth_deg >= 240.0 &&
        point.relative_azimuth_deg <= 300.0 &&
        point.zenith_deg >= 70.0 &&
        point.zenith_deg <= 88.0;
}

bool inAntisolarMidzen(const PointDiagnostics &point)
{
    return (point.relative_azimuth_deg <= 30.0 || point.relative_azimuth_deg >= 330.0) &&
        point.zenith_deg >= 30.0 &&
        point.zenith_deg <= 70.0;
}

bool inBrightHorizonArc(const PointDiagnostics &point)
{
    return point.relative_azimuth_deg >= 320.0 &&
        point.relative_azimuth_deg <= 350.0 &&
        point.zenith_deg >= 65.0 &&
        point.zenith_deg <= 88.0;
}

bool inNearZenith(const PointDiagnostics &point)
{
    return point.zenith_deg <= 25.0;
}

bool inAopFlipRelAz45Sector(const PointDiagnostics &point)
{
    return point.relative_azimuth_deg >= 30.0 &&
        point.relative_azimuth_deg <= 60.0 &&
        point.zenith_deg >= 5.0 &&
        point.zenith_deg <= 85.0;
}

bool inAopFlipRelAz215Sector(const PointDiagnostics &point)
{
    return point.relative_azimuth_deg >= 200.0 &&
        point.relative_azimuth_deg <= 230.0 &&
        point.zenith_deg >= 5.0 &&
        point.zenith_deg <= 85.0;
}

RegionSummary summarizeRegion(
    const std::string &name,
    const std::vector<PointDiagnostics> &points,
    bool (*predicate)(const PointDiagnostics &)
)
{
    RegionSummary summary {};
    summary.name = name;
    std::vector<double> dopErrors;
    std::vector<double> aopErrors;
    double sumSignedBias = 0.0;
    double sumFirst = 0.0;
    double sumSecond = 0.0;
    double sumHigher = 0.0;
    double sumSecondRr = 0.0;
    double sumSecondAr = 0.0;
    double sumSecondRa = 0.0;
    double sumSecondAa = 0.0;
    double sumNormRef = 0.0;
    double sumNormModel = 0.0;
    double worstDop = -1.0;

    for (const PointDiagnostics &point : points) {
        if (!predicate(point)) {
            continue;
        }
        summary.count += 1;
        sumFirst += point.first_frac;
        sumSecond += point.second_frac;
        sumHigher += point.higher_frac;
        sumSecondRr += point.second_rr_frac;
        sumSecondAr += point.second_ar_frac;
        sumSecondRa += point.second_ra_frac;
        sumSecondAa += point.second_aa_frac;
        sumNormRef += point.normalized_reference;
        sumNormModel += point.normalized_model;
        if (point.has_reference_dop) {
            dopErrors.push_back(point.dop_abs_error);
            sumSignedBias += point.signed_dolp_bias;
            if (point.dop_abs_error > worstDop) {
                worstDop = point.dop_abs_error;
                summary.max_dolp_abs = point.dop_abs_error;
                summary.max_dolp_zenith_deg = point.zenith_deg;
                summary.max_dolp_relative_azimuth_deg = point.relative_azimuth_deg;
            }
        }
        if (point.has_reference_aop && point.reference_dop >= 0.15) {
            aopErrors.push_back(point.aop_abs_error_deg);
        }
    }

    if (summary.count > 0) {
        summary.mean_first_frac = sumFirst / static_cast<double>(summary.count);
        summary.mean_second_frac = sumSecond / static_cast<double>(summary.count);
        summary.mean_higher_frac = sumHigher / static_cast<double>(summary.count);
        summary.mean_second_rr_frac = sumSecondRr / static_cast<double>(summary.count);
        summary.mean_second_ar_frac = sumSecondAr / static_cast<double>(summary.count);
        summary.mean_second_ra_frac = sumSecondRa / static_cast<double>(summary.count);
        summary.mean_second_aa_frac = sumSecondAa / static_cast<double>(summary.count);
        summary.mean_normalized_reference = sumNormRef / static_cast<double>(summary.count);
        summary.mean_normalized_model = sumNormModel / static_cast<double>(summary.count);
    }
    summary.dop_count = dopErrors.size();
    if (!dopErrors.empty()) {
        summary.median_dolp_abs = percentile(dopErrors, 0.5);
        summary.p95_dolp_abs = percentile(dopErrors, 0.95);
        summary.mean_signed_dolp_bias = sumSignedBias / static_cast<double>(dopErrors.size());
    }
    summary.aop_count = aopErrors.size();
    if (!aopErrors.empty()) {
        summary.median_aop_deg = percentile(aopErrors, 0.5);
        summary.p95_aop_deg = percentile(aopErrors, 0.95);
    }
    return summary;
}
}

int main(int argc, char **argv)
{
    try {
        const auto overallStart = std::chrono::steady_clock::now();
        std::cout.setf(std::ios::unitbuf);
        std::cerr.setf(std::ios::unitbuf);
        const RunnerOptions options = parseRunnerOptions(argc, argv);

        const auto configLoadStart = std::chrono::steady_clock::now();
        const SimulationConfig config = loadSimulationConfig(options.config_path.string());
        const double configLoadSeconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - configLoadStart).count();
        if (config.output.measurement_reference_csv.empty()) {
            throw std::runtime_error("measurement_reference_csv is not set in the supplied config.");
        }

        const auto referenceLoadStart = std::chrono::steady_clock::now();
        const std::vector<ReferencePoint> reference = loadReferenceCsv(config.output.measurement_reference_csv);
        const std::vector<SkyDirection> directions = absoluteDirections(reference, config.solar.azimuth_deg);
        std::vector<std::size_t> originalDirectionIndices;
        originalDirectionIndices.reserve(reference.size());
        for (const ReferencePoint &point : reference) {
            originalDirectionIndices.push_back(point.original_index);
        }
        const double referenceLoadSeconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - referenceLoadStart).count();

        if (!options.checkpoint_state_path.empty()) {
            if (reference.size() != 1 || directions.size() != 1 || originalDirectionIndices.size() != 1) {
                throw std::runtime_error("Checkpoint mode requires exactly one measurement reference direction.");
            }

            const std::filesystem::path reportDir =
                std::filesystem::path(config.output.output_dir) / "measurement_case_reports";
            std::filesystem::create_directories(reportDir);
            const std::filesystem::path comparisonCsvPath =
                reportDir / (config.output.case_id + "_comparison.csv");
            const std::filesystem::path reportPath =
                reportDir / (config.output.case_id + ".txt");
            const std::filesystem::path regionSummaryCsvPath =
                reportDir / (config.output.case_id + "_region_summary.csv");

            DirectionCheckpointState checkpointState {};
            loadDirectionCheckpointState(
                options.checkpoint_state_path,
                config.output.case_id,
                originalDirectionIndices.front(),
                checkpointState
            );

            const int totalTargetSamples = config.monte_carlo.photons_per_bin;
            const int currentCompletedSamples = checkpointState.higher_order.completed_samples;
            const int targetSamples =
                options.higher_order_block_size > 0
                    ? std::min(totalTargetSamples, currentCompletedSamples + options.higher_order_block_size)
                    : totalTargetSamples;

            std::cout << "[measurement-checkpoint] case_id=" << config.output.case_id
                      << " direction_index=" << originalDirectionIndices.front()
                      << " completed_samples=" << currentCompletedSamples
                      << " target_samples=" << targetSamples
                      << " total_target_samples=" << totalTargetSamples
                      << " checkpoint_state=" << options.checkpoint_state_path.string()
                      << "\n";

            DirectionTimingRecord timingRecord {};
            const SkyBinResult model = solveSkyDirection(
                config,
                directions.front(),
                originalDirectionIndices.front(),
                targetSamples,
                &checkpointState,
                [&](const SampleProgress &progress) {
                    timingRecord.direction_elapsed_seconds = progress.direction_elapsed_seconds;
                    timingRecord.first_order_seconds = progress.first_order_seconds;
                    timingRecord.first_order_view_samples = progress.first_order_view_samples;
                    timingRecord.first_order_steps = progress.first_order_steps;
                    timingRecord.solar_disk_nodes = progress.solar_disk_nodes;
                    timingRecord.second_order_seconds = progress.second_order_seconds;
                    timingRecord.second_order_incoming_single_scatter_seconds =
                        progress.second_order_incoming_single_scatter_seconds;
                    timingRecord.second_order_incoming_single_scatter_calls =
                        progress.second_order_incoming_single_scatter_calls;
                    timingRecord.second_order_nonzero_incoming_calls =
                        progress.second_order_nonzero_incoming_calls;
                    timingRecord.second_order_view_samples = progress.second_order_view_samples;
                    timingRecord.second_order_mu_phi_evaluations = progress.second_order_mu_phi_evaluations;
                    timingRecord.higher_order_seconds = progress.higher_order_seconds;
                    timingRecord.valid = true;

                    if (progress.stage == SampleProgress::Stage::higher_order_progress ||
                        progress.stage == SampleProgress::Stage::higher_order_complete ||
                        progress.stage == SampleProgress::Stage::direction_complete) {
                        const std::size_t completedSamples =
                            progress.stage == SampleProgress::Stage::direction_complete
                                ? static_cast<std::size_t>(checkpointState.higher_order.completed_samples)
                                : progress.stage_completed;
                        const std::size_t targetSamples =
                            progress.mc_sample_count > 0
                                ? static_cast<std::size_t>(progress.mc_sample_count)
                                : progress.stage_total;
                        std::cout << "[measurement-checkpoint-progress] stage="
                                  << sampleProgressStageName(progress.stage)
                                  << " direction_index=" << progress.direction_index
                                  << " completed_samples=" << completedSamples
                                  << " target_samples=" << targetSamples
                                  << " elapsed_s=" << progress.direction_elapsed_seconds
                                  << " first_s=" << progress.first_order_seconds
                                  << " second_s=" << progress.second_order_seconds
                                  << " higher_s=" << progress.higher_order_seconds
                                  << "\n";
                    }
                }
            );

            writeDirectionCheckpointState(
                options.checkpoint_state_path,
                config.output.case_id,
                originalDirectionIndices.front(),
                directions.front(),
                checkpointState
            );

            const bool complete = config.monte_carlo.single_scatter_only ||
                checkpointState.higher_order.completed_samples >= totalTargetSamples;
            if (!complete) {
                std::filesystem::remove(comparisonCsvPath);
                std::filesystem::remove(reportPath);
                std::filesystem::remove(regionSummaryCsvPath);
                std::cout << "[measurement-checkpoint] status=in_progress"
                          << " completed_samples=" << checkpointState.higher_order.completed_samples
                          << " remaining_samples=" << (totalTargetSamples - checkpointState.higher_order.completed_samples)
                          << "\n";
                return 0;
            }

            writeSingleDirectionComparisonOutputs(
                config,
                reference.front(),
                directions.front(),
                model,
                checkpointState,
                reportDir
            );
            std::filesystem::remove(regionSummaryCsvPath);
            std::cout << "[measurement-checkpoint] status=complete"
                      << " completed_samples=" << checkpointState.higher_order.completed_samples
                      << " comparison_csv=" << comparisonCsvPath.string()
                      << "\n";
            return 0;
        }

        std::cout << "[measurement-start] case_id=" << config.output.case_id
                  << " reference_points=" << reference.size()
                  << " photons_per_bin=" << config.monte_carlo.photons_per_bin
                  << " deterministic_second_scatter=" << (config.monte_carlo.deterministic_second_scatter ? 1 : 0)
                  << " wavelength_nm=[" << config.spectral.min_wavelength_nm << "," << config.spectral.max_wavelength_nm
                  << "] step=" << config.spectral.wavelength_step_nm
                  << " solar_disk_nodes=" << config.solar.solar_disk_quadrature_nodes
                  << "\n";

        double cumulativeDirectionSeconds = 0.0;
        double cumulativeFirstOrderSeconds = 0.0;
        std::size_t cumulativeFirstOrderViewSamples = 0;
        double cumulativeSecondOrderSeconds = 0.0;
        double cumulativeSecondOrderIncomingSeconds = 0.0;
        double cumulativeHigherOrderSeconds = 0.0;
        double maxDirectionSeconds = 0.0;
        double maxSecondOrderSeconds = 0.0;
        std::size_t slowestDirectionIndex = 0;
        std::size_t slowestSecondOrderDirectionIndex = 0;
        double slowestDirectionZenithDeg = 0.0;
        double slowestDirectionAzimuthDeg = 0.0;
        double slowestSecondOrderZenithDeg = 0.0;
        double slowestSecondOrderAzimuthDeg = 0.0;
        std::size_t cumulativeSecondOrderIncomingCalls = 0;
        std::size_t cumulativeSecondOrderNonzeroCalls = 0;
        std::size_t cumulativeSecondOrderViewSamples = 0;
        std::size_t cumulativeSecondOrderMuPhiEvaluations = 0;
        const std::size_t progressEvery =
            reference.size() <= 16 ? 1 : (reference.size() <= 64 ? 4 : 16);

        const auto samplingStart = std::chrono::steady_clock::now();
        const std::vector<SkyBinResult> model = sampleSkyDirections(
            config,
            directions,
            -1,
            [&](const SampleProgress &progress) {
                if (progress.stage == SampleProgress::Stage::first_order_band_complete) {
                    const bool shouldPrintBand =
                        progress.stage_completed == 1 ||
                        progress.stage_completed == progress.stage_total ||
                        progress.stage_completed % (progress.stage_total <= 12 ? 1 : 4) == 0;
                    if (shouldPrintBand) {
                        std::cout << "[measurement-stage] stage=" << sampleProgressStageName(progress.stage)
                                  << " direction_index=" << progress.direction_index
                                  << " band=" << progress.stage_completed << "/" << progress.stage_total
                                  << " direction_elapsed_s=" << progress.direction_elapsed_seconds
                                  << " first_order_s=" << progress.first_order_seconds
                                  << " first_order_view_samples=" << progress.first_order_view_samples
                                  << " first_order_steps=" << progress.first_order_steps
                                  << " solar_disk_nodes=" << progress.solar_disk_nodes
                                  << " zenith_deg=" << progress.zenith_deg
                                  << " azimuth_deg=" << progress.azimuth_deg
                                  << "\n";
                    }
                    return;
                }

                if (progress.stage == SampleProgress::Stage::second_order_band_complete) {
                    const bool shouldPrintBand =
                        progress.stage_completed == 1 ||
                        progress.stage_completed == progress.stage_total ||
                        progress.stage_completed % (progress.stage_total <= 12 ? 1 : 4) == 0;
                    if (shouldPrintBand) {
                        std::cout << "[measurement-stage] stage=" << sampleProgressStageName(progress.stage)
                                  << " direction_index=" << progress.direction_index
                                  << " band=" << progress.stage_completed << "/" << progress.stage_total
                                  << " direction_elapsed_s=" << progress.direction_elapsed_seconds
                                  << " second_order_s=" << progress.second_order_seconds
                                  << " second_inner_s=" << progress.second_order_incoming_single_scatter_seconds
                                  << " second_calls=" << progress.second_order_incoming_single_scatter_calls
                                  << " second_nonzero=" << progress.second_order_nonzero_incoming_calls
                                  << " second_view_samples=" << progress.second_order_view_samples
                                  << " second_mu_phi=" << progress.second_order_mu_phi_evaluations
                                  << " quadrature=" << progress.second_order_view_steps
                                  << "x" << progress.second_order_ray_steps
                                  << "x" << progress.second_order_mu_nodes
                                  << "x" << progress.second_order_phi_nodes
                                  << " zenith_deg=" << progress.zenith_deg
                                  << " azimuth_deg=" << progress.azimuth_deg
                                  << "\n";
                    }
                    return;
                }

                if (progress.stage == SampleProgress::Stage::higher_order_progress) {
                    std::cout << "[measurement-stage] stage=" << sampleProgressStageName(progress.stage)
                              << " direction_index=" << progress.direction_index
                              << " samples=" << progress.stage_completed << "/" << progress.stage_total
                              << " direction_elapsed_s=" << progress.direction_elapsed_seconds
                              << " higher_s=" << progress.higher_order_seconds
                              << " zenith_deg=" << progress.zenith_deg
                              << " azimuth_deg=" << progress.azimuth_deg
                              << "\n";
                    return;
                }

                if (progress.stage == SampleProgress::Stage::direction_complete) {
                    cumulativeDirectionSeconds += progress.direction_elapsed_seconds;
                    cumulativeFirstOrderSeconds += progress.first_order_seconds;
                    cumulativeFirstOrderViewSamples += progress.first_order_view_samples;
                    cumulativeSecondOrderSeconds += progress.second_order_seconds;
                    cumulativeSecondOrderIncomingSeconds += progress.second_order_incoming_single_scatter_seconds;
                    cumulativeHigherOrderSeconds += progress.higher_order_seconds;
                    cumulativeSecondOrderIncomingCalls += progress.second_order_incoming_single_scatter_calls;
                    cumulativeSecondOrderNonzeroCalls += progress.second_order_nonzero_incoming_calls;
                    cumulativeSecondOrderViewSamples += progress.second_order_view_samples;
                    cumulativeSecondOrderMuPhiEvaluations += progress.second_order_mu_phi_evaluations;

                    if (progress.direction_elapsed_seconds > maxDirectionSeconds) {
                        maxDirectionSeconds = progress.direction_elapsed_seconds;
                        slowestDirectionIndex = progress.direction_index;
                        slowestDirectionZenithDeg = progress.zenith_deg;
                        slowestDirectionAzimuthDeg = progress.azimuth_deg;
                    }
                    if (progress.second_order_seconds > maxSecondOrderSeconds) {
                        maxSecondOrderSeconds = progress.second_order_seconds;
                        slowestSecondOrderDirectionIndex = progress.direction_index;
                        slowestSecondOrderZenithDeg = progress.zenith_deg;
                        slowestSecondOrderAzimuthDeg = progress.azimuth_deg;
                    }

                    const bool shouldPrint =
                        progress.completed_count == 1 ||
                        progress.completed_count == progress.total_count ||
                        (progress.completed_count % progressEvery) == 0;
                    if (!shouldPrint) {
                        return;
                    }

                    const double meanDirectionSeconds =
                        cumulativeDirectionSeconds / static_cast<double>(progress.completed_count);
                    const double etaSeconds =
                        meanDirectionSeconds * static_cast<double>(progress.total_count - progress.completed_count);
                    std::cout << "[measurement-progress] completed="
                              << progress.completed_count << "/" << progress.total_count
                              << " last_direction_s=" << progress.direction_elapsed_seconds
                              << " mean_direction_s=" << meanDirectionSeconds
                              << " first_s=" << progress.first_order_seconds
                              << " first_view_samples=" << progress.first_order_view_samples
                              << " first_steps=" << progress.first_order_steps
                              << " second_s=" << progress.second_order_seconds
                              << " second_inner_s=" << progress.second_order_incoming_single_scatter_seconds
                              << " higher_s=" << progress.higher_order_seconds
                              << " second_calls=" << progress.second_order_incoming_single_scatter_calls
                              << " second_nonzero=" << progress.second_order_nonzero_incoming_calls
                              << " elapsed_s=" << progress.total_elapsed_seconds
                              << " eta_s=" << etaSeconds
                              << " zenith_deg=" << progress.zenith_deg
                              << " azimuth_deg=" << progress.azimuth_deg
                              << "\n";
                    return;
                }

                if (progress.stage == SampleProgress::Stage::first_order_complete ||
                    progress.stage == SampleProgress::Stage::second_order_complete ||
                    progress.stage == SampleProgress::Stage::higher_order_complete) {
                    if (progress.total_count <= 16) {
                        std::cout << "[measurement-stage] stage=" << sampleProgressStageName(progress.stage)
                                  << " direction_index=" << progress.direction_index
                                  << " direction_elapsed_s=" << progress.direction_elapsed_seconds
                                  << " first_s=" << progress.first_order_seconds
                                  << " first_view_samples=" << progress.first_order_view_samples
                                  << " first_steps=" << progress.first_order_steps
                                  << " second_s=" << progress.second_order_seconds
                                  << " second_inner_s=" << progress.second_order_incoming_single_scatter_seconds
                                  << " higher_s=" << progress.higher_order_seconds
                                  << " second_calls=" << progress.second_order_incoming_single_scatter_calls
                                  << " second_nonzero=" << progress.second_order_nonzero_incoming_calls
                                  << " zenith_deg=" << progress.zenith_deg
                                  << " azimuth_deg=" << progress.azimuth_deg
                                  << "\n";
                    }
                }
            },
            &originalDirectionIndices
        );
        const double samplingSeconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - samplingStart).count();
        const std::filesystem::path reportDir =
            std::filesystem::path(config.output.output_dir) / "measurement_case_reports";
        std::filesystem::create_directories(reportDir);

        const bool hasIntensityReference = std::any_of(
            reference.begin(),
            reference.end(),
            [](const ReferencePoint &point) { return point.has_intensity; }
        );

        std::vector<double> dopErrors;
        std::vector<double> aopErrorsDeg;
        std::vector<double> solarVerticalSignedDolpBias;
        double worstDopError = -1.0;
        std::size_t worstDopIndex = 0;
        double referencePeak = 0.0;
        double modelPeak = 0.0;
        double sumSquared = 0.0;
        double count = 0.0;
        std::size_t brightestReferenceIndex = 0;
        std::size_t brightestModelIndex = 0;
        double brightestReferenceValue = -1.0;
        double brightestModelValue = -1.0;
        const auto metricsStart = std::chrono::steady_clock::now();

        if (hasIntensityReference) {
            for (std::size_t index = 0; index < reference.size(); ++index) {
                if (!reference[index].has_intensity) {
                    continue;
                }
                referencePeak = std::max(referencePeak, reference[index].intensity);
                modelPeak = std::max(modelPeak, model[index].mean.I);
            }
        }

        const std::filesystem::path comparisonCsvPath =
            reportDir / (config.output.case_id + "_comparison.csv");
        const std::filesystem::path regionSummaryCsvPath =
            reportDir / (config.output.case_id + "_region_summary.csv");
        std::ofstream comparisonCsv(comparisonCsvPath);
        comparisonCsv << "index,zenith_deg,relative_azimuth_deg,absolute_azimuth_deg,"
                      << "reference_intensity,model_intensity,normalized_reference,normalized_model,"
                      << "reference_dop,model_dop,dop_abs_error,reference_aop_deg,model_aop_deg,aop_abs_error_deg,"
                      << "signed_dolp_bias,first_frac,second_frac,higher_frac,"
                      << "second_rr_frac,second_ar_frac,second_ra_frac,second_aa_frac,"
                      << "model_q,model_u,model_v\n";
        comparisonCsv << std::scientific << std::setprecision(10);
        std::vector<PointDiagnostics> diagnostics;
        diagnostics.reserve(reference.size());

        for (std::size_t index = 0; index < reference.size(); ++index) {
            const double normalizedReference = hasIntensityReference && reference[index].has_intensity
                ? reference[index].intensity / std::max(1.0e-12, referencePeak)
                : 0.0;
            const double normalizedModel = hasIntensityReference
                ? model[index].mean.I / std::max(1.0e-12, modelPeak)
                : 0.0;
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
            const double relativeAzimuthDeg = reference[index].uses_relative_azimuth
                ? reference[index].relative_azimuth_deg
                : relativeAzimuthFromAbsolute(directions[index].azimuth_deg, config.solar.azimuth_deg);
            const double totalModelI = std::max(1.0e-12, model[index].mean.I);
            const double firstFrac = model[index].first_order.I / totalModelI;
            const double secondFrac = model[index].second_order.I / totalModelI;
            const double higherFrac = model[index].higher_order.I / totalModelI;
            const double secondTotalI = std::max(1.0e-12, model[index].second_order.I);
            const double secondRrFrac = model[index].second_rr.I / secondTotalI;
            const double secondArFrac = model[index].second_ar.I / secondTotalI;
            const double secondRaFrac = model[index].second_ra.I / secondTotalI;
            const double secondAaFrac = model[index].second_aa.I / secondTotalI;

            diagnostics.push_back({
                reference[index].original_index,
                reference[index].zenith_deg,
                relativeAzimuthDeg,
                directions[index].azimuth_deg,
                normalizedReference,
                normalizedModel,
                reference[index].dop,
                reference[index].has_dop,
                modelDolp,
                dopAbsError,
                reference[index].aop_deg,
                reference[index].has_aop,
                modelAopDeg,
                aopAbsErrorDeg,
                signedDolpBias,
                firstFrac,
                secondFrac,
                higherFrac,
                secondRrFrac,
                secondArFrac,
                secondRaFrac,
                secondAaFrac,
            });

            if (reference[index].has_dop && dopAbsError > worstDopError) {
                worstDopError = dopAbsError;
                worstDopIndex = index;
            }

            if (hasIntensityReference && reference[index].has_intensity) {
                if (normalizedReference > brightestReferenceValue) {
                    brightestReferenceValue = normalizedReference;
                    brightestReferenceIndex = index;
                }
                if (normalizedModel > brightestModelValue) {
                    brightestModelValue = normalizedModel;
                    brightestModelIndex = index;
                }

                if (normalizedReference >= config.validation.measurement_mask_fraction_of_peak) {
                    const double diff = normalizedModel - normalizedReference;
                    sumSquared += diff * diff;
                    count += 1.0;
                    if (reference[index].has_dop) {
                        dopErrors.push_back(dopAbsError);
                        if (reference[index].zenith_deg >= 30.0 &&
                            reference[index].zenith_deg <= 70.0 &&
                            relativeAzimuthDeg >= 240.0 &&
                            relativeAzimuthDeg <= 300.0) {
                            solarVerticalSignedDolpBias.push_back(signedDolpBias);
                        }
                    }
                    if (reference[index].has_aop && reference[index].has_dop && reference[index].dop >= 0.15) {
                        aopErrorsDeg.push_back(aopAbsErrorDeg);
                    }
                }
            } else if (!hasIntensityReference && reference[index].has_dop) {
                dopErrors.push_back(dopAbsError);
                if (reference[index].zenith_deg >= 30.0 &&
                    reference[index].zenith_deg <= 70.0 &&
                    relativeAzimuthDeg >= 240.0 &&
                    relativeAzimuthDeg <= 300.0) {
                    solarVerticalSignedDolpBias.push_back(signedDolpBias);
                }
                if (reference[index].has_aop && reference[index].dop >= 0.15) {
                    aopErrorsDeg.push_back(aopAbsErrorDeg);
                }
            }

            comparisonCsv << reference[index].original_index << ","
                          << reference[index].zenith_deg << ","
                          << relativeAzimuthDeg << ","
                          << directions[index].azimuth_deg << ","
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
                          << firstFrac << ","
                          << secondFrac << ","
                          << higherFrac << ","
                          << secondRrFrac << ","
                          << secondArFrac << ","
                          << secondRaFrac << ","
                          << secondAaFrac << ","
                          << model[index].mean.Q << ","
                          << model[index].mean.U << ","
                          << model[index].mean.V << "\n";
        }

        const std::vector<RegionSummary> regionSummaries = {
            summarizeRegion("solar_vertical_midzen", diagnostics, inSolarVerticalMidzen),
            summarizeRegion("solar_vertical_lowelev", diagnostics, inSolarVerticalLowElevation),
            summarizeRegion("antisolar_midzen", diagnostics, inAntisolarMidzen),
            summarizeRegion("bright_horizon_arc", diagnostics, inBrightHorizonArc),
            summarizeRegion("near_zenith", diagnostics, inNearZenith),
            summarizeRegion("aop_flip_relaz_45_sector", diagnostics, inAopFlipRelAz45Sector),
            summarizeRegion("aop_flip_relaz_215_sector", diagnostics, inAopFlipRelAz215Sector),
        };

        std::ofstream regionSummaryCsv(regionSummaryCsvPath);
        regionSummaryCsv << "region_name,count,dop_count,aop_count,"
                         << "median_dolp_abs,p95_dolp_abs,median_aop_deg,p95_aop_deg,"
                         << "mean_signed_dolp_bias,mean_first_frac,mean_second_frac,mean_higher_frac,"
                         << "mean_second_rr_frac,mean_second_ar_frac,mean_second_ra_frac,mean_second_aa_frac,"
                         << "mean_normalized_reference,mean_normalized_model,"
                         << "max_dolp_abs,max_dolp_zenith_deg,max_dolp_relative_azimuth_deg\n";
        regionSummaryCsv << std::scientific << std::setprecision(10);
        for (const RegionSummary &summary : regionSummaries) {
            regionSummaryCsv << summary.name << ","
                             << summary.count << ","
                             << summary.dop_count << ","
                             << summary.aop_count << ","
                             << summary.median_dolp_abs << ","
                             << summary.p95_dolp_abs << ","
                             << summary.median_aop_deg << ","
                             << summary.p95_aop_deg << ","
                             << summary.mean_signed_dolp_bias << ","
                             << summary.mean_first_frac << ","
                             << summary.mean_second_frac << ","
                             << summary.mean_higher_frac << ","
                             << summary.mean_second_rr_frac << ","
                             << summary.mean_second_ar_frac << ","
                             << summary.mean_second_ra_frac << ","
                             << summary.mean_second_aa_frac << ","
                             << summary.mean_normalized_reference << ","
                             << summary.mean_normalized_model << ","
                             << summary.max_dolp_abs << ","
                             << summary.max_dolp_zenith_deg << ","
                             << summary.max_dolp_relative_azimuth_deg << "\n";
        }

        std::ostringstream report;
        report << "case_id=" << config.output.case_id << "\n";
        report << "reference_csv=" << config.output.measurement_reference_csv << "\n";
        report << "reference_points=" << reference.size() << "\n";
        report << "timing_config_load_seconds=" << configLoadSeconds << "\n";
        report << "timing_reference_load_seconds=" << referenceLoadSeconds << "\n";
        report << "timing_direction_sampling_seconds=" << samplingSeconds << "\n";
        report << "timing_progress_interval=" << progressEvery << "\n";
        report << "timing_mean_direction_seconds="
               << (reference.empty() ? 0.0 : cumulativeDirectionSeconds / static_cast<double>(reference.size()))
               << "\n";
        report << "timing_mean_first_order_seconds="
               << (reference.empty() ? 0.0 : cumulativeFirstOrderSeconds / static_cast<double>(reference.size()))
               << "\n";
        report << "timing_mean_first_order_view_samples="
               << (reference.empty() ? 0.0 : static_cast<double>(cumulativeFirstOrderViewSamples) / static_cast<double>(reference.size()))
               << "\n";
        report << "timing_mean_second_order_seconds="
               << (reference.empty() ? 0.0 : cumulativeSecondOrderSeconds / static_cast<double>(reference.size()))
               << "\n";
        report << "timing_mean_second_order_incoming_single_scatter_seconds="
               << (reference.empty() ? 0.0 : cumulativeSecondOrderIncomingSeconds / static_cast<double>(reference.size()))
               << "\n";
        report << "timing_mean_higher_order_seconds="
               << (reference.empty() ? 0.0 : cumulativeHigherOrderSeconds / static_cast<double>(reference.size()))
               << "\n";
        report << "timing_mean_second_order_incoming_single_scatter_calls="
               << (reference.empty() ? 0.0 : static_cast<double>(cumulativeSecondOrderIncomingCalls) / static_cast<double>(reference.size()))
               << "\n";
        report << "timing_mean_second_order_nonzero_incoming_calls="
               << (reference.empty() ? 0.0 : static_cast<double>(cumulativeSecondOrderNonzeroCalls) / static_cast<double>(reference.size()))
               << "\n";
        report << "timing_mean_second_order_view_samples="
               << (reference.empty() ? 0.0 : static_cast<double>(cumulativeSecondOrderViewSamples) / static_cast<double>(reference.size()))
               << "\n";
        report << "timing_mean_second_order_mu_phi_evaluations="
               << (reference.empty() ? 0.0 : static_cast<double>(cumulativeSecondOrderMuPhiEvaluations) / static_cast<double>(reference.size()))
               << "\n";
        report << "timing_max_direction_seconds=" << maxDirectionSeconds << "\n";
        report << "timing_slowest_direction_index=" << slowestDirectionIndex << "\n";
        report << "timing_slowest_direction_zenith_deg=" << slowestDirectionZenithDeg << "\n";
        report << "timing_slowest_direction_azimuth_deg=" << slowestDirectionAzimuthDeg << "\n";
        report << "timing_max_second_order_seconds=" << maxSecondOrderSeconds << "\n";
        report << "timing_slowest_second_order_direction_index=" << slowestSecondOrderDirectionIndex << "\n";
        report << "timing_slowest_second_order_zenith_deg=" << slowestSecondOrderZenithDeg << "\n";
        report << "timing_slowest_second_order_azimuth_deg=" << slowestSecondOrderAzimuthDeg << "\n";
        report << "timing_metrics_seconds="
               << std::chrono::duration<double>(std::chrono::steady_clock::now() - metricsStart).count()
               << "\n";
        if (hasIntensityReference && count > 0.0) {
            const double rmse = std::sqrt(sumSquared / count);
            const double brightestLocationError = angularSeparationDeg(
                directions[brightestReferenceIndex].zenith_deg,
                directions[brightestReferenceIndex].azimuth_deg,
                directions[brightestModelIndex].zenith_deg,
                directions[brightestModelIndex].azimuth_deg
            );
            report << "normalized_rmse=" << rmse << "\n";
            report << "brightest_location_deg=" << brightestLocationError << "\n";
        }
        if (!dopErrors.empty()) {
            report << "median_dolp_abs=" << percentile(dopErrors, 0.5) << "\n";
            report << "p95_dolp_abs=" << percentile(dopErrors, 0.95) << "\n";
            report << "max_dolp_abs=" << worstDopError << "\n";
            report << "max_dolp_zenith_deg=" << reference[worstDopIndex].zenith_deg << "\n";
            report << "max_dolp_relative_azimuth_deg="
                   << (reference[worstDopIndex].uses_relative_azimuth ? reference[worstDopIndex].relative_azimuth_deg : 0.0)
                   << "\n";
            report << "max_dolp_absolute_azimuth_deg=" << directions[worstDopIndex].azimuth_deg << "\n";
            report << "max_dolp_reference=" << reference[worstDopIndex].dop << "\n";
            report << "max_dolp_model=" << degreeOfLinearPolarization(model[worstDopIndex].mean) << "\n";
        } else {
            report << "median_dolp_abs=nan\n";
            report << "p95_dolp_abs=nan\n";
        }
        if (!aopErrorsDeg.empty()) {
            report << "median_aop_deg=" << percentile(aopErrorsDeg, 0.5) << "\n";
            report << "p95_aop_deg=" << percentile(aopErrorsDeg, 0.95) << "\n";
        } else {
            report << "median_aop_deg=nan\n";
            report << "p95_aop_deg=nan\n";
        }
        if (!solarVerticalSignedDolpBias.empty()) {
            const double meanBias = std::accumulate(
                solarVerticalSignedDolpBias.begin(),
                solarVerticalSignedDolpBias.end(),
                0.0
            ) / static_cast<double>(solarVerticalSignedDolpBias.size());
            report << "solar_vertical_signed_dolp_bias=" << meanBias << "\n";
        } else {
            report << "solar_vertical_signed_dolp_bias=nan\n";
        }

        report << "region_summary_csv=" << regionSummaryCsvPath.string() << "\n";
        for (const RegionSummary &summary : regionSummaries) {
            report << "region_" << summary.name << "_count=" << summary.count << "\n";
            report << "region_" << summary.name << "_dop_count=" << summary.dop_count << "\n";
            report << "region_" << summary.name << "_aop_count=" << summary.aop_count << "\n";
            report << "region_" << summary.name << "_median_dolp_abs=" << summary.median_dolp_abs << "\n";
            report << "region_" << summary.name << "_p95_dolp_abs=" << summary.p95_dolp_abs << "\n";
            report << "region_" << summary.name << "_median_aop_deg=" << summary.median_aop_deg << "\n";
            report << "region_" << summary.name << "_p95_aop_deg=" << summary.p95_aop_deg << "\n";
            report << "region_" << summary.name << "_mean_signed_dolp_bias=" << summary.mean_signed_dolp_bias << "\n";
            report << "region_" << summary.name << "_mean_first_frac=" << summary.mean_first_frac << "\n";
            report << "region_" << summary.name << "_mean_second_frac=" << summary.mean_second_frac << "\n";
            report << "region_" << summary.name << "_mean_higher_frac=" << summary.mean_higher_frac << "\n";
            report << "region_" << summary.name << "_mean_second_rr_frac=" << summary.mean_second_rr_frac << "\n";
            report << "region_" << summary.name << "_mean_second_ar_frac=" << summary.mean_second_ar_frac << "\n";
            report << "region_" << summary.name << "_mean_second_ra_frac=" << summary.mean_second_ra_frac << "\n";
            report << "region_" << summary.name << "_mean_second_aa_frac=" << summary.mean_second_aa_frac << "\n";
            report << "region_" << summary.name << "_mean_normalized_reference=" << summary.mean_normalized_reference << "\n";
            report << "region_" << summary.name << "_mean_normalized_model=" << summary.mean_normalized_model << "\n";
            report << "region_" << summary.name << "_max_dolp_abs=" << summary.max_dolp_abs << "\n";
            report << "region_" << summary.name << "_max_dolp_zenith_deg=" << summary.max_dolp_zenith_deg << "\n";
            report << "region_" << summary.name << "_max_dolp_relative_azimuth_deg=" << summary.max_dolp_relative_azimuth_deg << "\n";
        }
        report << "comparison_csv=" << comparisonCsvPath.string() << "\n";
        report << "timing_total_runtime_seconds="
               << std::chrono::duration<double>(std::chrono::steady_clock::now() - overallStart).count()
               << "\n";
        const std::filesystem::path reportPath = reportDir / (config.output.case_id + ".txt");
        std::ofstream reportFile(reportPath);
        reportFile << report.str();

        std::cout << report.str();
        std::cout << "report_txt=" << reportPath.string() << "\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Measurement case evaluation failed: " << error.what() << "\n";
        return 1;
    }
}
