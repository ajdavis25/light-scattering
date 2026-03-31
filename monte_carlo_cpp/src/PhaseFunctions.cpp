#include "PhaseFunctions.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double PI = 3.14159265358979323846;

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

double clamp(double value, double lower, double upper)
{
    if (value < lower) {
        return lower;
    }
    if (value > upper) {
        return upper;
    }
    return value;
}

PhaseMatrixCoefficients interpolateCoefficients(
    const PhaseMatrixCoefficients &lhs,
    const PhaseMatrixCoefficients &rhs,
    double t
)
{
    return {
        lhs.f11 + t * (rhs.f11 - lhs.f11),
        lhs.f12 + t * (rhs.f12 - lhs.f12),
        lhs.f22 + t * (rhs.f22 - lhs.f22),
        lhs.f33 + t * (rhs.f33 - lhs.f33),
        lhs.f34 + t * (rhs.f34 - lhs.f34),
        lhs.f44 + t * (rhs.f44 - lhs.f44),
    };
}

std::vector<AerosolPhaseMatrixTable::AngleEntry> normalizeTable(
    const std::vector<AerosolPhaseMatrixTable::AngleEntry> &table
)
{
    if (table.empty()) {
        return table;
    }

    double integral = 0.0;
    for (std::size_t index = 1; index < table.size(); ++index) {
        const double theta0 = table[index - 1].angle_deg * PI / 180.0;
        const double theta1 = table[index].angle_deg * PI / 180.0;
        const double dTheta = theta1 - theta0;
        const double p0 = table[index - 1].coeffs.f11 * std::sin(theta0);
        const double p1 = table[index].coeffs.f11 * std::sin(theta1);
        integral += 2.0 * PI * 0.5 * (p0 + p1) * dTheta;
    }

    const double normalization = integral > 0.0 ? integral : 1.0;
    std::vector<AerosolPhaseMatrixTable::AngleEntry> normalized = table;
    double cumulative = 0.0;
    normalized.front().cdf = 0.0;
    normalized.front().coeffs = interpolateCoefficients(normalized.front().coeffs, normalized.front().coeffs, 0.0);
    normalized.front().coeffs.f11 /= normalization;
    normalized.front().coeffs.f12 /= normalization;
    normalized.front().coeffs.f22 /= normalization;
    normalized.front().coeffs.f33 /= normalization;
    normalized.front().coeffs.f34 /= normalization;
    normalized.front().coeffs.f44 /= normalization;

    for (std::size_t index = 1; index < normalized.size(); ++index) {
        auto &entry = normalized[index];
        entry.coeffs.f11 /= normalization;
        entry.coeffs.f12 /= normalization;
        entry.coeffs.f22 /= normalization;
        entry.coeffs.f33 /= normalization;
        entry.coeffs.f34 /= normalization;
        entry.coeffs.f44 /= normalization;

        const double theta0 = normalized[index - 1].angle_deg * PI / 180.0;
        const double theta1 = entry.angle_deg * PI / 180.0;
        const double dTheta = theta1 - theta0;
        const double p0 = normalized[index - 1].coeffs.f11 * std::sin(theta0);
        const double p1 = entry.coeffs.f11 * std::sin(theta1);
        cumulative += 2.0 * PI * 0.5 * (p0 + p1) * dTheta;
        entry.cdf = cumulative;
    }

    if (normalized.back().cdf > 0.0) {
        for (auto &entry : normalized) {
            entry.cdf /= normalized.back().cdf;
        }
        normalized.back().cdf = 1.0;
    }

    return normalized;
}

PhaseMatrixCoefficients interpolateInAngle(
    const std::vector<AerosolPhaseMatrixTable::AngleEntry> &table,
    double angle_deg
)
{
    if (table.empty()) {
        return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    }
    if (angle_deg <= table.front().angle_deg) {
        return table.front().coeffs;
    }
    if (angle_deg >= table.back().angle_deg) {
        return table.back().coeffs;
    }

    for (std::size_t index = 1; index < table.size(); ++index) {
        if (angle_deg <= table[index].angle_deg) {
            const double angle0 = table[index - 1].angle_deg;
            const double angle1 = table[index].angle_deg;
            const double t = (angle_deg - angle0) / std::max(1.0e-12, angle1 - angle0);
            return interpolateCoefficients(table[index - 1].coeffs, table[index].coeffs, t);
        }
    }

    return table.back().coeffs;
}
}

void AerosolPhaseMatrixTable::loadCsv(const std::string &filename)
{
    wavelengths_nm_.clear();
    tables_.clear();

    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("Unable to open aerosol phase table: " + filename);
    }

    std::string line;
    bool headerSkipped = false;
    std::vector<double> rawWavelengths;
    std::vector<std::vector<AngleEntry>> rawTables;
    while (std::getline(input, line)) {
        const std::string stripped = trim(line);
        if (stripped.empty() || stripped[0] == '#') {
            continue;
        }

        if (!headerSkipped) {
            headerSkipped = true;
            continue;
        }

        const std::vector<std::string> values = splitCsvLine(stripped);
        if (values.size() < 8) {
            throw std::runtime_error("Aerosol phase row must have 8 columns: " + stripped);
        }

        const double wavelength = std::stod(values[0]);
        const double angle = std::stod(values[1]);
        PhaseMatrixCoefficients coeffs {
            std::stod(values[2]),
            std::stod(values[3]),
            std::stod(values[4]),
            std::stod(values[5]),
            std::stod(values[6]),
            std::stod(values[7]),
        };

        auto wavelengthIt = std::find(rawWavelengths.begin(), rawWavelengths.end(), wavelength);
        if (wavelengthIt == rawWavelengths.end()) {
            rawWavelengths.push_back(wavelength);
            rawTables.push_back({});
            wavelengthIt = std::prev(rawWavelengths.end());
        }

        const std::size_t index = static_cast<std::size_t>(std::distance(rawWavelengths.begin(), wavelengthIt));
        rawTables[index].push_back({angle, coeffs, 0.0});
    }

    for (std::size_t index = 0; index < rawWavelengths.size(); ++index) {
        auto &table = rawTables[index];
        std::sort(table.begin(), table.end(), [](const AngleEntry &lhs, const AngleEntry &rhs) {
            return lhs.angle_deg < rhs.angle_deg;
        });
        wavelengths_nm_.push_back(rawWavelengths[index]);
        tables_.push_back(normalizeTable(table));
    }

    if (wavelengths_nm_.empty()) {
        throw std::runtime_error("Aerosol phase table contained no usable rows: " + filename);
    }
}

std::vector<AerosolPhaseMatrixTable::AngleEntry> AerosolPhaseMatrixTable::interpolatedTable(double wavelength_nm) const
{
    if (wavelengths_nm_.empty()) {
        throw std::runtime_error("Aerosol phase table requested before loading a lookup table.");
    }

    if (wavelength_nm <= wavelengths_nm_.front()) {
        return tables_.front();
    }
    if (wavelength_nm >= wavelengths_nm_.back()) {
        return tables_.back();
    }

    for (std::size_t index = 1; index < wavelengths_nm_.size(); ++index) {
        if (wavelength_nm <= wavelengths_nm_[index]) {
            const double wavelength0 = wavelengths_nm_[index - 1];
            const double wavelength1 = wavelengths_nm_[index];
            const double t = (wavelength_nm - wavelength0) / std::max(1.0e-12, wavelength1 - wavelength0);
            const auto &lower = tables_[index - 1];
            const auto &upper = tables_[index];

            std::vector<AngleEntry> interpolated;
            interpolated.reserve(lower.size());
            for (const AngleEntry &entry : lower) {
                interpolated.push_back({
                    entry.angle_deg,
                    interpolateCoefficients(entry.coeffs, interpolateInAngle(upper, entry.angle_deg), t),
                    0.0,
                });
            }
            return normalizeTable(interpolated);
        }
    }

    return tables_.back();
}

PhaseMatrixCoefficients AerosolPhaseMatrixTable::coefficients(double wavelength_nm, double cos_theta) const
{
    if (wavelengths_nm_.empty()) {
        throw std::runtime_error("Aerosol phase table requested before loading a lookup table.");
    }

    const double angle_deg = std::acos(clamp(cos_theta, -1.0, 1.0)) * 180.0 / PI;
    if (wavelength_nm <= wavelengths_nm_.front()) {
        return interpolateInAngle(tables_.front(), angle_deg);
    }
    if (wavelength_nm >= wavelengths_nm_.back()) {
        return interpolateInAngle(tables_.back(), angle_deg);
    }

    for (std::size_t index = 1; index < wavelengths_nm_.size(); ++index) {
        if (wavelength_nm <= wavelengths_nm_[index]) {
            const double wavelength0 = wavelengths_nm_[index - 1];
            const double wavelength1 = wavelengths_nm_[index];
            const double t = (wavelength_nm - wavelength0) / std::max(1.0e-12, wavelength1 - wavelength0);
            const PhaseMatrixCoefficients lower = interpolateInAngle(tables_[index - 1], angle_deg);
            const PhaseMatrixCoefficients upper = interpolateInAngle(tables_[index], angle_deg);
            return interpolateCoefficients(lower, upper, t);
        }
    }

    return interpolateInAngle(tables_.back(), angle_deg);
}

double AerosolPhaseMatrixTable::phasePdf(double wavelength_nm, double cos_theta) const
{
    return std::max(0.0, coefficients(wavelength_nm, cos_theta).f11);
}

ScatteringSample AerosolPhaseMatrixTable::sampleDirection(double wavelength_nm, std::mt19937 &rng) const
{
    const std::vector<AngleEntry> table = interpolatedTable(wavelength_nm);
    std::uniform_real_distribution<double> uniform01(0.0, 1.0);
    const double xi = uniform01(rng);

    for (std::size_t index = 1; index < table.size(); ++index) {
        if (xi <= table[index].cdf) {
            const double cdf0 = table[index - 1].cdf;
            const double cdf1 = table[index].cdf;
            const double t = (xi - cdf0) / std::max(1.0e-12, cdf1 - cdf0);
            const double angle_deg = table[index - 1].angle_deg + t * (table[index].angle_deg - table[index - 1].angle_deg);
            const double cosTheta = std::cos(angle_deg * PI / 180.0);
            return {cosTheta, 2.0 * PI * uniform01(rng), phasePdf(wavelength_nm, cosTheta)};
        }
    }

    const double angle_deg = table.back().angle_deg;
    const double cosTheta = std::cos(angle_deg * PI / 180.0);
    return {cosTheta, 2.0 * PI * uniform01(rng), phasePdf(wavelength_nm, cosTheta)};
}

double rayleighPhase(double cosTheta)
{
    return (3.0 / (16.0 * PI)) * (1.0 + cosTheta * cosTheta);
}

double rayleighPolarizationFraction(double cosTheta)
{
    const double cos2 = cosTheta * cosTheta;
    const double denominator = 1.0 + cos2;
    if (denominator <= 0.0) {
        return 0.0;
    }
    return std::max(0.0, (1.0 - cos2) / denominator);
}

PhaseMatrixCoefficients rayleighPhaseMatrix(double cosTheta)
{
    const double cos2 = cosTheta * cosTheta;
    const double factor = 3.0 / (16.0 * PI);
    return {
        factor * (1.0 + cos2),
        -factor * (1.0 - cos2),
        factor * (1.0 + cos2),
        factor * (2.0 * cosTheta),
        0.0,
        factor * (2.0 * cosTheta),
    };
}

ScatteringSample sampleRayleighDirection(std::mt19937 &rng)
{
    std::uniform_real_distribution<double> uniform01(0.0, 1.0);
    while (true) {
        const double cosTheta = 2.0 * uniform01(rng) - 1.0;
        const double accept = 0.75 * uniform01(rng);
        const double pdfMu = 0.375 * (1.0 + cosTheta * cosTheta);
        if (accept <= pdfMu) {
            return {cosTheta, 2.0 * PI * uniform01(rng), rayleighPhase(cosTheta)};
        }
    }
}

double henyeyGreenstein(double cosTheta, double g)
{
    const double g2 = g * g;
    const double denom = 1.0 + g2 - 2.0 * g * cosTheta;
    return ((1.0 - g2) / (4.0 * PI)) / std::pow(denom, 1.5);
}
