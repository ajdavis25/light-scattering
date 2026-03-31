#include "Atmosphere.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double DEFAULT_TOA_ALTITUDE_M = 1.0e5;
constexpr double DEFAULT_AEROSOL_SCATTERING_ANGSTROM_EXPONENT = 1.6;
constexpr double DEFAULT_AEROSOL_ABSORPTION_ANGSTROM_EXPONENT = 1.0;

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

double linearInterpolate(double x0, double x1, double y0, double y1, double x)
{
    if (std::abs(x1 - x0) < 1.0e-12) {
        return y0;
    }
    const double t = (x - x0) / (x1 - x0);
    return y0 + t * (y1 - y0);
}

double positiveInterpolate(double x0, double x1, double y0, double y1, double x)
{
    if (y0 > 0.0 && y1 > 0.0) {
        const double logY = linearInterpolate(x0, x1, std::log(y0), std::log(y1), x);
        return std::exp(logY);
    }
    return linearInterpolate(x0, x1, y0, y1, x);
}

double getNamedOrPositionalValue(
    const std::vector<std::string> &values,
    const std::map<std::string, std::size_t> &headerIndex,
    const std::string &name,
    std::size_t position,
    double defaultValue,
    bool required = false
)
{
    const auto iterator = headerIndex.find(name);
    if (iterator != headerIndex.end()) {
        if (iterator->second >= values.size() || values[iterator->second].empty()) {
            if (required) {
                throw std::runtime_error("Missing required atmosphere column value: " + name);
            }
            return defaultValue;
        }
        return std::stod(values[iterator->second]);
    }

    if (position < values.size() && !values[position].empty()) {
        return std::stod(values[position]);
    }

    if (required) {
        throw std::runtime_error("Atmosphere profile is missing required column: " + name);
    }
    return defaultValue;
}

AtmosphereState makeState(const AtmosphereLayer &layer, double altitude_m)
{
    return {
        altitude_m,
        layer.pressure_pa,
        layer.temperature_k,
        layer.molecular_number_density_m3,
        layer.ozone_number_density_m3,
        layer.h2o_number_density_m3,
        layer.no2_number_density_m3,
        layer.aerosol_extinction_550_m_inv,
        layer.aerosol_single_scattering_albedo,
        layer.aerosol_asymmetry,
        layer.aerosol_scattering_angstrom_exponent,
        layer.aerosol_absorption_angstrom_exponent,
    };
}
}

void Atmosphere::loadProfileCsv(const std::string &filename)
{
    layers_.clear();

    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("Unable to open atmosphere profile: " + filename);
    }

    std::string line;
    std::vector<std::string> header;
    std::map<std::string, std::size_t> headerIndex;
    while (std::getline(input, line)) {
        const std::string stripped = trim(line);
        if (stripped.empty() || stripped[0] == '#') {
            continue;
        }

        if (header.empty()) {
            header = splitCsvLine(stripped);
            for (std::size_t index = 0; index < header.size(); ++index) {
                headerIndex[header[index]] = index;
            }
            continue;
        }

        const std::vector<std::string> values = splitCsvLine(stripped);
        if (values.size() < 8) {
            throw std::runtime_error("Atmosphere profile row must have 8 columns: " + stripped);
        }

        AtmosphereLayer layer {};
        layer.altitude_m = getNamedOrPositionalValue(values, headerIndex, "altitude_m", 0, 0.0, true);
        layer.pressure_pa = getNamedOrPositionalValue(values, headerIndex, "pressure_pa", 1, 0.0, true);
        layer.temperature_k = getNamedOrPositionalValue(values, headerIndex, "temperature_k", 2, 0.0, true);
        layer.molecular_number_density_m3 = getNamedOrPositionalValue(
            values,
            headerIndex,
            "molecular_number_density_m3",
            3,
            0.0,
            true
        );
        layer.ozone_number_density_m3 = getNamedOrPositionalValue(
            values,
            headerIndex,
            "ozone_number_density_m3",
            4,
            0.0,
            true
        );
        layer.h2o_number_density_m3 = getNamedOrPositionalValue(
            values,
            headerIndex,
            "h2o_number_density_m3",
            10,
            0.0
        );
        layer.no2_number_density_m3 = getNamedOrPositionalValue(
            values,
            headerIndex,
            "no2_number_density_m3",
            11,
            0.0
        );
        layer.aerosol_extinction_550_m_inv = getNamedOrPositionalValue(
            values,
            headerIndex,
            "aerosol_extinction_550_m_inv",
            5,
            0.0,
            true
        );
        layer.aerosol_single_scattering_albedo = getNamedOrPositionalValue(
            values,
            headerIndex,
            "aerosol_single_scattering_albedo",
            6,
            1.0,
            true
        );
        layer.aerosol_asymmetry = getNamedOrPositionalValue(
            values,
            headerIndex,
            "aerosol_asymmetry",
            7,
            0.6,
            true
        );
        layer.aerosol_scattering_angstrom_exponent = getNamedOrPositionalValue(
            values,
            headerIndex,
            "aerosol_scattering_angstrom_exponent",
            8,
            DEFAULT_AEROSOL_SCATTERING_ANGSTROM_EXPONENT
        );
        layer.aerosol_absorption_angstrom_exponent = getNamedOrPositionalValue(
            values,
            headerIndex,
            "aerosol_absorption_angstrom_exponent",
            9,
            DEFAULT_AEROSOL_ABSORPTION_ANGSTROM_EXPONENT
        );
        layers_.push_back(layer);
    }

    if (layers_.empty()) {
        throw std::runtime_error("Atmosphere profile contained no layers: " + filename);
    }

    std::sort(layers_.begin(), layers_.end(), [](const AtmosphereLayer &lhs, const AtmosphereLayer &rhs) {
        return lhs.altitude_m < rhs.altitude_m;
    });
}

void Atmosphere::loadLayerData(const std::string &filename)
{
    loadProfileCsv(filename);
}

double Atmosphere::topOfAtmosphereAltitudeM() const
{
    return layers_.empty() ? DEFAULT_TOA_ALTITUDE_M : layers_.back().altitude_m;
}

AtmosphereState Atmosphere::stateAtAltitude(double altitude_m) const
{
    if (layers_.empty()) {
        throw std::runtime_error("Atmosphere state requested before profile was loaded.");
    }

    if (altitude_m <= layers_.front().altitude_m) {
        return makeState(layers_.front(), altitude_m);
    }

    if (altitude_m >= layers_.back().altitude_m) {
        return makeState(layers_.back(), altitude_m);
    }

    for (std::size_t index = 1; index < layers_.size(); ++index) {
        const AtmosphereLayer &lower = layers_[index - 1];
        const AtmosphereLayer &upper = layers_[index];
        if (altitude_m <= upper.altitude_m) {
            return {
                altitude_m,
                positiveInterpolate(lower.altitude_m, upper.altitude_m, lower.pressure_pa, upper.pressure_pa, altitude_m),
                linearInterpolate(lower.altitude_m, upper.altitude_m, lower.temperature_k, upper.temperature_k, altitude_m),
                positiveInterpolate(
                    lower.altitude_m,
                    upper.altitude_m,
                    lower.molecular_number_density_m3,
                    upper.molecular_number_density_m3,
                    altitude_m
                ),
                positiveInterpolate(
                    lower.altitude_m,
                    upper.altitude_m,
                    lower.ozone_number_density_m3,
                    upper.ozone_number_density_m3,
                    altitude_m
                ),
                positiveInterpolate(
                    lower.altitude_m,
                    upper.altitude_m,
                    lower.h2o_number_density_m3,
                    upper.h2o_number_density_m3,
                    altitude_m
                ),
                positiveInterpolate(
                    lower.altitude_m,
                    upper.altitude_m,
                    lower.no2_number_density_m3,
                    upper.no2_number_density_m3,
                    altitude_m
                ),
                positiveInterpolate(
                    lower.altitude_m,
                    upper.altitude_m,
                    lower.aerosol_extinction_550_m_inv,
                    upper.aerosol_extinction_550_m_inv,
                    altitude_m
                ),
                linearInterpolate(
                    lower.altitude_m,
                    upper.altitude_m,
                    lower.aerosol_single_scattering_albedo,
                    upper.aerosol_single_scattering_albedo,
                    altitude_m
                ),
                linearInterpolate(
                    lower.altitude_m,
                    upper.altitude_m,
                    lower.aerosol_asymmetry,
                    upper.aerosol_asymmetry,
                    altitude_m
                ),
                linearInterpolate(
                    lower.altitude_m,
                    upper.altitude_m,
                    lower.aerosol_scattering_angstrom_exponent,
                    upper.aerosol_scattering_angstrom_exponent,
                    altitude_m
                ),
                linearInterpolate(
                    lower.altitude_m,
                    upper.altitude_m,
                    lower.aerosol_absorption_angstrom_exponent,
                    upper.aerosol_absorption_angstrom_exponent,
                    altitude_m
                ),
            };
        }
    }

    return makeState(layers_.back(), altitude_m);
}

double Atmosphere::pressurePa(double altitude_m) const
{
    return stateAtAltitude(altitude_m).pressure_pa;
}

double Atmosphere::temperatureK(double altitude_m) const
{
    return stateAtAltitude(altitude_m).temperature_k;
}

double Atmosphere::molecularNumberDensity(double altitude_m) const
{
    return stateAtAltitude(altitude_m).molecular_number_density_m3;
}

double Atmosphere::ozoneNumberDensity(double altitude_m) const
{
    return stateAtAltitude(altitude_m).ozone_number_density_m3;
}

double Atmosphere::h2oNumberDensity(double altitude_m) const
{
    return stateAtAltitude(altitude_m).h2o_number_density_m3;
}

double Atmosphere::no2NumberDensity(double altitude_m) const
{
    return stateAtAltitude(altitude_m).no2_number_density_m3;
}

double Atmosphere::aerosolExtinction550(double altitude_m) const
{
    return stateAtAltitude(altitude_m).aerosol_extinction_550_m_inv;
}

double Atmosphere::aerosolSingleScatteringAlbedo(double altitude_m) const
{
    return stateAtAltitude(altitude_m).aerosol_single_scattering_albedo;
}

double Atmosphere::aerosolAsymmetry(double altitude_m) const
{
    return stateAtAltitude(altitude_m).aerosol_asymmetry;
}

double Atmosphere::aerosolScatteringAngstromExponent(double altitude_m) const
{
    return stateAtAltitude(altitude_m).aerosol_scattering_angstrom_exponent;
}

double Atmosphere::aerosolAbsorptionAngstromExponent(double altitude_m) const
{
    return stateAtAltitude(altitude_m).aerosol_absorption_angstrom_exponent;
}

double Atmosphere::rayleighDensity(double altitude_m) const
{
    const double groundValue = layers_.empty() ? 1.0 : layers_.front().molecular_number_density_m3;
    if (groundValue <= 0.0) {
        return 0.0;
    }
    return molecularNumberDensity(altitude_m) / groundValue;
}

double Atmosphere::aerosolDensity(double altitude_m) const
{
    const double groundValue = layers_.empty() ? 1.0 : layers_.front().aerosol_extinction_550_m_inv;
    if (groundValue <= 0.0) {
        return 0.0;
    }
    return aerosolExtinction550(altitude_m) / groundValue;
}

double Atmosphere::absorptionCoeff(double altitude_m) const
{
    const double groundValue = layers_.empty() ? 1.0 : layers_.front().ozone_number_density_m3;
    if (groundValue <= 0.0) {
        return 0.0;
    }
    return ozoneNumberDensity(altitude_m) / groundValue;
}

double Atmosphere::atmosphericDensity(double altitude_m) const
{
    return rayleighDensity(altitude_m) + aerosolDensity(altitude_m);
}
