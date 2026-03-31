#include "WavelengthHandling.hpp"

#include "Atmosphere.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double DRY_AIR_O2_MIXING_RATIO = 0.20946;

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

std::vector<std::pair<double, double>> loadTwoColumnTable(const std::string &filename)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("Unable to open spectral table: " + filename);
    }

    std::vector<std::pair<double, double>> table;
    std::string line;
    bool headerSkipped = false;
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
        if (values.size() < 2) {
            throw std::runtime_error("Spectral table row must have 2 columns: " + stripped);
        }
        table.emplace_back(std::stod(values[0]), std::stod(values[1]));
    }

    if (table.empty()) {
        throw std::runtime_error("Spectral table was empty: " + filename);
    }

    std::sort(table.begin(), table.end(), [](const auto &lhs, const auto &rhs) { return lhs.first < rhs.first; });
    return table;
}

double interpolateTableValue(const std::vector<std::pair<double, double>> &table, double x)
{
    if (table.empty()) {
        return 0.0;
    }
    if (x <= table.front().first) {
        return table.front().second;
    }
    if (x >= table.back().first) {
        return table.back().second;
    }

    for (std::size_t index = 1; index < table.size(); ++index) {
        if (x <= table[index].first) {
            const double x0 = table[index - 1].first;
            const double x1 = table[index].first;
            const double y0 = table[index - 1].second;
            const double y1 = table[index].second;
            const double t = (x - x0) / std::max(1.0e-12, x1 - x0);
            return y0 + t * (y1 - y0);
        }
    }

    return table.back().second;
}
}

void WavelengthManager::configureGrid(double min_wavelength_nm, double max_wavelength_nm, double step_nm)
{
    min_wavelength_nm_ = min_wavelength_nm;
    max_wavelength_nm_ = max_wavelength_nm;
    step_nm_ = step_nm;
}

void WavelengthManager::loadSolarSpectrumCsv(const std::string &filename)
{
    solar_irradiance_table_ = loadTwoColumnTable(filename);
}

void WavelengthManager::loadInstrumentResponseCsv(const std::string &filename)
{
    instrument_response_table_ = loadTwoColumnTable(filename);
}

void WavelengthManager::loadRayleighCrossSectionCsv(const std::string &filename)
{
    rayleigh_cross_section_table_ = loadTwoColumnTable(filename);
}

void WavelengthManager::loadOzoneCrossSectionCsv(const std::string &filename)
{
    ozone_cross_section_table_ = loadTwoColumnTable(filename);
}

void WavelengthManager::loadO2CrossSectionCsv(const std::string &filename)
{
    o2_cross_section_table_ = loadTwoColumnTable(filename);
}

void WavelengthManager::loadO4CrossSectionCsv(const std::string &filename)
{
    o4_cross_section_table_ = loadTwoColumnTable(filename);
}

void WavelengthManager::loadH2OCrossSectionCsv(const std::string &filename)
{
    h2o_cross_section_table_ = loadTwoColumnTable(filename);
}

void WavelengthManager::loadNO2CrossSectionCsv(const std::string &filename)
{
    no2_cross_section_table_ = loadTwoColumnTable(filename);
}

void WavelengthManager::buildBands()
{
    bands_.clear();
    double totalWeight = 0.0;
    for (double wavelength = min_wavelength_nm_; wavelength <= max_wavelength_nm_ + 1.0e-9; wavelength += step_nm_) {
        SpectralBand band {};
        band.wavelength_nm = wavelength;
        band.instrument_response = instrumentResponse(wavelength);
        band.solar_irradiance_w_m2_nm = solarIrradiance(wavelength) * band.instrument_response;
        band.rayleigh_cross_section_m2 = rayleighCrossSection(wavelength);
        band.ozone_cross_section_m2 = ozoneCrossSection(wavelength);
        band.o2_cross_section_m2 = o2CrossSection(wavelength);
        band.o4_cross_section_m5 = o4CrossSection(wavelength);
        band.h2o_cross_section_m2 = h2oCrossSection(wavelength);
        band.no2_cross_section_m2 = no2CrossSection(wavelength);
        band.normalized_weight = band.solar_irradiance_w_m2_nm;
        totalWeight += band.solar_irradiance_w_m2_nm;
        bands_.push_back(band);
    }

    if (totalWeight <= 0.0) {
        throw std::runtime_error("Spectral weighting produced a zero-weight spectral grid.");
    }

    for (auto &band : bands_) {
        band.normalized_weight /= totalWeight;
    }
}

double WavelengthManager::solarIrradiance(double wavelength_nm) const
{
    return interpolateTableValue(solar_irradiance_table_, wavelength_nm);
}

double WavelengthManager::instrumentResponse(double wavelength_nm) const
{
    if (instrument_response_table_.empty()) {
        return 1.0;
    }
    return interpolateTableValue(instrument_response_table_, wavelength_nm);
}

double WavelengthManager::rayleighCrossSection(double wavelength_nm) const
{
    if (!rayleigh_cross_section_table_.empty()) {
        return interpolateTableValue(rayleigh_cross_section_table_, wavelength_nm);
    }
    return rayleighCrossSectionM2(wavelength_nm);
}

double WavelengthManager::ozoneCrossSection(double wavelength_nm) const
{
    return interpolateTableValue(ozone_cross_section_table_, wavelength_nm);
}

double WavelengthManager::o2CrossSection(double wavelength_nm) const
{
    if (o2_cross_section_table_.empty()) {
        return 0.0;
    }
    return interpolateTableValue(o2_cross_section_table_, wavelength_nm);
}

double WavelengthManager::o4CrossSection(double wavelength_nm) const
{
    if (o4_cross_section_table_.empty()) {
        return 0.0;
    }
    return interpolateTableValue(o4_cross_section_table_, wavelength_nm);
}

double WavelengthManager::h2oCrossSection(double wavelength_nm) const
{
    if (h2o_cross_section_table_.empty()) {
        return 0.0;
    }
    return interpolateTableValue(h2o_cross_section_table_, wavelength_nm);
}

double WavelengthManager::no2CrossSection(double wavelength_nm) const
{
    if (no2_cross_section_table_.empty()) {
        return 0.0;
    }
    return interpolateTableValue(no2_cross_section_table_, wavelength_nm);
}

double rayleighCrossSectionM2(double wavelength_nm)
{
    const double wavelengthRatio = 550.0 / wavelength_nm;
    return 5.45e-31 * std::pow(wavelengthRatio, 4.09);
}

LocalOpticalProperties computeOpticalProperties(
    const Atmosphere &atmosphere,
    const WavelengthManager &wavelength_manager,
    double altitude_m,
    double wavelength_nm
)
{
    const AtmosphereState state = atmosphere.stateAtAltitude(altitude_m);

    const double rayleighScattering = state.molecular_number_density_m3 *
        wavelength_manager.rayleighCrossSection(wavelength_nm);
    const double o2NumberDensity = state.molecular_number_density_m3 * DRY_AIR_O2_MIXING_RATIO;
    const double aerosolScattering550 = state.aerosol_extinction_550_m_inv * state.aerosol_single_scattering_albedo;
    const double aerosolAbsorption550 = state.aerosol_extinction_550_m_inv * (1.0 - state.aerosol_single_scattering_albedo);
    const double aerosolScattering = aerosolScattering550 *
        std::pow(550.0 / wavelength_nm, state.aerosol_scattering_angstrom_exponent);
    const double aerosolAbsorption = aerosolAbsorption550 *
        std::pow(550.0 / wavelength_nm, state.aerosol_absorption_angstrom_exponent);
    const double ozoneAbsorption = state.ozone_number_density_m3 * wavelength_manager.ozoneCrossSection(wavelength_nm);
    const double o2Absorption = o2NumberDensity * wavelength_manager.o2CrossSection(wavelength_nm);
    const double o4Absorption = o2NumberDensity * o2NumberDensity * wavelength_manager.o4CrossSection(wavelength_nm);
    const double h2oAbsorption = state.h2o_number_density_m3 * wavelength_manager.h2oCrossSection(wavelength_nm);
    const double no2Absorption = state.no2_number_density_m3 * wavelength_manager.no2CrossSection(wavelength_nm);
    const double extinction =
        rayleighScattering +
        aerosolScattering +
        aerosolAbsorption +
        ozoneAbsorption +
        o2Absorption +
        o4Absorption +
        h2oAbsorption +
        no2Absorption;

    LocalOpticalProperties properties {};
    properties.rayleigh_scattering_m_inv = rayleighScattering;
    properties.aerosol_scattering_m_inv = aerosolScattering;
    properties.aerosol_absorption_m_inv = aerosolAbsorption;
    properties.ozone_absorption_m_inv = ozoneAbsorption;
    properties.o2_absorption_m_inv = o2Absorption;
    properties.o4_absorption_m_inv = o4Absorption;
    properties.h2o_absorption_m_inv = h2oAbsorption;
    properties.no2_absorption_m_inv = no2Absorption;
    properties.extinction_m_inv = extinction;
    properties.single_scattering_albedo =
        extinction > 0.0 ? (rayleighScattering + aerosolScattering) / extinction : 0.0;
    properties.aerosol_asymmetry = state.aerosol_asymmetry;
    return properties;
}
