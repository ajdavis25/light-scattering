#include "Atmosphere.hpp"
#include "WavelengthHandling.hpp"

#include <iostream>

int main()
{
    Atmosphere atmosphere;
    atmosphere.loadProfileCsv("../data/atmosphere/clear_sky_midlatitude.csv");

    WavelengthManager wavelengthManager;
    wavelengthManager.configureGrid(400.0, 500.0, 50.0);
    wavelengthManager.loadSolarSpectrumCsv("../data/optics/solar_irradiance_reference.csv");
    wavelengthManager.loadInstrumentResponseCsv("../data/paper_cases/interim_gal_public_case/instrument_response.csv");
    wavelengthManager.loadOzoneCrossSectionCsv("../data/optics/ozone_cross_section_reference.csv");
    wavelengthManager.loadO2CrossSectionCsv("../data/optics/o2_cross_section_reference.csv");
    wavelengthManager.loadO4CrossSectionCsv("../data/optics/o4_cross_section_reference.csv");
    wavelengthManager.loadH2OCrossSectionCsv("../data/optics/h2o_cross_section_reference.csv");
    wavelengthManager.loadNO2CrossSectionCsv("../data/optics/no2_cross_section_reference.csv");
    wavelengthManager.buildBands();

    const auto &bands = wavelengthManager.getBands();
    if (bands.empty()) {
        std::cerr << "Spectral grid was not built.\n";
        return 1;
    }

    double weightSum = 0.0;
    for (const SpectralBand &band : bands) {
        weightSum += band.normalized_weight;
    }
    if (std::abs(weightSum - 1.0) > 1.0e-9) {
        std::cerr << "Instrument-weighted spectral grid is not normalized.\n";
        return 1;
    }

    const LocalOpticalProperties opticalProperties = computeOpticalProperties(
        atmosphere,
        wavelengthManager,
        0.0,
        450.0
    );
    if (opticalProperties.ozone_absorption_m_inv < 0.0 ||
        opticalProperties.o2_absorption_m_inv < 0.0 ||
        opticalProperties.o4_absorption_m_inv < 0.0 ||
        opticalProperties.h2o_absorption_m_inv < 0.0 ||
        opticalProperties.no2_absorption_m_inv < 0.0) {
        std::cerr << "Gas absorption terms must remain non-negative.\n";
        return 1;
    }

    if (opticalProperties.o4_absorption_m_inv > 1.0e-3) {
        std::cerr << "O4 absorption is implausibly large at 450 nm; check cm^5 to m^5 conversion.\n";
        return 1;
    }

    if (bands.front().instrument_response <= 0.0) {
        std::cerr << "Instrument response should be loaded for the paper-case test grid.\n";
        return 1;
    }

    return 0;
}
