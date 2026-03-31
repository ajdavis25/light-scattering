#include "Atmosphere.hpp"

#include <iostream>

int main()
{
    Atmosphere atmosphere;
    atmosphere.loadProfileCsv("../data/atmosphere/clear_sky_midlatitude.csv");

    if (atmosphere.empty()) {
        std::cerr << "Atmosphere profile did not load any layers.\n";
        return 1;
    }

    const AtmosphereState ground = atmosphere.stateAtAltitude(0.0);
    const AtmosphereState upper = atmosphere.stateAtAltitude(20000.0);
    const AtmosphereState mid = atmosphere.stateAtAltitude(3500.0);
    if (upper.molecular_number_density_m3 >= ground.molecular_number_density_m3) {
        std::cerr << "Molecular number density should decrease with altitude.\n";
        return 1;
    }

    if (upper.aerosol_extinction_550_m_inv >= ground.aerosol_extinction_550_m_inv) {
        std::cerr << "Aerosol extinction should decrease with altitude.\n";
        return 1;
    }

    if (ground.h2o_number_density_m3 < 0.0 || upper.no2_number_density_m3 < 0.0) {
        std::cerr << "Optional gas number densities should remain non-negative.\n";
        return 1;
    }

    if (ground.aerosol_scattering_angstrom_exponent <= 0.0 || mid.aerosol_absorption_angstrom_exponent < 0.0) {
        std::cerr << "Aerosol spectral exponents should load from the profile.\n";
        return 1;
    }

    if (mid.aerosol_extinction_550_m_inv <= 0.0) {
        std::cerr << "Interpolated aerosol extinction should remain positive.\n";
        return 1;
    }

    return 0;
}
