// WavelengthHandling.cpp
#include "WavelengthHandling.hpp"
#include "Atmosphere.hpp"   // so we can call atm.rayleighDensity(...), atm.aerosolDensity(...), etc.
#include <cmath>
#include <algorithm>

WavelengthManager::WavelengthManager()
{
    // Example bands (R/G/B or four bands)
    SpectralBand b1{400.0, 0.15}; // 400 nm
    SpectralBand b2{500.0, 0.40}; // 500 nm
    SpectralBand b3{600.0, 0.30}; // 600 nm
    SpectralBand b4{700.0, 0.15}; // 700 nm
    bands = {b1,b2,b3,b4};
    // You can add more or load from a file
}

/**
 * A refined scatteringCoefficient function that sums:
 *   - Rayleigh scattering ( ~ wave^-4 * rayleighDensity(alt) )
 *   - Aerosol scattering ( ~ wave^-p * aerosolDensity(alt) )
 */
double scatteringCoefficient(const Atmosphere& atm, double altitude_m, double wavelength_nm)
{
    // 1) Rayleigh: typical exponent ~4
    double waveFactorRay = std::pow((550.0 / wavelength_nm), 4.0);
    double densRay = atm.rayleighDensity(altitude_m);
    static const double RAY_CONST = 1e-5;  // or whatever reference you want
    double sigmaRay = densRay * waveFactorRay * RAY_CONST;

    // 2) Aerosols: exponent might be ~1.3 (or 1–2 depending on aerosol type)
    double waveFactorAer = std::pow((550.0 / wavelength_nm), 1.3);
    double densAero = atm.aerosolDensity(altitude_m);
    static const double AERO_CONST = 5e-5; // adjust to fit your aerosol intensity
    double sigmaAer = densAero * waveFactorAer * AERO_CONST;

    // Sum them
    double sigmaScat = sigmaRay + sigmaAer;
    return sigmaScat;
}
