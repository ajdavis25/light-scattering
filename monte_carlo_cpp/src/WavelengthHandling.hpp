// WavelengthHandling.hpp
#ifndef WAVELENGTHHANDLING_HPP
#define WAVELENGTHHANDLING_HPP

#include <vector>

// Forward-declare Atmosphere if you want to pass an object
class Atmosphere;

/**
 * A set of wavelength (nm) and flux fraction.
 * We'll keep the WavelengthManager class for storing multiple bands.
 */
struct SpectralBand
{
    double wavelength_nm;
    double fluxFraction;
};

class WavelengthManager
{
public:
    WavelengthManager();

    /**
     * Return the list of spectral bands
     */
    const std::vector<SpectralBand>& getBands() const { return bands; }

private:
    std::vector<SpectralBand> bands;
};

/**
 * Computes an approximate scattering coefficient 
 * at altitude (m) for a given wavelength (nm).
 * We'll pass a reference to an Atmosphere object 
 * if we want to query atmosphericDensity(...) etc.
 */
double scatteringCoefficient(const Atmosphere& atm, double altitude_m, double wavelength_nm);

#endif
