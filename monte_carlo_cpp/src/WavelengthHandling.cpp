// WavelengthHandling.cpp
#include "WavelengthHandling.hpp"
#include "Atmosphere.hpp"
#include <cmath>
#include <algorithm>

/**
 * Example: Rayleigh ~ wave^-4 * density
 * plus some user-chosen reference. 
 * We ignore aerosols for brevity, but you can add them too.
 */
double scatteringCoefficient(double alt_m, double wavelength_nm)
{
    // reference ~ at 550 nm => 1.0
    // wave term: (550 / wavelength)^4
    double waveFactor = std::pow((550.0 / wavelength_nm), 4.0);

    // use atmosphericDensity from Atmosphere module
    double dens = atmosphericDensity(alt_m);
    // some reference scaling
    static const double SCAT_REF = 1e-5; 

    return SCAT_REF * waveFactor * dens;
}
