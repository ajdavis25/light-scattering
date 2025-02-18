/*******************************************************
 * MonteCarloDriver.hpp
 *
 * Updated header with:
 *   - Photon struct
 *   - function prototypes for launchPhoton, tracePhoton, runMonteCarloSimulation
 *   - <vector> for std::vector usage
 *******************************************************/
#ifndef MONTECARLODRIVER_HPP
#define MONTECARLODRIVER_HPP

#include <vector>
#include "Vec3.hpp"
#include "Polarization.hpp"
#include "SurfaceReflection.hpp"

// Our Photon struct
struct Photon
{
    Vec3 pos;
    Vec3 dir;
    double weight;
    double wavelength_nm;
    StokesVector stokes;
};

// function prototypes
Photon launchPhoton(double wave_nm, double fluxFrac);
void tracePhoton(Photon &p, int bandIndex,
                 std::vector<std::vector<double>> &localRadiance);
void runMonteCarloSimulation();

#endif // MONTECARLODRIVER_HPP
