/**************************************************************
 * main.cpp
 *
 * Example main that:
 *   - Does snippet tests
 *   - Calls runMonteCarloSimulation() for a more physical twilight scenario
 **************************************************************/
#include <iostream>
#include "PhaseFunctions.hpp"
#include "SurfaceReflection.hpp"
#include "Polarization.hpp"
#include "Atmosphere.hpp"
#include "WavelengthHandling.hpp"
#include "MonteCarloDriver.hpp"  

int main()
{
    {
        // 1) snippet tests
        double val = rayleighPhase(0.5);
        std::cout << "Rayleigh( cosTheta=0.5 ) => " << val << "\n";

        auto refl = reflectLambertian(0.2);
        std::cout << "Lambert reflection dir z=" << refl.dir.z
                  << ", weight mult=" << refl.weightMultiplier << "\n";

        auto unpol = initUnpolarized(1.0);
        auto outPol = applyIdentityMueller(unpol);
        std::cout << "Unpolarized stokes => I=" << outPol.I << "\n";

        Atmosphere atm;
        atm.loadLayerData("");
        double alt = 5000.0; 
        double rho = atm.atmosphericDensity(alt);
        double absorb = atm.absorptionCoeff(alt);
        std::cout << "At alt=" << alt
                  << "m => density=" << rho
                  << ", absorption=" << absorb << "\n";

        double wave_nm = 550.0;
        double scat = scatteringCoefficient(atm, alt, wave_nm);
        std::cout << "Scattering at " << wave_nm
                  << " nm => " << scat << "\n";
    }

    // 2) MAIN multi-scattering simulation
    std::cout << "\n=== Running Full Monte Carlo Simulation (Twilight) ===\n";
    runMonteCarloSimulation();
    std::cout << "=== Simulation complete. ===\n";

    return 0;
}
