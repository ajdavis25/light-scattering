// main.cpp
#include <iostream>
#include "PhaseFunctions.hpp"
#include "SurfaceReflection.hpp"
#include "Polarization.hpp"
#include "Atmosphere.hpp"
#include "WavelengthHandling.hpp"

int main()
{
    // 1) Basic test: PhaseFunctions
    double val = rayleighPhase(0.5);
    std::cout << "Rayleigh( cosTheta=0.5 ) => " << val << "\n";

    // 2) Surface reflection test
    auto refl = reflectLambertian(0.2);
    std::cout << "Lambert reflection dir z=" << refl.dir.z 
              << ", weight mult=" << refl.weightMultiplier << "\n";

    // 3) Polarization test
    auto unpol = initUnpolarized(1.0);
    auto outPol = applyIdentityMueller(unpol);
    std::cout << "Unpolarized stokes => I=" << outPol.I << "\n";

    // 4) Atmosphere test
    double alt = 5000.0; // 5 km
    double rho = atmosphericDensity(alt);
    double absorb = absorptionCoefficient(alt);
    std::cout << "At alt=" << alt 
              << "m => density=" << rho 
              << ", absorption=" << absorb << "\n";

    // 5) Wavelength test
    double wave_nm = 550.0; // green
    double scat = scatteringCoefficient(alt, wave_nm);
    std::cout << "Scattering at " << wave_nm << " nm => " << scat << "\n";

    return 0;
}
