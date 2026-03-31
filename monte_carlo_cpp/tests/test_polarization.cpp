#include "Polarization.hpp"

#include <cmath>
#include <iostream>

int main()
{
    const StokesVector input = initUnpolarized(1.0);
    const StokesVector scattered = applyRayleighMueller(input, 0.0);
    if (!isPhysicallyValid(scattered)) {
        std::cerr << "Rayleigh Mueller output is not physically valid.\n";
        return 1;
    }

    const StokesVector rotated = rotateStokes(scattered, 0.5);
    const StokesVector unrotated = rotateStokes(rotated, -0.5);
    if (std::abs(unrotated.I - scattered.I) > 1.0e-9 ||
        std::abs(unrotated.Q - scattered.Q) > 1.0e-9 ||
        std::abs(unrotated.U - scattered.U) > 1.0e-9 ||
        std::abs(unrotated.V - scattered.V) > 1.0e-9) {
        std::cerr << "Stokes rotation is not reversible.\n";
        return 1;
    }

    return 0;
}
