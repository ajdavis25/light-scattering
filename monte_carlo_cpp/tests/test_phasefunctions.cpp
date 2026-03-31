#include "PhaseFunctions.hpp"

#include <cmath>
#include <iostream>
#include <random>

int main()
{
    constexpr double pi = 3.14159265358979323846;
    const int samples = 4000;
    double integral = 0.0;
    for (int index = 0; index < samples; ++index) {
        const double mu = -1.0 + (2.0 * index + 1.0) / samples;
        integral += 2.0 * pi * rayleighPhase(mu) * (2.0 / samples);
    }

    if (std::abs(integral - 1.0) > 1.0e-3) {
        std::cerr << "Rayleigh phase failed normalization: " << integral << "\n";
        return 1;
    }

    std::mt19937 rng(1234);
    for (int index = 0; index < 100; ++index) {
        const ScatteringSample sample = sampleRayleighDirection(rng);
        if (sample.cos_theta < -1.0 || sample.cos_theta > 1.0 || sample.pdf <= 0.0) {
            std::cerr << "Invalid Rayleigh sample returned.\n";
            return 1;
        }
    }

    AerosolPhaseMatrixTable table;
    table.loadCsv("../data/optics/aerosol_phase_matrix_reference.csv");
    const PhaseMatrixCoefficients blue = table.coefficients(350.0, 0.0);
    const PhaseMatrixCoefficients red = table.coefficients(800.0, 0.0);
    const PhaseMatrixCoefficients green = table.coefficients(525.0, 0.0);

    if (std::abs(blue.f11 - red.f11) < 1.0e-8) {
        std::cerr << "Aerosol phase table should vary with wavelength.\n";
        return 1;
    }

    if (green.f11 <= std::min(blue.f11, red.f11) || green.f11 >= std::max(blue.f11, red.f11)) {
        std::cerr << "Interpolated aerosol coefficients should lie between tabulated wavelengths.\n";
        return 1;
    }

    const ScatteringSample aerosolSample = table.sampleDirection(525.0, rng);
    if (aerosolSample.cos_theta < -1.0 || aerosolSample.cos_theta > 1.0 || aerosolSample.pdf <= 0.0) {
        std::cerr << "Invalid aerosol sample returned.\n";
        return 1;
    }

    return 0;
}
