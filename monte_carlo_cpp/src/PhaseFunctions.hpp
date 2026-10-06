#ifndef PHASEFUNCTIONS_HPP
#define PHASEFUNCTIONS_HPP

#include <random>
#include <string>
#include <vector>

struct PhaseMatrixCoefficients
{
    double f11;
    double f12;
    double f22;
    double f33;
    double f34;
    double f44;
};

struct ScatteringSample
{
    double cos_theta;
    double azimuth_rad;
    double pdf;
};

class AerosolPhaseMatrixTable
{
public:
    struct AngleEntry
    {
        double angle_deg;
        PhaseMatrixCoefficients coeffs;
        double cdf;
    };

    void loadCsv(const std::string &filename);
    PhaseMatrixCoefficients coefficients(double wavelength_nm, double cos_theta) const;
    double phasePdf(double wavelength_nm, double cos_theta) const;
    ScatteringSample sampleDirection(double wavelength_nm, std::mt19937 &rng) const;

private:
    std::vector<AngleEntry> interpolatedTable(double wavelength_nm) const;
    std::vector<double> wavelengths_nm_;
    std::vector<std::vector<AngleEntry>> tables_;
};

double rayleighPhase(double cosTheta);
double rayleighPolarizationFraction(double cosTheta);
PhaseMatrixCoefficients rayleighPhaseMatrix(double cosTheta);
// Depolarized (King-factor) Rayleigh phase matrix after Hansen & Travis
// (1974): delta = (1-rho)/(1+rho/2) scales the anisotropic part, an isotropic
// unpolarized term (1-delta)/(4 pi) keeps f11 normalized, and f44 carries the
// extra delta' = (1-2 rho)/(1-rho). rho = 0 reproduces the pure matrix
// bit-identically.
PhaseMatrixCoefficients rayleighPhaseMatrix(double cosTheta, double depolarizationFactor);
ScatteringSample sampleRayleighDirection(std::mt19937 &rng);
double henyeyGreenstein(double cosTheta, double g);

#endif // PHASEFUNCTIONS_HPP
