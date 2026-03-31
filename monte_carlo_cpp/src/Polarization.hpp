#ifndef POLARIZATION_HPP
#define POLARIZATION_HPP

#include "PhaseFunctions.hpp"

struct StokesVector
{
    double I;
    double Q;
    double U;
    double V;
};

struct MuellerMatrix
{
    double m[4][4];
};

StokesVector initUnpolarized(double intensity);
MuellerMatrix identityMueller();
MuellerMatrix rotationMueller(double angle_rad);
MuellerMatrix scaleMueller(const MuellerMatrix &matrix, double scalar);
MuellerMatrix multiply(const MuellerMatrix &lhs, const MuellerMatrix &rhs);
StokesVector multiply(const MuellerMatrix &matrix, const StokesVector &vector);
StokesVector rotateStokes(const StokesVector &vector, double angle_rad);

MuellerMatrix rayleighMuellerMatrix(double cosTheta);
MuellerMatrix aerosolMuellerMatrix(const PhaseMatrixCoefficients &coefficients);
StokesVector applyRayleighMueller(const StokesVector &in, double cosTheta);
StokesVector applyAerosolMueller(const StokesVector &in, const PhaseMatrixCoefficients &coefficients);

double degreeOfLinearPolarization(const StokesVector &vector);
double angleOfLinearPolarizationRad(const StokesVector &vector);
bool isPhysicallyValid(const StokesVector &vector);

#endif // POLARIZATION_HPP
