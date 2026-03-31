#include "Polarization.hpp"

#include <algorithm>
#include <cmath>

namespace {
MuellerMatrix makeZeroMatrix()
{
    MuellerMatrix matrix {};
    for (auto &row : matrix.m) {
        for (double &value : row) {
            value = 0.0;
        }
    }
    return matrix;
}
}

StokesVector initUnpolarized(double intensity)
{
    return {intensity, 0.0, 0.0, 0.0};
}

MuellerMatrix identityMueller()
{
    MuellerMatrix matrix = makeZeroMatrix();
    for (int index = 0; index < 4; ++index) {
        matrix.m[index][index] = 1.0;
    }
    return matrix;
}

MuellerMatrix rotationMueller(double angle_rad)
{
    MuellerMatrix matrix = identityMueller();
    const double cos2 = std::cos(2.0 * angle_rad);
    const double sin2 = std::sin(2.0 * angle_rad);

    matrix.m[1][1] = cos2;
    matrix.m[1][2] = sin2;
    matrix.m[2][1] = -sin2;
    matrix.m[2][2] = cos2;
    return matrix;
}

MuellerMatrix scaleMueller(const MuellerMatrix &matrix, double scalar)
{
    MuellerMatrix scaled = matrix;
    for (auto &row : scaled.m) {
        for (double &value : row) {
            value *= scalar;
        }
    }
    return scaled;
}

MuellerMatrix multiply(const MuellerMatrix &lhs, const MuellerMatrix &rhs)
{
    MuellerMatrix result = makeZeroMatrix();
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            for (int inner = 0; inner < 4; ++inner) {
                result.m[row][col] += lhs.m[row][inner] * rhs.m[inner][col];
            }
        }
    }
    return result;
}

StokesVector multiply(const MuellerMatrix &matrix, const StokesVector &vector)
{
    const double input[4] = {vector.I, vector.Q, vector.U, vector.V};
    double output[4] = {0.0, 0.0, 0.0, 0.0};
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            output[row] += matrix.m[row][col] * input[col];
        }
    }
    return {output[0], output[1], output[2], output[3]};
}

StokesVector rotateStokes(const StokesVector &vector, double angle_rad)
{
    return multiply(rotationMueller(angle_rad), vector);
}

MuellerMatrix rayleighMuellerMatrix(double cosTheta)
{
    const PhaseMatrixCoefficients coefficients = rayleighPhaseMatrix(cosTheta);
    return aerosolMuellerMatrix(coefficients);
}

MuellerMatrix aerosolMuellerMatrix(const PhaseMatrixCoefficients &coefficients)
{
    MuellerMatrix matrix = makeZeroMatrix();
    matrix.m[0][0] = coefficients.f11;
    matrix.m[0][1] = coefficients.f12;
    matrix.m[1][0] = coefficients.f12;
    matrix.m[1][1] = coefficients.f22;
    matrix.m[2][2] = coefficients.f33;
    matrix.m[2][3] = coefficients.f34;
    matrix.m[3][2] = -coefficients.f34;
    matrix.m[3][3] = coefficients.f44;
    return matrix;
}

StokesVector applyRayleighMueller(const StokesVector &in, double cosTheta)
{
    return multiply(rayleighMuellerMatrix(cosTheta), in);
}

StokesVector applyAerosolMueller(const StokesVector &in, const PhaseMatrixCoefficients &coefficients)
{
    return multiply(aerosolMuellerMatrix(coefficients), in);
}

double degreeOfLinearPolarization(const StokesVector &vector)
{
    if (vector.I <= 0.0) {
        return 0.0;
    }
    return std::sqrt(vector.Q * vector.Q + vector.U * vector.U) / vector.I;
}

double angleOfLinearPolarizationRad(const StokesVector &vector)
{
    return 0.5 * std::atan2(vector.U, vector.Q);
}

bool isPhysicallyValid(const StokesVector &vector)
{
    const double polarizedMagnitude = std::sqrt(
        vector.Q * vector.Q + vector.U * vector.U + vector.V * vector.V
    );
    return vector.I >= -1.0e-12 && polarizedMagnitude <= vector.I + 1.0e-12;
}
