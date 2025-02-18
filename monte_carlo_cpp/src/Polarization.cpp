/******************************************************
 *               Polarization.cpp
 ******************************************************/
#include "Polarization.hpp"
#include <cmath>
#include <algorithm>

/**
 * Simple unpolarized initialization: I = intensity, Q=U=V=0
 */
StokesVector initUnpolarized(double intensity)
{
    return {intensity, 0.0, 0.0, 0.0};
}

/**
 * Identity transform => no change
 */
StokesVector applyIdentityMueller(const StokesVector &in)
{
    return in;
}

/**
 * Rayleigh scattering Mueller matrix, 4×4, for scattering angle theta. 
 * This is a typical representation for atmospheric Rayleigh 
 * ignoring any rotation of reference plane if needed. 
 *
 * We do a direct partial matrix multiplication for an "incoming" 
 * Stokes vector. We assume the scattering plane is properly oriented 
 * so that U_in=0 if we define reference plane along Q.
 */
StokesVector applyRayleighMueller(const StokesVector &in, double cosTheta)
{
    // For a physically correct matrix, we often see something like:
    // factor = 3/4, 
    // M(0,0)= (1+cos^2θ), M(0,1)= -(1 - cos^2θ)/2, etc.
    // See references: e.g. Bohren & Huffman, Chandrasekhar, etc.

    double cos2 = cosTheta*cosTheta;
    double onePlusCos2 = (1.0 + cos2);
    double oneMinusCos2= (1.0 - cos2);
    // The standard factor for unpolarized incident is (3/4)
    double factor = 0.75;  // i.e. 3/4

    // We'll define:
    // I_out = factor[ (1+cos^2θ)*I_in  - (1-cos^2θ)/2 * Q_in ]
    // Q_out = factor[ -(1-cos^2θ)/2 * I_in + (1+cos^2θ)*Q_in ]
    // U_out = factor[  2cosθ * U_in ] 
    // V_out = factor[  2cosθ * V_in ] 
    // This approach is a common one for a scattering plane orientation 
    // where U_in=0 if the incoming beam was unpolarized.

    double Iout = factor * ( onePlusCos2* in.I - 0.5*oneMinusCos2* in.Q );
    double Qout = factor * (-0.5*oneMinusCos2* in.I + onePlusCos2* in.Q );
    double Uout = factor * (2.0*cosTheta * in.U);
    double Vout = factor * (2.0*cosTheta * in.V);

    return { Iout, Qout, Uout, Vout };
}

/**
 * Mie scattering or aerosol scattering often has a complicated 4×4 matrix
 * that depends on size distribution, angle, etc.
 * We'll do a trivial approach: partial depolarization => reduce (Q,U,V).
 */
StokesVector applyMieMueller(const StokesVector &in, double cosTheta)
{
    // If we had a LUT, we'd do M= MieMuellerLUT(cosTheta, wave, etc.),
    // then out= M * in. For demonstration, let's do a partial factor:
    double factor = 0.7; // more strongly depolarizing than your old 0.8

    StokesVector out;
    out.I = in.I; // keep same intensity for now
    out.Q = in.Q * factor;
    out.U = in.U * factor;
    out.V = in.V * factor;

    // More advanced approach => table or expansions
    return out;
}
