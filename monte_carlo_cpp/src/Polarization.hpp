#ifndef POLARIZATION_HPP
#define POLARIZATION_HPP

/**
 * StokesVector: (I, Q, U, V)
 * 
 * - I : total intensity
 * - Q, U: linear polarization components
 * - V : circular polarization (often negligible in atmospheric scattering)
 */
struct StokesVector
{
    double I;  
    double Q;  
    double U;
    double V;
};

/**
 * Initialize unpolarized light with total intensity => Q=U=V=0.
 */
StokesVector initUnpolarized(double intensity);

/**
 * A trivial identity Mueller transform (for testing).
 */
StokesVector applyIdentityMueller(const StokesVector &in);

/**
 * Apply a more realistic Rayleigh Mueller transform for scattering angle theta
 * for an incoming Stokes vector "in".
 * 
 * This code uses a standard 4×4 matrix for Rayleigh scattering of unpolarized or partially polarized light.
 * Reference matrix form can be found in many atmospheric optics resources.
 */
StokesVector applyRayleighMueller(const StokesVector &in, double cosTheta);

/**
 * Mie or aerosol Mueller. 
 * Typically a 4×4 matrix from a LUT or from an approximate formula. 
 * We keep a placeholder with partial depolarization for demonstration.
 */
StokesVector applyMieMueller(const StokesVector &in, double cosTheta);

#endif // POLARIZATION_HPP
