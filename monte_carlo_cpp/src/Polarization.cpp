// Polarization.cpp
#include "Polarization.hpp"

StokesVector initUnpolarized(double intensity)
{
    StokesVector sv;
    sv.I = intensity;
    sv.Q = 0.0;
    sv.U = 0.0;
    sv.V = 0.0;
    return sv;
}

StokesVector applyIdentityMueller(const StokesVector &in)
{
    // No change for identity
    return in;
}

// TODO: Implement real Rayleigh or Mie Mueller matrix
// to transform (I,Q,U,V) at scattering event
