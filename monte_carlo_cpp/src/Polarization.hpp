// Polarization.hpp
#ifndef POLARIZATION_HPP
#define POLARIZATION_HPP

struct StokesVector
{
    double I, Q, U, V;
};

StokesVector initUnpolarized(double intensity);

// Example applying identity
StokesVector applyIdentityMueller(const StokesVector &in);

#endif
