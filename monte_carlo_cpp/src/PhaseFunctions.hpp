// PhaseFunctions.hpp
#ifndef PHASEFUNCTIONS_HPP
#define PHASEFUNCTIONS_HPP

// Rayleigh phase function (unpolarized)
double rayleighPhase(double cosTheta);

// Henyey-Greenstein aerosol
double henyeyGreenstein(double cosTheta, double g);

// Combined, if you do fractionRay for Rayleigh vs. aerosol
double combinedPhase(double cosTheta, double fracRay, double g);

#endif
