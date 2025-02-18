// PhaseFunctions.hpp
#ifndef PHASEFUNCTIONS_HPP
#define PHASEFUNCTIONS_HPP

double rayleighPhase(double cosTheta);
double henyeyGreenstein(double cosTheta, double g);

/**
 * For Mie scattering, we might do table lookups 
 * from angle->phaseFunction(angle,wavelength). 
 * We'll store or read from a 2D LUT if we want advanced logic.
 */
double miePhase(double cosTheta, double wavelength_nm);

/**
 * Combine Rayleigh & Mie or Rayleigh & HG with fraction 
 * based on local ratio of densities.
 */
double combinedPhase(double cosTheta, double fracRay, double fracMie);

#endif
