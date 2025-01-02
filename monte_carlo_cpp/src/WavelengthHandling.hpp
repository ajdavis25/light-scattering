// WavelengthHandling.hpp
#ifndef WAVELENGTHHANDLING_HPP
#define WAVELENGTHHANDLING_HPP

/**
 * Returns scattering coefficient at altitude & wavelength
 * E.g. Rayleigh ~ wave^-4, plus an altitude factor
 */
double scatteringCoefficient(double altitude_m, double wavelength_nm);

#endif
