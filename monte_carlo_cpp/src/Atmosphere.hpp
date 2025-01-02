// Atmosphere.hpp
#ifndef ATMOSPHERE_HPP
#define ATMOSPHERE_HPP

/**
 * Returns a dimensionless "density" at altitude (m).
 * For example: exponential atmosphere.
 */
double atmosphericDensity(double altitude_m);

/**
 * Returns an approximate absorption coefficient at altitude (m).
 * Could represent O2, O3, H2O, etc.
 */
double absorptionCoefficient(double altitude_m);

#endif
