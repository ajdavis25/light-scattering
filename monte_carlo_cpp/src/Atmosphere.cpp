// Atmosphere.cpp
#include "Atmosphere.hpp"
#include <cmath>

// Example scale height for density
static const double H_SCALE = 8000.0; 
static const double ABS_0   = 1e-5;   // baseline absorption at ground
static const double H_SCALE_ABS = 7000.0; // scale for absorption

double atmosphericDensity(double alt_m)
{
    if(alt_m < 0.0) return 0.0;
    // simple exponential
    return std::exp(-alt_m / H_SCALE);
}

double absorptionCoefficient(double alt_m)
{
    // If below ground
    if(alt_m < 0.0) return 999999.0; 
    // simple exponential, you can add more advanced formula
    return ABS_0 * std::exp(-alt_m / H_SCALE_ABS);
}
