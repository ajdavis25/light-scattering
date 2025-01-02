// SurfaceReflection.cpp
#include "SurfaceReflection.hpp"
#include <cmath>
#include <random>
#include <algorithm>

static std::mt19937 rngSurf(1234);
static std::uniform_real_distribution<double> uni(0.0,1.0);

ReflectionResult reflectLambertian(double albedo)
{
    ReflectionResult result;
    double mu = uni(rngSurf); // cos(theta) in [0..1]
    double phi= 2.0*M_PI*uni(rngSurf);
    double sinTheta = std::sqrt(std::max(0.0,1.0 - mu*mu));

    // z is "up"
    result.dir.x = sinTheta*std::cos(phi);
    result.dir.y = sinTheta*std::sin(phi);
    result.dir.z = mu;

    result.weightMultiplier = albedo;
    return result;
}
