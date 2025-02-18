// SurfaceReflection.cpp
#include "SurfaceReflection.hpp"
#include <random>
#include <cmath>
#include <algorithm>

static std::mt19937 rngSurf(1234);
static std::uniform_real_distribution<double> uni(0.0,1.0);

ReflectionResult reflectLambertian(double albedo)
{
    ReflectionResult rr;
    double mu  = uni(rngSurf);
    double phi = 2.0*M_PI*uni(rngSurf);
    double st  = std::sqrt(std::max(0.0, 1.0 - mu*mu));
    rr.dir = {st*std::cos(phi), st*std::sin(phi), mu};
    rr.weightMultiplier = albedo;
    return rr;
}

ReflectionResult reflectFresnelWater(const Vec3 &incoming)
{
    // Very rough approach:
    // 1) find cos(incAngle) = -incoming.z if z is normal 
    double cosi = -incoming.z;
    // water ref index
    double n1=1.0, n2=1.33;
    // TODO: compute Fresnel reflect or refract
    // For demonstration:
    double R = 0.02; // placeholder
    ReflectionResult rr;
    if(uni(rngSurf) < R)
    {
        // reflect specularly
        rr.dir = {incoming.x, incoming.y, -incoming.z};
        rr.weightMultiplier = R;
    }
    else
    {
        // transmit => might go below surface, or partial
        // if we do a real ocean, that means photon is lost or scattered below
        rr.dir = {0,0,0}; 
        rr.weightMultiplier = 0.0; // kill or partial
    }
    return rr;
}
