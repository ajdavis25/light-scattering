// SurfaceReflection.hpp
#ifndef SURFACEREFLECTION_HPP
#define SURFACEREFLECTION_HPP

#include "Vec3.hpp" // Use the single definition

struct ReflectionResult
{
    Vec3 dir;
    double weightMultiplier;
};

ReflectionResult reflectLambertian(double albedo);

/**
 * Attempt a simple Fresnel reflection if we have 
 * water index of refraction ~1.33 or so.
 */
ReflectionResult reflectFresnelWater(const Vec3 &incoming);

#endif
