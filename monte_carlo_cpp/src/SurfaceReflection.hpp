// SurfaceReflection.hpp
#ifndef SURFACEREFLECTION_HPP
#define SURFACEREFLECTION_HPP

struct Vec3
{
    double x,y,z;
};

struct ReflectionResult
{
    Vec3 dir;
    double weightMultiplier;
};

// Simple lambertian reflection
ReflectionResult reflectLambertian(double albedo);

// TODO: function for BRDF-based reflection if needed

#endif
