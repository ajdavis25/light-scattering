#include "SurfaceReflection.hpp"
#include "Vec3.hpp"

#include <iostream>
#include <random>

int main()
{
    std::mt19937 rng(42);
    const Vec3 normal {0.0, 0.0, 1.0};
    for (int index = 0; index < 128; ++index) {
        const ReflectionResult result = sampleLambertianReflection(normal, 0.2, rng);
        if (dot(result.dir, normal) < -1.0e-12) {
            std::cerr << "Lambertian reflection sampled below the surface.\n";
            return 1;
        }
        if (result.weightMultiplier < 0.0 || result.weightMultiplier > 1.0) {
            std::cerr << "Lambertian reflection returned invalid albedo weight.\n";
            return 1;
        }
    }

    SurfaceModel surface;
    surface.setModel("coxmunk_ocean");
    surface.setOceanWindSpeed(5.0);
    const double oceanRadiance = surface.directSolarRadiance(
        normal,
        normalize(Vec3 {0.0, 0.0, 1.0}),
        normalize(Vec3 {0.0, 0.0, 1.0}),
        550.0,
        1.0
    );
    if (oceanRadiance < 0.0) {
        std::cerr << "Cox-Munk direct solar radiance should remain non-negative.\n";
        return 1;
    }

    return 0;
}
