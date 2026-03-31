#ifndef SURFACEREFLECTION_HPP
#define SURFACEREFLECTION_HPP

#include "Vec3.hpp"

#include <random>
#include <string>
#include <utility>
#include <vector>

struct ReflectionResult
{
    Vec3 dir;
    double weightMultiplier;
};

class SurfaceModel
{
public:
    void loadAlbedoCsv(const std::string &filename);
    void setModel(const std::string &model_name);
    void setDefaultAlbedo(double value);
    void setOceanWindSpeed(double value);
    void loadParameterJson(const std::string &filename);
    double albedo(double wavelength_nm) const;
    double directSolarRadiance(
        const Vec3 &surfaceNormal,
        const Vec3 &sunDirection,
        const Vec3 &viewDirection,
        double wavelength_nm,
        double directSolarIrradiance
    ) const;

private:
    std::string model_name_ = "lambertian_land";
    double default_albedo_ = 0.15;
    double ocean_wind_speed_m_s_ = 5.0;
    std::vector<std::pair<double, double>> table_;
};

ReflectionResult sampleLambertianReflection(
    const Vec3 &surfaceNormal,
    double spectralAlbedo,
    std::mt19937 &rng
);

#endif // SURFACEREFLECTION_HPP
