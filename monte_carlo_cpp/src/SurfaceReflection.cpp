#include "SurfaceReflection.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {
constexpr double PI = 3.14159265358979323846;

std::string trim(const std::string &value)
{
    const std::string whitespace = " \t\r\n";
    const std::size_t begin = value.find_first_not_of(whitespace);
    if (begin == std::string::npos) {
        return "";
    }
    const std::size_t end = value.find_last_not_of(whitespace);
    return value.substr(begin, end - begin + 1);
}

std::vector<std::string> splitCsvLine(const std::string &line)
{
    std::vector<std::string> tokens;
    std::stringstream stream(line);
    std::string token;
    while (std::getline(stream, token, ',')) {
        tokens.push_back(trim(token));
    }
    return tokens;
}

double interpolateTableValue(const std::vector<std::pair<double, double>> &table, double wavelength_nm)
{
    if (table.empty()) {
        return 0.15;
    }
    if (wavelength_nm <= table.front().first) {
        return table.front().second;
    }
    if (wavelength_nm >= table.back().first) {
        return table.back().second;
    }

    for (std::size_t index = 1; index < table.size(); ++index) {
        if (wavelength_nm <= table[index].first) {
            const double wave0 = table[index - 1].first;
            const double wave1 = table[index].first;
            const double value0 = table[index - 1].second;
            const double value1 = table[index].second;
            const double t = (wavelength_nm - wave0) / std::max(1.0e-12, wave1 - wave0);
            return value0 + t * (value1 - value0);
        }
    }

    return table.back().second;
}

std::string lowercase(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

double coxMunkSlopeVariance(double windSpeedMps)
{
    return std::max(1.0e-4, 0.003 + 0.00512 * std::max(0.0, windSpeedMps));
}

double fresnelUnpolarizedReflectance(double cosThetaI, double refractiveIndex)
{
    const double clampedCosThetaI = std::clamp(cosThetaI, 0.0, 1.0);
    const double sinThetaTSquared = (1.0 - clampedCosThetaI * clampedCosThetaI) / (refractiveIndex * refractiveIndex);
    if (sinThetaTSquared >= 1.0) {
        return 1.0;
    }

    const double cosThetaT = std::sqrt(std::max(0.0, 1.0 - sinThetaTSquared));
    const double rs = (clampedCosThetaI - refractiveIndex * cosThetaT) /
        std::max(1.0e-12, clampedCosThetaI + refractiveIndex * cosThetaT);
    const double rp = (refractiveIndex * clampedCosThetaI - cosThetaT) /
        std::max(1.0e-12, refractiveIndex * clampedCosThetaI + cosThetaT);
    return 0.5 * (rs * rs + rp * rp);
}

double coxMunkBrdf(
    const Vec3 &surfaceNormal,
    const Vec3 &sunDirection,
    const Vec3 &viewDirection,
    double windSpeedMps
)
{
    const double muIn = std::max(0.0, dot(surfaceNormal, sunDirection));
    const double muOut = std::max(0.0, dot(surfaceNormal, viewDirection));
    if (muIn <= 0.0 || muOut <= 0.0) {
        return 0.0;
    }

    const Vec3 halfVector = normalize(sunDirection + viewDirection);
    const double nh = std::max(0.0, dot(surfaceNormal, halfVector));
    const double slopeVariance = coxMunkSlopeVariance(windSpeedMps);
    const double tanThetaHSquared = (1.0 - nh * nh) / std::max(1.0e-12, nh * nh);
    const double slopePdf = std::exp(-tanThetaHSquared / slopeVariance) /
        (PI * slopeVariance * std::max(1.0e-12, nh * nh * nh * nh));
    const double fresnel = fresnelUnpolarizedReflectance(std::max(0.0, dot(halfVector, sunDirection)), 1.333);
    return fresnel * slopePdf / std::max(1.0e-12, 4.0 * muIn * muOut);
}
}

void SurfaceModel::loadAlbedoCsv(const std::string &filename)
{
    table_.clear();

    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("Unable to open surface albedo table: " + filename);
    }

    std::string line;
    bool headerSkipped = false;
    while (std::getline(input, line)) {
        const std::string stripped = trim(line);
        if (stripped.empty() || stripped[0] == '#') {
            continue;
        }
        if (!headerSkipped) {
            headerSkipped = true;
            continue;
        }

        const std::vector<std::string> values = splitCsvLine(stripped);
        if (values.size() < 2) {
            throw std::runtime_error("Surface albedo row must have 2 columns: " + stripped);
        }
        table_.emplace_back(std::stod(values[0]), std::stod(values[1]));
    }

    if (table_.empty()) {
        throw std::runtime_error("Surface albedo table contained no rows: " + filename);
    }

    std::sort(table_.begin(), table_.end(), [](const auto &lhs, const auto &rhs) { return lhs.first < rhs.first; });
}

void SurfaceModel::setModel(const std::string &model_name)
{
    model_name_ = lowercase(trim(model_name));
}

void SurfaceModel::setDefaultAlbedo(double value)
{
    default_albedo_ = std::clamp(value, 0.0, 1.0);
}

void SurfaceModel::setOceanWindSpeed(double value)
{
    ocean_wind_speed_m_s_ = std::max(0.0, value);
}

void SurfaceModel::loadParameterJson(const std::string &filename)
{
    if (filename.empty()) {
        return;
    }

    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("Unable to open surface parameter JSON: " + filename);
    }

    const std::string content(
        (std::istreambuf_iterator<char>(input)),
        std::istreambuf_iterator<char>()
    );

    const auto readNumber = [&](const std::string &key, double fallback) {
        const std::string token = "\"" + key + "\"";
        const std::size_t keyPos = content.find(token);
        if (keyPos == std::string::npos) {
            return fallback;
        }
        const std::size_t colonPos = content.find(':', keyPos + token.size());
        if (colonPos == std::string::npos) {
            return fallback;
        }
        const std::size_t start = content.find_first_of("-0123456789.", colonPos + 1);
        if (start == std::string::npos) {
            return fallback;
        }
        const std::size_t end = content.find_first_not_of("0123456789+-.eE", start);
        return std::stod(content.substr(start, end - start));
    };

    ocean_wind_speed_m_s_ = readNumber("ocean_wind_speed_m_s", ocean_wind_speed_m_s_);
    default_albedo_ = std::clamp(readNumber("default_albedo", default_albedo_), 0.0, 1.0);
}

double SurfaceModel::albedo(double wavelength_nm) const
{
    if (table_.empty()) {
        return default_albedo_;
    }
    return interpolateTableValue(table_, wavelength_nm);
}

double SurfaceModel::directSolarRadiance(
    const Vec3 &surfaceNormal,
    const Vec3 &sunDirection,
    const Vec3 &viewDirection,
    double wavelength_nm,
    double directSolarIrradiance
) const
{
    const double muIn = std::max(0.0, dot(surfaceNormal, sunDirection));
    const double muOut = std::max(0.0, dot(surfaceNormal, viewDirection));
    if (muIn <= 0.0 || muOut <= 0.0) {
        return 0.0;
    }

    if (model_name_ == "coxmunk_ocean") {
        return directSolarIrradiance * muIn *
            coxMunkBrdf(surfaceNormal, sunDirection, viewDirection, ocean_wind_speed_m_s_);
    }

    return albedo(wavelength_nm) * directSolarIrradiance * muIn / PI;
}

ReflectionResult sampleLambertianReflection(
    const Vec3 &surfaceNormal,
    double spectralAlbedo,
    std::mt19937 &rng
)
{
    std::uniform_real_distribution<double> uniform01(0.0, 1.0);
    const double mu = std::sqrt(uniform01(rng));
    const double phi = 2.0 * PI * uniform01(rng);
    const double sinTheta = std::sqrt(std::max(0.0, 1.0 - mu * mu));

    Vec3 tangent = normalize(cross(surfaceNormal, {0.0, 0.0, 1.0}));
    if (norm(tangent) < 1.0e-12) {
        tangent = normalize(cross(surfaceNormal, {0.0, 1.0, 0.0}));
    }
    const Vec3 bitangent = normalize(cross(surfaceNormal, tangent));

    ReflectionResult result {};
    result.dir = normalize(
        tangent * (sinTheta * std::cos(phi)) +
        bitangent * (sinTheta * std::sin(phi)) +
        surfaceNormal * mu
    );
    result.weightMultiplier = std::clamp(spectralAlbedo, 0.0, 1.0);
    return result;
}
