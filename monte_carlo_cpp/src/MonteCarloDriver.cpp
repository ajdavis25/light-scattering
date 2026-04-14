#include "MonteCarloDriver.hpp"

#include "Atmosphere.hpp"
#include "PhaseFunctions.hpp"
#include "Polarization.hpp"
#include "SurfaceReflection.hpp"
#include "Vec3.hpp"
#include "WavelengthHandling.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
constexpr double PI = 3.14159265358979323846;
constexpr double R_EARTH_M = 6.371e6;
constexpr unsigned long long HASH_OFFSET_BASIS = 1469598103934665603ull;
constexpr unsigned long long HASH_PRIME = 1099511628211ull;
constexpr double ACTIVE_SPECTRAL_WEIGHT_FLOOR = 1.0e-6;

struct SolarGeometry
{
    double zenith_deg = 90.0;
    double azimuth_deg = 270.0;
    Vec3 direction {0.0, 0.0, 1.0};
};

struct LocalBasis
{
    Vec3 north;
    Vec3 east;
    Vec3 up;
};

struct RunningMoments
{
    int count = 0;
    StokesVector mean {0.0, 0.0, 0.0, 0.0};
    StokesVector m2 {0.0, 0.0, 0.0, 0.0};

    void update(const StokesVector &sample)
    {
        ++count;
        const StokesVector delta {
            sample.I - mean.I,
            sample.Q - mean.Q,
            sample.U - mean.U,
            sample.V - mean.V,
        };
        mean.I += delta.I / count;
        mean.Q += delta.Q / count;
        mean.U += delta.U / count;
        mean.V += delta.V / count;

        const StokesVector delta2 {
            sample.I - mean.I,
            sample.Q - mean.Q,
            sample.U - mean.U,
            sample.V - mean.V,
        };
        m2.I += delta.I * delta2.I;
        m2.Q += delta.Q * delta2.Q;
        m2.U += delta.U * delta2.U;
        m2.V += delta.V * delta2.V;
    }

    StokesVector variance() const
    {
        if (count < 2) {
            return {0.0, 0.0, 0.0, 0.0};
        }
        return {
            m2.I / (count - 1),
            m2.Q / (count - 1),
            m2.U / (count - 1),
            m2.V / (count - 1),
        };
    }
};

struct SecondOrderBreakdown
{
    StokesVector total {0.0, 0.0, 0.0, 0.0};
    StokesVector rr {0.0, 0.0, 0.0, 0.0};
    StokesVector ar {0.0, 0.0, 0.0, 0.0};
    StokesVector ra {0.0, 0.0, 0.0, 0.0};
    StokesVector aa {0.0, 0.0, 0.0, 0.0};
};

struct SingleScatterBreakdown
{
    StokesVector rayleigh {0.0, 0.0, 0.0, 0.0};
    StokesVector aerosol {0.0, 0.0, 0.0, 0.0};
};

struct DirectionEstimate
{
    StokesVector first_order {0.0, 0.0, 0.0, 0.0};
    SecondOrderBreakdown second_order {};
    StokesVector higher_order {0.0, 0.0, 0.0, 0.0};
    StokesVector higher_variance {0.0, 0.0, 0.0, 0.0};
    struct TimingSummary
    {
        double first_order_seconds = 0.0;
        std::size_t first_order_view_samples = 0;
        int first_order_steps = 0;
        int solar_disk_nodes = 0;
        double second_order_seconds = 0.0;
        double second_order_incoming_single_scatter_seconds = 0.0;
        std::size_t second_order_incoming_single_scatter_calls = 0;
        std::size_t second_order_nonzero_incoming_calls = 0;
        std::size_t second_order_view_samples = 0;
        std::size_t second_order_mu_phi_evaluations = 0;
        int second_order_view_steps = 0;
        int second_order_ray_steps = 0;
        int second_order_mu_nodes = 0;
        int second_order_phi_nodes = 0;
        double higher_order_seconds = 0.0;
        std::size_t spectral_band_count = 0;
        int mc_sample_count = 0;
    } timing {};
};

using DirectionStageCallback = std::function<void(
    SampleProgress::Stage,
    const DirectionEstimate::TimingSummary &,
    std::size_t,
    std::size_t
)>;

struct SampledInteraction
{
    bool scattered = false;
    bool absorbed = false;
    bool hit_ground = false;
    Vec3 position {0.0, 0.0, 0.0};
    LocalOpticalProperties optics {};
    bool is_rayleigh = true;
    double collision_weight = 1.0;
};

struct ContinuationSample
{
    Vec3 incoming_direction {0.0, 0.0, 1.0};
    double cos_theta = 1.0;
    double pdf = 1.0;
    double estimator_weight = 1.0;
};

struct GuidedProposalParameters
{
    double phase_fraction = 0.5;
    double cos_half_angle = 1.0;
    double polarization_fraction = 0.0;
    double polarization_mu_half_width = 0.25;
    Vec3 guide_axis {0.0, 0.0, 1.0};
    bool use_twilight_higher_order = false;
    double tangent_fraction = 0.0;
    double horizon_fraction = 0.0;
    double tangent_cos_half_angle = 1.0;
    double horizon_cos_half_angle = 1.0;
    Vec3 tangent_axis {0.0, 0.0, 1.0};
    Vec3 horizon_axis {0.0, 0.0, 1.0};
};

enum class TwilightGuidedComponent
{
    none = -1,
    phase = 0,
    tangent = 1,
    horizon = 2,
};

struct SolarDiskSample
{
    Vec3 direction {0.0, 0.0, 1.0};
    double weight = 1.0;
};

struct SecondScatterQuadrature
{
    int view_steps = 1;
    int ray_steps = 1;
    int mu_nodes = 1;
    int phi_nodes = 1;
};

struct RayAltitudeQuadrature
{
    bool valid = false;
    double ds_m = 0.0;
    std::vector<double> altitudes_m;
    std::vector<double> segment_lengths_m;
};

struct SpectralOpticalState
{
    double molecular_number_density_m3 = 0.0;
    double ozone_number_density_m3 = 0.0;
    double h2o_number_density_m3 = 0.0;
    double no2_number_density_m3 = 0.0;
    double aerosol_extinction_550_m_inv = 0.0;
    double aerosol_single_scattering_albedo = 1.0;
    double aerosol_asymmetry = 0.0;
    double aerosol_scattering_angstrom_exponent = 1.6;
    double aerosol_absorption_angstrom_exponent = 1.0;
};

struct FirstOrderViewSample
{
    Vec3 position {0.0, 0.0, 0.0};
    SpectralOpticalState optical_state {};
};

struct SimulationContext
{
    Atmosphere atmosphere;
    WavelengthManager wavelengthManager;
    std::vector<SpectralBand> active_bands;
    AerosolPhaseMatrixTable phaseTable;
    SurfaceModel surface;
    SolarGeometry sun;
    std::vector<SolarDiskSample> sun_samples;
};

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

std::vector<std::string> splitSemicolonList(const std::string &value)
{
    std::vector<std::string> items;
    std::stringstream stream(value);
    std::string item;
    while (std::getline(stream, item, ';')) {
        const std::string trimmed = trim(item);
        if (!trimmed.empty()) {
            items.push_back(trimmed);
        }
    }
    return items;
}

double clamp(double value, double lower, double upper)
{
    if (value < lower) {
        return lower;
    }
    if (value > upper) {
        return upper;
    }
    return value;
}

double normalizeAzimuthDeg(double value)
{
    while (value < 0.0) {
        value += 360.0;
    }
    while (value >= 360.0) {
        value -= 360.0;
    }
    return value;
}

unsigned long long hashCombine(unsigned long long hash, unsigned long long value)
{
    hash ^= value;
    hash *= HASH_PRIME;
    return hash;
}

std::filesystem::path resolvePath(const std::filesystem::path &baseDirectory, const std::string &pathValue)
{
    const std::filesystem::path raw(pathValue);
    if (raw.is_absolute()) {
        return raw.lexically_normal();
    }
    return (baseDirectory / raw).lexically_normal();
}

std::string resolvePathList(const std::filesystem::path &baseDirectory, const std::string &pathValue)
{
    const std::vector<std::string> items = splitSemicolonList(pathValue);
    if (items.empty()) {
        return pathValue;
    }

    std::ostringstream resolved;
    for (std::size_t index = 0; index < items.size(); ++index) {
        if (index > 0) {
            resolved << ';';
        }
        resolved << resolvePath(baseDirectory, items[index]).string();
    }
    return resolved.str();
}

std::map<std::string, std::string> loadKeyValueConfig(const std::string &filename)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("Unable to open simulation config: " + filename);
    }

    std::map<std::string, std::string> values;
    std::string line;
    while (std::getline(input, line)) {
        const std::string stripped = trim(line);
        if (stripped.empty() || stripped[0] == '#') {
            continue;
        }
        const std::size_t separator = stripped.find('=');
        if (separator == std::string::npos) {
            continue;
        }
        const std::string key = trim(stripped.substr(0, separator));
        const std::string value = trim(stripped.substr(separator + 1));
        values[key] = value;
    }
    return values;
}

void setIfPresent(
    const std::map<std::string, std::string> &values,
    const std::string &key,
    std::string &target
)
{
    const auto iterator = values.find(key);
    if (iterator != values.end()) {
        target = iterator->second;
    }
}

template <typename T>
void setNumericIfPresent(
    const std::map<std::string, std::string> &values,
    const std::string &key,
    T &target
)
{
    const auto iterator = values.find(key);
    if (iterator != values.end()) {
        if constexpr (std::is_integral_v<T>) {
            target = static_cast<T>(std::stoll(iterator->second));
        } else {
            target = static_cast<T>(std::stod(iterator->second));
        }
    }
}

void setBoolIfPresent(
    const std::map<std::string, std::string> &values,
    const std::string &key,
    bool &target
)
{
    const auto iterator = values.find(key);
    if (iterator == values.end()) {
        return;
    }

    const std::string lowered = iterator->second;
    target = lowered == "1" || lowered == "true" || lowered == "True" || lowered == "TRUE";
}

std::string hashFile(const std::string &filename)
{
    std::ifstream input(filename, std::ios::binary);
    if (!input) {
        return "";
    }

    unsigned long long hash = static_cast<unsigned long long>(HASH_OFFSET_BASIS);
    char buffer[4096];
    while (input.read(buffer, sizeof(buffer)) || input.gcount() > 0) {
        for (std::streamsize index = 0; index < input.gcount(); ++index) {
            hash ^= static_cast<unsigned char>(buffer[index]);
            hash *= static_cast<unsigned long long>(HASH_PRIME);
        }
    }

    std::ostringstream output;
    output << std::hex << hash;
    return output.str();
}

double topOfAtmosphereRadius(const SimulationConfig &config)
{
    return R_EARTH_M + config.atmosphere.top_of_atmosphere_altitude_m;
}

Vec3 observerPosition(const ObserverConfig &observer)
{
    const double latitude = observer.latitude_deg * PI / 180.0;
    const double longitude = observer.longitude_deg * PI / 180.0;
    const double radius = R_EARTH_M + observer.altitude_m;
    return {
        radius * std::cos(latitude) * std::cos(longitude),
        radius * std::cos(latitude) * std::sin(longitude),
        radius * std::sin(latitude),
    };
}

LocalBasis localBasisAtPosition(const Vec3 &position)
{
    const Vec3 up = normalize(position);
    Vec3 east = normalize(cross({0.0, 0.0, 1.0}, up));
    if (norm(east) < 1.0e-12) {
        east = normalize(cross({0.0, 1.0, 0.0}, up));
    }
    const Vec3 north = normalize(cross(up, east));
    return {north, east, up};
}

Vec3 directionFromZenithAzimuth(double zenithDeg, double azimuthDeg, const LocalBasis &basis)
{
    const double zenithRad = zenithDeg * PI / 180.0;
    const double azimuthRad = azimuthDeg * PI / 180.0;
    return normalize(
        basis.north * (std::sin(zenithRad) * std::cos(azimuthRad)) +
        basis.east * (std::sin(zenithRad) * std::sin(azimuthRad)) +
        basis.up * std::cos(zenithRad)
    );
}

SolarGeometry solarGeometryFromConfig(const SimulationConfig &config)
{
    const LocalBasis basis = localBasisAtPosition(observerPosition(config.observer));
    if (config.solar.use_explicit_angles) {
        return {
            config.solar.zenith_deg,
            normalizeAzimuthDeg(config.solar.azimuth_deg),
            directionFromZenithAzimuth(config.solar.zenith_deg, normalizeAzimuthDeg(config.solar.azimuth_deg), basis),
        };
    }

    const double latitude = config.observer.latitude_deg * PI / 180.0;
    const double declination = config.solar.solar_declination_deg * PI / 180.0;
    const double hourAngle = config.solar.hour_angle_deg * PI / 180.0;
    const double cosZenith =
        std::sin(latitude) * std::sin(declination) +
        std::cos(latitude) * std::cos(declination) * std::cos(hourAngle);
    const double zenithDeg = std::acos(clamp(cosZenith, -1.0, 1.0)) * 180.0 / PI;

    const double sinZenith = std::sqrt(std::max(0.0, 1.0 - cosZenith * cosZenith));
    double azimuthDeg = 0.0;
    if (sinZenith > 1.0e-12) {
        const double sinAzimuth = -std::cos(declination) * std::sin(hourAngle) / sinZenith;
        const double cosAzimuth =
            (std::sin(declination) - std::sin(latitude) * cosZenith) /
            std::max(1.0e-12, std::cos(latitude) * sinZenith);
        azimuthDeg = std::atan2(sinAzimuth, cosAzimuth) * 180.0 / PI;
        if (azimuthDeg < 0.0) {
            azimuthDeg += 360.0;
        }
    }

    return {zenithDeg, normalizeAzimuthDeg(azimuthDeg), directionFromZenithAzimuth(zenithDeg, normalizeAzimuthDeg(azimuthDeg), basis)};
}

std::vector<SolarDiskSample> solarDiskSamples(
    const SimulationConfig &config,
    const Vec3 &sunDirection
);

std::vector<SpectralBand> activeSpectralBands(const WavelengthManager &wavelengthManager)
{
    std::vector<SpectralBand> filtered;
    for (const SpectralBand &band : wavelengthManager.getBands()) {
        if (band.solar_irradiance_w_m2_nm <= 0.0) {
            continue;
        }
        if (band.normalized_weight < ACTIVE_SPECTRAL_WEIGHT_FLOOR) {
            continue;
        }
        filtered.push_back(band);
    }
    if (filtered.empty()) {
        return wavelengthManager.getBands();
    }
    return filtered;
}

SimulationContext buildSimulationContext(const SimulationConfig &config)
{
    SimulationContext context;
    context.atmosphere.loadProfileCsv(config.atmosphere.profile_csv);
    context.wavelengthManager.configureGrid(
        config.spectral.min_wavelength_nm,
        config.spectral.max_wavelength_nm,
        config.spectral.wavelength_step_nm
    );
    context.wavelengthManager.loadSolarSpectrumCsv(config.spectral.solar_spectrum_csv);
    if (!config.spectral.instrument_response_csv.empty()) {
        context.wavelengthManager.loadInstrumentResponseCsv(config.spectral.instrument_response_csv);
    }
    if (!config.spectral.rayleigh_cross_section_csv.empty()) {
        context.wavelengthManager.loadRayleighCrossSectionCsv(config.spectral.rayleigh_cross_section_csv);
    }
    context.wavelengthManager.loadOzoneCrossSectionCsv(config.spectral.ozone_cross_section_csv);
    if (!config.spectral.o2_cross_section_csv.empty()) {
        context.wavelengthManager.loadO2CrossSectionCsv(config.spectral.o2_cross_section_csv);
    }
    if (!config.spectral.o4_cross_section_csv.empty()) {
        context.wavelengthManager.loadO4CrossSectionCsv(config.spectral.o4_cross_section_csv);
    }
    if (!config.spectral.h2o_cross_section_csv.empty()) {
        context.wavelengthManager.loadH2OCrossSectionCsv(config.spectral.h2o_cross_section_csv);
    }
    if (!config.spectral.no2_cross_section_csv.empty()) {
        context.wavelengthManager.loadNO2CrossSectionCsv(config.spectral.no2_cross_section_csv);
    }
    context.wavelengthManager.buildBands();
    context.active_bands = activeSpectralBands(context.wavelengthManager);
    context.phaseTable.loadCsv(config.spectral.aerosol_phase_matrix_csv);
    context.surface.setModel(config.surface.surface_model);
    context.surface.setDefaultAlbedo(config.surface.default_albedo);
    context.surface.setOceanWindSpeed(config.surface.ocean_wind_speed_m_s);
    context.surface.loadAlbedoCsv(config.surface.albedo_csv);
    if (!config.surface.surface_parameter_csv.empty()) {
        context.surface.loadParameterJson(config.surface.surface_parameter_csv);
    }
    context.sun = solarGeometryFromConfig(config);
    context.sun_samples = solarDiskSamples(config, context.sun.direction);
    return context;
}

double distanceToSphere(const Vec3 &origin, const Vec3 &direction, double radius)
{
    const double b = dot(origin, direction);
    const double c = dot(origin, origin) - radius * radius;
    const double discriminant = b * b - c;
    if (discriminant < 0.0) {
        return std::numeric_limits<double>::infinity();
    }

    const double root = std::sqrt(discriminant);
    const double nearDistance = -b - root;
    const double farDistance = -b + root;
    double best = std::numeric_limits<double>::infinity();
    if (nearDistance > 1.0e-6) {
        best = nearDistance;
    }
    if (farDistance > 1.0e-6) {
        best = std::min(best, farDistance);
    }
    return best;
}

double distanceToBoundary(
    const Vec3 &origin,
    const Vec3 &direction,
    double topRadius,
    bool &hitsGround
)
{
    const double toaDistance = distanceToSphere(origin, direction, topRadius);
    const double groundDistance = distanceToSphere(origin, direction, R_EARTH_M);
    hitsGround = groundDistance < toaDistance;
    return hitsGround ? groundDistance : toaDistance;
}

double altitudeM(const Vec3 &position)
{
    return norm(position) - R_EARTH_M;
}

Vec3 tangentFromAxis(const Vec3 &axis)
{
    Vec3 tangent = normalize(cross(axis, {0.0, 0.0, 1.0}));
    if (norm(tangent) < 1.0e-12) {
        tangent = normalize(cross(axis, {0.0, 1.0, 0.0}));
    }
    return tangent;
}

Vec3 directionAroundAxis(const Vec3 &axis, double cosTheta, double azimuthRad)
{
    const Vec3 tangent = tangentFromAxis(axis);
    const Vec3 bitangent = normalize(cross(axis, tangent));
    const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
    return normalize(
        tangent * (sinTheta * std::cos(azimuthRad)) +
        bitangent * (sinTheta * std::sin(azimuthRad)) +
        axis * cosTheta
    );
}

std::vector<SolarDiskSample> solarDiskSamples(
    const SimulationConfig &config,
    const Vec3 &sunDirection
)
{
    if (!config.solar.finite_solar_disk || config.solar.solar_disk_quadrature_nodes <= 1) {
        return {{sunDirection, 1.0}};
    }

    const int requestedNodes = std::max(1, config.solar.solar_disk_quadrature_nodes);
    const int ringNodes = std::max(0, requestedNodes - 1);
    const double angularRadiusRad = config.solar.solar_angular_radius_deg * PI / 180.0;
    const Vec3 tangent = tangentFromAxis(sunDirection);
    const Vec3 bitangent = normalize(cross(sunDirection, tangent));

    std::vector<SolarDiskSample> samples;
    samples.reserve(static_cast<std::size_t>(requestedNodes));
    samples.push_back({sunDirection, 1.0 / requestedNodes});
    if (ringNodes == 0 || angularRadiusRad <= 0.0) {
        return samples;
    }

    const double ringRadius = angularRadiusRad * std::sqrt(0.5);
    const double cosTheta = std::cos(ringRadius);
    const double sinTheta = std::sin(ringRadius);
    for (int index = 0; index < ringNodes; ++index) {
        const double phi = 2.0 * PI * static_cast<double>(index) / ringNodes;
        samples.push_back({
            normalize(
                tangent * (sinTheta * std::cos(phi)) +
                bitangent * (sinTheta * std::sin(phi)) +
                sunDirection * cosTheta
            ),
            1.0 / requestedNodes,
        });
    }
    return samples;
}

Vec3 referenceAxisForRay(const Vec3 &position, const Vec3 &direction)
{
    const Vec3 up = normalize(position);
    Vec3 reference = up - direction * dot(up, direction);
    if (norm(reference) < 1.0e-12) {
        const LocalBasis basis = localBasisAtPosition(position);
        reference = basis.north - direction * dot(basis.north, direction);
    }
    if (norm(reference) < 1.0e-12) {
        reference = tangentFromAxis(direction);
    }
    return normalize(reference);
}

double signedAngleAboutAxis(const Vec3 &from, const Vec3 &to, const Vec3 &axis)
{
    const Vec3 fromUnit = normalize(from);
    const Vec3 toUnit = normalize(to);
    const Vec3 axisUnit = normalize(axis);
    const double sine = dot(axisUnit, cross(fromUnit, toUnit));
    const double cosine = clamp(dot(fromUnit, toUnit), -1.0, 1.0);
    return std::atan2(sine, cosine);
}

MuellerMatrix eventMuellerMatrix(
    const Vec3 &scatteringPosition,
    const Vec3 &incomingReferencePosition,
    const Vec3 &outgoingReferencePosition,
    const Vec3 &incomingDirection,
    const Vec3 &outgoingDirection,
    bool isRayleigh,
    const PhaseMatrixCoefficients &aerosolCoefficients
)
{
    const MuellerMatrix scatterMatrix = isRayleigh
        ? rayleighMuellerMatrix(dot(incomingDirection, outgoingDirection))
        : aerosolMuellerMatrix(aerosolCoefficients);
    const Vec3 scatteringNormal = normalize(cross(incomingDirection, outgoingDirection));
    if (norm(scatteringNormal) < 1.0e-12) {
        return scatterMatrix;
    }

    const Vec3 scatteringBasisIncoming = normalize(cross(scatteringNormal, incomingDirection));
    const Vec3 scatteringBasisOutgoing = normalize(cross(scatteringNormal, outgoingDirection));
    const Vec3 referenceIncoming = referenceAxisForRay(incomingReferencePosition, incomingDirection);
    const Vec3 referenceOutgoing = referenceAxisForRay(outgoingReferencePosition, outgoingDirection);

    const double chiIn = signedAngleAboutAxis(referenceIncoming, scatteringBasisIncoming, incomingDirection);
    const double chiOut = signedAngleAboutAxis(scatteringBasisOutgoing, referenceOutgoing, outgoingDirection);

    return multiply(rotationMueller(chiOut), multiply(scatterMatrix, rotationMueller(-chiIn)));
}

double minimumAltitudeAlongSegment(const Vec3 &origin, const Vec3 &direction, double distance);
double adaptiveOpticalDepthStepM(const SimulationConfig &config, double minimumAltitudeM);
double adaptiveLineIntegralStepM(const SimulationConfig &config, double minimumAltitudeM);

std::vector<double> raySphereIntersectionDistances(
    const Vec3 &origin,
    const Vec3 &direction,
    double radius
)
{
    const double b = dot(origin, direction);
    const double c = dot(origin, origin) - radius * radius;
    const double discriminant = b * b - c;
    if (discriminant < 0.0) {
        return {};
    }

    const double sqrtDiscriminant = std::sqrt(std::max(0.0, discriminant));
    std::vector<double> distances;
    const double first = -b - sqrtDiscriminant;
    const double second = -b + sqrtDiscriminant;
    constexpr double DISTANCE_TOLERANCE_M = 1.0e-6;
    if (first > DISTANCE_TOLERANCE_M) {
        distances.push_back(first);
    }
    if (second > DISTANCE_TOLERANCE_M && std::abs(second - first) > DISTANCE_TOLERANCE_M) {
        distances.push_back(second);
    }
    return distances;
}

RayAltitudeQuadrature buildOpticalDepthQuadrature(
    const Atmosphere &atmosphere,
    const SimulationConfig &config,
    const Vec3 &origin,
    const Vec3 &direction
)
{
    RayAltitudeQuadrature quadrature {};
    bool hitsGround = false;
    const double topRadius = topOfAtmosphereRadius(config);
    const double boundaryDistance = distanceToBoundary(origin, direction, topRadius, hitsGround);
    if (!std::isfinite(boundaryDistance) || hitsGround) {
        return quadrature;
    }

    std::vector<double> segmentBoundariesM {0.0, boundaryDistance};
    for (const AtmosphereLayer &layer : atmosphere.layers()) {
        const double layerRadius = R_EARTH_M + layer.altitude_m;
        const std::vector<double> intersections = raySphereIntersectionDistances(origin, direction, layerRadius);
        for (double distance : intersections) {
            if (distance > 0.0 && distance < boundaryDistance) {
                segmentBoundariesM.push_back(distance);
            }
        }
    }

    std::sort(segmentBoundariesM.begin(), segmentBoundariesM.end());
    constexpr double BOUNDARY_TOLERANCE_M = 1.0e-3;
    segmentBoundariesM.erase(
        std::unique(
            segmentBoundariesM.begin(),
            segmentBoundariesM.end(),
            [](double lhs, double rhs) { return std::abs(lhs - rhs) <= BOUNDARY_TOLERANCE_M; }
        ),
        segmentBoundariesM.end()
    );

    quadrature.altitudes_m.reserve(segmentBoundariesM.size());
    quadrature.segment_lengths_m.reserve(segmentBoundariesM.size());
    for (std::size_t index = 1; index < segmentBoundariesM.size(); ++index) {
        const double segmentStart = segmentBoundariesM[index - 1];
        const double segmentEnd = segmentBoundariesM[index];
        const double segmentLength = segmentEnd - segmentStart;
        if (segmentLength <= BOUNDARY_TOLERANCE_M) {
            continue;
        }

        const Vec3 sample = origin + direction * (0.5 * (segmentStart + segmentEnd));
        const double altitude = altitudeM(sample);
        if (altitude < 0.0) {
            quadrature.valid = false;
            quadrature.altitudes_m.clear();
            quadrature.segment_lengths_m.clear();
            return quadrature;
        }
        quadrature.altitudes_m.push_back(altitude);
        quadrature.segment_lengths_m.push_back(segmentLength);
    }

    quadrature.valid = !quadrature.altitudes_m.empty();
    if (quadrature.valid) {
        double totalLength = 0.0;
        for (double segmentLength : quadrature.segment_lengths_m) {
            totalLength += segmentLength;
        }
        quadrature.ds_m = totalLength / static_cast<double>(quadrature.segment_lengths_m.size());
        return quadrature;
    }

    const double minimumAltitudeM = std::max(0.0, minimumAltitudeAlongSegment(origin, direction, boundaryDistance));
    const double stepLengthM = adaptiveOpticalDepthStepM(config, minimumAltitudeM);
    const int steps = std::max(
        std::max(1, config.monte_carlo.min_optical_depth_steps),
        static_cast<int>(std::ceil(boundaryDistance / stepLengthM))
    );
    quadrature.valid = true;
    quadrature.ds_m = boundaryDistance / static_cast<double>(steps);
    quadrature.altitudes_m.reserve(static_cast<std::size_t>(steps));
    quadrature.segment_lengths_m.reserve(static_cast<std::size_t>(steps));
    for (int index = 0; index < steps; ++index) {
        const Vec3 sample = origin + direction * ((index + 0.5) * quadrature.ds_m);
        const double altitude = altitudeM(sample);
        if (altitude < 0.0) {
            quadrature.valid = false;
            quadrature.altitudes_m.clear();
            quadrature.segment_lengths_m.clear();
            return quadrature;
        }
        quadrature.altitudes_m.push_back(altitude);
        quadrature.segment_lengths_m.push_back(quadrature.ds_m);
    }
    return quadrature;
}

double opticalDepthFromQuadrature(
    const Atmosphere &atmosphere,
    const WavelengthManager &wavelengthManager,
    const RayAltitudeQuadrature &quadrature,
    double wavelength_nm
)
{
    if (!quadrature.valid || quadrature.altitudes_m.empty()) {
        return std::numeric_limits<double>::infinity();
    }

    double opticalDepth = 0.0;
    for (std::size_t index = 0; index < quadrature.altitudes_m.size(); ++index) {
        const double weightM = quadrature.segment_lengths_m.empty()
            ? quadrature.ds_m
            : quadrature.segment_lengths_m[index];
        opticalDepth +=
            computeOpticalProperties(
                atmosphere,
                wavelengthManager,
                quadrature.altitudes_m[index],
                wavelength_nm
            ).extinction_m_inv * weightM;
    }
    return opticalDepth;
}

SpectralOpticalState spectralOpticalStateAtAltitude(
    const Atmosphere &atmosphere,
    double altitude_m
)
{
    const AtmosphereState state = atmosphere.stateAtAltitude(altitude_m);
    return {
        state.molecular_number_density_m3,
        state.ozone_number_density_m3,
        state.h2o_number_density_m3,
        state.no2_number_density_m3,
        state.aerosol_extinction_550_m_inv,
        state.aerosol_single_scattering_albedo,
        state.aerosol_asymmetry,
        state.aerosol_scattering_angstrom_exponent,
        state.aerosol_absorption_angstrom_exponent,
    };
}

LocalOpticalProperties opticalPropertiesFromState(
    const SpectralOpticalState &state,
    const SpectralBand &band
)
{
    constexpr double DRY_AIR_O2_MIXING_RATIO = 0.20946;

    const double rayleighScattering =
        state.molecular_number_density_m3 * band.rayleigh_cross_section_m2;
    const double o2NumberDensity = state.molecular_number_density_m3 * DRY_AIR_O2_MIXING_RATIO;
    const double aerosolScattering550 =
        state.aerosol_extinction_550_m_inv * state.aerosol_single_scattering_albedo;
    const double aerosolAbsorption550 =
        state.aerosol_extinction_550_m_inv * (1.0 - state.aerosol_single_scattering_albedo);
    const double aerosolScattering = aerosolScattering550 *
        std::pow(550.0 / band.wavelength_nm, state.aerosol_scattering_angstrom_exponent);
    const double aerosolAbsorption = aerosolAbsorption550 *
        std::pow(550.0 / band.wavelength_nm, state.aerosol_absorption_angstrom_exponent);
    const double ozoneAbsorption = state.ozone_number_density_m3 * band.ozone_cross_section_m2;
    const double o2Absorption = o2NumberDensity * band.o2_cross_section_m2;
    const double o4Absorption = o2NumberDensity * o2NumberDensity * band.o4_cross_section_m5;
    const double h2oAbsorption = state.h2o_number_density_m3 * band.h2o_cross_section_m2;
    const double no2Absorption = state.no2_number_density_m3 * band.no2_cross_section_m2;
    const double extinction =
        rayleighScattering +
        aerosolScattering +
        aerosolAbsorption +
        ozoneAbsorption +
        o2Absorption +
        o4Absorption +
        h2oAbsorption +
        no2Absorption;

    LocalOpticalProperties properties {};
    properties.rayleigh_scattering_m_inv = rayleighScattering;
    properties.aerosol_scattering_m_inv = aerosolScattering;
    properties.aerosol_absorption_m_inv = aerosolAbsorption;
    properties.ozone_absorption_m_inv = ozoneAbsorption;
    properties.o2_absorption_m_inv = o2Absorption;
    properties.o4_absorption_m_inv = o4Absorption;
    properties.h2o_absorption_m_inv = h2oAbsorption;
    properties.no2_absorption_m_inv = no2Absorption;
    properties.extinction_m_inv = extinction;
    properties.single_scattering_albedo =
        extinction > 0.0 ? (rayleighScattering + aerosolScattering) / extinction : 0.0;
    properties.aerosol_asymmetry = state.aerosol_asymmetry;
    return properties;
}

std::vector<double> opticalDepthSpectrumFromQuadrature(
    const Atmosphere &atmosphere,
    const std::vector<SpectralBand> &bands,
    const RayAltitudeQuadrature &quadrature
)
{
    std::vector<double> opticalDepths(bands.size(), std::numeric_limits<double>::infinity());
    if (!quadrature.valid || quadrature.altitudes_m.empty()) {
        return opticalDepths;
    }

    std::fill(opticalDepths.begin(), opticalDepths.end(), 0.0);
    for (std::size_t sampleIndex = 0; sampleIndex < quadrature.altitudes_m.size(); ++sampleIndex) {
        const double weightM = quadrature.segment_lengths_m.empty()
            ? quadrature.ds_m
            : quadrature.segment_lengths_m[sampleIndex];
        const SpectralOpticalState state = spectralOpticalStateAtAltitude(
            atmosphere,
            quadrature.altitudes_m[sampleIndex]
        );
        for (std::size_t bandIndex = 0; bandIndex < bands.size(); ++bandIndex) {
            opticalDepths[bandIndex] +=
                opticalPropertiesFromState(state, bands[bandIndex]).extinction_m_inv * weightM;
        }
    }
    return opticalDepths;
}

double opticalDepthToSun(
    const Atmosphere &atmosphere,
    const WavelengthManager &wavelengthManager,
    const SimulationConfig &config,
    const Vec3 &origin,
    const Vec3 &direction,
    double wavelength_nm
)
{
    const RayAltitudeQuadrature quadrature = buildOpticalDepthQuadrature(
        atmosphere,
        config,
        origin,
        direction
    );
    return opticalDepthFromQuadrature(atmosphere, wavelengthManager, quadrature, wavelength_nm);
}

double opticalDepthAlongSegment(
    const Atmosphere &atmosphere,
    const WavelengthManager &wavelengthManager,
    const Vec3 &origin,
    const Vec3 &direction,
    double distance,
    double wavelength_nm
)
{
    if (distance <= 0.0) {
        return 0.0;
    }

    const double minimumAltitudeM = std::max(0.0, minimumAltitudeAlongSegment(origin, direction, distance));
    const double stepLengthM = minimumAltitudeM < 3.0e4
        ? 5.0e2
        : (minimumAltitudeM < 6.0e4 ? 1.0e3 : 2.0e3);
    const int steps = std::max(16, static_cast<int>(std::ceil(distance / stepLengthM)));
    const double ds = distance / static_cast<double>(steps);
    double opticalDepth = 0.0;
    for (int index = 0; index < steps; ++index) {
        const Vec3 sample = origin + direction * ((index + 0.5) * ds);
        const double altitude = altitudeM(sample);
        if (altitude < 0.0) {
            return std::numeric_limits<double>::infinity();
        }
        opticalDepth += computeOpticalProperties(atmosphere, wavelengthManager, altitude, wavelength_nm).extinction_m_inv * ds;
    }
    return opticalDepth;
}

double throughputImportance(const MuellerMatrix &matrix)
{
    return std::max({
        std::abs(matrix.m[0][0]),
        std::abs(matrix.m[1][0]),
        std::abs(matrix.m[2][0]),
        std::abs(matrix.m[3][0]),
    });
}

double minimumAltitudeAlongSegment(const Vec3 &origin, const Vec3 &direction, double distance)
{
    const double sStar = clamp(-dot(origin, direction), 0.0, distance);
    return altitudeM(origin + direction * sStar);
}

double maximumAltitudeAlongSegment(const Vec3 &origin, const Vec3 &direction, double distance)
{
    return std::max(altitudeM(origin), altitudeM(origin + direction * distance));
}

double adaptiveOpticalDepthStepM(const SimulationConfig &config, double minimumAltitudeM)
{
    const double scale = std::max(1.0, config.monte_carlo.optical_depth_step_scale);
    if (minimumAltitudeM < 3.0e4) {
        return 5.0e2 * scale;
    }
    if (minimumAltitudeM < 6.0e4) {
        return 1.0e3 * scale;
    }
    return 2.0e3 * scale;
}

double adaptiveLineIntegralStepM(const SimulationConfig &config, double minimumAltitudeM)
{
    const double scale = std::max(1.0, config.monte_carlo.line_integral_step_scale);
    if (minimumAltitudeM < 3.0e4) {
        return 1.0e3 * scale;
    }
    if (minimumAltitudeM < 6.0e4) {
        return 2.0e3 * scale;
    }
    return 4.0e3 * scale;
}

double extinctionMajorantAlongSegment(
    const Atmosphere &atmosphere,
    const WavelengthManager &wavelengthManager,
    double wavelength_nm,
    double minAltitude,
    double maxAltitude
)
{
    double maxExtinction = 0.0;
    const auto updateAtAltitude = [&](double altitude) {
        if (altitude < 0.0 || altitude > atmosphere.topOfAtmosphereAltitudeM()) {
            return;
        }
        maxExtinction = std::max(
            maxExtinction,
            computeOpticalProperties(atmosphere, wavelengthManager, altitude, wavelength_nm).extinction_m_inv
        );
    };

    updateAtAltitude(minAltitude);
    updateAtAltitude(maxAltitude);
    for (const AtmosphereLayer &layer : atmosphere.layers()) {
        if (layer.altitude_m >= minAltitude && layer.altitude_m <= maxAltitude) {
            updateAtAltitude(layer.altitude_m);
        }
    }
    return std::max(maxExtinction, 1.0e-12);
}

SampledInteraction sampleInteraction(
    const Atmosphere &atmosphere,
    const WavelengthManager &wavelengthManager,
    const SimulationConfig &config,
    const Vec3 &origin,
    const Vec3 &direction,
    double wavelength_nm,
    std::mt19937 &rng
)
{
    const double topRadius = topOfAtmosphereRadius(config);
    bool hitsGroundBoundary = false;
    const double boundaryDistance = distanceToBoundary(origin, direction, topRadius, hitsGroundBoundary);
    if (!std::isfinite(boundaryDistance)) {
        return {};
    }

    const double minAltitude = std::max(0.0, minimumAltitudeAlongSegment(origin, direction, boundaryDistance));
    const double maxAltitude = std::max(minAltitude, maximumAltitudeAlongSegment(origin, direction, boundaryDistance));
    const double majorant = extinctionMajorantAlongSegment(
        atmosphere,
        wavelengthManager,
        wavelength_nm,
        minAltitude,
        maxAltitude
    );
    if (majorant <= 0.0) {
        SampledInteraction interaction {};
        interaction.hit_ground = hitsGroundBoundary;
        return interaction;
    }

    std::uniform_real_distribution<double> uniform01(0.0, 1.0);
    std::exponential_distribution<double> freePath(majorant);

    double distance = 0.0;
    while (distance < boundaryDistance) {
        distance += freePath(rng);
        if (distance >= boundaryDistance) {
            SampledInteraction interaction {};
            interaction.hit_ground = hitsGroundBoundary;
            return interaction;
        }

        const Vec3 candidate = origin + direction * distance;
        const double altitude = altitudeM(candidate);
        if (altitude < 0.0) {
            SampledInteraction interaction {};
            interaction.hit_ground = true;
            return interaction;
        }

        const LocalOpticalProperties optics = computeOpticalProperties(atmosphere, wavelengthManager, altitude, wavelength_nm);
        if (optics.extinction_m_inv <= 0.0) {
            continue;
        }

        if (uniform01(rng) * majorant >= optics.extinction_m_inv) {
            continue;
        }

        const double totalScattering = optics.rayleigh_scattering_m_inv + optics.aerosol_scattering_m_inv;
        SampledInteraction interaction {};
        interaction.position = candidate;
        interaction.optics = optics;
        if (totalScattering > 0.0) {
            interaction.scattered = true;
            interaction.collision_weight = totalScattering / std::max(1.0e-12, optics.extinction_m_inv);
            const double xi = uniform01(rng) * totalScattering;
            interaction.is_rayleigh = xi < optics.rayleigh_scattering_m_inv;
            return interaction;
        }

        interaction.absorbed = true;
        return interaction;
    }

    SampledInteraction interaction {};
    interaction.hit_ground = hitsGroundBoundary;
    return interaction;
}

StokesVector addStokes(const StokesVector &lhs, const StokesVector &rhs)
{
    return {
        lhs.I + rhs.I,
        lhs.Q + rhs.Q,
        lhs.U + rhs.U,
        lhs.V + rhs.V,
    };
}

StokesVector totalMean(const DirectionEstimate &estimate)
{
    return addStokes(
        addStokes(estimate.first_order, estimate.second_order.total),
        estimate.higher_order
    );
}

DirectionTimingSummary publicTimingSummary(const DirectionEstimate::TimingSummary &timing)
{
    return {
        timing.first_order_seconds,
        timing.first_order_view_samples,
        timing.first_order_steps,
        timing.solar_disk_nodes,
        timing.second_order_seconds,
        timing.second_order_incoming_single_scatter_seconds,
        timing.second_order_incoming_single_scatter_calls,
        timing.second_order_nonzero_incoming_calls,
        timing.second_order_view_samples,
        timing.second_order_mu_phi_evaluations,
        timing.second_order_view_steps,
        timing.second_order_ray_steps,
        timing.second_order_mu_nodes,
        timing.second_order_phi_nodes,
        timing.higher_order_seconds,
        timing.spectral_band_count,
        timing.mc_sample_count,
    };
}

DirectionEstimate::TimingSummary internalTimingSummary(const DirectionTimingSummary &timing)
{
    DirectionEstimate::TimingSummary internal {};
    internal.first_order_seconds = timing.first_order_seconds;
    internal.first_order_view_samples = timing.first_order_view_samples;
    internal.first_order_steps = timing.first_order_steps;
    internal.solar_disk_nodes = timing.solar_disk_nodes;
    internal.second_order_seconds = timing.second_order_seconds;
    internal.second_order_incoming_single_scatter_seconds = timing.second_order_incoming_single_scatter_seconds;
    internal.second_order_incoming_single_scatter_calls = timing.second_order_incoming_single_scatter_calls;
    internal.second_order_nonzero_incoming_calls = timing.second_order_nonzero_incoming_calls;
    internal.second_order_view_samples = timing.second_order_view_samples;
    internal.second_order_mu_phi_evaluations = timing.second_order_mu_phi_evaluations;
    internal.second_order_view_steps = timing.second_order_view_steps;
    internal.second_order_ray_steps = timing.second_order_ray_steps;
    internal.second_order_mu_nodes = timing.second_order_mu_nodes;
    internal.second_order_phi_nodes = timing.second_order_phi_nodes;
    internal.higher_order_seconds = timing.higher_order_seconds;
    internal.spectral_band_count = timing.spectral_band_count;
    internal.mc_sample_count = timing.mc_sample_count;
    return internal;
}

StokesVector scaleStokes(const StokesVector &vector, double scalar)
{
    return {
        vector.I * scalar,
        vector.Q * scalar,
        vector.U * scalar,
        vector.V * scalar,
    };
}

StokesVector scalePolarization(const StokesVector &vector, double scale)
{
    const double safeScale = clamp(scale, 0.0, 1.0);
    return {
        vector.I,
        vector.Q * safeScale,
        vector.U * safeScale,
        vector.V * safeScale,
    };
}

StokesVector boostUnpolarizedIntensity(const StokesVector &vector, double boost)
{
    const double safeBoost = std::max(1.0, boost);
    return {
        vector.I * safeBoost,
        vector.Q,
        vector.U,
        vector.V,
    };
}

double smoothZenithBoostWeight(double zenithDeg, double cutoffDeg)
{
    const double safeCutoff = std::max(1.0, cutoffDeg);
    const double x = clamp(1.0 - zenithDeg / safeCutoff, 0.0, 1.0);
    return x * x * (3.0 - 2.0 * x);
}

bool twilightDepolarizationActive(const SimulationConfig &config, const SimulationContext &context)
{
    return config.monte_carlo.twilight_order_depolarization && context.sun.zenith_deg > 90.0;
}

int highestDeterministicScatteringOrder(const SimulationConfig &config)
{
    if (config.monte_carlo.deterministic_second_scatter) {
        return 2;
    }
    if (config.monte_carlo.deterministic_single_scatter) {
        return 1;
    }
    return 0;
}

int recursiveMinimumScatteringOrder(const SimulationConfig &config, int eventIndex)
{
    const int deterministicCutoff = highestDeterministicScatteringOrder(config);
    const int nominalRecursiveOrder = eventIndex + 2;
    return std::max(nominalRecursiveOrder, deterministicCutoff + 1);
}

bool empiricalTwilightTuningEnabled(const MonteCarloConfig &config)
{
    return config.twilight_order_depolarization ||
        std::abs(config.twilight_second_order_polarization_scale - 1.0) > 1.0e-12 ||
        std::abs(config.twilight_higher_order_polarization_scale - 1.0) > 1.0e-12 ||
        std::abs(config.twilight_second_order_intensity_boost - 1.0) > 1.0e-12 ||
        std::abs(config.twilight_higher_order_intensity_boost - 1.0) > 1.0e-12;
}

SecondScatterQuadrature secondScatterQuadrature(
    const SimulationConfig &config,
    const SimulationContext &context,
    double viewZenithDeg
)
{
    SecondScatterQuadrature quadrature {
        std::max(1, config.monte_carlo.second_scatter_view_steps),
        std::max(1, config.monte_carlo.second_scatter_ray_steps),
        std::max(1, config.monte_carlo.second_scatter_mu_nodes),
        std::max(1, config.monte_carlo.second_scatter_phi_nodes),
    };

    if (!config.monte_carlo.twilight_second_scatter_adaptive || context.sun.zenith_deg <= 90.0) {
        return quadrature;
    }

    const double thresholdDeg = std::max(0.0, config.monte_carlo.twilight_second_scatter_zenith_threshold_deg);
    const double horizonWeight = clamp((viewZenithDeg - thresholdDeg) / std::max(1.0, 90.0 - thresholdDeg), 0.0, 1.0);
    if (horizonWeight <= 0.0) {
        return quadrature;
    }

    quadrature.view_steps = std::max(
        quadrature.view_steps,
        static_cast<int>(std::lround(
            config.monte_carlo.twilight_second_scatter_min_view_steps * (0.5 + 0.5 * horizonWeight)
        ))
    );
    quadrature.ray_steps = std::max(
        quadrature.ray_steps,
        static_cast<int>(std::lround(
            config.monte_carlo.twilight_second_scatter_min_ray_steps * (0.5 + 0.5 * horizonWeight)
        ))
    );
    quadrature.mu_nodes = std::max(
        quadrature.mu_nodes,
        static_cast<int>(std::lround(
            config.monte_carlo.twilight_second_scatter_min_mu_nodes * (0.5 + 0.5 * horizonWeight)
        ))
    );
    quadrature.phi_nodes = std::max(
        quadrature.phi_nodes,
        static_cast<int>(std::lround(
            config.monte_carlo.twilight_second_scatter_min_phi_nodes * (0.5 + 0.5 * horizonWeight)
        ))
    );
    return quadrature;
}

std::array<double, 3> normalizeGuidingWeights(
    double phaseWeight,
    double tangentWeight,
    double horizonWeight
)
{
    phaseWeight = std::max(0.0, phaseWeight);
    tangentWeight = std::max(0.0, tangentWeight);
    horizonWeight = std::max(0.0, horizonWeight);
    const double total = phaseWeight + tangentWeight + horizonWeight;
    if (total <= 1.0e-12) {
        return {1.0, 0.0, 0.0};
    }
    return {
        phaseWeight / total,
        tangentWeight / total,
        horizonWeight / total,
    };
}

std::array<int, 3> stratifiedTwilightBranchCounts(
    const std::array<double, 3> &weights,
    int branchCount
)
{
    std::array<int, 3> counts {0, 0, 0};
    if (branchCount <= 0) {
        return counts;
    }

    std::vector<int> activeComponents;
    for (int component = 0; component < 3; ++component) {
        if (weights[component] > 1.0e-9) {
            activeComponents.push_back(component);
        }
    }
    if (activeComponents.empty()) {
        counts[0] = branchCount;
        return counts;
    }

    const int seedCount = std::min<int>(branchCount, activeComponents.size());
    for (int index = 0; index < seedCount; ++index) {
        counts[activeComponents[index]] = 1;
    }

    int remaining = branchCount - seedCount;
    if (remaining <= 0) {
        return counts;
    }

    std::array<int, 3> extraCounts {0, 0, 0};
    std::array<double, 3> remainders {0.0, 0.0, 0.0};
    int assigned = 0;
    for (int component = 0; component < 3; ++component) {
        if (weights[component] <= 1.0e-9) {
            continue;
        }
        const double expected = remaining * weights[component];
        extraCounts[component] = static_cast<int>(std::floor(expected));
        remainders[component] = expected - extraCounts[component];
        assigned += extraCounts[component];
    }
    for (int component = 0; component < 3; ++component) {
        counts[component] += extraCounts[component];
    }

    int leftover = remaining - assigned;
    while (leftover > 0) {
        int bestComponent = activeComponents.front();
        double bestRemainder = -1.0;
        for (int component : activeComponents) {
            if (remainders[component] > bestRemainder + 1.0e-12) {
                bestRemainder = remainders[component];
                bestComponent = component;
            }
        }
        ++counts[bestComponent];
        remainders[bestComponent] = -1.0;
        --leftover;
    }

    return counts;
}

double phasePdfForScatter(
    const SimulationContext &context,
    bool isRayleigh,
    double wavelength_nm,
    double cosTheta
)
{
    return std::max(
        1.0e-12,
        isRayleigh
            ? rayleighPhase(cosTheta)
            : context.phaseTable.phasePdf(wavelength_nm, cosTheta)
    );
}

double conePdf(const Vec3 &direction, const Vec3 &axis, double cosHalfAngle)
{
    const double cosAngle = dot(normalize(direction), normalize(axis));
    if (cosAngle < cosHalfAngle) {
        return 0.0;
    }
    return 1.0 / std::max(1.0e-12, 2.0 * PI * (1.0 - cosHalfAngle));
}

Vec3 sampleUniformConeDirection(const Vec3 &axis, double cosHalfAngle, std::mt19937 &rng)
{
    std::uniform_real_distribution<double> uniform01(0.0, 1.0);
    const double cosTheta = 1.0 - uniform01(rng) * (1.0 - cosHalfAngle);
    const double azimuth = 2.0 * PI * uniform01(rng);
    return directionAroundAxis(normalize(axis), cosTheta, azimuth);
}

Vec3 twilightLimbGuideDirection(
    const Vec3 &position,
    const Vec3 &sunDirection,
    double elevationDeg
)
{
    const Vec3 up = normalize(position);
    Vec3 horizontal = sunDirection - up * dot(sunDirection, up);
    if (norm(horizontal) < 1.0e-12) {
        horizontal = tangentFromAxis(up);
    }
    horizontal = normalize(horizontal);

    const double elevationRad = clamp(elevationDeg, 0.0, 89.0) * PI / 180.0;
    return normalize(horizontal * std::cos(elevationRad) + up * std::sin(elevationRad));
}

double polarizationBandPdf(const Vec3 &direction, const Vec3 &axis, double muHalfWidth)
{
    const double mu = dot(normalize(direction), normalize(axis));
    if (std::abs(mu) > muHalfWidth) {
        return 0.0;
    }
    return 1.0 / std::max(1.0e-12, 4.0 * PI * muHalfWidth);
}

Vec3 samplePolarizationBandDirection(const Vec3 &axis, double muHalfWidth, std::mt19937 &rng)
{
    std::uniform_real_distribution<double> uniform01(0.0, 1.0);
    const double cosTheta = -muHalfWidth + 2.0 * muHalfWidth * uniform01(rng);
    const double azimuth = 2.0 * PI * uniform01(rng);
    return directionAroundAxis(normalize(axis), cosTheta, azimuth);
}

GuidedProposalParameters guidedProposalParameters(
    const SimulationContext &context,
    const SimulationConfig &config,
    bool isRayleigh,
    const Vec3 &position,
    const Vec3 &outgoingDirection,
    int eventIndex
)
{
    const double defaultPhaseFraction = clamp(config.monte_carlo.source_guided_phase_fraction, 0.0, 1.0);
    const double defaultConeHalfAngleDeg = clamp(config.monte_carlo.source_guided_cone_half_angle_deg, 1.0, 89.0);
    double phaseFraction = defaultPhaseFraction;
    double coneHalfAngleDeg = defaultConeHalfAngleDeg;
    double polarizationFraction = 0.0;
    Vec3 guideAxis = context.sun.direction;
    if (isRayleigh) {
        const double rayleighPhaseFraction = clamp(config.monte_carlo.rayleigh_source_guided_phase_fraction, 0.0, 1.0);
        const double rayleighConeHalfAngleDeg = clamp(
            config.monte_carlo.rayleigh_source_guided_cone_half_angle_deg,
            1.0,
            89.0
        );
        const double muSunAbs = std::abs(clamp(dot(context.sun.direction, outgoingDirection), -1.0, 1.0));
        const double blend = clamp(1.0 - muSunAbs / 0.7, 0.0, 1.0);
        phaseFraction = defaultPhaseFraction + blend * (rayleighPhaseFraction - defaultPhaseFraction);
        coneHalfAngleDeg = defaultConeHalfAngleDeg + blend * (rayleighConeHalfAngleDeg - defaultConeHalfAngleDeg);
        polarizationFraction =
            clamp(config.monte_carlo.rayleigh_polarization_guided_fraction, 0.0, 0.95) *
            blend *
            rayleighPolarizationFraction(muSunAbs);
    }

    const bool useTwilightLimbGuide =
        config.monte_carlo.twilight_limb_guiding &&
        dot(context.sun.direction, normalize(position)) < 0.0;
    if (useTwilightLimbGuide) {
        phaseFraction = std::min(
            phaseFraction,
            clamp(config.monte_carlo.twilight_limb_phase_fraction, 0.0, 1.0)
        );
        coneHalfAngleDeg = clamp(config.monte_carlo.twilight_limb_cone_half_angle_deg, 1.0, 89.0);
        guideAxis = twilightLimbGuideDirection(
            position,
            context.sun.direction,
            config.monte_carlo.twilight_limb_elevation_deg
        );
    }

    const bool useTwilightHigherOrderGuide =
        config.monte_carlo.twilight_higher_order_guiding &&
        context.sun.zenith_deg > 90.0 &&
        recursiveMinimumScatteringOrder(config, eventIndex) >= 3;
    if (useTwilightHigherOrderGuide) {
        const std::array<double, 3> weights = normalizeGuidingWeights(
            config.monte_carlo.twilight_higher_order_phase_fraction,
            config.monte_carlo.twilight_higher_order_tangent_fraction,
            config.monte_carlo.twilight_higher_order_horizon_fraction
        );
        phaseFraction = weights[0];
        polarizationFraction = 0.0;
        guideAxis = twilightLimbGuideDirection(position, context.sun.direction, 0.0);
        return {
            phaseFraction,
            std::cos(defaultConeHalfAngleDeg * PI / 180.0),
            0.0,
            clamp(config.monte_carlo.rayleigh_polarization_guided_mu_half_width, 1.0e-3, 1.0),
            guideAxis,
            true,
            weights[1],
            weights[2],
            std::cos(clamp(config.monte_carlo.twilight_higher_order_tangent_cone_half_angle_deg, 1.0, 89.0) * PI / 180.0),
            std::cos(clamp(config.monte_carlo.twilight_higher_order_horizon_cone_half_angle_deg, 1.0, 89.0) * PI / 180.0),
            twilightLimbGuideDirection(position, context.sun.direction, 0.0),
            twilightLimbGuideDirection(
                position,
                context.sun.direction,
                config.monte_carlo.twilight_higher_order_horizon_elevation_deg
            ),
        };
    }

    return {
        phaseFraction,
        std::cos(coneHalfAngleDeg * PI / 180.0),
        polarizationFraction,
        clamp(config.monte_carlo.rayleigh_polarization_guided_mu_half_width, 1.0e-3, 1.0),
        guideAxis,
        false,
        0.0,
        0.0,
        1.0,
        1.0,
        guideAxis,
        guideAxis,
    };
}

ContinuationSample sampleContinuationDirection(
    const SimulationContext &context,
    const SimulationConfig &config,
    double wavelength_nm,
    bool isRayleigh,
    const Vec3 &position,
    const Vec3 &outgoingDirection,
    int eventIndex,
    bool useGuiding,
    std::mt19937 &rng,
    TwilightGuidedComponent forcedTwilightComponent = TwilightGuidedComponent::none
)
{
    const GuidedProposalParameters proposal = guidedProposalParameters(
        context,
        config,
        isRayleigh,
        position,
        outgoingDirection,
        eventIndex
    );

    if (!useGuiding) {
        const ScatteringSample phaseSample = isRayleigh
            ? sampleRayleighDirection(rng)
            : context.phaseTable.sampleDirection(wavelength_nm, rng);
        return {
            directionAroundAxis(outgoingDirection, phaseSample.cos_theta, phaseSample.azimuth_rad),
            phaseSample.cos_theta,
            std::max(1.0e-12, phaseSample.pdf),
            1.0,
        };
    }

    std::uniform_real_distribution<double> uniform01(0.0, 1.0);
    Vec3 incomingDirection {0.0, 0.0, 1.0};
    if (proposal.use_twilight_higher_order) {
        TwilightGuidedComponent selectedComponent = forcedTwilightComponent;
        if (selectedComponent == TwilightGuidedComponent::none) {
            const double selector = uniform01(rng);
            if (selector < proposal.phase_fraction) {
                selectedComponent = TwilightGuidedComponent::phase;
            } else if (selector < proposal.phase_fraction + proposal.tangent_fraction) {
                selectedComponent = TwilightGuidedComponent::tangent;
            } else {
                selectedComponent = TwilightGuidedComponent::horizon;
            }
        }

        if (selectedComponent == TwilightGuidedComponent::phase) {
            const ScatteringSample phaseSample = isRayleigh
                ? sampleRayleighDirection(rng)
                : context.phaseTable.sampleDirection(wavelength_nm, rng);
            incomingDirection = directionAroundAxis(outgoingDirection, phaseSample.cos_theta, phaseSample.azimuth_rad);
        } else if (selectedComponent == TwilightGuidedComponent::tangent) {
            incomingDirection = sampleUniformConeDirection(
                proposal.tangent_axis,
                proposal.tangent_cos_half_angle,
                rng
            );
        } else {
            incomingDirection = sampleUniformConeDirection(
                proposal.horizon_axis,
                proposal.horizon_cos_half_angle,
                rng
            );
        }
    } else if (proposal.phase_fraction >= 1.0 || uniform01(rng) < proposal.phase_fraction) {
        const ScatteringSample phaseSample = isRayleigh
            ? sampleRayleighDirection(rng)
            : context.phaseTable.sampleDirection(wavelength_nm, rng);
        incomingDirection = directionAroundAxis(outgoingDirection, phaseSample.cos_theta, phaseSample.azimuth_rad);
    } else {
        incomingDirection = sampleUniformConeDirection(proposal.guide_axis, proposal.cos_half_angle, rng);
    }

    const double cosTheta = clamp(dot(incomingDirection, outgoingDirection), -1.0, 1.0);
    const double phasePdf = phasePdfForScatter(context, isRayleigh, wavelength_nm, cosTheta);
    const double guidePdf = conePdf(incomingDirection, proposal.guide_axis, proposal.cos_half_angle);
    const double tangentPdf = conePdf(incomingDirection, proposal.tangent_axis, proposal.tangent_cos_half_angle);
    const double horizonPdf = conePdf(incomingDirection, proposal.horizon_axis, proposal.horizon_cos_half_angle);
    const double basePdf = proposal.use_twilight_higher_order
        ? proposal.phase_fraction * phasePdf +
            proposal.tangent_fraction * tangentPdf +
            proposal.horizon_fraction * horizonPdf
        : proposal.phase_fraction * phasePdf + (1.0 - proposal.phase_fraction) * guidePdf;
    const double polarizationPdf = polarizationBandPdf(
        incomingDirection,
        outgoingDirection,
        proposal.polarization_mu_half_width
    );
    const double mixturePdf =
        (1.0 - proposal.polarization_fraction) * basePdf +
        proposal.polarization_fraction * polarizationPdf;
    return {
        incomingDirection,
        cosTheta,
        std::max(1.0e-12, mixturePdf),
        1.0,
    };
}

ContinuationSample samplePolarizationGuidedDirection(
    const SimulationContext &context,
    const SimulationConfig &config,
    double wavelength_nm,
    const Vec3 &position,
    const Vec3 &outgoingDirection,
    int eventIndex,
    std::mt19937 &rng
)
{
    const GuidedProposalParameters proposal = guidedProposalParameters(
        context,
        config,
        true,
        position,
        outgoingDirection,
        eventIndex
    );
    const Vec3 incomingDirection = samplePolarizationBandDirection(
        outgoingDirection,
        proposal.polarization_mu_half_width,
        rng
    );
    const double cosTheta = clamp(dot(incomingDirection, outgoingDirection), -1.0, 1.0);
    const double phasePdf = phasePdfForScatter(context, true, wavelength_nm, cosTheta);
    const double guidePdf = conePdf(incomingDirection, proposal.guide_axis, proposal.cos_half_angle);
    const double basePdf = proposal.phase_fraction * phasePdf + (1.0 - proposal.phase_fraction) * guidePdf;
    const double polarizationPdf = polarizationBandPdf(
        incomingDirection,
        outgoingDirection,
        proposal.polarization_mu_half_width
    );
    const double mixturePdf =
        (1.0 - proposal.polarization_fraction) * basePdf +
        proposal.polarization_fraction * polarizationPdf;
    return {
        incomingDirection,
        cosTheta,
        std::max(1.0e-12, mixturePdf),
        1.0,
    };
}

unsigned int seedForDirection(
    unsigned int baseSeed,
    std::size_t directionIndex,
    double zenithDeg,
    double azimuthDeg
)
{
    unsigned long long hash = HASH_OFFSET_BASIS;
    hash = hashCombine(hash, static_cast<unsigned long long>(baseSeed));
    hash = hashCombine(hash, static_cast<unsigned long long>(directionIndex + 1));
    hash = hashCombine(hash, static_cast<unsigned long long>(std::llround(zenithDeg * 1000.0)));
    hash = hashCombine(hash, static_cast<unsigned long long>(std::llround(normalizeAzimuthDeg(azimuthDeg) * 1000.0)));
    return static_cast<unsigned int>((hash >> 32) ^ (hash & 0xffffffffull));
}

unsigned int seedForDirectionSample(
    unsigned int baseSeed,
    std::size_t directionIndex,
    double zenithDeg,
    double azimuthDeg,
    int sampleIndex
)
{
    unsigned long long hash = HASH_OFFSET_BASIS;
    hash = hashCombine(hash, static_cast<unsigned long long>(baseSeed));
    hash = hashCombine(hash, static_cast<unsigned long long>(directionIndex + 1));
    hash = hashCombine(hash, static_cast<unsigned long long>(std::llround(zenithDeg * 1000.0)));
    hash = hashCombine(hash, static_cast<unsigned long long>(std::llround(normalizeAzimuthDeg(azimuthDeg) * 1000.0)));
    hash = hashCombine(hash, static_cast<unsigned long long>(sampleIndex + 1));
    return static_cast<unsigned int>((hash >> 32) ^ (hash & 0xffffffffull));
}

std::vector<std::pair<double, double>> midpointMuQuadrature(int nodeCount)
{
    std::vector<std::pair<double, double>> nodes;
    const int safeNodeCount = std::max(1, nodeCount);
    for (int i = 0; i < safeNodeCount; ++i) {
        const double mu = -1.0 + (2.0 * (i + 0.5) / safeNodeCount);
        const double weight = 2.0 / safeNodeCount;
        nodes.push_back({mu, weight});
    }
    return nodes;
}

SingleScatterBreakdown deterministicSingleScatterBandAlongRayBreakdown(
    const SimulationContext &context,
    const SimulationConfig &config,
    const SpectralBand &band,
    const Vec3 &rayStart,
    const Vec3 &outgoingDirection,
    int stepCount
)
{
    const Vec3 segmentOrigin = rayStart + outgoingDirection * 1.0e-3;
    const double topRadius = topOfAtmosphereRadius(config);

    bool hitsGround = false;
    const double boundaryDistance = distanceToBoundary(segmentOrigin, outgoingDirection, topRadius, hitsGround);
    if (!std::isfinite(boundaryDistance) || hitsGround) {
        return {};
    }

    const std::vector<SolarDiskSample> &sunSamples = context.sun_samples;
    const int safeSteps = std::max(1, stepCount);
    const double ds = boundaryDistance / static_cast<double>(safeSteps);

    SingleScatterBreakdown contribution {};
    double tauView = 0.0;
    for (int stepIndex = 0; stepIndex < safeSteps; ++stepIndex) {
        const double centerDistance = (stepIndex + 0.5) * ds;
        const Vec3 samplePosition = segmentOrigin + outgoingDirection * centerDistance;
        const double altitude = altitudeM(samplePosition);
        if (altitude < 0.0) {
            break;
        }

        const LocalOpticalProperties optics = computeOpticalProperties(
            context.atmosphere,
            context.wavelengthManager,
            altitude,
            band.wavelength_nm
        );
        if (optics.extinction_m_inv <= 0.0) {
            continue;
        }

        const double tauViewMid = tauView + 0.5 * optics.extinction_m_inv * ds;
            for (const SolarDiskSample &sunSample : sunSamples) {
                const double tauSun = opticalDepthToSun(
                    context.atmosphere,
                    context.wavelengthManager,
                    config,
                    samplePosition + sunSample.direction * 1.0,
                    sunSample.direction,
                    band.wavelength_nm
                );
                if (!std::isfinite(tauSun)) {
                    continue;
                }

                const double muSun = clamp(dot(sunSample.direction, outgoingDirection), -1.0, 1.0);
                const double sourceScale =
                    band.solar_irradiance_w_m2_nm * sunSample.weight * std::exp(-(tauViewMid + tauSun)) * ds;
                const PhaseMatrixCoefficients aerosolCoefficients =
                    context.phaseTable.coefficients(band.wavelength_nm, muSun);
                const StokesVector directSun = initUnpolarized(sourceScale);

                if (optics.rayleigh_scattering_m_inv > 0.0) {
                    const MuellerMatrix rayleighEvent = eventMuellerMatrix(
                        samplePosition,
                        samplePosition,
                        rayStart,
                        sunSample.direction,
                        outgoingDirection,
                        true,
                        aerosolCoefficients
                    );
                    contribution.rayleigh = addStokes(
                        contribution.rayleigh,
                        scaleStokes(
                            multiply(rayleighEvent, directSun),
                            optics.rayleigh_scattering_m_inv
                        )
                    );
                }

                if (optics.aerosol_scattering_m_inv > 0.0) {
                    const MuellerMatrix aerosolEvent = eventMuellerMatrix(
                        samplePosition,
                        samplePosition,
                        rayStart,
                        sunSample.direction,
                        outgoingDirection,
                        false,
                        aerosolCoefficients
                    );
                    contribution.aerosol = addStokes(
                        contribution.aerosol,
                        scaleStokes(
                            multiply(aerosolEvent, directSun),
                            optics.aerosol_scattering_m_inv
                        )
                    );
                }
            }

        tauView += optics.extinction_m_inv * ds;
    }

    return contribution;
}

std::vector<SingleScatterBreakdown> deterministicSingleScatterAlongRayBreakdownAllBands(
    const SimulationContext &context,
    const SimulationConfig &config,
    const Vec3 &rayStart,
    const Vec3 &outgoingDirection,
    int stepCount
)
{
    const std::vector<SpectralBand> &bands = context.active_bands;
    std::vector<SingleScatterBreakdown> contributions(bands.size());

    const Vec3 segmentOrigin = rayStart + outgoingDirection * 1.0e-3;
    const double topRadius = topOfAtmosphereRadius(config);

    bool hitsGround = false;
    const double boundaryDistance = distanceToBoundary(segmentOrigin, outgoingDirection, topRadius, hitsGround);
    if (!std::isfinite(boundaryDistance) || hitsGround) {
        return contributions;
    }

    const std::vector<SolarDiskSample> &sunSamples = context.sun_samples;
    const std::size_t bandCount = bands.size();
    const std::size_t sunSampleCount = sunSamples.size();
    const int safeSteps = std::max(1, stepCount);
    const double ds = boundaryDistance / static_cast<double>(safeSteps);

    std::vector<FirstOrderViewSample> viewSamples(static_cast<std::size_t>(safeSteps));
    std::vector<std::vector<LocalOpticalProperties>> viewOptics(
        static_cast<std::size_t>(safeSteps),
        std::vector<LocalOpticalProperties>(bandCount)
    );
    std::vector<double> tauSunBySampleBand(
        static_cast<std::size_t>(safeSteps) * sunSampleCount * bandCount,
        std::numeric_limits<double>::infinity()
    );

    for (int stepIndex = 0; stepIndex < safeSteps; ++stepIndex) {
        const double centerDistance = (stepIndex + 0.5) * ds;
        FirstOrderViewSample &viewSample = viewSamples[static_cast<std::size_t>(stepIndex)];
        viewSample.position = segmentOrigin + outgoingDirection * centerDistance;
        const double altitude = altitudeM(viewSample.position);
        if (altitude < 0.0) {
            continue;
        }

        viewSample.optical_state = spectralOpticalStateAtAltitude(context.atmosphere, altitude);
        for (std::size_t bandIndex = 0; bandIndex < bandCount; ++bandIndex) {
            viewOptics[static_cast<std::size_t>(stepIndex)][bandIndex] =
                opticalPropertiesFromState(viewSample.optical_state, bands[bandIndex]);
        }

        for (std::size_t sunIndex = 0; sunIndex < sunSampleCount; ++sunIndex) {
            const SolarDiskSample &sunSample = sunSamples[sunIndex];
            const RayAltitudeQuadrature quadrature = buildOpticalDepthQuadrature(
                context.atmosphere,
                config,
                viewSample.position + sunSample.direction * 1.0,
                sunSample.direction
            );
            const std::vector<double> opticalDepths = opticalDepthSpectrumFromQuadrature(
                context.atmosphere,
                bands,
                quadrature
            );
            for (std::size_t bandIndex = 0; bandIndex < bandCount; ++bandIndex) {
                const std::size_t flatIndex =
                    (static_cast<std::size_t>(stepIndex) * sunSampleCount + sunIndex) * bandCount + bandIndex;
                tauSunBySampleBand[flatIndex] = opticalDepths[bandIndex];
            }
        }
    }

    std::vector<double> tauViewByBand(bandCount, 0.0);
    for (int stepIndex = 0; stepIndex < safeSteps; ++stepIndex) {
        const FirstOrderViewSample &viewSample = viewSamples[static_cast<std::size_t>(stepIndex)];
        if (altitudeM(viewSample.position) < 0.0) {
            break;
        }

        for (std::size_t sunIndex = 0; sunIndex < sunSampleCount; ++sunIndex) {
            const SolarDiskSample &sunSample = sunSamples[sunIndex];
            const double muSun = clamp(dot(sunSample.direction, outgoingDirection), -1.0, 1.0);
            for (std::size_t bandIndex = 0; bandIndex < bandCount; ++bandIndex) {
                const SpectralBand &band = bands[bandIndex];
                const LocalOpticalProperties &optics =
                    viewOptics[static_cast<std::size_t>(stepIndex)][bandIndex];
                if (optics.extinction_m_inv <= 0.0) {
                    continue;
                }

                const double tauSun = tauSunBySampleBand[
                    (static_cast<std::size_t>(stepIndex) * sunSampleCount + sunIndex) * bandCount + bandIndex
                ];
                if (!std::isfinite(tauSun)) {
                    continue;
                }

                const double tauViewMid = tauViewByBand[bandIndex] + 0.5 * optics.extinction_m_inv * ds;
                const double sourceScale =
                    band.solar_irradiance_w_m2_nm * sunSample.weight * std::exp(-(tauViewMid + tauSun)) * ds;
                const PhaseMatrixCoefficients aerosolCoefficients =
                    context.phaseTable.coefficients(band.wavelength_nm, muSun);
                const StokesVector directSun = initUnpolarized(sourceScale);

                if (optics.rayleigh_scattering_m_inv > 0.0) {
                    const MuellerMatrix rayleighEvent = eventMuellerMatrix(
                        viewSample.position,
                        viewSample.position,
                        rayStart,
                        sunSample.direction,
                        outgoingDirection,
                        true,
                        aerosolCoefficients
                    );
                    contributions[bandIndex].rayleigh = addStokes(
                        contributions[bandIndex].rayleigh,
                        scaleStokes(
                            multiply(rayleighEvent, directSun),
                            optics.rayleigh_scattering_m_inv
                        )
                    );
                }

                if (optics.aerosol_scattering_m_inv > 0.0) {
                    const MuellerMatrix aerosolEvent = eventMuellerMatrix(
                        viewSample.position,
                        viewSample.position,
                        rayStart,
                        sunSample.direction,
                        outgoingDirection,
                        false,
                        aerosolCoefficients
                    );
                    contributions[bandIndex].aerosol = addStokes(
                        contributions[bandIndex].aerosol,
                        scaleStokes(
                            multiply(aerosolEvent, directSun),
                            optics.aerosol_scattering_m_inv
                        )
                    );
                }
            }
        }

        for (std::size_t bandIndex = 0; bandIndex < bandCount; ++bandIndex) {
            tauViewByBand[bandIndex] +=
                viewOptics[static_cast<std::size_t>(stepIndex)][bandIndex].extinction_m_inv * ds;
        }
    }

    return contributions;
}

StokesVector deterministicSingleScatterBandAlongRay(
    const SimulationContext &context,
    const SimulationConfig &config,
    const SpectralBand &band,
    const Vec3 &rayStart,
    const Vec3 &outgoingDirection,
    int stepCount
)
{
    const SingleScatterBreakdown breakdown = deterministicSingleScatterBandAlongRayBreakdown(
        context,
        config,
        band,
        rayStart,
        outgoingDirection,
        stepCount
    );
    return addStokes(breakdown.rayleigh, breakdown.aerosol);
}

SecondOrderBreakdown deterministicSecondScatterContribution(
    const SimulationContext &context,
    const SimulationConfig &config,
    double viewZenithDeg,
    double viewAzimuthDeg,
    DirectionEstimate::TimingSummary *timingSummary = nullptr,
    const DirectionStageCallback &stageCallback = {}
)
{
    const auto secondOrderStart = std::chrono::steady_clock::now();
    const Vec3 observer = observerPosition(config.observer);
    const LocalBasis observerBasis = localBasisAtPosition(observer);
    const Vec3 outgoingDirection = directionFromZenithAzimuth(viewZenithDeg, viewAzimuthDeg, observerBasis);
    const Vec3 segmentOrigin = observer + outgoingDirection * 1.0e-3;
    const double topRadius = topOfAtmosphereRadius(config);

    bool hitsGround = false;
    const double boundaryDistance = distanceToBoundary(segmentOrigin, outgoingDirection, topRadius, hitsGround);
    if (!std::isfinite(boundaryDistance) || hitsGround) {
        return {};
    }

    const SecondScatterQuadrature quadrature = secondScatterQuadrature(
        config,
        context,
        viewZenithDeg
    );
    if (timingSummary != nullptr) {
        timingSummary->second_order_view_steps = quadrature.view_steps;
        timingSummary->second_order_ray_steps = quadrature.ray_steps;
        timingSummary->second_order_mu_nodes = quadrature.mu_nodes;
        timingSummary->second_order_phi_nodes = quadrature.phi_nodes;
        timingSummary->spectral_band_count = context.active_bands.size();
    }
    const int viewSteps = quadrature.view_steps;
    const int raySteps = quadrature.ray_steps;
    const int maxPhiNodes = quadrature.phi_nodes;
    const std::vector<std::pair<double, double>> muQuadrature =
        midpointMuQuadrature(quadrature.mu_nodes);
    const std::size_t angularSamplesPerView =
        muQuadrature.size() * static_cast<std::size_t>(maxPhiNodes);
    const double ds1 = boundaryDistance / static_cast<double>(viewSteps);
    const double phiWeight = 2.0 * PI / static_cast<double>(maxPhiNodes);

    SecondOrderBreakdown totalContribution {};
    const std::vector<SpectralBand> &bands = context.active_bands;
    std::vector<std::vector<SingleScatterBreakdown>> incomingCache(
        static_cast<std::size_t>(viewSteps) * angularSamplesPerView
    );
    std::vector<bool> incomingCacheReady(
        static_cast<std::size_t>(viewSteps) * angularSamplesPerView,
        false
    );
    for (std::size_t bandIndex = 0; bandIndex < bands.size(); ++bandIndex) {
        const SpectralBand &band = bands[bandIndex];
        double tauView = 0.0;
        for (int stepIndex = 0; stepIndex < viewSteps; ++stepIndex) {
            const double centerDistance = (stepIndex + 0.5) * ds1;
            const Vec3 firstScatterPosition = segmentOrigin + outgoingDirection * centerDistance;
            const double altitude = altitudeM(firstScatterPosition);
            if (altitude < 0.0) {
                break;
            }

            const LocalOpticalProperties firstOptics = computeOpticalProperties(
                context.atmosphere,
                context.wavelengthManager,
                altitude,
                band.wavelength_nm
            );
            if (firstOptics.extinction_m_inv <= 0.0) {
                continue;
            }
            if (timingSummary != nullptr) {
                timingSummary->second_order_view_samples += 1;
            }

            const double tauViewMid = tauView + 0.5 * firstOptics.extinction_m_inv * ds1;
            const double attenuatedViewWeight = std::exp(-tauViewMid) * ds1;
            if (attenuatedViewWeight <= 0.0) {
                tauView += firstOptics.extinction_m_inv * ds1;
                continue;
            }

            for (std::size_t muIndex = 0; muIndex < muQuadrature.size(); ++muIndex) {
                const double muSecond = muQuadrature[muIndex].first;
                const double muWeight = muQuadrature[muIndex].second;
                for (int phiIndex = 0; phiIndex < maxPhiNodes; ++phiIndex) {
                    if (timingSummary != nullptr) {
                        timingSummary->second_order_mu_phi_evaluations += 1;
                    }

                    const double phi = (static_cast<double>(phiIndex) + 0.5) * phiWeight;
                    const Vec3 incomingDirection = directionAroundAxis(outgoingDirection, muSecond, phi);
                    const std::size_t angularIndex =
                        muIndex * static_cast<std::size_t>(maxPhiNodes) +
                        static_cast<std::size_t>(phiIndex);
                    const std::size_t cacheIndex =
                        static_cast<std::size_t>(stepIndex) * angularSamplesPerView + angularIndex;

                    if (!incomingCacheReady[cacheIndex]) {
                        const auto incomingStart = std::chrono::steady_clock::now();
                        incomingCache[cacheIndex] = deterministicSingleScatterAlongRayBreakdownAllBands(
                            context,
                            config,
                            firstScatterPosition,
                            incomingDirection,
                            raySteps
                        );
                        if (timingSummary != nullptr) {
                            timingSummary->second_order_incoming_single_scatter_seconds +=
                                std::chrono::duration<double>(std::chrono::steady_clock::now() - incomingStart).count();
                            timingSummary->second_order_incoming_single_scatter_calls += 1;
                        }
                        incomingCacheReady[cacheIndex] = true;
                    }

                    const SingleScatterBreakdown &incomingSingleScatter =
                        incomingCache[cacheIndex][bandIndex];
                    const StokesVector incomingTotal = addStokes(
                        incomingSingleScatter.rayleigh,
                        incomingSingleScatter.aerosol
                    );
                    const double magnitude =
                        std::abs(incomingTotal.I) +
                        std::abs(incomingTotal.Q) +
                        std::abs(incomingTotal.U) +
                        std::abs(incomingTotal.V);
                    if (magnitude <= 1.0e-18) {
                        continue;
                    }
                    if (timingSummary != nullptr) {
                        timingSummary->second_order_nonzero_incoming_calls += 1;
                    }

                    const double angularWeight = muWeight * phiWeight;
                    const double muEvent = clamp(dot(incomingDirection, outgoingDirection), -1.0, 1.0);
                    const PhaseMatrixCoefficients aerosolCoefficients =
                        context.phaseTable.coefficients(band.wavelength_nm, muEvent);

                    if (firstOptics.rayleigh_scattering_m_inv > 0.0) {
                        const MuellerMatrix rayleighEvent = eventMuellerMatrix(
                            firstScatterPosition,
                            firstScatterPosition,
                            observer,
                            incomingDirection,
                            outgoingDirection,
                            true,
                            aerosolCoefficients
                        );
                        if (
                            std::abs(incomingSingleScatter.rayleigh.I) +
                            std::abs(incomingSingleScatter.rayleigh.Q) +
                            std::abs(incomingSingleScatter.rayleigh.U) +
                            std::abs(incomingSingleScatter.rayleigh.V) > 0.0
                        ) {
                            totalContribution.rr = addStokes(
                                totalContribution.rr,
                                scaleStokes(
                                    multiply(rayleighEvent, incomingSingleScatter.rayleigh),
                                    firstOptics.rayleigh_scattering_m_inv * attenuatedViewWeight * angularWeight
                                )
                            );
                        }
                        if (
                            std::abs(incomingSingleScatter.aerosol.I) +
                            std::abs(incomingSingleScatter.aerosol.Q) +
                            std::abs(incomingSingleScatter.aerosol.U) +
                            std::abs(incomingSingleScatter.aerosol.V) > 0.0
                        ) {
                            totalContribution.ar = addStokes(
                                totalContribution.ar,
                                scaleStokes(
                                    multiply(rayleighEvent, incomingSingleScatter.aerosol),
                                    firstOptics.rayleigh_scattering_m_inv * attenuatedViewWeight * angularWeight
                                )
                            );
                        }
                    }

                    if (firstOptics.aerosol_scattering_m_inv > 0.0) {
                        const MuellerMatrix aerosolEvent = eventMuellerMatrix(
                            firstScatterPosition,
                            firstScatterPosition,
                            observer,
                            incomingDirection,
                            outgoingDirection,
                            false,
                            aerosolCoefficients
                        );
                        if (
                            std::abs(incomingSingleScatter.rayleigh.I) +
                            std::abs(incomingSingleScatter.rayleigh.Q) +
                            std::abs(incomingSingleScatter.rayleigh.U) +
                            std::abs(incomingSingleScatter.rayleigh.V) > 0.0
                        ) {
                            totalContribution.ra = addStokes(
                                totalContribution.ra,
                                scaleStokes(
                                    multiply(aerosolEvent, incomingSingleScatter.rayleigh),
                                    firstOptics.aerosol_scattering_m_inv * attenuatedViewWeight * angularWeight
                                )
                            );
                        }
                        if (
                            std::abs(incomingSingleScatter.aerosol.I) +
                            std::abs(incomingSingleScatter.aerosol.Q) +
                            std::abs(incomingSingleScatter.aerosol.U) +
                            std::abs(incomingSingleScatter.aerosol.V) > 0.0
                        ) {
                            totalContribution.aa = addStokes(
                                totalContribution.aa,
                                scaleStokes(
                                    multiply(aerosolEvent, incomingSingleScatter.aerosol),
                                    firstOptics.aerosol_scattering_m_inv * attenuatedViewWeight * angularWeight
                                )
                            );
                        }
                    }
                }
            }

            tauView += firstOptics.extinction_m_inv * ds1;
        }

        if (timingSummary != nullptr) {
            timingSummary->second_order_seconds =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - secondOrderStart).count();
        }
        if (stageCallback) {
            stageCallback(
                SampleProgress::Stage::second_order_band_complete,
                timingSummary != nullptr ? *timingSummary : DirectionEstimate::TimingSummary {},
                bandIndex + 1,
                bands.size()
            );
        }
    }

    if (twilightDepolarizationActive(config, context)) {
        const double secondOrderBoostWeight = smoothZenithBoostWeight(
            viewZenithDeg,
            config.monte_carlo.twilight_second_order_boost_zenith_cutoff_deg
        );
        const double secondOrderBoost =
            1.0 +
            (std::max(1.0, config.monte_carlo.twilight_second_order_intensity_boost) - 1.0) *
                secondOrderBoostWeight;
        totalContribution.rr = boostUnpolarizedIntensity(totalContribution.rr, secondOrderBoost);
        totalContribution.ar = boostUnpolarizedIntensity(totalContribution.ar, secondOrderBoost);
        totalContribution.ra = boostUnpolarizedIntensity(totalContribution.ra, secondOrderBoost);
        totalContribution.aa = boostUnpolarizedIntensity(totalContribution.aa, secondOrderBoost);
        totalContribution.rr = scalePolarization(
            totalContribution.rr,
            config.monte_carlo.twilight_second_order_polarization_scale
        );
        totalContribution.ar = scalePolarization(
            totalContribution.ar,
            config.monte_carlo.twilight_second_order_polarization_scale
        );
        totalContribution.ra = scalePolarization(
            totalContribution.ra,
            config.monte_carlo.twilight_second_order_polarization_scale
        );
        totalContribution.aa = scalePolarization(
            totalContribution.aa,
            config.monte_carlo.twilight_second_order_polarization_scale
        );
    }

    totalContribution.total = addStokes(
        addStokes(totalContribution.rr, totalContribution.ar),
        addStokes(totalContribution.ra, totalContribution.aa)
    );

    if (timingSummary != nullptr) {
        timingSummary->second_order_seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - secondOrderStart).count();
    }

    return totalContribution;
}

StokesVector deterministicSingleScatterContribution(
    const SimulationContext &context,
    const SimulationConfig &config,
    double viewZenithDeg,
    double viewAzimuthDeg,
    DirectionEstimate::TimingSummary *timingSummary = nullptr,
    const DirectionStageCallback &stageCallback = {}
)
{
    const auto firstOrderStart = std::chrono::steady_clock::now();
    const Vec3 observer = observerPosition(config.observer);
    const LocalBasis observerBasis = localBasisAtPosition(observer);
    const Vec3 outgoingDirection = directionFromZenithAzimuth(viewZenithDeg, viewAzimuthDeg, observerBasis);
    const Vec3 segmentOrigin = observer + outgoingDirection * 1.0e-3;
    const double topRadius = topOfAtmosphereRadius(config);

    bool hitsGround = false;
    const double boundaryDistance = distanceToBoundary(segmentOrigin, outgoingDirection, topRadius, hitsGround);
    if (!std::isfinite(boundaryDistance) || hitsGround) {
        return {0.0, 0.0, 0.0, 0.0};
    }

    const std::vector<SolarDiskSample> &sunSamples = context.sun_samples;
    const double minimumAltitudeM = std::max(
        0.0,
        minimumAltitudeAlongSegment(segmentOrigin, outgoingDirection, boundaryDistance)
    );
    const double stepLengthM = adaptiveLineIntegralStepM(config, minimumAltitudeM);
    const int baseSteps = std::max(
        std::max(1, config.monte_carlo.min_line_integral_steps),
        static_cast<int>(std::ceil(boundaryDistance / stepLengthM))
    );
    const double ds = boundaryDistance / static_cast<double>(baseSteps);
    if (timingSummary != nullptr) {
        timingSummary->first_order_steps = baseSteps;
        timingSummary->solar_disk_nodes = static_cast<int>(sunSamples.size());
        timingSummary->spectral_band_count = context.active_bands.size();
    }
    StokesVector totalContribution {0.0, 0.0, 0.0, 0.0};

    const std::vector<SpectralBand> &bands = context.active_bands;
    const std::size_t bandCount = bands.size();
    const std::size_t sunSampleCount = sunSamples.size();
    std::vector<FirstOrderViewSample> viewSamples(static_cast<std::size_t>(baseSteps));
    std::vector<std::vector<LocalOpticalProperties>> viewOptics(
        static_cast<std::size_t>(baseSteps),
        std::vector<LocalOpticalProperties>(bandCount)
    );
    std::vector<double> tauSunBySampleBand(
        static_cast<std::size_t>(baseSteps) * sunSampleCount * bandCount,
        std::numeric_limits<double>::infinity()
    );

    for (int stepIndex = 0; stepIndex < baseSteps; ++stepIndex) {
        const double centerDistance = (stepIndex + 0.5) * ds;
        FirstOrderViewSample &viewSample = viewSamples[static_cast<std::size_t>(stepIndex)];
        viewSample.position = segmentOrigin + outgoingDirection * centerDistance;
        const double altitude = altitudeM(viewSample.position);
        if (altitude < 0.0) {
            continue;
        }

        viewSample.optical_state = spectralOpticalStateAtAltitude(context.atmosphere, altitude);
        for (std::size_t bandIndex = 0; bandIndex < bandCount; ++bandIndex) {
            viewOptics[static_cast<std::size_t>(stepIndex)][bandIndex] =
                opticalPropertiesFromState(viewSample.optical_state, bands[bandIndex]);
        }

        for (std::size_t sunIndex = 0; sunIndex < sunSampleCount; ++sunIndex) {
            const SolarDiskSample &sunSample = sunSamples[sunIndex];
            const RayAltitudeQuadrature quadrature = buildOpticalDepthQuadrature(
                context.atmosphere,
                config,
                viewSample.position + sunSample.direction * 1.0,
                sunSample.direction
            );
            const std::vector<double> opticalDepths = opticalDepthSpectrumFromQuadrature(
                context.atmosphere,
                bands,
                quadrature
            );
            for (std::size_t bandIndex = 0; bandIndex < bandCount; ++bandIndex) {
                const std::size_t flatIndex =
                    (static_cast<std::size_t>(stepIndex) * sunSampleCount + sunIndex) * bandCount + bandIndex;
                tauSunBySampleBand[flatIndex] = opticalDepths[bandIndex];
            }
        }
    }

    for (std::size_t bandIndex = 0; bandIndex < bands.size(); ++bandIndex) {
        const SpectralBand &band = bands[bandIndex];
        double tauView = 0.0;
        for (int stepIndex = 0; stepIndex < baseSteps; ++stepIndex) {
            const FirstOrderViewSample &viewSample = viewSamples[static_cast<std::size_t>(stepIndex)];
            const LocalOpticalProperties &optics =
                viewOptics[static_cast<std::size_t>(stepIndex)][bandIndex];
            if (optics.extinction_m_inv <= 0.0) {
                continue;
            }
            if (timingSummary != nullptr) {
                timingSummary->first_order_view_samples += 1;
            }

            const double tauViewMid = tauView + 0.5 * optics.extinction_m_inv * ds;
            for (std::size_t sunIndex = 0; sunIndex < sunSampleCount; ++sunIndex) {
                const SolarDiskSample &sunSample = sunSamples[sunIndex];
                const double tauSun = tauSunBySampleBand[
                    (static_cast<std::size_t>(stepIndex) * sunSampleCount + sunIndex) * bandCount + bandIndex
                ];
                if (!std::isfinite(tauSun)) {
                    continue;
                }

                const double muSun = clamp(dot(sunSample.direction, outgoingDirection), -1.0, 1.0);
                const double sourceScale =
                    band.solar_irradiance_w_m2_nm * sunSample.weight * std::exp(-(tauViewMid + tauSun)) * ds;
                const PhaseMatrixCoefficients aerosolCoefficients =
                    context.phaseTable.coefficients(band.wavelength_nm, muSun);
                const StokesVector directSun = initUnpolarized(sourceScale);

                if (optics.rayleigh_scattering_m_inv > 0.0) {
                    const MuellerMatrix rayleighEvent = eventMuellerMatrix(
                        viewSample.position,
                        viewSample.position,
                        observer,
                        sunSample.direction,
                        outgoingDirection,
                        true,
                        aerosolCoefficients
                    );
                    totalContribution = addStokes(
                        totalContribution,
                        scaleStokes(
                            multiply(rayleighEvent, directSun),
                            optics.rayleigh_scattering_m_inv
                        )
                    );
                }

                if (optics.aerosol_scattering_m_inv > 0.0) {
                    const MuellerMatrix aerosolEvent = eventMuellerMatrix(
                        viewSample.position,
                        viewSample.position,
                        observer,
                        sunSample.direction,
                        outgoingDirection,
                        false,
                        aerosolCoefficients
                    );
                    totalContribution = addStokes(
                        totalContribution,
                        scaleStokes(
                            multiply(aerosolEvent, directSun),
                            optics.aerosol_scattering_m_inv
                        )
                    );
                }
            }

            tauView += optics.extinction_m_inv * ds;
        }

        if (timingSummary != nullptr) {
            timingSummary->first_order_seconds =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - firstOrderStart).count();
        }
        if (stageCallback) {
            stageCallback(
                SampleProgress::Stage::first_order_band_complete,
                timingSummary != nullptr ? *timingSummary : DirectionEstimate::TimingSummary {},
                bandIndex + 1,
                bands.size()
            );
        }
    }

    if (timingSummary != nullptr) {
        timingSummary->first_order_seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - firstOrderStart).count();
    }

    return totalContribution;
}

StokesVector traceBandPath(
    const SimulationContext &context,
    const SimulationConfig &config,
    const SpectralBand &band,
    const Vec3 &position,
    const Vec3 &outgoingDirection,
    const MuellerMatrix &throughput,
    int eventIndex,
    std::mt19937 &rng
)
{
    if (eventIndex >= std::max(8, config.monte_carlo.max_events_guard)) {
        return {0.0, 0.0, 0.0, 0.0};
    }

    std::uniform_real_distribution<double> uniform01(0.0, 1.0);
    StokesVector totalContribution {0.0, 0.0, 0.0, 0.0};
    const SampledInteraction interaction = sampleInteraction(
        context.atmosphere,
        context.wavelengthManager,
        config,
        position + outgoingDirection * 1.0e-3,
        outgoingDirection,
        band.wavelength_nm,
        rng
    );
    if (interaction.absorbed) {
        return totalContribution;
    }

    if (!interaction.scattered) {
        if (interaction.hit_ground) {
            const double groundDistance = distanceToSphere(
                position + outgoingDirection * 1.0e-3,
                outgoingDirection,
                R_EARTH_M
            );
            if (std::isfinite(groundDistance)) {
                const Vec3 groundPoint = position + outgoingDirection * groundDistance;
                const Vec3 groundNormal = normalize(groundPoint);
                const std::vector<SolarDiskSample> &sunSamples = context.sun_samples;
                for (const SolarDiskSample &sunSample : sunSamples) {
                    const double tauSun = opticalDepthToSun(
                        context.atmosphere,
                        context.wavelengthManager,
                        config,
                        groundPoint + groundNormal * 1.0,
                        sunSample.direction,
                        band.wavelength_nm
                    );
                    if (!std::isfinite(tauSun)) {
                        continue;
                    }

                    const double directIrradiance =
                        band.solar_irradiance_w_m2_nm * sunSample.weight * std::exp(-tauSun);
                    const double radiance = context.surface.directSolarRadiance(
                        groundNormal,
                        sunSample.direction,
                        outgoingDirection * -1.0,
                        band.wavelength_nm,
                        directIrradiance
                    );
                    if (radiance > 0.0) {
                        totalContribution = addStokes(
                            totalContribution,
                            multiply(throughput, initUnpolarized(radiance))
                        );
                    }
                }
            }
        }
        return totalContribution;
    }

    const bool suppressDirectSource =
        (config.monte_carlo.deterministic_single_scatter && eventIndex == 0) ||
        (config.monte_carlo.deterministic_second_scatter && eventIndex == 1);
    if (!suppressDirectSource) {
        const std::vector<SolarDiskSample> &sunSamples = context.sun_samples;
        for (const SolarDiskSample &sunSample : sunSamples) {
            const double tauSun = opticalDepthToSun(
                context.atmosphere,
                context.wavelengthManager,
                config,
                interaction.position + sunSample.direction * 1.0,
                sunSample.direction,
                band.wavelength_nm
            );
            if (!std::isfinite(tauSun)) {
                continue;
            }

            const double muSun = clamp(dot(sunSample.direction, outgoingDirection), -1.0, 1.0);
            const PhaseMatrixCoefficients aerosolCoefficients =
                context.phaseTable.coefficients(band.wavelength_nm, muSun);
            const MuellerMatrix sourceEvent = eventMuellerMatrix(
                interaction.position,
                interaction.position,
                position,
                sunSample.direction,
                outgoingDirection,
                interaction.is_rayleigh,
                aerosolCoefficients
            );
            const StokesVector directSun = initUnpolarized(
                band.solar_irradiance_w_m2_nm * sunSample.weight * std::exp(-tauSun)
            );
            totalContribution = addStokes(
                totalContribution,
                multiply(
                    throughput,
                    multiply(
                        scaleMueller(sourceEvent, interaction.collision_weight),
                        directSun
                    )
                )
            );
        }
    }

    if (config.monte_carlo.single_scatter_only) {
        return totalContribution;
    }

    const double muSun = clamp(dot(context.sun.direction, outgoingDirection), -1.0, 1.0);
    const bool useSourceGuiding =
        config.monte_carlo.source_guided_first_scatter &&
        eventIndex <= std::max(0, config.monte_carlo.source_guided_max_event_index);
    const bool useTwilightHigherOrderGuiding =
        config.monte_carlo.twilight_higher_order_guiding &&
        context.sun.zenith_deg > 90.0 &&
        recursiveMinimumScatteringOrder(config, eventIndex) >= 3;
    const bool useGuidedContinuation = useSourceGuiding || useTwilightHigherOrderGuiding;

    const GuidedProposalParameters proposal = guidedProposalParameters(
        context,
        config,
        interaction.is_rayleigh,
        interaction.position,
        outgoingDirection,
        eventIndex
    );

    int branchCount = 1;
    if (useGuidedContinuation) {
        branchCount = std::max(1, config.monte_carlo.source_guided_first_scatter_branches);
        const double muSunAbs = std::abs(muSun);
        if (useSourceGuiding && interaction.is_rayleigh && muSunAbs < 0.7) {
            branchCount = std::max(branchCount, config.monte_carlo.rayleigh_source_guided_first_scatter_branches);
        }
        if (proposal.use_twilight_higher_order) {
            int activeTwilightComponents = 0;
            activeTwilightComponents += proposal.phase_fraction > 1.0e-9 ? 1 : 0;
            activeTwilightComponents += proposal.tangent_fraction > 1.0e-9 ? 1 : 0;
            activeTwilightComponents += proposal.horizon_fraction > 1.0e-9 ? 1 : 0;
            branchCount = std::max(branchCount, config.monte_carlo.twilight_higher_order_branches);
            branchCount = std::max(branchCount, activeTwilightComponents);
        }
    }

    int polarizationBranches = 0;
    int baseBranches = branchCount;
    double polarizationFraction = proposal.polarization_fraction;
    if (useSourceGuiding && interaction.is_rayleigh && branchCount > 1 && !proposal.use_twilight_higher_order) {
        if (polarizationFraction > 0.0) {
            const int requestedBranches = static_cast<int>(std::llround(branchCount * polarizationFraction));
            if (requestedBranches > 0) {
                polarizationBranches = std::clamp(
                    std::min(config.monte_carlo.rayleigh_polarization_guided_branches, requestedBranches),
                    1,
                    branchCount - 1
                );
                baseBranches = branchCount - polarizationBranches;
            }
        }
    }

    std::array<int, 3> twilightComponentBranches {0, 0, 0};
    if (proposal.use_twilight_higher_order && branchCount > 1 && polarizationBranches == 0) {
        twilightComponentBranches = stratifiedTwilightBranchCounts(
            {
                proposal.phase_fraction,
                proposal.tangent_fraction,
                proposal.horizon_fraction,
            },
            branchCount
        );
    }

    for (int branchIndex = 0; branchIndex < branchCount; ++branchIndex) {
        ContinuationSample sample {};
        double componentWeight = 1.0;
        int componentBranchCount = branchCount;
        if (branchIndex < polarizationBranches) {
            sample = samplePolarizationGuidedDirection(
                context,
                config,
                band.wavelength_nm,
                interaction.position,
                outgoingDirection,
                eventIndex,
                rng
            );
            componentWeight = polarizationFraction;
            componentBranchCount = polarizationBranches;
        } else if (proposal.use_twilight_higher_order && polarizationBranches == 0) {
            TwilightGuidedComponent forcedComponent = TwilightGuidedComponent::phase;
            int localIndex = branchIndex;
            if (localIndex < twilightComponentBranches[0]) {
                forcedComponent = TwilightGuidedComponent::phase;
                componentWeight = proposal.phase_fraction;
                componentBranchCount = twilightComponentBranches[0];
            } else {
                localIndex -= twilightComponentBranches[0];
                if (localIndex < twilightComponentBranches[1]) {
                    forcedComponent = TwilightGuidedComponent::tangent;
                    componentWeight = proposal.tangent_fraction;
                    componentBranchCount = twilightComponentBranches[1];
                } else {
                    forcedComponent = TwilightGuidedComponent::horizon;
                    componentWeight = proposal.horizon_fraction;
                    componentBranchCount = twilightComponentBranches[2];
                }
            }
            sample = sampleContinuationDirection(
                context,
                config,
                band.wavelength_nm,
                interaction.is_rayleigh,
                interaction.position,
                outgoingDirection,
                eventIndex,
                useGuidedContinuation,
                rng,
                forcedComponent
            );
        } else {
            sample = sampleContinuationDirection(
                context,
                config,
                band.wavelength_nm,
                interaction.is_rayleigh,
                interaction.position,
                outgoingDirection,
                eventIndex,
                useGuidedContinuation,
                rng
            );
            if (polarizationBranches > 0) {
                componentWeight = 1.0 - polarizationFraction;
                componentBranchCount = baseBranches;
            }
        }

        sample.estimator_weight = componentWeight / std::max(1, componentBranchCount);
        const PhaseMatrixCoefficients aerosolCoefficients =
            context.phaseTable.coefficients(band.wavelength_nm, sample.cos_theta);
        const MuellerMatrix scatterMatrix = eventMuellerMatrix(
            interaction.position,
            interaction.position,
            position,
            sample.incoming_direction,
            outgoingDirection,
            interaction.is_rayleigh,
            aerosolCoefficients
        );

        MuellerMatrix branchThroughput = multiply(
            throughput,
            scaleMueller(
                scatterMatrix,
                interaction.collision_weight *
                sample.estimator_weight /
                std::max(1.0e-12, sample.pdf)
            )
        );

        const double importance = throughputImportance(branchThroughput);
        if (importance < config.monte_carlo.russian_roulette_threshold) {
            const double survivalProbability = 0.5;
            if (uniform01(rng) > survivalProbability) {
                continue;
            }
            branchThroughput = scaleMueller(branchThroughput, 1.0 / survivalProbability);
        }

        StokesVector recursiveContribution = traceBandPath(
            context,
            config,
            band,
            interaction.position,
            sample.incoming_direction,
            branchThroughput,
            eventIndex + 1,
            rng
        );
        if (
            twilightDepolarizationActive(config, context) &&
            recursiveMinimumScatteringOrder(config, eventIndex) >= 3
        ) {
            recursiveContribution = boostUnpolarizedIntensity(
                recursiveContribution,
                config.monte_carlo.twilight_higher_order_intensity_boost
            );
            recursiveContribution = scalePolarization(
                recursiveContribution,
                config.monte_carlo.twilight_higher_order_polarization_scale
            );
        }
        totalContribution = addStokes(totalContribution, recursiveContribution);
    }

    return totalContribution;
}

StokesVector traceOnePath(
    const SimulationContext &context,
    const SimulationConfig &config,
    double viewZenithDeg,
    double viewAzimuthDeg,
    std::mt19937 &rng
)
{
    const Vec3 observer = observerPosition(config.observer);
    const LocalBasis observerBasis = localBasisAtPosition(observer);
    const Vec3 initialDirection = directionFromZenithAzimuth(viewZenithDeg, viewAzimuthDeg, observerBasis);
    StokesVector totalContribution {0.0, 0.0, 0.0, 0.0};

    for (const SpectralBand &band : context.active_bands) {
        totalContribution = addStokes(
            totalContribution,
            traceBandPath(
                context,
                config,
                band,
                observer,
                initialDirection,
                identityMueller(),
                0,
                rng
            )
        );
    }

    return totalContribution;
}

DirectionEstimate estimateDirectionMoments(
    const SimulationContext &context,
    const SimulationConfig &config,
    const SkyDirection &direction,
    int sampleCount,
    std::size_t directionIndex,
    const DirectionStageCallback &stageCallback = {},
    DirectionCheckpointState *checkpointState = nullptr
)
{
    DirectionEstimate estimate;
    if (checkpointState != nullptr) {
        estimate.timing = internalTimingSummary(checkpointState->timing);
    }
    estimate.timing.spectral_band_count = context.active_bands.size();
    estimate.timing.mc_sample_count = sampleCount;

    if (stageCallback) {
        stageCallback(
            SampleProgress::Stage::direction_started,
            estimate.timing,
            0,
            0
        );
    }

    StokesVector deterministicSingleScatter {0.0, 0.0, 0.0, 0.0};
    if (config.monte_carlo.deterministic_single_scatter) {
        if (checkpointState != nullptr && checkpointState->has_first_order) {
            deterministicSingleScatter = checkpointState->first_order;
        } else {
            deterministicSingleScatter = deterministicSingleScatterContribution(
                context,
                config,
                direction.zenith_deg,
                direction.azimuth_deg,
                &estimate.timing,
                stageCallback
            );
            if (checkpointState != nullptr) {
                checkpointState->has_first_order = true;
                checkpointState->first_order = deterministicSingleScatter;
                checkpointState->timing = publicTimingSummary(estimate.timing);
            }
        }
    }
    if (stageCallback) {
        stageCallback(
            SampleProgress::Stage::first_order_complete,
            estimate.timing,
            1,
            1
        );
    }

    SecondOrderBreakdown deterministicSecondScatter {};
    if (config.monte_carlo.deterministic_second_scatter &&
        config.monte_carlo.deterministic_single_scatter &&
        !config.monte_carlo.single_scatter_only) {
        if (checkpointState != nullptr && checkpointState->has_second_order) {
            deterministicSecondScatter.total = checkpointState->second_total;
            deterministicSecondScatter.rr = checkpointState->second_rr;
            deterministicSecondScatter.ar = checkpointState->second_ar;
            deterministicSecondScatter.ra = checkpointState->second_ra;
            deterministicSecondScatter.aa = checkpointState->second_aa;
        } else {
            deterministicSecondScatter = deterministicSecondScatterContribution(
                context,
                config,
                direction.zenith_deg,
                direction.azimuth_deg,
                &estimate.timing,
                stageCallback
            );
            if (checkpointState != nullptr) {
                checkpointState->has_second_order = true;
                checkpointState->second_total = deterministicSecondScatter.total;
                checkpointState->second_rr = deterministicSecondScatter.rr;
                checkpointState->second_ar = deterministicSecondScatter.ar;
                checkpointState->second_ra = deterministicSecondScatter.ra;
                checkpointState->second_aa = deterministicSecondScatter.aa;
                checkpointState->timing = publicTimingSummary(estimate.timing);
            }
        }
    }

    estimate.first_order = deterministicSingleScatter;
    estimate.second_order = deterministicSecondScatter;

    if (config.monte_carlo.deterministic_second_scatter && stageCallback) {
        stageCallback(
            SampleProgress::Stage::second_order_complete,
            estimate.timing,
            1,
            1
        );
    }

    if (config.monte_carlo.single_scatter_only && config.monte_carlo.deterministic_single_scatter) {
        return estimate;
    }

    RunningMoments moments;
    if (checkpointState != nullptr) {
        moments.count = checkpointState->higher_order.completed_samples;
        moments.mean = checkpointState->higher_order.mean;
        moments.m2 = checkpointState->higher_order.m2;
    }
    const int higherOrderInitialCompletedSamples = std::max(0, moments.count);
    const int higherOrderProgressInterval =
        sampleCount <= 16 ? 1 : (sampleCount <= 64 ? 4 : (sampleCount <= 256 ? 16 : 32));
    const double higherOrderBaseSeconds = estimate.timing.higher_order_seconds;
    const auto higherOrderStart = std::chrono::steady_clock::now();
    const auto reportHigherOrderProgress = [&](int completedSamples) {
        if (!stageCallback || completedSamples <= higherOrderInitialCompletedSamples) {
            return;
        }
        const bool shouldReport =
            completedSamples == higherOrderInitialCompletedSamples + 1 ||
            completedSamples == sampleCount ||
            (completedSamples % higherOrderProgressInterval) == 0;
        if (!shouldReport) {
            return;
        }

        DirectionEstimate::TimingSummary timing = estimate.timing;
        timing.higher_order_seconds = higherOrderBaseSeconds +
            std::chrono::duration<double>(std::chrono::steady_clock::now() - higherOrderStart).count();
        stageCallback(
            SampleProgress::Stage::higher_order_progress,
            timing,
            static_cast<std::size_t>(completedSamples),
            static_cast<std::size_t>(sampleCount)
        );
    };
    if (checkpointState != nullptr) {
        const int startSample = std::max(0, moments.count);
        for (int sampleIndex = startSample; sampleIndex < sampleCount; ++sampleIndex) {
            std::mt19937 sampleRng(seedForDirectionSample(
                config.monte_carlo.random_seed,
                directionIndex,
                direction.zenith_deg,
                direction.azimuth_deg,
                sampleIndex
            ));
            moments.update(traceOnePath(
                context,
                config,
                direction.zenith_deg,
                direction.azimuth_deg,
                sampleRng
            ));
            if (checkpointState != nullptr) {
                estimate.timing.higher_order_seconds = higherOrderBaseSeconds +
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - higherOrderStart).count();
                checkpointState->higher_order.completed_samples = moments.count;
                checkpointState->higher_order.mean = moments.mean;
                checkpointState->higher_order.m2 = moments.m2;
                checkpointState->timing = publicTimingSummary(estimate.timing);
            }
            reportHigherOrderProgress(moments.count);
        }
    } else {
        std::mt19937 rng(seedForDirection(
            config.monte_carlo.random_seed,
            directionIndex,
            direction.zenith_deg,
            direction.azimuth_deg
        ));
        for (int sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex) {
            moments.update(traceOnePath(
                context,
                config,
                direction.zenith_deg,
                direction.azimuth_deg,
                rng
            ));
            if (checkpointState != nullptr) {
                estimate.timing.higher_order_seconds = higherOrderBaseSeconds +
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - higherOrderStart).count();
                checkpointState->higher_order.completed_samples = moments.count;
                checkpointState->higher_order.mean = moments.mean;
                checkpointState->higher_order.m2 = moments.m2;
                checkpointState->timing = publicTimingSummary(estimate.timing);
            }
            reportHigherOrderProgress(moments.count);
        }
    }
    estimate.timing.higher_order_seconds = higherOrderBaseSeconds +
        std::chrono::duration<double>(std::chrono::steady_clock::now() - higherOrderStart).count();
    estimate.higher_order = moments.mean;
    estimate.higher_variance = moments.variance();
    if (checkpointState != nullptr) {
        checkpointState->higher_order.completed_samples = moments.count;
        checkpointState->higher_order.mean = moments.mean;
        checkpointState->higher_order.m2 = moments.m2;
        checkpointState->timing = publicTimingSummary(estimate.timing);
    }
    if (stageCallback) {
        stageCallback(
            SampleProgress::Stage::higher_order_complete,
            estimate.timing,
            static_cast<std::size_t>(moments.count),
            static_cast<std::size_t>(sampleCount)
        );
    }
    return estimate;
}

double estimateHemisphericFlux(const SkyResult &result)
{
    const double zenithStep = 90.0 / static_cast<double>(result.config.monte_carlo.zenith_bins);
    const double azimuthStep = 360.0 / static_cast<double>(result.config.monte_carlo.azimuth_bins);
    const double dTheta = zenithStep * PI / 180.0;
    const double dPhi = azimuthStep * PI / 180.0;

    double flux = 0.0;
    for (const SkyBinResult &bin : result.bins) {
        const double zenithRad = bin.zenith_deg * PI / 180.0;
        const double dOmega = std::sin(zenithRad) * dTheta * dPhi;
        flux += bin.mean.I * std::cos(zenithRad) * dOmega;
    }
    return flux;
}
}

SimulationConfig defaultSimulationConfig()
{
    SimulationConfig config;
    config.config_path = "monte_carlo_cpp/config/default_clear_sky.cfg";
    config.atmosphere.profile_csv = "monte_carlo_cpp/data/atmosphere/clear_sky_midlatitude.csv";
    config.spectral.solar_spectrum_csv = "monte_carlo_cpp/data/optics/solar_irradiance_reference.csv";
    config.spectral.instrument_response_csv = "";
    config.spectral.rayleigh_cross_section_csv = "";
    config.spectral.ozone_cross_section_csv = "monte_carlo_cpp/data/optics/ozone_cross_section_reference.csv";
    config.spectral.o2_cross_section_csv = "";
    config.spectral.o4_cross_section_csv = "";
    config.spectral.h2o_cross_section_csv = "";
    config.spectral.no2_cross_section_csv = "";
    config.spectral.aerosol_phase_matrix_csv = "monte_carlo_cpp/data/optics/aerosol_phase_matrix_reference.csv";
    config.surface.albedo_csv = "monte_carlo_cpp/data/optics/lambertian_land_albedo.csv";
    return config;
}

SimulationConfig loadSimulationConfig(const std::string &config_path)
{
    SimulationConfig config = defaultSimulationConfig();
    config.config_path = config_path;

    const std::filesystem::path configFile(config_path);
    const std::filesystem::path configDir = configFile.parent_path();
    const std::map<std::string, std::string> values = loadKeyValueConfig(config_path);

    setIfPresent(values, "case_id", config.output.case_id);
    setIfPresent(values, "output_dir", config.output.output_dir);
    setIfPresent(values, "profile_csv", config.atmosphere.profile_csv);
    setIfPresent(values, "solar_spectrum_csv", config.spectral.solar_spectrum_csv);
    setIfPresent(values, "instrument_response_csv", config.spectral.instrument_response_csv);
    setIfPresent(values, "rayleigh_cross_section_csv", config.spectral.rayleigh_cross_section_csv);
    setIfPresent(values, "ozone_cross_section_csv", config.spectral.ozone_cross_section_csv);
    setIfPresent(values, "o2_cross_section_csv", config.spectral.o2_cross_section_csv);
    setIfPresent(values, "o4_cross_section_csv", config.spectral.o4_cross_section_csv);
    setIfPresent(values, "h2o_cross_section_csv", config.spectral.h2o_cross_section_csv);
    setIfPresent(values, "no2_cross_section_csv", config.spectral.no2_cross_section_csv);
    setIfPresent(values, "aerosol_phase_matrix_csv", config.spectral.aerosol_phase_matrix_csv);
    setIfPresent(values, "surface_albedo_csv", config.surface.albedo_csv);
    setIfPresent(values, "surface_model", config.surface.surface_model);
    setIfPresent(values, "surface_parameter_csv", config.surface.surface_parameter_csv);
    setIfPresent(values, "benchmark_case_config", config.output.benchmark_case_config);
    setIfPresent(values, "benchmark_reference_csv", config.output.benchmark_reference_csv);
    setIfPresent(values, "benchmark_metadata_json", config.output.benchmark_metadata_json);
    setIfPresent(values, "measurement_case_config", config.output.measurement_case_config);
    setIfPresent(values, "measurement_reference_csv", config.output.measurement_reference_csv);
    setIfPresent(values, "measurement_metadata_json", config.output.measurement_metadata_json);
    setIfPresent(values, "paper_primary_measurement_case_config", config.output.paper_primary_measurement_case_config);
    setIfPresent(values, "paper_case_provenance_json", config.output.paper_case_provenance_json);

    setNumericIfPresent(values, "observer_latitude_deg", config.observer.latitude_deg);
    setNumericIfPresent(values, "observer_longitude_deg", config.observer.longitude_deg);
    setNumericIfPresent(values, "observer_altitude_m", config.observer.altitude_m);
    setNumericIfPresent(values, "solar_declination_deg", config.solar.solar_declination_deg);
    setNumericIfPresent(values, "hour_angle_deg", config.solar.hour_angle_deg);
    setNumericIfPresent(values, "solar_zenith_deg", config.solar.zenith_deg);
    setNumericIfPresent(values, "solar_azimuth_deg", config.solar.azimuth_deg);
    setBoolIfPresent(values, "use_explicit_solar_angles", config.solar.use_explicit_angles);
    setBoolIfPresent(values, "finite_solar_disk", config.solar.finite_solar_disk);
    setNumericIfPresent(values, "solar_angular_radius_deg", config.solar.solar_angular_radius_deg);
    setNumericIfPresent(values, "solar_disk_quadrature_nodes", config.solar.solar_disk_quadrature_nodes);
    setNumericIfPresent(values, "top_of_atmosphere_altitude_m", config.atmosphere.top_of_atmosphere_altitude_m);
    setNumericIfPresent(values, "min_wavelength_nm", config.spectral.min_wavelength_nm);
    setNumericIfPresent(values, "max_wavelength_nm", config.spectral.max_wavelength_nm);
    setNumericIfPresent(values, "wavelength_step_nm", config.spectral.wavelength_step_nm);
    setNumericIfPresent(values, "zenith_bins", config.monte_carlo.zenith_bins);
    setNumericIfPresent(values, "azimuth_bins", config.monte_carlo.azimuth_bins);
    setNumericIfPresent(values, "photons_per_bin", config.monte_carlo.photons_per_bin);
    setNumericIfPresent(values, "russian_roulette_threshold", config.monte_carlo.russian_roulette_threshold);
    setNumericIfPresent(values, "random_seed", config.monte_carlo.random_seed);
    setNumericIfPresent(values, "max_events_guard", config.monte_carlo.max_events_guard);
    setNumericIfPresent(values, "optical_depth_step_scale", config.monte_carlo.optical_depth_step_scale);
    setNumericIfPresent(values, "min_optical_depth_steps", config.monte_carlo.min_optical_depth_steps);
    setNumericIfPresent(values, "line_integral_step_scale", config.monte_carlo.line_integral_step_scale);
    setNumericIfPresent(values, "min_line_integral_steps", config.monte_carlo.min_line_integral_steps);
    setBoolIfPresent(values, "single_scatter_only", config.monte_carlo.single_scatter_only);
    setBoolIfPresent(values, "deterministic_single_scatter", config.monte_carlo.deterministic_single_scatter);
    setBoolIfPresent(values, "deterministic_second_scatter", config.monte_carlo.deterministic_second_scatter);
    setNumericIfPresent(values, "second_scatter_view_steps", config.monte_carlo.second_scatter_view_steps);
    setNumericIfPresent(values, "second_scatter_ray_steps", config.monte_carlo.second_scatter_ray_steps);
    setNumericIfPresent(values, "second_scatter_mu_nodes", config.monte_carlo.second_scatter_mu_nodes);
    setNumericIfPresent(values, "second_scatter_phi_nodes", config.monte_carlo.second_scatter_phi_nodes);
    setBoolIfPresent(values, "twilight_second_scatter_adaptive", config.monte_carlo.twilight_second_scatter_adaptive);
    setNumericIfPresent(
        values,
        "twilight_second_scatter_zenith_threshold_deg",
        config.monte_carlo.twilight_second_scatter_zenith_threshold_deg
    );
    setNumericIfPresent(
        values,
        "twilight_second_scatter_min_view_steps",
        config.monte_carlo.twilight_second_scatter_min_view_steps
    );
    setNumericIfPresent(
        values,
        "twilight_second_scatter_min_ray_steps",
        config.monte_carlo.twilight_second_scatter_min_ray_steps
    );
    setNumericIfPresent(
        values,
        "twilight_second_scatter_min_mu_nodes",
        config.monte_carlo.twilight_second_scatter_min_mu_nodes
    );
    setNumericIfPresent(
        values,
        "twilight_second_scatter_min_phi_nodes",
        config.monte_carlo.twilight_second_scatter_min_phi_nodes
    );
    setBoolIfPresent(values, "source_guided_first_scatter", config.monte_carlo.source_guided_first_scatter);
    setNumericIfPresent(
        values,
        "source_guided_max_event_index",
        config.monte_carlo.source_guided_max_event_index
    );
    setNumericIfPresent(
        values,
        "source_guided_first_scatter_branches",
        config.monte_carlo.source_guided_first_scatter_branches
    );
    setNumericIfPresent(values, "source_guided_phase_fraction", config.monte_carlo.source_guided_phase_fraction);
    setNumericIfPresent(
        values,
        "source_guided_cone_half_angle_deg",
        config.monte_carlo.source_guided_cone_half_angle_deg
    );
    setNumericIfPresent(
        values,
        "rayleigh_source_guided_phase_fraction",
        config.monte_carlo.rayleigh_source_guided_phase_fraction
    );
    setNumericIfPresent(
        values,
        "rayleigh_source_guided_cone_half_angle_deg",
        config.monte_carlo.rayleigh_source_guided_cone_half_angle_deg
    );
    setNumericIfPresent(
        values,
        "rayleigh_source_guided_first_scatter_branches",
        config.monte_carlo.rayleigh_source_guided_first_scatter_branches
    );
    setNumericIfPresent(
        values,
        "rayleigh_polarization_guided_fraction",
        config.monte_carlo.rayleigh_polarization_guided_fraction
    );
    setNumericIfPresent(
        values,
        "rayleigh_polarization_guided_mu_half_width",
        config.monte_carlo.rayleigh_polarization_guided_mu_half_width
    );
    setNumericIfPresent(
        values,
        "rayleigh_polarization_guided_branches",
        config.monte_carlo.rayleigh_polarization_guided_branches
    );
    setBoolIfPresent(values, "twilight_limb_guiding", config.monte_carlo.twilight_limb_guiding);
    setNumericIfPresent(
        values,
        "twilight_limb_phase_fraction",
        config.monte_carlo.twilight_limb_phase_fraction
    );
    setNumericIfPresent(
        values,
        "twilight_limb_cone_half_angle_deg",
        config.monte_carlo.twilight_limb_cone_half_angle_deg
    );
    setNumericIfPresent(
        values,
        "twilight_limb_elevation_deg",
        config.monte_carlo.twilight_limb_elevation_deg
    );
    setBoolIfPresent(
        values,
        "twilight_higher_order_guiding",
        config.monte_carlo.twilight_higher_order_guiding
    );
    setNumericIfPresent(
        values,
        "twilight_higher_order_branches",
        config.monte_carlo.twilight_higher_order_branches
    );
    setNumericIfPresent(
        values,
        "twilight_higher_order_phase_fraction",
        config.monte_carlo.twilight_higher_order_phase_fraction
    );
    setNumericIfPresent(
        values,
        "twilight_higher_order_tangent_fraction",
        config.monte_carlo.twilight_higher_order_tangent_fraction
    );
    setNumericIfPresent(
        values,
        "twilight_higher_order_horizon_fraction",
        config.monte_carlo.twilight_higher_order_horizon_fraction
    );
    setNumericIfPresent(
        values,
        "twilight_higher_order_tangent_cone_half_angle_deg",
        config.monte_carlo.twilight_higher_order_tangent_cone_half_angle_deg
    );
    setNumericIfPresent(
        values,
        "twilight_higher_order_horizon_cone_half_angle_deg",
        config.monte_carlo.twilight_higher_order_horizon_cone_half_angle_deg
    );
    setNumericIfPresent(
        values,
        "twilight_higher_order_horizon_elevation_deg",
        config.monte_carlo.twilight_higher_order_horizon_elevation_deg
    );
    setBoolIfPresent(
        values,
        "twilight_order_depolarization",
        config.monte_carlo.twilight_order_depolarization
    );
    setNumericIfPresent(
        values,
        "twilight_second_order_polarization_scale",
        config.monte_carlo.twilight_second_order_polarization_scale
    );
    setNumericIfPresent(
        values,
        "twilight_higher_order_polarization_scale",
        config.monte_carlo.twilight_higher_order_polarization_scale
    );
    setNumericIfPresent(
        values,
        "twilight_second_order_intensity_boost",
        config.monte_carlo.twilight_second_order_intensity_boost
    );
    setNumericIfPresent(
        values,
        "twilight_higher_order_intensity_boost",
        config.monte_carlo.twilight_higher_order_intensity_boost
    );
    setNumericIfPresent(
        values,
        "twilight_second_order_boost_zenith_cutoff_deg",
        config.monte_carlo.twilight_second_order_boost_zenith_cutoff_deg
    );
    setNumericIfPresent(values, "default_surface_albedo", config.surface.default_albedo);
    setNumericIfPresent(values, "ocean_wind_speed_m_s", config.surface.ocean_wind_speed_m_s);
    setBoolIfPresent(values, "paper_primary_measurement_frozen", config.output.paper_primary_measurement_frozen);
    setBoolIfPresent(values, "strict_paper_mode", config.output.strict_paper_mode);
    setNumericIfPresent(values, "benchmark_mask_fraction_of_peak", config.validation.benchmark_mask_fraction_of_peak);
    setNumericIfPresent(values, "measurement_mask_fraction_of_peak", config.validation.measurement_mask_fraction_of_peak);
    setNumericIfPresent(values, "median_intensity_error_limit", config.validation.median_intensity_error_limit);
    setNumericIfPresent(values, "p95_intensity_error_limit", config.validation.p95_intensity_error_limit);
    setNumericIfPresent(values, "median_dolp_abs_error_limit", config.validation.median_dolp_abs_error_limit);
    setNumericIfPresent(values, "p95_dolp_abs_error_limit", config.validation.p95_dolp_abs_error_limit);
    setNumericIfPresent(values, "median_aop_error_deg_limit", config.validation.median_aop_error_deg_limit);
    setNumericIfPresent(values, "p95_aop_error_deg_limit", config.validation.p95_aop_error_deg_limit);
    setNumericIfPresent(
        values,
        "solar_vertical_signed_dolp_bias_limit",
        config.validation.solar_vertical_signed_dolp_bias_limit
    );
    setNumericIfPresent(values, "normalized_rmse_limit", config.validation.normalized_rmse_limit);
    setNumericIfPresent(values, "brightest_region_deg_limit", config.validation.brightest_region_deg_limit);
    setNumericIfPresent(
        values,
        "neutral_point_location_deg_limit",
        config.validation.neutral_point_location_deg_limit
    );

    if (config.monte_carlo.deterministic_second_scatter) {
        config.monte_carlo.deterministic_single_scatter = true;
    }

    config.output.output_dir = resolvePath(configDir, config.output.output_dir).string();
    config.atmosphere.profile_csv = resolvePath(configDir, config.atmosphere.profile_csv).string();
    config.spectral.solar_spectrum_csv = resolvePath(configDir, config.spectral.solar_spectrum_csv).string();
    if (!config.spectral.instrument_response_csv.empty()) {
        config.spectral.instrument_response_csv =
            resolvePath(configDir, config.spectral.instrument_response_csv).string();
    }
    if (!config.spectral.rayleigh_cross_section_csv.empty()) {
        config.spectral.rayleigh_cross_section_csv =
            resolvePath(configDir, config.spectral.rayleigh_cross_section_csv).string();
    }
    config.spectral.ozone_cross_section_csv = resolvePath(configDir, config.spectral.ozone_cross_section_csv).string();
    if (!config.spectral.o2_cross_section_csv.empty()) {
        config.spectral.o2_cross_section_csv = resolvePath(configDir, config.spectral.o2_cross_section_csv).string();
    }
    if (!config.spectral.o4_cross_section_csv.empty()) {
        config.spectral.o4_cross_section_csv = resolvePath(configDir, config.spectral.o4_cross_section_csv).string();
    }
    if (!config.spectral.h2o_cross_section_csv.empty()) {
        config.spectral.h2o_cross_section_csv = resolvePath(configDir, config.spectral.h2o_cross_section_csv).string();
    }
    if (!config.spectral.no2_cross_section_csv.empty()) {
        config.spectral.no2_cross_section_csv = resolvePath(configDir, config.spectral.no2_cross_section_csv).string();
    }
    config.spectral.aerosol_phase_matrix_csv = resolvePath(configDir, config.spectral.aerosol_phase_matrix_csv).string();
    config.surface.albedo_csv = resolvePath(configDir, config.surface.albedo_csv).string();
    if (!config.surface.surface_parameter_csv.empty()) {
        config.surface.surface_parameter_csv = resolvePath(configDir, config.surface.surface_parameter_csv).string();
    }
    if (!config.output.benchmark_reference_csv.empty()) {
        config.output.benchmark_reference_csv = resolvePath(configDir, config.output.benchmark_reference_csv).string();
    }
    if (!config.output.benchmark_case_config.empty()) {
        config.output.benchmark_case_config = resolvePathList(configDir, config.output.benchmark_case_config);
    }
    if (!config.output.benchmark_metadata_json.empty()) {
        config.output.benchmark_metadata_json = resolvePath(configDir, config.output.benchmark_metadata_json).string();
    }
    if (!config.output.measurement_reference_csv.empty()) {
        config.output.measurement_reference_csv = resolvePath(configDir, config.output.measurement_reference_csv).string();
    }
    if (!config.output.measurement_case_config.empty()) {
        config.output.measurement_case_config = resolvePathList(configDir, config.output.measurement_case_config);
    }
    if (!config.output.measurement_metadata_json.empty()) {
        config.output.measurement_metadata_json = resolvePath(configDir, config.output.measurement_metadata_json).string();
    }
    if (!config.output.paper_primary_measurement_case_config.empty()) {
        config.output.paper_primary_measurement_case_config = resolvePath(
            configDir,
            config.output.paper_primary_measurement_case_config
        ).string();
    }
    if (!config.output.paper_case_provenance_json.empty()) {
        config.output.paper_case_provenance_json = resolvePath(
            configDir,
            config.output.paper_case_provenance_json
        ).string();
    }

    if (config.output.strict_paper_mode && empiricalTwilightTuningEnabled(config.monte_carlo)) {
        throw std::runtime_error(
            "strict_paper_mode=true rejects empirical twilight tuning fields. "
            "Disable twilight_order_depolarization and reset twilight_*_polarization_scale and "
            "twilight_*_intensity_boost to 1.0 in paper configs."
        );
    }

    config.config_hash = hashFile(config_path);
    return config;
}

std::vector<SkyBinResult> sampleSkyDirections(
    const SimulationConfig &config,
    const std::vector<SkyDirection> &directions,
    int samples_override,
    const std::function<void(const SampleProgress &)> &progress_callback,
    const std::vector<std::size_t> *original_direction_indices
)
{
    const SimulationContext context = buildSimulationContext(config);
    const int sampleCount = samples_override > 0 ? samples_override : config.monte_carlo.photons_per_bin;
    std::vector<SkyBinResult> bins(directions.size());
    const auto startTime = std::chrono::steady_clock::now();
    std::atomic<std::size_t> completedCount {0};
    std::mutex progressMutex;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int directionIndex = 0; directionIndex < static_cast<int>(directions.size()); ++directionIndex) {
        const std::size_t reportedDirectionIndex =
            original_direction_indices != nullptr && directionIndex < static_cast<int>(original_direction_indices->size())
                ? (*original_direction_indices)[static_cast<std::size_t>(directionIndex)]
                : static_cast<std::size_t>(directionIndex);
        const auto directionStart = std::chrono::steady_clock::now();
        if (progress_callback) {
            const SampleProgress progress {
                reportedDirectionIndex,
                completedCount.load(),
                directions.size(),
                directions[directionIndex].zenith_deg,
                normalizeAzimuthDeg(directions[directionIndex].azimuth_deg),
                0.0,
                std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count(),
                SampleProgress::Stage::direction_started,
                0,
                0,
                0.0,
                0,
                0,
                0,
                };
            std::lock_guard<std::mutex> lock(progressMutex);
            progress_callback(progress);
        }

        const DirectionEstimate estimate = estimateDirectionMoments(
            context,
            config,
            directions[directionIndex],
            sampleCount,
            reportedDirectionIndex,
            [&](SampleProgress::Stage stage,
                const DirectionEstimate::TimingSummary &timing,
                std::size_t stageCompleted,
                std::size_t stageTotal) {
                if (!progress_callback || stage == SampleProgress::Stage::direction_started) {
                    return;
                }
                const auto now = std::chrono::steady_clock::now();
                const SampleProgress progress {
                    reportedDirectionIndex,
                    completedCount.load(),
                    directions.size(),
                    directions[directionIndex].zenith_deg,
                    normalizeAzimuthDeg(directions[directionIndex].azimuth_deg),
                    std::chrono::duration<double>(now - directionStart).count(),
                    std::chrono::duration<double>(now - startTime).count(),
                    stage,
                    stageCompleted,
                    stageTotal,
                    timing.first_order_seconds,
                    timing.first_order_view_samples,
                    timing.first_order_steps,
                    timing.solar_disk_nodes,
                    timing.second_order_seconds,
                    timing.second_order_incoming_single_scatter_seconds,
                    timing.second_order_incoming_single_scatter_calls,
                    timing.second_order_nonzero_incoming_calls,
                    timing.second_order_view_samples,
                    timing.second_order_mu_phi_evaluations,
                    timing.second_order_view_steps,
                    timing.second_order_ray_steps,
                    timing.second_order_mu_nodes,
                    timing.second_order_phi_nodes,
                    timing.higher_order_seconds,
                    timing.spectral_band_count,
                    timing.mc_sample_count,
                };
                std::lock_guard<std::mutex> lock(progressMutex);
                progress_callback(progress);
            }
        );
        bins[directionIndex] = {
            directions[directionIndex].zenith_deg,
            normalizeAzimuthDeg(directions[directionIndex].azimuth_deg),
            totalMean(estimate),
            estimate.higher_variance,
            estimate.first_order,
            estimate.second_order.total,
            estimate.higher_order,
            estimate.higher_variance,
            estimate.second_order.rr,
            estimate.second_order.ar,
            estimate.second_order.ra,
            estimate.second_order.aa,
        };

        if (progress_callback) {
            const auto directionEnd = std::chrono::steady_clock::now();
            const SampleProgress progress {
                reportedDirectionIndex,
                completedCount.fetch_add(1) + 1,
                directions.size(),
                directions[directionIndex].zenith_deg,
                normalizeAzimuthDeg(directions[directionIndex].azimuth_deg),
                std::chrono::duration<double>(directionEnd - directionStart).count(),
                std::chrono::duration<double>(directionEnd - startTime).count(),
                SampleProgress::Stage::direction_complete,
                directions.size(),
                directions.size(),
                estimate.timing.first_order_seconds,
                estimate.timing.first_order_view_samples,
                estimate.timing.first_order_steps,
                estimate.timing.solar_disk_nodes,
                estimate.timing.second_order_seconds,
                estimate.timing.second_order_incoming_single_scatter_seconds,
                estimate.timing.second_order_incoming_single_scatter_calls,
                estimate.timing.second_order_nonzero_incoming_calls,
                estimate.timing.second_order_view_samples,
                estimate.timing.second_order_mu_phi_evaluations,
                estimate.timing.second_order_view_steps,
                estimate.timing.second_order_ray_steps,
                estimate.timing.second_order_mu_nodes,
                estimate.timing.second_order_phi_nodes,
                estimate.timing.higher_order_seconds,
                estimate.timing.spectral_band_count,
                estimate.timing.mc_sample_count,
            };
            std::lock_guard<std::mutex> lock(progressMutex);
            progress_callback(progress);
        }
    }
    return bins;
}

SkyBinResult solveSkyDirection(
    const SimulationConfig &config,
    const SkyDirection &direction,
    std::size_t direction_index,
    int higher_order_sample_target,
    DirectionCheckpointState *checkpoint_state,
    const std::function<void(const SampleProgress &)> &progress_callback
)
{
    const SimulationContext context = buildSimulationContext(config);
    const int sampleTarget =
        higher_order_sample_target > 0 ? higher_order_sample_target : config.monte_carlo.photons_per_bin;
    const auto startTime = std::chrono::steady_clock::now();
    const auto directionStart = startTime;

    const DirectionEstimate estimate = estimateDirectionMoments(
        context,
        config,
        direction,
        sampleTarget,
        direction_index,
        [&](SampleProgress::Stage stage,
            const DirectionEstimate::TimingSummary &timing,
            std::size_t stageCompleted,
            std::size_t stageTotal) {
            if (!progress_callback) {
                return;
            }
            const auto now = std::chrono::steady_clock::now();
            const SampleProgress progress {
                direction_index,
                0,
                1,
                direction.zenith_deg,
                normalizeAzimuthDeg(direction.azimuth_deg),
                std::chrono::duration<double>(now - directionStart).count(),
                std::chrono::duration<double>(now - startTime).count(),
                stage,
                stageCompleted,
                stageTotal,
                timing.first_order_seconds,
                timing.first_order_view_samples,
                timing.first_order_steps,
                timing.solar_disk_nodes,
                timing.second_order_seconds,
                timing.second_order_incoming_single_scatter_seconds,
                timing.second_order_incoming_single_scatter_calls,
                timing.second_order_nonzero_incoming_calls,
                timing.second_order_view_samples,
                timing.second_order_mu_phi_evaluations,
                timing.second_order_view_steps,
                timing.second_order_ray_steps,
                timing.second_order_mu_nodes,
                timing.second_order_phi_nodes,
                timing.higher_order_seconds,
                timing.spectral_band_count,
                timing.mc_sample_count,
            };
            progress_callback(progress);
        },
        checkpoint_state
    );

    if (progress_callback) {
        const auto directionEnd = std::chrono::steady_clock::now();
        const SampleProgress progress {
            direction_index,
            1,
            1,
            direction.zenith_deg,
            normalizeAzimuthDeg(direction.azimuth_deg),
            std::chrono::duration<double>(directionEnd - directionStart).count(),
            std::chrono::duration<double>(directionEnd - startTime).count(),
            SampleProgress::Stage::direction_complete,
            1,
            1,
            estimate.timing.first_order_seconds,
            estimate.timing.first_order_view_samples,
            estimate.timing.first_order_steps,
            estimate.timing.solar_disk_nodes,
            estimate.timing.second_order_seconds,
            estimate.timing.second_order_incoming_single_scatter_seconds,
            estimate.timing.second_order_incoming_single_scatter_calls,
            estimate.timing.second_order_nonzero_incoming_calls,
            estimate.timing.second_order_view_samples,
            estimate.timing.second_order_mu_phi_evaluations,
            estimate.timing.second_order_view_steps,
            estimate.timing.second_order_ray_steps,
            estimate.timing.second_order_mu_nodes,
            estimate.timing.second_order_phi_nodes,
            estimate.timing.higher_order_seconds,
            estimate.timing.spectral_band_count,
            estimate.timing.mc_sample_count,
        };
        progress_callback(progress);
    }

    return {
        direction.zenith_deg,
        normalizeAzimuthDeg(direction.azimuth_deg),
        totalMean(estimate),
        estimate.higher_variance,
        estimate.first_order,
        estimate.second_order.total,
        estimate.higher_order,
        estimate.higher_variance,
        estimate.second_order.rr,
        estimate.second_order.ar,
        estimate.second_order.ra,
        estimate.second_order.aa,
    };
}

SkyResult runMonteCarloSimulation(const SimulationConfig &config)
{
    const SolarGeometry sun = solarGeometryFromConfig(config);
    const double zenithStep = 90.0 / static_cast<double>(config.monte_carlo.zenith_bins);
    const double azimuthStep = 360.0 / static_cast<double>(config.monte_carlo.azimuth_bins);
    std::vector<SkyDirection> directions;
    directions.reserve(static_cast<std::size_t>(config.monte_carlo.zenith_bins * config.monte_carlo.azimuth_bins));

    for (int zenithIndex = 0; zenithIndex < config.monte_carlo.zenith_bins; ++zenithIndex) {
        const double viewZenithDeg = (zenithIndex + 0.5) * zenithStep;
        std::cout << "Evaluating zenith row " << (zenithIndex + 1)
                  << "/" << config.monte_carlo.zenith_bins << "\n";
        for (int azimuthIndex = 0; azimuthIndex < config.monte_carlo.azimuth_bins; ++azimuthIndex) {
            directions.push_back({
                viewZenithDeg,
                (azimuthIndex + 0.5) * azimuthStep,
            });
        }
    }

    const std::vector<SkyBinResult> bins = sampleSkyDirections(config, directions);

    SkyResult result;
    result.config = config;
    result.sun_zenith_deg = sun.zenith_deg;
    result.sun_azimuth_deg = sun.azimuth_deg;
    result.atmosphere_hash = hashFile(config.atmosphere.profile_csv);
    result.bins = bins;

    for (const SkyBinResult &bin : result.bins) {
        result.peak_intensity = std::max(result.peak_intensity, bin.mean.I);
    }

    const double peakMask = std::max(1.0e-12, result.peak_intensity * 0.01);
    for (const SkyBinResult &bin : result.bins) {
        if (bin.mean.I >= peakMask) {
            result.peak_dolp = std::max(result.peak_dolp, degreeOfLinearPolarization(bin.mean));
        }
    }
    result.hemispheric_flux_estimate = estimateHemisphericFlux(result);
    return result;
}

void writeSkyResult(const SkyResult &result)
{
    const auto jsonEscape = [](std::string value) {
        std::string escaped;
        escaped.reserve(value.size());
        for (char ch : value) {
            if (ch == '\\' || ch == '"') {
                escaped.push_back('\\');
            }
            escaped.push_back(ch);
        }
        return escaped;
    };

    const std::filesystem::path outputDir =
        std::filesystem::path(result.config.output.output_dir) / result.config.output.case_id;
    std::filesystem::create_directories(outputDir);

    const std::filesystem::path csvPath = outputDir / "sky_result.csv";
    std::ofstream csv(csvPath);
    csv << "zenith_deg,azimuth_deg,"
        << "I,Q,U,V,var_I,var_Q,var_U,var_V,dolp,aop_rad,"
        << "first_I,first_Q,first_U,first_V,"
        << "second_I,second_Q,second_U,second_V,"
        << "higher_I,higher_Q,higher_U,higher_V,"
        << "higher_var_I,higher_var_Q,higher_var_U,higher_var_V,"
        << "second_rr_I,second_rr_Q,second_rr_U,second_rr_V,"
        << "second_ar_I,second_ar_Q,second_ar_U,second_ar_V,"
        << "second_ra_I,second_ra_Q,second_ra_U,second_ra_V,"
        << "second_aa_I,second_aa_Q,second_aa_U,second_aa_V\n";
    csv << std::scientific << std::setprecision(10);
    for (const SkyBinResult &bin : result.bins) {
        csv << bin.zenith_deg << ","
            << bin.azimuth_deg << ","
            << bin.mean.I << ","
            << bin.mean.Q << ","
            << bin.mean.U << ","
            << bin.mean.V << ","
            << bin.variance.I << ","
            << bin.variance.Q << ","
            << bin.variance.U << ","
            << bin.variance.V << ","
            << degreeOfLinearPolarization(bin.mean) << ","
            << angleOfLinearPolarizationRad(bin.mean) << ","
            << bin.first_order.I << ","
            << bin.first_order.Q << ","
            << bin.first_order.U << ","
            << bin.first_order.V << ","
            << bin.second_order.I << ","
            << bin.second_order.Q << ","
            << bin.second_order.U << ","
            << bin.second_order.V << ","
            << bin.higher_order.I << ","
            << bin.higher_order.Q << ","
            << bin.higher_order.U << ","
            << bin.higher_order.V << ","
            << bin.higher_variance.I << ","
            << bin.higher_variance.Q << ","
            << bin.higher_variance.U << ","
            << bin.higher_variance.V << ","
            << bin.second_rr.I << ","
            << bin.second_rr.Q << ","
            << bin.second_rr.U << ","
            << bin.second_rr.V << ","
            << bin.second_ar.I << ","
            << bin.second_ar.Q << ","
            << bin.second_ar.U << ","
            << bin.second_ar.V << ","
            << bin.second_ra.I << ","
            << bin.second_ra.Q << ","
            << bin.second_ra.U << ","
            << bin.second_ra.V << ","
            << bin.second_aa.I << ","
            << bin.second_aa.Q << ","
            << bin.second_aa.U << ","
            << bin.second_aa.V << "\n";
    }

    const std::filesystem::path jsonPath = outputDir / "sky_result_metadata.json";
    std::ofstream json(jsonPath);
    json << "{\n";
    json << "  \"case_id\": \"" << jsonEscape(result.config.output.case_id) << "\",\n";
    json << "  \"config_path\": \"" << jsonEscape(result.config.config_path) << "\",\n";
    json << "  \"config_hash\": \"" << jsonEscape(result.config.config_hash) << "\",\n";
    json << "  \"atmosphere_hash\": \"" << jsonEscape(result.atmosphere_hash) << "\",\n";
    json << "  \"sun_zenith_deg\": " << result.sun_zenith_deg << ",\n";
    json << "  \"sun_azimuth_deg\": " << result.sun_azimuth_deg << ",\n";
    json << "  \"finite_solar_disk\": " << (result.config.solar.finite_solar_disk ? "true" : "false") << ",\n";
    json << "  \"strict_paper_mode\": " << (result.config.output.strict_paper_mode ? "true" : "false") << ",\n";
    json << "  \"solar_angular_radius_deg\": " << result.config.solar.solar_angular_radius_deg << ",\n";
    json << "  \"solar_disk_quadrature_nodes\": " << result.config.solar.solar_disk_quadrature_nodes << ",\n";
    json << "  \"zenith_bins\": " << result.config.monte_carlo.zenith_bins << ",\n";
    json << "  \"azimuth_bins\": " << result.config.monte_carlo.azimuth_bins << ",\n";
    json << "  \"photons_per_bin\": " << result.config.monte_carlo.photons_per_bin << ",\n";
    json << "  \"optical_depth_step_scale\": " << result.config.monte_carlo.optical_depth_step_scale << ",\n";
    json << "  \"min_optical_depth_steps\": " << result.config.monte_carlo.min_optical_depth_steps << ",\n";
    json << "  \"line_integral_step_scale\": " << result.config.monte_carlo.line_integral_step_scale << ",\n";
    json << "  \"min_line_integral_steps\": " << result.config.monte_carlo.min_line_integral_steps << ",\n";
    json << "  \"single_scatter_only\": " << (result.config.monte_carlo.single_scatter_only ? "true" : "false") << ",\n";
    json << "  \"deterministic_single_scatter\": " << (result.config.monte_carlo.deterministic_single_scatter ? "true" : "false") << ",\n";
    json << "  \"deterministic_second_scatter\": " << (result.config.monte_carlo.deterministic_second_scatter ? "true" : "false") << ",\n";
    json << "  \"second_scatter_view_steps\": " << result.config.monte_carlo.second_scatter_view_steps << ",\n";
    json << "  \"second_scatter_ray_steps\": " << result.config.monte_carlo.second_scatter_ray_steps << ",\n";
    json << "  \"second_scatter_mu_nodes\": " << result.config.monte_carlo.second_scatter_mu_nodes << ",\n";
    json << "  \"second_scatter_phi_nodes\": " << result.config.monte_carlo.second_scatter_phi_nodes << ",\n";
    json << "  \"twilight_second_scatter_adaptive\": " << (result.config.monte_carlo.twilight_second_scatter_adaptive ? "true" : "false") << ",\n";
    json << "  \"twilight_second_scatter_zenith_threshold_deg\": " << result.config.monte_carlo.twilight_second_scatter_zenith_threshold_deg << ",\n";
    json << "  \"twilight_second_scatter_min_view_steps\": " << result.config.monte_carlo.twilight_second_scatter_min_view_steps << ",\n";
    json << "  \"twilight_second_scatter_min_ray_steps\": " << result.config.monte_carlo.twilight_second_scatter_min_ray_steps << ",\n";
    json << "  \"twilight_second_scatter_min_mu_nodes\": " << result.config.monte_carlo.twilight_second_scatter_min_mu_nodes << ",\n";
    json << "  \"twilight_second_scatter_min_phi_nodes\": " << result.config.monte_carlo.twilight_second_scatter_min_phi_nodes << ",\n";
    json << "  \"source_guided_first_scatter\": " << (result.config.monte_carlo.source_guided_first_scatter ? "true" : "false") << ",\n";
    json << "  \"source_guided_first_scatter_branches\": " << result.config.monte_carlo.source_guided_first_scatter_branches << ",\n";
    json << "  \"source_guided_phase_fraction\": " << result.config.monte_carlo.source_guided_phase_fraction << ",\n";
    json << "  \"source_guided_cone_half_angle_deg\": " << result.config.monte_carlo.source_guided_cone_half_angle_deg << ",\n";
    json << "  \"rayleigh_source_guided_phase_fraction\": " << result.config.monte_carlo.rayleigh_source_guided_phase_fraction << ",\n";
    json << "  \"rayleigh_source_guided_cone_half_angle_deg\": " << result.config.monte_carlo.rayleigh_source_guided_cone_half_angle_deg << ",\n";
    json << "  \"rayleigh_source_guided_first_scatter_branches\": " << result.config.monte_carlo.rayleigh_source_guided_first_scatter_branches << ",\n";
    json << "  \"rayleigh_polarization_guided_fraction\": " << result.config.monte_carlo.rayleigh_polarization_guided_fraction << ",\n";
    json << "  \"rayleigh_polarization_guided_mu_half_width\": " << result.config.monte_carlo.rayleigh_polarization_guided_mu_half_width << ",\n";
    json << "  \"rayleigh_polarization_guided_branches\": " << result.config.monte_carlo.rayleigh_polarization_guided_branches << ",\n";
    json << "  \"twilight_limb_guiding\": " << (result.config.monte_carlo.twilight_limb_guiding ? "true" : "false") << ",\n";
    json << "  \"twilight_limb_phase_fraction\": " << result.config.monte_carlo.twilight_limb_phase_fraction << ",\n";
    json << "  \"twilight_limb_cone_half_angle_deg\": " << result.config.monte_carlo.twilight_limb_cone_half_angle_deg << ",\n";
    json << "  \"twilight_limb_elevation_deg\": " << result.config.monte_carlo.twilight_limb_elevation_deg << ",\n";
    json << "  \"twilight_higher_order_guiding\": " << (result.config.monte_carlo.twilight_higher_order_guiding ? "true" : "false") << ",\n";
    json << "  \"twilight_higher_order_branches\": " << result.config.monte_carlo.twilight_higher_order_branches << ",\n";
    json << "  \"twilight_higher_order_phase_fraction\": " << result.config.monte_carlo.twilight_higher_order_phase_fraction << ",\n";
    json << "  \"twilight_higher_order_tangent_fraction\": " << result.config.monte_carlo.twilight_higher_order_tangent_fraction << ",\n";
    json << "  \"twilight_higher_order_horizon_fraction\": " << result.config.monte_carlo.twilight_higher_order_horizon_fraction << ",\n";
    json << "  \"twilight_higher_order_tangent_cone_half_angle_deg\": " << result.config.monte_carlo.twilight_higher_order_tangent_cone_half_angle_deg << ",\n";
    json << "  \"twilight_higher_order_horizon_cone_half_angle_deg\": " << result.config.monte_carlo.twilight_higher_order_horizon_cone_half_angle_deg << ",\n";
    json << "  \"twilight_higher_order_horizon_elevation_deg\": " << result.config.monte_carlo.twilight_higher_order_horizon_elevation_deg << ",\n";
    json << "  \"twilight_order_depolarization\": " << (result.config.monte_carlo.twilight_order_depolarization ? "true" : "false") << ",\n";
    json << "  \"twilight_second_order_polarization_scale\": " << result.config.monte_carlo.twilight_second_order_polarization_scale << ",\n";
    json << "  \"twilight_higher_order_polarization_scale\": " << result.config.monte_carlo.twilight_higher_order_polarization_scale << ",\n";
    json << "  \"twilight_second_order_intensity_boost\": " << result.config.monte_carlo.twilight_second_order_intensity_boost << ",\n";
    json << "  \"twilight_higher_order_intensity_boost\": " << result.config.monte_carlo.twilight_higher_order_intensity_boost << ",\n";
    json << "  \"twilight_second_order_boost_zenith_cutoff_deg\": " << result.config.monte_carlo.twilight_second_order_boost_zenith_cutoff_deg << ",\n";
    json << "  \"random_seed\": " << result.config.monte_carlo.random_seed << ",\n";
    json << "  \"peak_intensity\": " << result.peak_intensity << ",\n";
    json << "  \"peak_dolp\": " << result.peak_dolp << ",\n";
    json << "  \"hemispheric_flux_estimate\": " << result.hemispheric_flux_estimate << ",\n";
    json << "  \"solar_spectrum_csv\": \"" << jsonEscape(result.config.spectral.solar_spectrum_csv) << "\",\n";
    json << "  \"instrument_response_csv\": \"" << jsonEscape(result.config.spectral.instrument_response_csv) << "\",\n";
    json << "  \"rayleigh_cross_section_csv\": \"" << jsonEscape(result.config.spectral.rayleigh_cross_section_csv) << "\",\n";
    json << "  \"ozone_cross_section_csv\": \"" << jsonEscape(result.config.spectral.ozone_cross_section_csv) << "\",\n";
    json << "  \"o2_cross_section_csv\": \"" << jsonEscape(result.config.spectral.o2_cross_section_csv) << "\",\n";
    json << "  \"o4_cross_section_csv\": \"" << jsonEscape(result.config.spectral.o4_cross_section_csv) << "\",\n";
    json << "  \"h2o_cross_section_csv\": \"" << jsonEscape(result.config.spectral.h2o_cross_section_csv) << "\",\n";
    json << "  \"no2_cross_section_csv\": \"" << jsonEscape(result.config.spectral.no2_cross_section_csv) << "\",\n";
    json << "  \"aerosol_phase_matrix_csv\": \"" << jsonEscape(result.config.spectral.aerosol_phase_matrix_csv) << "\",\n";
    json << "  \"surface_model\": \"" << jsonEscape(result.config.surface.surface_model) << "\",\n";
    json << "  \"surface_parameter_csv\": \"" << jsonEscape(result.config.surface.surface_parameter_csv) << "\",\n";
    json << "  \"paper_case_provenance_json\": \"" << jsonEscape(result.config.output.paper_case_provenance_json) << "\"\n";
    json << "}\n";

    std::cout << "Sky result written to " << csvPath.string() << "\n";
    std::cout << "Metadata written to " << jsonPath.string() << "\n";
}
