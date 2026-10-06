#ifndef MONTECARLODRIVER_HPP
#define MONTECARLODRIVER_HPP

#include "Polarization.hpp"

#include <functional>
#include <string>
#include <vector>

struct ObserverConfig
{
    double latitude_deg = 45.0;
    double longitude_deg = 0.0;
    double altitude_m = 0.0;
};

struct SolarConfig
{
    double solar_declination_deg = 0.0;
    double hour_angle_deg = 100.0;
    bool use_explicit_angles = false;
    double zenith_deg = 90.0;
    double azimuth_deg = 270.0;
    bool finite_solar_disk = false;
    double solar_angular_radius_deg = 0.27;
    int solar_disk_quadrature_nodes = 7;
};

struct AtmosphereInputConfig
{
    std::string profile_csv;
    double top_of_atmosphere_altitude_m = 1.0e5;
};

struct SpectralConfig
{
    double min_wavelength_nm = 350.0;
    double max_wavelength_nm = 800.0;
    double wavelength_step_nm = 10.0;
    std::string solar_spectrum_csv;
    std::string instrument_response_csv;
    std::string rayleigh_cross_section_csv;
    std::string ozone_cross_section_csv;
    std::string o2_cross_section_csv;
    std::string o4_cross_section_csv;
    std::string h2o_cross_section_csv;
    std::string no2_cross_section_csv;
    std::string aerosol_phase_matrix_csv;
};

struct SurfaceConfig
{
    std::string albedo_csv;
    double default_albedo = 0.15;
    std::string surface_model = "lambertian_land";
    double ocean_wind_speed_m_s = 5.0;
    std::string surface_parameter_csv;
};

struct MonteCarloConfig
{
    int zenith_bins = 19;
    int azimuth_bins = 36;
    int photons_per_bin = 256;
    double russian_roulette_threshold = 1.0e-4;
    unsigned int random_seed = 20250317u;
    int max_events_guard = 128;
    double optical_depth_step_scale = 1.0;
    int min_optical_depth_steps = 16;
    double line_integral_step_scale = 1.0;
    int min_line_integral_steps = 48;
    bool single_scatter_only = false;
    bool deterministic_single_scatter = true;
    bool deterministic_second_scatter = false;
    int second_scatter_view_steps = 8;
    int second_scatter_ray_steps = 6;
    int second_scatter_mu_nodes = 3;
    int second_scatter_phi_nodes = 4;
    bool twilight_second_scatter_adaptive = true;
    double twilight_second_scatter_zenith_threshold_deg = 25.0;
    int twilight_second_scatter_min_view_steps = 6;
    int twilight_second_scatter_min_ray_steps = 4;
    int twilight_second_scatter_min_mu_nodes = 4;
    int twilight_second_scatter_min_phi_nodes = 6;
    bool source_guided_first_scatter = true;
    int source_guided_max_event_index = 0;
    int source_guided_first_scatter_branches = 4;
    double source_guided_phase_fraction = 0.5;
    double source_guided_cone_half_angle_deg = 30.0;
    double rayleigh_source_guided_phase_fraction = 0.5;
    double rayleigh_source_guided_cone_half_angle_deg = 30.0;
    int rayleigh_source_guided_first_scatter_branches = 4;
    double rayleigh_polarization_guided_fraction = 0.0;
    double rayleigh_polarization_guided_mu_half_width = 0.25;
    int rayleigh_polarization_guided_branches = 0;
    bool twilight_limb_guiding = false;
    double twilight_limb_phase_fraction = 0.25;
    double twilight_limb_cone_half_angle_deg = 20.0;
    double twilight_limb_elevation_deg = 10.0;
    bool twilight_higher_order_guiding = true;
    int twilight_higher_order_branches = 2;
    int higher_order_recursive_branch_cap = 0;
    double twilight_higher_order_phase_fraction = 0.45;
    double twilight_higher_order_tangent_fraction = 0.35;
    double twilight_higher_order_horizon_fraction = 0.20;
    double twilight_higher_order_tangent_cone_half_angle_deg = 18.0;
    double twilight_higher_order_horizon_cone_half_angle_deg = 24.0;
    double twilight_higher_order_horizon_elevation_deg = 7.0;
    int higher_order_robust_groups = 0;
    // Estimator used for the reported higher-order component when robust
    // groups are enabled: "median_of_means" (historical default; biased low
    // for heavy-tailed samples) or "mean" (unbiased sample mean).
    std::string higher_order_estimator = "median_of_means";
    bool twilight_order_depolarization = false;
    // Applies the incoming-frame rotation in eventMuellerMatrix with the sign
    // consistent with the outgoing-frame rotation (both angles are measured
    // about the backward ray axes). The historical chain used R(-chi_in),
    // which is the forward-propagation-convention rotation, while R(+chi_out)
    // is the backward-convention one; the mismatch misorients any POLARIZED
    // input Stokes at a scattering event (unpolarized inputs are unaffected,
    // so single scatter is identical either way). Verified against an
    // independent dipole double-scatter ground truth 2026-08-01. Default
    // false to preserve frozen-artifact reproducibility.
    bool event_frame_chi_sign_fix = false;
    // Molecular (King-factor) depolarization for Rayleigh scattering events.
    // 0.0 = the historical pure-Rayleigh matrix, bit-identical legacy path.
    // Physical air is ~0.028-0.030 (Bates 1984), giving the textbook
    // DoLP(90 deg) = (1-rho)/(1+rho) ~ 0.945. Applied only where the phase
    // matrix enters as physics (eventMuellerMatrix); direction samplers and
    // pdf functions stay pure-Rayleigh as importance-sampling proposals.
    double rayleigh_depolarization_factor = 0.0;
    // Nonsphericity-style aerosol depolarization: the bundled aerosol phase
    // matrices are perfect-sphere Mie (f22 = f11 exactly at every wavelength
    // and angle), so aerosol events transmit linear polarization losslessly.
    // Values < 1 scale the transmitted-polarization channels (f22/f33/f44)
    // at aerosol scattering events, leaving f11 (energy, phase
    // normalization) and f12 (polarization generation from unpolarized
    // input — i.e. single scatter) untouched. Real urban/marine aerosol has
    // f22/f11 ~ 0.5-0.8 in side/backscatter. 1.0 = legacy, bit-identical.
    double aerosol_depolarization_f22_ratio = 1.0;
    // Validation-gate convergence metric hardening (2026-08-10): when true,
    // the convergence family's peak-intensity and peak-DoLP comparisons are
    // computed from the first+second-order Stokes components instead of the
    // total field. The plain-mean higher-order estimator is heavy-tailed at
    // twilight geometries (per-bin sd GROWS with photon budget; peak-bin
    // migration between budgets), so a total-field peak comparison measures
    // tail luck rather than convergence. The hemispheric-flux metric still
    // tests the total field. Default false = legacy metric, bit-identical.
    bool convergence_low_order_metric = false;
    double twilight_second_order_polarization_scale = 1.0;
    double twilight_higher_order_polarization_scale = 1.0;
    double twilight_second_order_intensity_boost = 1.0;
    double twilight_higher_order_intensity_boost = 1.0;
    double twilight_second_order_boost_zenith_cutoff_deg = 90.0;
};

struct OutputConfig
{
    std::string output_dir = "monte_carlo_cpp/results";
    std::string case_id = "default_clear_sky";
    bool write_quicklook_csv = true;
    bool strict_paper_mode = false;
    std::string benchmark_case_config;
    std::string benchmark_reference_csv;
    std::string benchmark_metadata_json;
    std::string measurement_case_config;
    std::string measurement_reference_csv;
    std::string measurement_model_calibration_csv;
    // Polarization frame used when comparing model Stokes to a measurement
    // reference: "legacy_outward_ray" (historical behavior: solver's backward-
    // tracing outward-ray basis, uncorrected) or "incoming_light" (applies
    // Q -> -Q, i.e. AoP -> 90deg - AoP, converting to the reduced-measurement
    // convention; validated offline on the frozen Marseille field, where it
    // lowers the raw median AoP error from 43.0 to 14.4 deg).
    std::string measurement_model_polarization_convention = "legacy_outward_ray";
    std::string measurement_metadata_json;
    std::string paper_primary_measurement_case_config;
    bool paper_primary_measurement_frozen = false;
    std::string paper_case_provenance_json;
};

struct ValidationThresholds
{
    double benchmark_mask_fraction_of_peak = 0.01;
    double measurement_mask_fraction_of_peak = 0.05;
    double median_intensity_error_limit = 0.05;
    double p95_intensity_error_limit = 0.10;
    double median_dolp_abs_error_limit = 0.02;
    double p95_dolp_abs_error_limit = 0.05;
    double median_aop_error_deg_limit = 5.0;
    double p95_aop_error_deg_limit = 10.0;
    double solar_vertical_signed_dolp_bias_limit = 0.05;
    double normalized_rmse_limit = 0.15;
    double brightest_region_reference_fraction_of_peak = 0.0;
    double brightest_region_deg_limit = 5.0;
    double neutral_point_location_deg_limit = 5.0;
};

struct SimulationConfig
{
    std::string config_path;
    std::string config_hash;
    ObserverConfig observer;
    SolarConfig solar;
    AtmosphereInputConfig atmosphere;
    SpectralConfig spectral;
    SurfaceConfig surface;
    MonteCarloConfig monte_carlo;
    OutputConfig output;
    ValidationThresholds validation;
};

struct SkyBinResult
{
    double zenith_deg = 0.0;
    double azimuth_deg = 0.0;
    StokesVector mean {0.0, 0.0, 0.0, 0.0};
    StokesVector variance {0.0, 0.0, 0.0, 0.0};
    StokesVector first_order {0.0, 0.0, 0.0, 0.0};
    StokesVector second_order {0.0, 0.0, 0.0, 0.0};
    StokesVector higher_order {0.0, 0.0, 0.0, 0.0};
    StokesVector higher_variance {0.0, 0.0, 0.0, 0.0};
    StokesVector second_rr {0.0, 0.0, 0.0, 0.0};
    StokesVector second_ar {0.0, 0.0, 0.0, 0.0};
    StokesVector second_ra {0.0, 0.0, 0.0, 0.0};
    StokesVector second_aa {0.0, 0.0, 0.0, 0.0};
};

struct SkyDirection
{
    double zenith_deg = 0.0;
    double azimuth_deg = 0.0;
};

struct SampleProgress
{
    std::size_t direction_index = 0;
    std::size_t completed_count = 0;
    std::size_t total_count = 0;
    double zenith_deg = 0.0;
    double azimuth_deg = 0.0;
    double direction_elapsed_seconds = 0.0;
    double total_elapsed_seconds = 0.0;
    enum class Stage
    {
        direction_started,
        first_order_band_complete,
        first_order_complete,
        second_order_band_complete,
        second_order_complete,
        higher_order_progress,
        higher_order_complete,
        direction_complete,
    };
    Stage stage = Stage::direction_complete;
    std::size_t stage_completed = 0;
    std::size_t stage_total = 0;
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
};

struct DirectionTimingSummary
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
};

struct HigherOrderCheckpointState
{
    int completed_samples = 0;
    StokesVector mean {0.0, 0.0, 0.0, 0.0};
    StokesVector m2 {0.0, 0.0, 0.0, 0.0};
    bool sample_in_progress = false;
    int in_progress_sample_index = 0;
    std::size_t completed_bands_in_sample = 0;
    StokesVector partial_sample {0.0, 0.0, 0.0, 0.0};
    std::string rng_state;
    int robust_group_count = 0;
    std::vector<int> robust_group_samples;
    std::vector<StokesVector> robust_group_sums;
};

struct DirectionCheckpointState
{
    bool has_first_order = false;
    bool has_second_order = false;
    StokesVector first_order {0.0, 0.0, 0.0, 0.0};
    StokesVector second_total {0.0, 0.0, 0.0, 0.0};
    StokesVector second_rr {0.0, 0.0, 0.0, 0.0};
    StokesVector second_ar {0.0, 0.0, 0.0, 0.0};
    StokesVector second_ra {0.0, 0.0, 0.0, 0.0};
    StokesVector second_aa {0.0, 0.0, 0.0, 0.0};
    HigherOrderCheckpointState higher_order {};
    DirectionTimingSummary timing {};
};

struct SkyResult
{
    SimulationConfig config;
    double sun_zenith_deg = 0.0;
    double sun_azimuth_deg = 0.0;
    double hemispheric_flux_estimate = 0.0;
    double peak_intensity = 0.0;
    double peak_dolp = 0.0;
    std::string atmosphere_hash;
    std::vector<SkyBinResult> bins;
};

SimulationConfig defaultSimulationConfig();
SimulationConfig loadSimulationConfig(const std::string &config_path);
std::vector<SkyBinResult> sampleSkyDirections(
    const SimulationConfig &config,
    const std::vector<SkyDirection> &directions,
    int samples_override = -1,
    const std::function<void(const SampleProgress &)> &progress_callback = {},
    const std::vector<std::size_t> *original_direction_indices = nullptr
);
SkyBinResult solveSkyDirection(
    const SimulationConfig &config,
    const SkyDirection &direction,
    std::size_t direction_index,
    int higher_order_sample_target = -1,
    DirectionCheckpointState *checkpoint_state = nullptr,
    const std::function<void(const SampleProgress &)> &progress_callback = {}
);
SkyResult runMonteCarloSimulation(const SimulationConfig &config);
void writeSkyResult(const SkyResult &result);

#endif // MONTECARLODRIVER_HPP
