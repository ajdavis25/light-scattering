# this file makes the scripts directory a package
"""
light-scattering.scripts package
=======================================
this package provides scripts for analyzing and plotting time constants.
"""


# import from src package    
from spherical import(
    generate_sky_angles,
    calculate_intensity,
    generate_unpolarized_fisheye_data,
    generate_polarized_fisheye_data,
    plot_fisheye,
    plot_unpolarized_fisheye,
    plot_polarized_fisheye,
    intensity_at_ground_spherical,
    generate_intensity_spherical,
    plot_intensity_vs_theta_spherical,
    calculate_scattering_angle_for_polarization,
    polarized_intensity_at_ground_spherical,
    calculate_polarization_stokes,
    photon_unit_vector_spherical,
    sun_position_vector,
    scattering_angle_spherical,
    rayleigh_phase_function,
    mie_phase_function,
    calculate_sun_position,
    calculate_scattering_angle,
    save_plot,
    air_mass,
    refraction_correction,
    create_twilight_colormap,
    gaussian_kernel,
    apply_gaussian_smoothing
)


"""
# import test utilities (if necessary)
from tests import (
    # Add imports from test modules if needed, like test utilities or helper functions
)
"""


# define the __all__ to control what's exposed when importing from the root package
__all__ = [
    'generate_sky_angles',
    'calculate_intensity',
    'generate_unpolarized_fisheye_data',
    'generate_polarized_fisheye_data',
    'plot_fisheye',
    'plot_unpolarized_fisheye',
    'plot_polarized_fisheye',
    'intensity_at_ground_spherical',
    'generate_intensity_spherical',
    'plot_intensity_vs_theta_spherical',
    'calculate_scattering_angle_for_polarization',
    'polarized_intensity_at_ground_spherical',
    'calculate_polarization_stokes',
    'photon_unit_vector_spherical',
    'sun_position_vector',
    'scattering_angle_spherical',
    'rayleigh_phase_function',
    'mie_phase_function',
    'calculate_sun_position',
    'calculate_scattering_angle',
    'save_plot',
    'air_mass',
    'refraction_correction',
    'create_twilight_colormap',
    'gaussian_kernel',
    'apply_gaussian_smoothing'
]
