# this file makes the scripts directory a package
"""
light-scattering.scripts package
=======================================
this package provides scripts for analyzing and plotting time constants.
"""


# initialize logging settings for the entire package
from logging_setup import setup_logging


# imports for easy access
from plane_parallel import (
    rayleigh_phase_function,
    scattering_angle,
    photon_unit_vector,
    scattered_intensity_reflected,
    scattered_intensity_transmitted,
    intensity_at_ground,
    intensity_at_top,
    scattered_intensity_at_tau,
    intensity_at_tau
)


# import the spherical model functions
from spherical_model import (
    photon_unit_vector_spherical,
    sun_position_vector,
    scattering_angle_spherical
)


# import vector calculation utilities
from vector_operations import (
    direction_vector,
    dot_product,
    normalize,
    rayleigh_matrix,
    rotation_matrix,
    rotation_angles
)


# import plotting functions
from intensity_vs_theta import (
    generate_intensity,
    plot_intensity_vs_theta
)
from polarization_vs_theta import (
    calculate_polarization,
    generate_polarization_data,
    plot_polarization_vs_theta
)
from single_scatter_plot import (
    calculate_quantity,
    make_contour_plot
)
from intensity_vs_theta_spherical import (
    air_mass,
    refraction_correction,
    generate_intensity_spherical,
    plot_intensity_vs_theta_spherical
)
from single_scatter_plot_spherical import (
    calculate_sun_position,
    generate_fisheye_data,
    create_twilight_colormap,
    plot_fisheye_twilight_sky
)
    

# import utils if there are any general utility functions
from utils import format_sig_figs


# expose certain names to package level for easier access
__all__ = [
    'setup_logging',
    'rayleigh_phase_function',
    'scattering_angle',
    'photon_unit_vector',
    'scattered_intensity_reflected',
    'scattered_intensity_transmitted',
    'intensity_at_ground',
    'intensity_at_top',
    'scattered_intensity_at_tau',
    'intensity_at_tau',
    'photon_unit_vector_spherical',
    'sun_position_vector',
    'scattering_angle_spherical',
    'direction_vector',
    'dot_product',
    'normalize',
    'rayleigh_matrix',
    'rotation_matrix',
    'rotation_angles',
    'generate_intensity',
    'plot_intensity_vs_theta',
    'calculate_polarization',
    'generate_polarization_data',
    'plot_polarization_vs_theta',
    'calculate_quantity',
    'make_contour_plot',
    'air_mass',
    'refraction_correction',
    'generate_intensity_spherical',
    'plot_intensity_vs_theta_spherical',
    'calculate_sun_position',
    'generate_fisheye_data',
    'create_twilight_colormap',
    'plot_fisheye_twilight_sky',
    'format_sig_figs'
]
