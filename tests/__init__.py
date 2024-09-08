# this file makes the scripts directory a package
"""
light-scattering.scripts package
=======================================
this package provides scripts for analyzing and plotting time constants.
"""


# initialize logging settings for the entire package
from src.logging_setup import setup_logging


# imports for easy access
from src.plane_parallel import (
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


# import vector calculation utilities
from src.vector_operations import (
    direction_vector,
    dot_product,
    normalize,
    rayleigh_matrix,
    rotation_matrix,
    rotation_angles
)


# import plotting functions
from src.intensity_vs_theta import (
    generate_intensity,
    plot_intensity_vs_theta
)
from src.polarization_vs_theta import (
    calculate_polarization,
    generate_polarization_data,
    plot_polarization_vs_theta
)
from src.single_scatter_plot import (
    calculate_quantity,
    make_contour_plot
)


# import utils if there are any general utility functions
from src.utils import format_sig_figs


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
    'format_sig_figs'
]
