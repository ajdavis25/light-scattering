# This file makes the scripts directory a package
"""
light-scattering.scripts package
=======================================
This package provides scripts for analyzing and plotting time constants.
"""


# Initialize logging settings for the entire package
from .logging_setup import setup_logging
setup_logging()


# Imports for easy access
from .plane_parallel import (
    rayleigh_phase_function,
    scattering_angle,
    photon_unit_vector,
    scattered_intensity_reflected,
    scattered_intensity_transmitted,
    intensity_at_ground,
    intensity_at_top,
    scattered_intensity_at_tau,
    intensity_at_tau,
)


# Import vector calculation utilities
from .vector_operations import (
    direction_vector,
    dot_product,
    normalize,
    rayleigh_matrix,
    rotation_matrix,
    rotation_angles
)


# Import plotting functions
from intensity_vs_theta import plot_intensity_vs_theta, generate_intensity
from polarization_vs_theta import plot_polarization_vs_theta


# Import utils if there are any general utility functions
from .utils import format_sig_figs


# Expose certain names to package level for easier access
__all__ = [
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
    'plot_polarization_vs_theta',
    'format_sig_figs'
]
