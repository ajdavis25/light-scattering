# this file makes the scripts directory a package
"""
light-scattering.scripts package
=======================================
this package provides scripts for analyzing and plotting time constants.
"""


# import from src package
from src import(
    setup_logging,
    rayleigh_phase_function,
    scattering_angle,
    photon_unit_vector,
    scattered_intensity_reflected,
    scattered_intensity_transmitted,
    intensity_at_ground,
    intensity_at_top,
    scattered_intensity_at_tau,
    intensity_at_tau,
    direction_vector,
    dot_product,
    normalize,
    rayleigh_matrix,
    rotation_matrix,
    rotation_angles,
    generate_intensity,
    plot_intensity_vs_theta,
    calculate_polarization,
    generate_polarization_data,
    plot_polarization_vs_theta,
    calculate_quantity,
    make_contour_plot,
    format_sig_figs
)


"""
# import test utilities (if necessary)
from tests import (
    # Add imports from test modules if needed, like test utilities or helper functions
)
"""


# define the __all__ to control what's exposed when importing from the root package
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
