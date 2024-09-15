# this file makes the scripts directory a package
"""
light-scattering.scripts package
=======================================
this package provides scripts for analyzing and plotting time constants.
"""


# import atmospheric model functions
from spherical.v1.atmospheric_model import (
    air_mass,
    refraction_correction,
    rayleigh_phase_function,
    mie_phase_function
)


# import energy loss functions
from spherical.v1.energy import (
    calculate_energy_loss
)


#import fisheye projection functions
from spherical.v1.fisheye_projection import(
    generate_sky_angles,
    calculate_intensity,
    generate_unpolarized_fisheye_data,
    generate_polarized_fisheye_data
)


# import intensity vs. theta functions
from spherical.v1.intensity_vs_theta_spherical import (
    intensity_at_ground_spherical,
    generate_intensity_spherical
)


# import plot utility functions
from spherical.v1.plot_utils import (
    plot_intensity_vs_theta_spherical,
    plot_fisheye,
    plot_unpolarized_fisheye,
    plot_polarized_fisheye
)


# import single scatter polarization functions
from spherical.v1.polarization_model import (
    polarized_intensity,
    calculate_polarization_stokes
)


# import signal processing functions
from spherical.v1.signal_processing import (
    gaussian_kernel,
    apply_gaussian_smoothing
)


# import the spherical model functions
from spherical_model import (
    photon_unit_vector_spherical,
    sun_position_vector,
    scattering_angle_spherical,
    calculate_sun_position,
    calculate_scattering_angle
)


# import utility functions
from utils import (
    save_plot,
    create_twilight_colormap
)


# expose certain names to package level for easier access
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
    'polarized_intensity',
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
    'apply_gaussian_smoothing',
    'calculate_energy_loss'
]
