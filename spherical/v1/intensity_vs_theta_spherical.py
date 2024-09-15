""" single scattering of intensity vs. theta using a spherical model of the earth and a sun that isn't an infinitely small point """

import math, numpy as np
from typing import List, Tuple
from spherical.v1.signal_processing import (
    apply_gaussian_smoothing
)
from spherical.v1.atmospheric_model import (
    air_mass,
    refraction_correction,
    rayleigh_phase_function
)
from spherical_model import (
    photon_unit_vector_spherical,
    sun_position_vector
)


def intensity_at_ground_spherical(
        latitude: float,
        longitude: float,
        solar_declination: float,
        hour_angle: float,
        tau_max: float,
        num_layers: int = 100,
        scattering_angle: float = None # scattering angle parameter
) -> float:
    """
    calculate the intensity at the ground level based on the spherical earth model and atmospheric scattering effects
    by integrating across multiple atmospheric layers, considering the scattering angle

    args:
        latitude: the observer's latitude
        longitude: the observer's longitude
        solar_declination: the sun's declination angle
        hour_angle: the sun's hour angle
        tau_max: maximum optical depth
        num_layers: number of layers in the atmosphere for integration
        scattering_angle: the angle between the sun and the observer's location for scattering

    returns:
        the scattered intensity at the ground level
    """

    # calculate the earth's surface normal vector at the observer's location
    earth_surface_vector = photon_unit_vector_spherical(latitude, longitude)

    # calculate the sun's position vector
    sun_vector = sun_position_vector(solar_declination, hour_angle)

    # compute the cosine of the zenith angle (angle between the observer and the sun)
    cos_zenith_angle = np.dot(earth_surface_vector, sun_vector) / (np.linalg.norm(earth_surface_vector) * np.linalg.norm(sun_vector))
    cos_zenith_angle = np.clip(cos_zenith_angle, -1.0, 1.0) # avoid precision errors outside [-1, 1]

    # calculate the actual zenith angle (in degrees)
    zenith_angle = math.degrees(math.acos(cos_zenith_angle))

    # apply refraction correction to the zenith angle near the horizon
    zenith_angle -= refraction_correction(zenith_angle)

    # air mass factor to account for increased optical depth near the horizon
    airmass_factor = air_mass(zenith_angle)

    # integrate the scattering effects over multiple atmospheric layers
    intensity = 0.0
    delta_tau = tau_max / num_layers # optical depth per layer

    for i in range(num_layers):
        # calculate the optical depth at the current layer
        tau_layer = (i + 0.5) * delta_tau

        # avoid division by zero when zenith angle is near the horizon (mu near 0)
        mu = math.cos(math.radians(zenith_angle))
        if mu < 1e-5: # small threshold to prevent division by very small numbers
            continue

        # incorporate the scattering angle into the phase function (rayleigh scattering)
        if scattering_angle is not None:
            phase_function_value = rayleigh_phase_function(scattering_angle)
        else:
            phase_function_value = rayleigh_phase_function(zenith_angle) # fallback to zenith angle if scattering angle is not provided

        # calculate the intensity contribution from this atmospheric layer
        layer_intensity = phase_function_value * math.exp(-tau_layer / (mu * airmass_factor)) * delta_tau

        # accumulate the intensity from this layer
        intensity += layer_intensity

    return intensity


def generate_intensity_spherical(
        latitude: float,
        longitude: float,
        solar_declination: float,
        tau_atm: float,
        num_layers: int
) -> Tuple[List[float], List[float]]:
    """
    generate the intensity of light at ground level for a range of theta angles, given a latitude and longitude

    args:
        latitude: observer's latitude
        longitude: observer's longitude
        solar_declination: sun's declination angle
        tau_atm: atmospheric tau value
        num_layers: number of layers in the atmosphere for integration

    returns:
        a tuple of lists containing theta observation angles corresponding intensity values
    """
    theta_obs = []
    intensities = []

    # simulate for hour angles from -180 to 180 degrees (i.e., a full day)
    for hour_angle in range(-180, 180, 10):
        theta_obs.append(hour_angle)
        intensities.append(intensity_at_ground_spherical(latitude, longitude, solar_declination, hour_angle, tau_atm, num_layers))

    # apply a gaussian smoothing filter to reduce sharp transitions
    smoothed_intensities = apply_gaussian_smoothing(intensities, sigma=2)

    return theta_obs, smoothed_intensities
