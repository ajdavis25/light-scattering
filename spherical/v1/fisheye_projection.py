import numpy as np
from typing import Tuple
from spherical_model import (
    calculate_sun_position,
    calculate_scattering_angle
)
from spherical.v1.polarization_model import (
    calculate_polarization_stokes
)
from spherical.v1.intensity_vs_theta_spherical import (
    intensity_at_ground_spherical
)
from spherical.v1.atmospheric_model import (
    rayleigh_phase_function
)


def generate_sky_angles(num_zenith_points: int, num_azimuth_points: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    generates arrays of zenith and azimuth angles for a specified number of points

    args:
        num_zenith_points: the number of points to generate along the zenith angle axis (0 to 90 degrees)
        num_azimuth_points: the number of points to generate along the azimuth angle axis (0 to 360 degrees)

    returns:
        a tuple containing two NumPy arrays:
        - zenith_angles: an array of zenith angles in degrees
        - azimuth_angles: an array of azimuth angles in degrees
    """
    zenith_angles = np.linspace(0, 90, num_zenith_points)
    azimuth_angles = np.linspace(0, 2 * np.pi, num_azimuth_points)

    return zenith_angles, azimuth_angles


def calculate_intensity(
        zenith_angles: np.ndarray,
        azimuth_angles: np.ndarray,
        sun_zenith: float,
        sun_azimuth: float,
        I_0: float,
        latitude: float,
        longitude: float,
        solar_declination: float,
        tau_atm: float,
        num_layers: int
) -> np.ndarray:
    """
    calculates the scattered intensity of light at various sky positions

    args:
        zenith_angles: an array of zenith angles in degrees
        azimuth_angles: an array of azimuth angles in degrees
        sun_zenith: the zenith angle of the sun in degrees
        sun_azimuth: the azimuth angle of the sun in degrees
        I_0: the intensity of the incoming solar radiation
        latitude: the latitude of the observation point in degrees
        longitude: the longitude of the observation point in degrees
        solar_declination: the solar declination in degrees
        tau_atm: the atmospheric optical depth
        num_layers: the number of atmospheric layers to consider in the calculation

    returns:
        a 2D NumPy array representing the scattered intensity at each combination of zenith and azimuth angles
    """
    intensities = np.zeros((len(zenith_angles), len(azimuth_angles)))

    for i, zenith in enumerate(zenith_angles):
        for j, azimuth in enumerate(azimuth_angles):
            scattering_angle = calculate_scattering_angle(zenith, np.degrees(azimuth), sun_zenith, np.degrees(sun_azimuth))
            phase_function_value = rayleigh_phase_function(scattering_angle)
            base_intensity = intensity_at_ground_spherical(
                latitude, longitude, solar_declination, np.degrees(azimuth) - 180, tau_atm, num_layers
                )
            intensities[i, j] = I_0 * base_intensity * phase_function_value

    return intensities


def generate_unpolarized_fisheye_data(
        latitude: float,
        longitude: float,
        solar_declination: float,
        tau_atm: float,
        num_layers: int,
        hour_angle: float,
        num_zenith_points: int = 90,
        num_azimuth_points: int = 180
) -> np.ndarray:
    """
    generates fisheye data (scattered intensity) for a given location and time

    args:
        latitude: the latitude of the observation point in degrees
        longitude: the longitude of the observation point in degrees
        solar_declination: the solar declination in degrees
        tau_atm: the atmospheric optical depth
        num_layers: the number of atmospheric layers to consider in the calculation
        hour_angle: the hour angle of the sun in degrees
        num_zenith_points: the number of points to generate along the zenith angle axis (default: 90)
        num_azimuth_points: the number of points to generate along the azimuth angle axis (default: 180)

    returns:
        a 2D NumPy array representing the scattered intensity at various sky positions
    """
    zenith_angles, azimuth_angles = generate_sky_angles(num_zenith_points, num_azimuth_points)
    sun_zenith, sun_azimuth = calculate_sun_position(latitude, solar_declination, hour_angle)
    print("--------------------")
    print(f"GENERATED DATA || Unpolarized Sun Zenith: {sun_zenith}, Unpolarized Sun Azimuth: {sun_azimuth}")
    print("--------------------")

    I_0 = 1.0 # scaling factor for total intensity

    return calculate_intensity(
        zenith_angles, azimuth_angles, sun_zenith, sun_azimuth, I_0, latitude,
        longitude, solar_declination, tau_atm, num_layers
    )


def generate_polarized_fisheye_data(
        latitude: float,
        longitude: float,
        solar_declination: float,
        tau_atm: float,
        num_layers: int,
        hour_angle: float,
        num_zenith_points: int = 90,
        num_azimuth_points: int = 180
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    generates fisheye data (scattered intensity and polarization) for a given location and time

    args:
        latitude: the latitude of the observation point in degrees
        longitude: the longitude of the observation point in degrees
        solar_declination: the solar declination in degrees
        tau_atm: the atmospheric optical depth
        num_layers: the number of atmospheric layers to consider in the calculation
        hour_angle: the hour angle of the sun in degrees
        num_zenith_points: the number of points to generate along the zenith angle axis (default: 90)
        num_azimuth_points: the number of points to generate along the azimuth angle axis (default: 180)

    returns:
        a tuple containing three 2D NumPy arrays:
            - I: the intensity component of the stokes parameters
            - Q: the Q component of the stokes parameters
            - U: the U component of the stokes parameters
    """
    zenith_angles, azimuth_angles = generate_sky_angles(num_zenith_points, num_azimuth_points)
    sun_zenith, sun_azimuth = calculate_sun_position(latitude, solar_declination, hour_angle)
    print("--------------------")
    print(f"GENERATED DATA Polarized Sun Zenith: {sun_zenith}, Polarized Sun Azimuth: {sun_azimuth}")
    print("--------------------")

    I_0 = 1.0 # scaling factor for total intensity

    return calculate_polarization_stokes(
        zenith_angles=zenith_angles,
        azimuth_angles=azimuth_angles,
        sun_zenith=sun_zenith,
        sun_azimuth=sun_azimuth,
        I_0=I_0,
        latitude=latitude,
        longitude=longitude,
        solar_declination=solar_declination,
        tau_atm=tau_atm,
        num_layers=num_layers
    )
