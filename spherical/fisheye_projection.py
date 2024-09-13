import numpy as np, matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from typing import Tuple
from utils import (
    save_plot,
    create_twilight_colormap
)
from spherical_model import (
    rayleigh_phase_function,
    calculate_sun_position,
    calculate_scattering_angle
)
from polarization_model import calculate_polarization_stokes
from intensity_vs_theta_spherical import intensity_at_ground_spherical


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
    I_0 = 1.0 # scaling factor for total intensity

    return calculate_polarization_stokes(
        zenith_angles, azimuth_angles, sun_zenith, sun_azimuth, I_0, latitude,
        longitude, solar_declination, tau_atm, num_layers
    )


def plot_fisheye(
    intensities: np.ndarray,
    sun_position: Tuple[float, float],
    title: str,
    save_path: str,
    save_name: str
) -> None:
    """
    plots a fisheye projection of the given intensity data

    args:
        intensities: a 2D NumPy array representing the intensity values at different zenith and azimuth angles
        sun_position: a tuple containing the sun's zenith and azimuth angles in degrees
        title: the title for the plot
        save_path: the path to the directory where the plot will be saved
        save_name: the filename for the saved plot
    """
    num_zenith_points, num_azimuth_points = intensities.shape
    zenith_angles, azimuth_angles = generate_sky_angles(num_zenith_points, num_azimuth_points)
    theta, r = np.meshgrid(azimuth_angles[:-1], zenith_angles[:-1])

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(6, 6))
    twilight_colormap = create_twilight_colormap()
 
    norm = Normalize(vmin=np.min(intensities), vmax=np.max(intensities))
    cax = ax.pcolormesh(theta, r, intensities[:-1, :-1], shading='auto', norm=norm, cmap=twilight_colormap)
    fig.colorbar(cax)

    # rotate the plot and adjust azimuth direction
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    ax.scatter(np.radians(sun_position[1]), sun_position[0], color='yellow', s=100, label='Sun')
    ax.set_ylim(0, 90)
    ax.set_title(title)

    save_plot(fig, save_path, save_name)


def plot_unpolarized_fisheye(
    intensities: np.ndarray,
    latitude: float,
    solar_declination: float,
    hour_angle: float,
    save_path: str = './plots',
    save_name: str = 'twilight_fisheye_unpolarized.png'
) -> None:
    """
    plots an unpolarized fisheye projection of the given intensity data

    args:
        intensities: a 2D NumPy array representing the intensity values at various sky positions
        latitude: the latitude of the observation point in degrees
        solar_declination: the solar declination in degrees
        hour_angle: the hour angle of the sun in degrees
        save_path: the path to the directory where the plot will be saved
        save_name: the filename for the saved plot
    """
    sun_zenith, sun_azimuth = calculate_sun_position(latitude, solar_declination, hour_angle)
    plot_fisheye(
        intensities, (sun_zenith, sun_azimuth), 'Unpolarized Twilight Sky Intensity', save_path, save_name
    )


def plot_polarized_fisheye(
    latitude: float,
    solar_declination: float,
    tau_atm: float,
    hour_angle: float,
    num_zenith_points: int = 90,
    num_azimuth_points: int = 180,
    save_path: str = './plots',
    save_name: str = 'twilight_fisheye_polarized.png',
    roll_amount: int = 0
) -> None:
    """
    plots a polarized fisheye projection of the twilight sky

    args:
        latitude: the latitude of the observation point in degrees
        solar_declination: the solar declination in degrees
        tau_atm: the atmospheric optical depth
        hour_angle: the hour angle of the sun in degrees
        num_zenith_points: the number of points to generate along the zenith angle axis (default: 90)
        num_azimuth_points: the number of points to generate along the azimuth angle axis (default: 180)
        save_path: the path to the directory where the plot will be saved
        save_name: the filename for the saved plot
        roll_amount: the amount to roll the data along the azimuth axis (default: 0)
    """
    I, Q, U = generate_polarized_fisheye_data(
        latitude, solar_declination, tau_atm, hour_angle, num_zenith_points, num_azimuth_points
        )
    
    # roll the data
    if roll_amount != 0:
        I = np.roll(I, roll_amount, axis=1)
        Q = np.roll(Q, roll_amount, axis=1)
        U = np.roll(U, roll_amount, axis=1)

    # calculate the degree of polarization
    P = np.sqrt(Q ** 2 + U ** 2) / I
    P = np.clip(P, 0, None) # clip values to avoid negatives

    sun_zenith, sun_azimuth = calculate_sun_position(latitude, solar_declination, hour_angle)

    plot_fisheye(
        P, (sun_zenith, sun_azimuth), 'Polarized Twilight Sky Intensity', save_path, save_name
    )
