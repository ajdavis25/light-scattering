import os, numpy as np, matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from typing import Tuple
from utils import (
    save_plot,
    create_twilight_colormap
)
from spherical.v1.intensity_vs_theta_spherical import (
    generate_intensity_spherical
)
from spherical_model import (
    calculate_sun_position
)
from spherical.v1.fisheye_projection import (
    generate_sky_angles,
    generate_polarized_fisheye_data
)


def plot_intensity_vs_theta_spherical(
        latitude: float = 30.0,
        longitude: float = 0.0,
        solar_declination: float = 23.5, # approx for summer solstice
        tau_atm: float = 0.5,
        num_layers: int = 20,
        save_path: str = './plots',
        save_name: str = 'intensity_vs_theta_spherical.png'
) -> None:
    """
    plot intensity vs. theta (hour angle) for a given latitude, longitude, and solar declination

    args:
        latitude: observer's latitude
        longitude: observer's longitude
        solar_declination: sun's declination angle
        tau_atm: atmospheric tau value (default: 0.5)
        num_layers: number of layers in the atmosphere for integration (default: 10)
        save_path: directory to save the plot
        save_name: file name for the saved plot
    """

    # create directory if it doesn't exist
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # generate intensity data for a range of hour angles (theta values)
    theta_obs, intensity = generate_intensity_spherical(latitude, longitude, solar_declination, tau_atm, num_layers)

    # plot the data
    plt.plot(theta_obs, intensity, label=f'Lat: {latitude}, Solar Dec: {solar_declination}')

    plt.legend()
    plt.title('Intensity vs. Hour Angle (Spherical Model)')
    plt.xlabel('Hour Angle (Degrees)')
    plt.ylabel('Intensity')

    # save and display the plot
    plt.savefig(os.path.join(save_path, save_name))
    plt.show()


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

    ax.scatter(np.radians(sun_position[1]), sun_position[0], color='yellow', s=100, label='Sun')
    print("--------------------")
    print(f"PLOTTING || Sun Zenith: {sun_position[0]}, Sun Azimuth: {sun_position[1]} (Radians)")
    print("--------------------")

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
    print("--------------------")
    print(f"PLOTTED VALUES || Unpolarized Sun Zenith: {sun_zenith}, Unpolarized Sun Azimuth: {sun_azimuth}")
    print("--------------------")

    plot_fisheye(
        intensities=intensities,
        sun_position=(sun_zenith, sun_azimuth),
        title="Unpolarized Twilight Sky Intensity",
        save_path=save_path,
        save_name=save_name
    )


def plot_polarized_fisheye(
        longitude: float,
        latitude: float,
        solar_declination: float,
        tau_atm: float,
        hour_angle: float,
        num_layers: int,
        num_zenith_points: int = 90,
        num_azimuth_points: int = 180,
        save_path: str = './plots',
        save_name: str = 'twilight_fisheye_polarized.png'
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

    # calculate the degree of polarization
    P = np.sqrt(Q ** 2 + U ** 2) / I
    P = np.clip(P, 0, None)  # clip values to avoid negatives

    sun_zenith, sun_azimuth = calculate_sun_position(latitude, solar_declination, hour_angle)
    print("--------------------")
    print(f"PLOTTING VALUES || Plotted Polarized Sun Zenith: {sun_zenith}, Plotted Polarized Sun Azimuth: {sun_azimuth}")
    print("--------------------")

    plot_fisheye(
        intensities=P,
        sun_position=(sun_zenith, sun_azimuth),
        title="Polarized Twilight Sky Intensity",
        save_path=save_path,
        save_name=save_name
    )
