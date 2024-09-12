import os, numpy as np, matplotlib.pyplot as plt
from typing import Tuple, List, Dict
from matplotlib.colors import Normalize, LinearSegmentedColormap
from intensity_vs_theta_spherical import intensity_at_ground_spherical


def calculate_sun_position(solar_declination: float, hour_angle: float) -> Tuple[float, float]:
    """
    calculate the sun's position in zenith and azimuth angles based on solar declination and hour angle

    args:
        solar_declination: angle of the sun relative to the equator in degrees
        hour_angle: the sun's position relative to the observer's local meridian in degrees

    returns:
        tuple of (zenith_angle, azimuth_angle) in degrees
    """
    zenith_angle = 90 - solar_declination # zenith angle (90 is horizon)
    azimuth_angle = hour_angle # azimuth angle directly corresponds to hour angle

    return zenith_angle, azimuth_angle


def generate_fisheye_data(
        latitude: float,
        longitude: float,
        solar_declination: float,
        tau_atm: float,
        num_layers: int,
        num_zenith_points: int = 90,
        num_azimuth_points: int = 180
) -> np.ndarray:
    """
    generate intensity values across the sky for a fisheye plot using the spherical earth model

    args:
        latitude: observer's latitude in degrees
        longitude: observer's longitude in degrees
        solar_declination: solar declination angle in degrees
        tau_atm: atmospheric optical depth
        num_layers: number of atmospheric layers to integrate over
        num_zenith_points: number of radial (zenith) points to calculate
        num_azimuth_points: number of azimuthal points to calculate

    returns:
        2D numpy array of intensity values across the sky for the fisheye plot
    """
    # arrays to hold the zenith angles and azimuth angles
    zenith_angles = np.linspace(0, 90, num_zenith_points) # 0 (zenith) to 90 (horizon)
    azimuth_angles = np.linspace(0, 360, num_azimuth_points) # 360 degrees around

    # create a 2D grid to store intensities for each (zenith, azimuth) pair
    intensities = np.zeros((num_zenith_points, num_azimuth_points))

    # compute intensity for each zenith and azimuth angle
    for i, zenith in enumerate(zenith_angles):
        for j, azimuth in enumerate(azimuth_angles):
            # adjust the hour angle based on azimuth
            hour_angle = azimuth - 180 # map azimuth to [-180, 180] degrees
            intensities[i, j] = intensity_at_ground_spherical(latitude, longitude, solar_declination, hour_angle, tau_atm, num_layers)

    return intensities


def create_twilight_colormap() -> LinearSegmentedColormap:
    """
    create a custom colormap to represent the twilight sky

    returns:
        a LinearSegmentedColormap object
    """
    # define key colors and positions in the gradient (0 = far from sun, 1 = near the sun)
    cdict: Dict[str, List[Tuple[float, float, float]]] = {
        'red':   [(0.0, 0.0, 0.0),   # black (night sky)
                  (0.2, 0.2, 0.2),   # dark blue
                  (0.5, 0.5, 0.5),   # purple
                  (0.7, 0.9, 0.9),   # red/orange (sunset colors)
                  (1.0, 1.0, 1.0)],  # yellow (bright sunlight near horizon)
        
        'green': [(0.0, 0.0, 0.0),   # black
                  (0.2, 0.1, 0.1),   # dark blue
                  (0.5, 0.0, 0.0),   # purple
                  (0.7, 0.5, 0.5),   # red/orange
                  (1.0, 1.0, 1.0)],  # yellow
        
        'blue':  [(0.0, 0.0, 0.0),   # black
                  (0.2, 0.5, 0.5),   # dark blue
                  (0.5, 0.5, 0.5),   # purple
                  (0.7, 0.0, 0.0),   # red/orange
                  (1.0, 0.0, 0.0)]   # yellow
    }

    # create the colormap object
    twilight_colormap = LinearSegmentedColormap('twilight_sky', cdict)

    return twilight_colormap


def plot_fisheye_twilight_sky(intensities: np.ndarray, solar_declination: float, hour_angle: float, save_path: str = './plots', save_name: str = 'twilight_fisheye_spherical.png'):
    """
    create a fisheye plot representing the twilight sky based on scattering intensities

    args:
        intensities: 2D array of scattering intensity values
        solar_declination: angle of the sun relative to the equator in degrees
        hour_angle: the sun's position relative to the observer's local meridian in degrees
        save_path: directory to save the plot
        save_name: name of the saved plot file
    """
    # create directory if it doesn't exist
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    num_zenith_points, num_azimuth_points = intensities.shape

    # create polar plot
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(6, 6))

    # create the meshgrid for zenith and azimuth
    zenith_angles = np.linspace(0, 90, num_zenith_points) # radius in degrees (0 to 90)
    azimuth_angles = np.linspace(0, 2 * np.pi, num_azimuth_points) # azimuth in radians

    # normalize the intensity data for better color mapping
    norm = Normalize(vmin=np.min(intensities), vmax=np.max(intensities))

    # use the custom twilight colormap
    twilight_colormap = create_twilight_colormap()

    # plot the twilight sky as a color-mapped polar plot
    theta, r = np.meshgrid(azimuth_angles, zenith_angles)
    ax.pcolormesh(theta, r, intensities, shading='auto', norm=norm, cmap=twilight_colormap)

    # set radial limits (zenith angle, 0 = zenith, 90 = horizon)
    ax.set_ylim(0, 90)
    ax.set_yticklabels([]) # hide radial ticks for cleaner look

    # plot the sun's position as a yellow circle
    sun_zenith, sun_azimuth = calculate_sun_position(solar_declination, hour_angle)
    ax.scatter(np.radians(sun_azimuth), sun_zenith, color='yellow', s=100, label='Sun')

    # customize the look of the plot
    ax.set_title('Twilight Sky Intensity (Spherical Model)', va='bottom')
    ax.legend(loc='upper right')

    # save the plot
    plt.savefig(os.path.join(save_path, save_name))
    plt.show()
