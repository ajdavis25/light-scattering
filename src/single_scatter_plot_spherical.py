import os, numpy as np, matplotlib.pyplot as plt
from typing import Tuple, List, Dict
from matplotlib.colors import Normalize, LinearSegmentedColormap
from intensity_vs_theta_spherical import intensity_at_ground_spherical
from spherical_model import rayleigh_phase_function


def calculate_sun_position(latitude: float, solar_declination: float, hour_angle: float) -> Tuple[float, float]:
    """
    calculate the sun's position (zenith angle and azimuth angle) based on the latitude, solar declination, and hour angle

    args:
        latitude: the observer's latitude in degrees
        solar_declination: the sun's declination in degrees
        hour_angle: the sun's hour angle relative to the observer's meridian in degrees

    returns:
        tuple of (zenith_angle, azimuth_angle) in degrees
    """
    # convert angles to radians for calculation
    latitude_rad = np.radians(latitude)
    declination_rad = np.radians(solar_declination)
    hour_angle_rad = np.radians(hour_angle)

    # compute the cosine of the zenith angle using the spherical trigonometry formula
    cos_zenith_angle = (np.sin(latitude_rad) * np.sin(declination_rad) +
                        np.cos(latitude_rad) * np.cos(declination_rad) * np.cos(hour_angle_rad))

    # ensure the value is within [-1, 1] to avoid numerical errors
    cos_zenith_angle = np.clip(cos_zenith_angle, -1.0, 1.0)

    # zenith angle is the arccos of the cosine value
    zenith_angle = np.degrees(np.arccos(cos_zenith_angle))

    # azimuth angle (can approximate directly from the hour angle)
    azimuth_angle = hour_angle

    return zenith_angle, azimuth_angle


def generate_fisheye_data(
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
    generate intensity values across the sky for a fisheye plot using the spherical earth model

    args:
        latitude: observer's latitude in degrees
        longitude: observer's longitude in degrees
        solar_declination: solar declination angle in degrees
        tau_atm: atmospheric optical depth
        num_layers: number of atmospheric layers to integrate over
        hour_angle: hour angle of the sun (affects azimuth position of the sun)
        num_zenith_points: number of radial (zenith) points to calculate
        num_azimuth_points: number of azimuthal points to calculate

    returns:
        2D numpy array of intensity values across the sky for the fisheye plot
    """
    # arrays to hold the zenith and azimuth angles
    zenith_angles = np.linspace(0, 90, num_zenith_points) # zenith: 0 (zenith) to 90 (horizon)
    azimuth_angles = np.linspace(0, 2 * np.pi, num_azimuth_points) # azimuth in radians (0, 2pi)

    # calculate the Sun's position (zenith and azimuth angles)
    sun_zenith, sun_azimuth = calculate_sun_position(latitude, solar_declination, hour_angle) # adjust hour_angle as needed

    # create a 2D grid to store intensities for each (zenith, azimuth) pair
    intensities = np.zeros((num_zenith_points, num_azimuth_points))

    # compute intensity for each zenith and azimuth angle based on the scattering angle
    for i, zenith in enumerate(zenith_angles):
        for j, azimuth in enumerate(azimuth_angles):
            # calculate the scattering angle between the sun's position and the point in the sky
            scattering_angle = calculate_scattering_angle(zenith, np.degrees(azimuth), sun_zenith, np.degrees(sun_azimuth))

            # compute the phase function based on the scattering angle
            phase_function_value = rayleigh_phase_function(scattering_angle)

            # calculate the base intensity using the spherical model, considering the sun's position and atmospheric layers
            base_intensity = intensity_at_ground_spherical(latitude, longitude, solar_declination, np.degrees(azimuth) - 180, tau_atm, num_layers)

            # multiply the base intensity by the phase function
            intensities[i, j] = base_intensity * phase_function_value

    return intensities


def calculate_scattering_angle(zenith1: float, azimuth1: float, zenith2: float, azimuth2: float) -> float:
    """
    calculate the scattering angle between two points on the celestial sphere

    args:
        zenith1: zenith angle of the first point (degrees)
        azimuth1: azimuth angle of the first point (degrees)
        zenith2: zenith angle of the second point (degrees)
        azimuth2: azimuth angle of the second point (degrees)

    returns:
        the scattering angle (in degrees) between the two points
    """
    # convert angles to radians for calculation
    zenith1_rad = np.radians(zenith1)
    azimuth1_rad = np.radians(azimuth1)
    zenith2_rad = np.radians(zenith2)
    azimuth2_rad = np.radians(azimuth2)

    # apply spherical trigonometry to calculate the scattering angle
    cos_scattering_angle = (np.sin(zenith1_rad) * np.sin(zenith2_rad) *
                            np.cos(azimuth1_rad - azimuth2_rad) +
                            np.cos(zenith1_rad) * np.cos(zenith2_rad))

    # clip the result to avoid numerical issues with arccos
    cos_scattering_angle = np.clip(cos_scattering_angle, -1.0, 1.0)

    # convert to degrees
    scattering_angle = np.degrees(np.arccos(cos_scattering_angle))

    return scattering_angle


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


def plot_fisheye_twilight_sky(
    intensities: np.ndarray,
    latitude: float,
    solar_declination: float,
    hour_angle: float,
    save_path: str = './plots',
    save_name: str = 'twilight_fisheye_spherical.png'
) -> None:
    """
    create a fisheye plot representing the twilight sky based on scattering intensities

    args:
        intensities: 2D array of scattering intensity values
        latitude: the observer's latitude
        solar_declination: the sun's declination angle
        hour_angle: the sun's hour angle relative to the observer's meridian
        save_path: directory to save the plot
        save_name: file name to save the plot
    """
    # create directory if it doesn't exist
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    num_zenith_points, num_azimuth_points = intensities.shape

    # create polar plot
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(6, 6))

    # create the meshgrid for zenith and azimuth
    zenith_angles = np.linspace(0, 90, num_zenith_points) # radial distance: 0 = horizon, 90 = zenith, 180 = opposite horizon
    azimuth_angles = np.linspace(0, 2 * np.pi, num_azimuth_points) # azimuth in radians

    # azimuth is theta, zenith is r
    theta, r = np.meshgrid(azimuth_angles, zenith_angles)

    # normalize the intensity data for better color mapping
    norm = Normalize(vmin=np.min(intensities), vmax=np.max(intensities))

    # use the custom twilight colormap
    twilight_colormap = create_twilight_colormap()

    # plot the twilight sky as a color-mapped polar plot
    ax.pcolormesh(theta, r, intensities, shading='auto', norm=norm, cmap=twilight_colormap)

    # range from one horizon (0) to the opposite horizon (90)
    ax.set_ylim(0, 90)
    ax.set_yticklabels([]) # hide radial ticks for cleaner look

    # plot the sun's position as a yellow circle, where sun_zenith is now directly the radial distance
    sun_zenith, sun_azimuth = calculate_sun_position(latitude, solar_declination, hour_angle)
    ax.scatter(np.radians(sun_azimuth), sun_zenith, color='yellow', s=100, label='Sun')

    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1) # counterclockwise azimuth
    
    # customize the look of the plot
    ax.set_title('Twilight Sky Intensity (Spherical Model)', va='bottom')
    ax.legend(loc='upper right')

    # save the plot
    plt.savefig(os.path.join(save_path, save_name))
    plt.show()
