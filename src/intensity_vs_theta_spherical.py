import os, math, matplotlib.pyplot as plt, numpy as np
from typing import List, Tuple
from spherical_model import photon_unit_vector_spherical, sun_position_vector, scattering_angle_spherical, rayleigh_phase_function


def intensity_at_ground_spherical(
    latitude: float,
    longitude: float,
    solar_declination: float,
    hour_angle: float,
    tau_max: float,
    num_layers: int = 10
) -> float:
    """
    calculate the intensity at the ground level based on the spherical Earth model and atmospheric scattering effects
    by integrating across multiple atmospheric layers

    args:
        latitude: the observer's latitude
        longitude: the observer's longitude
        solar_declination: the sun's hour angle
        hour_angle: the sun's hour angle
        tau_max: maximum optical depth
        num_layers: number of layers in the atmosphere for integration (default is 10)

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

        # calculate the phase function for scattering at this layer
        layer_intensity = rayleigh_phase_function(zenith_angle) * math.exp(-tau_layer / mu) * delta_tau

        # accumulate the intensity from this layer
        intensity += layer_intensity

    return intensity


def generate_intensity_spherical(
        latitude: float, longitude: float, solar_declination: float, tau_atm: float, num_layers: int
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

    return theta_obs, intensities


def plot_intensity_vs_theta_spherical(
        latitude: float = 30.0, # example latitude
        longitude: float = 0.0, # example longitude
        solar_declination: float = 23.5, # example solar declination (approx for summer solstice)
        tau_atm: float = 0.5,
        num_layers: int = 10,
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
