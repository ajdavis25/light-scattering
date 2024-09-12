""" single scattering of intensity vs. theta using a spherical model of the earth and a sun that isn't an infinitely small point """

import os, math, numpy as np, matplotlib.pyplot as plt
from typing import List, Tuple
from spherical_model import (
    photon_unit_vector_spherical,
    sun_position_vector,
    rayleigh_phase_function
)


def air_mass(zenith_angle: float) -> float:
    """
    calculate the air mass for a given zenith angle
    """
    if zenith_angle > 89.0:
        return 40.0 # approximate max air mass near the horizon
    
    return 1.0 / (math.cos(math.radians(zenith_angle)) + 0.15 * (93.885 - zenith_angle) ** -1.253) # kasten-young model


def refraction_correction(zenith_angle: float) -> float:
    """
    apply a basic atmospheric refraction correction for light near the horizon
    """
    if zenith_angle < 90.0:
        z = math.radians(zenith_angle)
        # approximate refraction formula for standard atmosphere
        return (1.02 / math.tan(z + (10.3 / (z + 5.11)))) / 60.0 # returns refraction in degrees
    
    return 0.0


def gaussian_kernel(size: int, sigma: float) -> np.ndarray:
    """
    create a 1D gaussian kernel using the specified size and standard deviation (sigma)

    args:
        size: size of the gaussian kernel (typically an odd number)
        sigma: standard deviation of the gaussian distribution

    returns:
        1D numpy array representing the gaussian kernel
    """
    # create an array of values [-size // 2, size // 2]
    x = np.arange(-size // 2 + 1, size // 2 + 1)

    # compute the 1D gaussian distribution for these values
    kernel = np.exp(-x ** 2 / (2 * sigma ** 2))

    # normalize the kernel to ensure the sum is 1
    return kernel / np.sum(kernel)


def apply_gaussian_smoothing(data: np.ndarray, sigma: float) -> np.ndarray:
    """
    apply a gaussian smoothing filter to the input data

    args:
        data: 1D numpy array of intensity values to smooth
        sigma: standard deviation of the gaussian kernel

    returns:
        smoothed 1D numpy array
    """
    # choose the size of the gaussian kernel (3 * sigma is a good rule of thumb)
    kernel_size = int(3 * sigma) * 2 + 1

    # generate the gaussian kernel
    kernel = gaussian_kernel(kernel_size, sigma)

    # convolve the data with the kernel using 'same' mode to keep the same size
    smoothed_data = np.convolve(data, kernel, mode='same')

    return smoothed_data


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
