import numpy as np
from typing import Tuple, List
from spherical.v1.atmospheric_model import (
    rayleigh_phase_function
)
from spherical_model import (
    calculate_scattering_angle
)


def polarized_intensity(
    latitude: float, 
    longitude: float, 
    solar_declination: float, 
    azimuth: float, 
    tau_atm: float, 
    num_layers: int,
    sun_azimuth: float,
    sun_zenith: float
) -> List[float]:
    """
    calculate the stokes parameters at ground level using the spherical earth model

    returns a list of four stokes parameters: [I, Q, U, V]
    """
    zenith_angle = 90 - latitude
    # calculate the scattering angle based on solar and observer geometry
    scattering_angle = calculate_scattering_angle(
        zenith_angle, azimuth, sun_azimuth, sun_zenith
    )

    # rayleigh scattering phase function for total intensity
    I_0 = 1.0 # base intensity (adjust based on optical depth)

    # calculate the phase function P(theta) for the scattering angle
    phase_function = rayleigh_phase_function(scattering_angle)

    # calculate total intensity I using the phase function
    I = I_0 * phase_function

    # degree of polarization for rayleigh scattering
    P_pol = (1 - np.cos(scattering_angle) ** 2) / (1 + np.cos(scattering_angle) ** 2)
    
    # linear polarization components Q and U
    Q = I * P_pol # linear polarization using the scattering plane
    
    # normalize Q relative to I to ensure it's not too large
    if Q > I:
        Q = I
        
    U = 0.0 # assume no U component for simplicity

    # circular polarization V is 0 for rayleigh scattering
    V = 0.0

    return [I, Q, U, V]


def calculate_polarization_stokes(
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
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    calculates the stokes parameters (I, Q, U) for scattered light at various sky positions

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
        a tuple containing three 2D NumPy arrays:
            - I: the intensity component of the stokes parameters
            - Q: the Q component of the stokes parameters
            - U: the U component of the stokes parameters
    """
    I = np.zeros((len(zenith_angles), len(azimuth_angles)))
    Q = np.zeros((len(zenith_angles), len(azimuth_angles)))
    U = np.zeros((len(zenith_angles), len(azimuth_angles)))

    for i, zenith in enumerate(zenith_angles):
        for j, azimuth in enumerate(azimuth_angles):
            scattering_angle = calculate_scattering_angle(zenith, np.degrees(azimuth), sun_zenith, np.degrees(sun_azimuth))
            phase_function_value = rayleigh_phase_function(scattering_angle)
            
            stokes_parameters = polarized_intensity(
                latitude, longitude, solar_declination, np.degrees(azimuth) - 180, tau_atm, num_layers, sun_azimuth, sun_zenith
            )
            
            I[i, j] = stokes_parameters[0] * phase_function_value
            Q[i, j] = stokes_parameters[1] * phase_function_value
            U[i, j] = stokes_parameters[2] * phase_function_value

    return I, Q, U
