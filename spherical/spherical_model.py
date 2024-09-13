import math, numpy as np
from typing import Tuple


# constants
EARTH_RADIUS = 6371.0 # in km
SUN_DISTANCE = 1.496e8 # average distance from earth to sun in km


def photon_unit_vector_spherical(latitude: float, longitude: float, radius: float = EARTH_RADIUS) -> np.ndarray:
    """
    create a direction vector for a point on the spherical earth

    args:
        latitude: the latitude angle in degrees (-90 to 90)
        longitude: the longitude angle in degrees (-180 to 180)
        radius: the radius of the earth in km (default to 6371.0)

    returns:
        a 3D unit vector as a NumPy array representing a point on the earth's surface
    """
    lat_rad = math.radians(latitude)
    lon_rad = math.radians(longitude)

    x = radius * math.cos(lat_rad) * math.cos(lon_rad)
    y = radius * math.cos(lat_rad) * math.sin(lon_rad)
    z = radius * math.sin(lat_rad)

    return np.array([x, y, z])


def sun_position_vector(solar_declination: float, hour_angle: float, sun_distance: float = SUN_DISTANCE) -> np.ndarray:
    """
    create a direction vector for the sun based on solar declination and hour angle

    args:
        solar_declination: the angle of the Sun relative to the equator (degrees)
        hour_angle: the sun's position relative to the observer's local meridian (degrees)
        sun_distance: the distance from the earth to the sun (default is 1.496e8 km)

    returns:
        a 3D vector representing the sun's position relative to the earth
    """
    dec_rad = math.radians(solar_declination)
    hour_rad = math.radians(hour_angle)

    x = sun_distance * math.cos(dec_rad) * math.cos(hour_rad)
    y = sun_distance * math.cos(dec_rad) * math.sin(hour_rad)
    z = sun_distance * math.sin(dec_rad)

    return np.array([x, y, z])


def scattering_angle_spherical(vector_1: np.ndarray, vector_2: np.ndarray) -> float:
    """
    calculate the scattering angle between two vectors on a spherical earth

    args:
        vector_1: first 3D vector (e.g., the earth's surface normal)
        vector_2: second 3D vector (e.g., the sun's direction)

    returns:
        the scattering angle in degrees between the two vectors
    """
    dot_product = np.dot(vector_1, vector_2)
    dot_product = np.clip(dot_product, -1.0, 1.0) # to avoid precision errors outside [-1, 1]
    psi_rad = math.acos(dot_product)
    
    return math.degrees(psi_rad)


def rayleigh_phase_function(psi: float) -> float:
    """
    calculate the rayleigh phase function for a given scattering angle

    args:
        psi: scattering angle in degrees

    returns:
        the value of the rayleigh phase function at the given scattering angle
    """
    psi_rad = math.radians(psi)
    phase_function = (3.0 / (16.0 * math.pi)) * (1.0 + math.cos(psi_rad) ** 2)

    return phase_function


def mie_phase_function(psi: float) -> float:
    """
    calculate the mie phase function for a given scattering angle

    args:
        psi: scattering angle in degrees

    returns:
        the value of the mie phase function at the given scattering angle
    """
    psi_rad = math.radians(psi)
    g = 0.75 # asymmetry parameter for forward scattering
    phase_function = (1 - g**2) / (4 * math.pi * (1 + g**2 - 2 * g * math.cos(psi_rad))**(3/2))

    return phase_function


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
