import numpy as np
from typing import Tuple


def calculate_energy_loss(
        intensities: np.ndarray,
        sun_zenith: float,
        zenith_angles: np.ndarray,
        azimuth_angles: np.ndarray,
        I_0: float = 1.0
) -> Tuple[float, float, float]:
    """
    calculates the energy loss due to scattering based on intensity data.

    args:
        intensities: a 2D NumPy array representing the scattered intensities at different zenith and azimuth angles
        sun_zenith: the zenith angle of the sun in degrees
        zenith_angles: an array of zenith angles in degrees
        azimuth_angles: an array of azimuth angles in degrees
        I_0: the intensity of the incoming solar radiation (default: 1.0)

    returns:
        a tuple containing:
            - E_in: the incoming energy
            - E_out: the outgoing energy
            - energy_loss_ratio: the ratio of outgoing energy to incoming energy
    """
    E_in = I_0 * max(0, np.cos(np.radians(sun_zenith)))
    print("--------------------")
    print(f"Energy in: {E_in}")

    # calculate E_out (numerical integration over zenith and azimuth angles)
    d_zenith = np.radians(np.diff(zenith_angles)) # small zenith angle intervals in radians
    d_azimuth = np.radians(np.diff(azimuth_angles)) # small azimuth angle intervals in radians

    cos_zenith = np.cos(np.radians(zenith_angles[:-1])) # exclude last point to match diff dimentions
    sin_zenith = np.sin(np.radians(zenith_angles[:-1]))

    # E_out is the sum of the intensity contributions over the entire sky
    E_out = np.sum(intensities[:-1, :-1] * sin_zenith[:, None] * cos_zenith[:, None] * d_zenith[:, None] * d_azimuth)
    print(f"Energy out: {E_out}")

    # optional: calculate and return energy loss ratio
    energy_loss_ratio = E_out / E_in if E_in != 0 else None
    print(f"Energy loss ratio (E_out / E_in): {energy_loss_ratio}")
    print("--------------------")

    return E_in, E_out, energy_loss_ratio
