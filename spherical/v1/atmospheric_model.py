import math


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


def rayleigh_phase_function(psi: float) -> float:
    """
    calculate the rayleigh phase function for a given scattering angle

    args:
        psi: scattering angle in degrees

    returns:
        the value of the rayleigh phase function at the given scattering angle
    """
    psi_rad = math.radians(psi)
    phase_function = (3.0 / (16.0 * math.pi)) * (1.0 + math.cos(psi_rad)**2)

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
