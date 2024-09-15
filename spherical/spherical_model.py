import math

def scattered_intensity_at_tau(separation_angle, theta_obs, tau, tau_max):
    ''' Rayleigh scattering model for a spherical atmosphere '''
    mu_obs = max(abs(math.cos(math.radians(theta_obs))), 1e-6)  # Avoid division by zero
    mu_0 = max(abs(math.cos(math.radians(separation_angle))), 1e-6)  # Avoid division by zero

    # Rayleigh scattering factor: I(theta) ~ (1 + cos^2(theta)) / r^2
    scattering_angle = math.radians(separation_angle)  # Scattering angle in radians
    rayleigh_factor = (1 + math.cos(scattering_angle)**2)

    # Distance attenuation: decay due to distance traveled through the atmosphere
    distance_factor = 1.0 / (1.0 + separation_angle**2 / 100.0)  # You can adjust the 100 for realism

    # Overall scattered intensity
    intensity = rayleigh_factor * distance_factor

    return intensity
