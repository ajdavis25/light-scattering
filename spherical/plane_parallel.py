import math
import numpy as np
from vector_operations import rayleigh_matrix, rotation_angles, rotation_matrix

# Constants
EARTH_RADIUS = 6371  # Earth radius in kilometers
SUN_ANGULAR_DIAMETER = 0.53  # approx. angular size of Sun in degrees

def rayleigh(psi):
    ''' Rayleigh phase function '''
    psi_rad = math.radians(psi)
    phase_function = (3.0/(16.0*math.pi))*(1.0 + math.pow(math.cos(psi_rad), 2))
    return phase_function

def photon_unit_vector(theta, phi, observer_height=0):
    ''' Modify direction vector for a photon, considering a spherical Earth '''
    theta_rad = math.radians(theta)
    phi_rad = math.radians(phi)
    # Adjust for observer's height above the Earth's surface
    observer_offset = (EARTH_RADIUS + observer_height) / EARTH_RADIUS

    u_z = math.cos(theta_rad) * observer_offset
    u_x = math.sin(theta_rad) * math.cos(phi_rad) * observer_offset
    u_y = math.sin(theta_rad) * math.sin(phi_rad) * observer_offset

    return np.array([u_x, u_y, u_z])

def sun_position(solar_zenith_angle):
    ''' Calculate Sun's apparent position relative to the horizon '''
    sun_lower_bound = solar_zenith_angle + SUN_ANGULAR_DIAMETER / 2.0
    sun_upper_bound = solar_zenith_angle - SUN_ANGULAR_DIAMETER / 2.0
    return sun_lower_bound, sun_upper_bound

def scattered_intensity_at_tau(theta_0, theta_obs, tau, tau_max):
    ''' Improved scattered intensity calculation for twilight and spherical Earth with twilight cutoff '''
    
    mu_obs = abs(math.cos(math.radians(theta_obs)))
    mu_0 = abs(math.cos(math.radians(theta_0)))

    epsilon = 1e-8  # Small value to prevent division by zero
    max_exp_arg = 709  # Maximum argument for safe exp() calculation

    # Check if we're within the twilight range
    if theta_obs > 120.0:
        return 0.0  # No scattered light beyond 120 degrees
    
    if mu_0 == 0 or mu_obs == 0:
        return 0.0  # Sun is below the horizon or observer is at the horizon, no scattered intensity

    # Factor to adjust for increased optical depth
    factor0 = math.exp(-(tau_max - tau) / mu_0)

    # Prevent division by zero by ensuring mu_0 and mu_obs don't become too close
    if abs(mu_0 - mu_obs) < epsilon:
        return 0.0  # Return zero intensity when they are too close

    if theta_obs < 90.0:
        factor1 = mu_0 / (mu_0 + mu_obs)
        arg_inter = tau * (1.0 / mu_obs + 1.0 / mu_0)
        
        # Clamp the argument to prevent overflow
        if arg_inter > max_exp_arg:
            factor2 = 1.0  # exp(-arg_inter) becomes effectively 0 for large arg_inter
        else:
            factor2 = 1.0 - math.exp(-arg_inter)
        
        return factor0 * factor1 * factor2

    if theta_0 != theta_obs:
        factor1 = mu_0 * math.exp(-tau / mu_0) / (mu_0 - mu_obs + epsilon)  # Add epsilon to prevent zero division
        arg_inter = (tau_max - tau) * (1.0 / mu_obs - 1.0 / mu_0)

        # Clamp the argument to prevent overflow
        if arg_inter > max_exp_arg:
            factor2 = 1.0  # exp(-arg_inter) becomes effectively 0 for large arg_inter
        else:
            factor2 = 1.0 - math.exp(-arg_inter)
        
        return factor0 * factor1 * factor2

    # Final case if no overflow, using clamped values
    final_arg = -(tau_max - tau) / mu_obs
    if final_arg > max_exp_arg:
        return 0.0  # Avoid overflow
    return (tau_max - tau) * math.exp(final_arg) / mu_obs

def intensity_at_tau(theta_0, theta_obs, phi_obs, tau, tau_max):
    ''' Scattered intensity with Sun's angular size and twilight effects '''
    u_0 = photon_unit_vector(theta_0, 0.0)
    u_1 = photon_unit_vector(theta_obs, phi_obs)
    
    if phi_obs == 0.0 and theta_0 == theta_obs:
        phi, psi, theta = 0.0, 0.0, 0.0
    else:
        phi, psi, theta = rotation_angles(u_0, u_1)

    intensity = scattered_intensity_at_tau(theta_0, theta_obs, tau, tau_max) * rayleigh(theta)
    mueller_scatter = rayleigh_matrix(theta)
    mueller_rotation = rotation_matrix(psi)
    stokes_init = np.array([1.0, 0.0, 0.0, 0.0])
    scattered_stokes = np.dot(mueller_rotation, np.dot(mueller_scatter, stokes_init))

    return intensity * scattered_stokes

def main():
    ''' Main function with twilight and spherical Earth modeling '''
    theta_0 = 130.0  # Solar zenith angle
    theta_1 = 120.0
    phi_0 = 60.0
    phi_1 = 70.0

    u_0 = photon_unit_vector(theta_0, phi_0)
    u_1 = photon_unit_vector(theta_1, phi_1)
    phi, psi, theta = rotation_angles(u_0, u_1)

    print(f"Angle between directions: {theta}")

    # Test scattered intensity with updated model
    print("Testing scattered intensity with twilight and spherical Earth:")
    print(intensity_at_tau(130.0, 131.0, 30.0, 0.0, 0.5))
    print(intensity_at_tau(170.0, 170.0, 180.0, 0.0, 1.0))
    print(intensity_at_tau(170.0, 157.0, 180.0, 0.0, 1.0))

if __name__ == '__main__':
    main()
