import math
import numpy as np
import matplotlib.pyplot as plt
from utils import create_twilight_colormap
from vector_operations import rayleigh_matrix, normalize, dot_product


# constants
R_EARTH = 6371e3 # radius of the Earth in meters
SUN_RADIUS = 6.957e8 # radius of the Sun in meters
AU = 1.496e11 # astronomical Unit in meters
SUN_ANGULAR_RADIUS = np.degrees(np.arctan(SUN_RADIUS / AU)) # angular radius of the sun in degrees (~0.27 degrees)


def rayleigh_phase_function(psi):
    psi_rad = math.radians(psi)
    return (3.0 / (16.0 * math.pi)) * (1.0 + math.pow(math.cos(psi_rad), 2))


def photon_unit_vector(observer_position, target_position):
    direction_vector = normalize(target_position - observer_position)
    return direction_vector


def calculate_optical_depth(theta_obs, tau, tau_max):
    mu_obs = abs(math.cos(math.radians(theta_obs)))
    optical_depth = (tau_max - tau) / mu_obs if mu_obs > 0 else tau_max
    return max(optical_depth, 0.01) # ensure non-negative optical depth


def observer_position_on_earth(zenith_angle, azimuth_angle):
    zenith_rad = math.radians(zenith_angle)
    azimuth_rad = math.radians(azimuth_angle)
    x = R_EARTH * math.sin(zenith_rad) * math.cos(azimuth_rad)
    y = R_EARTH * math.sin(zenith_rad) * math.sin(azimuth_rad)
    z = R_EARTH * math.cos(zenith_rad)
    return np.array([x, y, z])


def scattered_intensity_at_tau(theta_0, theta_obs, phi_obs, tau, tau_max, sun_position):
    observer_position = observer_position_on_earth(theta_obs, phi_obs)
    # compute the vector from the observer to the sun and normalize
    photon_direction = photon_unit_vector(observer_position, sun_position)
    
    # dot product between the photon direction and the observer's view direction
    dot_prod = np.clip(dot_product(photon_direction, observer_position), -1.0, 1.0)
    scattering_angle = math.degrees(math.acos(dot_prod))
    phase_function = rayleigh_phase_function(scattering_angle)
    optical_depth = calculate_optical_depth(theta_obs, tau, tau_max)
    
    # attenuate the intensity based on optical depth and phase function
    intensity = math.exp(-optical_depth) * phase_function
    if theta_obs >= 80:
        intensity *= 2

    return max(intensity * math.cos(math.radians(theta_0)), 0)


def calculate_sun_position(solar_zenith_deg, solar_azimuth_deg):
    theta_sun = np.radians(solar_zenith_deg)
    phi_sun = np.radians(solar_azimuth_deg)
    sin_theta = np.sin(theta_sun)
    
    # compute 3D position of the sun
    sun_pos = AU * np.array([sin_theta * np.cos(phi_sun), sin_theta * np.sin(phi_sun), np.cos(theta_sun)])
    return sun_pos


def compute_radiance_spherical(theta_0, tau_obs, tau_atm, sun_position):
    azimuths = np.radians(np.linspace(0, 360, 361)) # ensure full 360° coverage
    zeniths = np.linspace(0, 90, 91) # capture from zenith to horizon
    r, theta = np.meshgrid(zeniths, azimuths)

    values = np.zeros((azimuths.size, zeniths.size))

    # loop over azimuths and zeniths
    for i in range(len(azimuths)):
        for j in range(len(zeniths)):
            phi = math.degrees(azimuths[i])
            theta_obs = zeniths[j]

            # compute intensity at each angle
            intensity = scattered_intensity_at_tau(theta_0, theta_obs, phi, tau_obs, tau_atm, sun_position)
            values[i, j] = intensity

            # debug output for key values
            if i % 90 == 0 and j % 45 == 0: # print for select grid points
                print(f"At azimuth {phi}, zenith {theta_obs}: intensity={intensity}")

    return r, theta, values


def generate_fisheye_plot(theta_0, tau_obs, tau_atm, sun_zenith, sun_azimuth):
    sun_position = calculate_sun_position(sun_zenith, sun_azimuth)
    r, theta, values = compute_radiance_spherical(theta_0, tau_obs, tau_atm, sun_position)

    # normalize radiance values
    max_radiance = np.max(values)
    values = values / max_radiance if max_radiance != 0 else values

    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    twilight_colormap = create_twilight_colormap()

    # ensure full coverage of azimuthal angles (0 to 2*pi for 360°)
    contour = ax.contourf(theta, r, values, 100, cmap=twilight_colormap, extend='both')

    # plot the sun at the correct location (sun can now be below the horizon)
    sun_r = 90 - sun_zenith # allow sun to go below the horizon for twilight
    sun_azimuth_rad = np.radians(sun_azimuth)

    # plot the sun as a yellow circle
    sun_circle = plt.Circle((sun_azimuth_rad, max(sun_r, 0)), 4.0, transform=ax.transData._b, color='yellow', label="Sun")
    ax.add_artist(sun_circle)

    plt.colorbar(contour, label='Normalized Radiance')
    ax.set_theta_zero_location("N") # ensure north is upward
    ax.set_theta_direction(-1) # clockwise direction for azimuthal angles
    plt.legend(loc="upper right")

    plt.title('Twilight Sky Intensity (Spherical Model)')
    plt.show()


def my_main():
    tau_atm = 0.5 # adjusted for twilight
    tau_obs = 0.0 # observer optical depth
    solar_zenith_angle = 90.0 # sun is below the horizon for nautical twilight
    solar_azimuth_angle = 0.0 # set to any valid azimuth for now
    theta_0 = 180.0 - solar_zenith_angle # observer's zenith relative to the sun

    generate_fisheye_plot(theta_0, tau_obs, tau_atm, solar_zenith_angle, solar_azimuth_angle)

if __name__ == '__main__':
    my_main()
