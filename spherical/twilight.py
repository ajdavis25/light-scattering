#!/usr/bin/python3
''' This code produces a view of the sky at one wavelength for a spherical Earth and Sun-like star model, incorporating twilight scattering. '''

import math
import numpy as np
import matplotlib.pyplot as plt
from utils import create_twilight_colormap
from spherical_model import scattered_intensity_at_tau

# constants
EARTH_RADIUS = 6371 # in kilometers
SUN_ANGULAR_DIAMETER = 0.53 # approx degrees

"""
def sun_position(solar_zenith_angle, sun_altitude):
    # convert the altitude of the sun to account for twilight and horizon effects
    sun_lower_bound = solar_zenith_angle + SUN_ANGULAR_DIAMETER / 2.0
    sun_upper_bound = solar_zenith_angle - SUN_ANGULAR_DIAMETER / 2.0

    if sun_altitude < -SUN_ANGULAR_DIAMETER:
        # sun is completely below the horizon
        return False, 0
    return True, (sun_lower_bound, sun_upper_bound)
"""

def angular_separation(theta_0, theta_p, phi):
    # theta_0 is the solar zenith angle (position of the sun)
    # theta_p is the observation zenith angle (direction we're looking at)
    # phi is the azimuthal angle between the sun and the observation direction

    # convert angles to radians
    theta_0_rad = math.radians(theta_0)
    theta_p_rad = math.radians(theta_p)
    phi_rad = math.radians(phi)

    # calculate the angular separation between two directions in spherical coordinates
    separation = math.acos(math.sin(theta_0_rad) * math.sin(theta_p_rad) +
                           math.cos(theta_0_rad) * math.cos(theta_p_rad) * math.cos(phi_rad))
    
    # return separation in degrees
    return math.degrees(separation)


def calculate_energy(values, azimuths, zeniths):

    total_energy = 0.0
    for i in range(len(azimuths)):
        for j in range(len(zeniths)):
            # intensity value at the given point
            intensity = values[i, j]

            # calculate the area element in spherical coordinates
            dA = np.sin(np.radians(zeniths[j])) * (np.radians(azimuths[1] - azimuths[0])) * (np.radians(zeniths[1] - zeniths[0]))

            # add to total energy
            total_energy += intensity * dA

    return total_energy


def make_contour_plot(optical_quantity, solar_zenith_angle, tau_obs, tau_atm, direction):
    
    # generate data
    azimuths = np.radians(np.linspace(0, 360, 181)) # azimuth angles in radians
    zeniths = np.arange(0, 91, 2) # zenith angles from 0 to 90 degrees
    r, theta = np.meshgrid(zeniths, azimuths)
    values = np.zeros((azimuths.size, zeniths.size)) # initialize values to zero

    # calculate intensity or polarization based on the direction
    for i in range(len(azimuths)):
        for j in range(len(zeniths)):
            phi = math.degrees(azimuths[i])
            if direction == "downwelling":
                theta_p = 180.0 - zeniths[j]  # Downwelling for sunset
            else:
                theta_p = zeniths[j]
            
            # calculate angular separation between sun's position and the observation point
            separation_angle = angular_separation(solar_zenith_angle, theta_p, phi)

            # calculate intensity as a function of angular separation
            intensity = scattered_intensity_at_tau(separation_angle, theta_p, tau_obs, tau_atm)

            # emphasize scattering towards the sun's direction, especially at the horizon
            directional_factor = math.cos(math.radians(separation_angle / 2.0)) # adjust this for asymmetry
            
            # scale intensity with directional emphasis (removing the artificial decay)
            if optical_quantity == 'I':
                values[i, j] = intensity * directional_factor # direct intensity from scattering
            elif optical_quantity == 'Pol':
                # calculate polarization as a percentage of intensity
                values[i, j] = intensity * 0.5 * directional_factor # adjust polarization calculation if needed

    # calculate total energy in and out
    energy_in = 2 * np.pi # total energy from the sun (hemisphere)
    energy_out = calculate_energy(values, np.degrees(azimuths), zeniths)
    energy_ratio = energy_out / energy_in

    print(f"Total energy in: {energy_in:.5f}, total energy out: {energy_out:.5f}, energy ratio: {energy_ratio:.5f}")

    # plot the sky radiance using polar coordinates
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    twilight_colormap = create_twilight_colormap() # custom colormap for twilight effect
    contour = ax.contourf(theta, r, values, 70, cmap=twilight_colormap)
    plt.colorbar(contour, ax=ax, label='Normalized Radiance')
    ax.set_theta_zero_location("N")
    ax.set(xticks=np.arange(0, 2 * math.pi, math.pi / 6.0),
           yticks=np.arange(0, 91, 30))

    # corrected sun position plotting based on solar zenith
    if direction == "downwelling":
        if solar_zenith_angle <= 90.0:
            sun_zenith_pos = 180 - solar_zenith_angle # sun's zenith position on the plot
            sun = plt.Circle((0, sun_zenith_pos), 4.0, transform=ax.transData._b, color='red', label='Sun')
            ax.add_artist(sun)
        else:
            print(f"Sun is below the horizon at solar zenith angle = {solar_zenith_angle:.2f}, not plotting.")

    plt.title(f'{optical_quantity} - {direction.capitalize()}')
    plt.show()


def my_main():

    tau_atm = 0.4 # atmospheric optical depth (adjusted for twilight)
    tau_obs = 0.02 # observer's optical depth during twilight
    solar_zenith_angle = 96.0 # Sun's zenith angle, which controls the sun's position

    # Generate contour plots for intensity and polarization
    make_contour_plot('I', solar_zenith_angle, tau_obs, tau_atm, "downwelling")
    make_contour_plot('Pol', solar_zenith_angle, tau_obs, tau_atm, "downwelling")

if __name__ == '__main__':
    my_main()
