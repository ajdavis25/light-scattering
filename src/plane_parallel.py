#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 15 15:23:36 2022

@author: jadams

This code simulates the downwelling light field for an observer
standing at the base of a plane parallel atmosphere. This is a scalar
model using a standard Rayleigh phase function for scattering.

The ground is completely absorbing.
"""

import math, numpy as np
from vector_operations import rayleigh_matrix
from vector_operations import rotation_angles
from vector_operations import rotation_matrix


def rayleigh(psi):
    ''' Rayleigh phase function '''

    # Return the Rayleigh phase function for a scattering angle psi
    # expressed in degrees.

    psi_rad = math.radians(psi)
    phase_function = (3.0/(16.0*math.pi))*(1.0 + math.pow(math.cos(psi_rad), 2))

    return phase_function


def scattering_angle(vector_1, vector_2):
    ''' find angle between two vectors '''

    if np.dot(vector_1, vector_2) > 1.0:
        return 0.0
    elif np.dot(vector_1, vector_2) < -1.0:
        return 180.0
    else:
        psi = math.acos(np.dot(vector_1, vector_2))

    return math.degrees(psi)


def photon_unit_vector(theta, phi):
    ''' Create the direction vector for a photon heading in the
    theta, phi direction. '''

    theta_rad = math.radians(theta)
    phi_rad = math.radians(phi)

    u_z = math.cos(theta_rad)
    u_x = math.sin(theta_rad)*math.cos(phi_rad)
    u_y = math.sin(theta_rad)*math.sin(phi_rad)

    return np.array([u_x, u_y, u_z])


def scattered_intensity_reflected(theta_0, theta, tau_max):
    ''' Give the scattered intensity (without the phase function) of
        light at the top of the atmosphere. '''

    mu_obs = abs(math.cos(math.radians(theta)))
    mu_0 = abs(math.cos(math.radians(theta_0)))

    if theta_0 != theta:
        factor1 = mu_0/(mu_0 + mu_obs)
        arg_inter = tau_max*(1.0/mu_obs + 1.0/mu_0)
        factor2 = 1.0 - math.exp(-arg_inter)
        return factor1 * factor2

    return tau_max*math.exp(-tau_max/mu_0)/mu_obs


def scattered_intensity_transmitted(theta_0, theta, tau_max):
    ''' Give the scattered intensity (without the phase function) of
        light at the bottom of the atmosphere. '''

    mu_obs = abs(math.cos(math.radians(theta)))
    mu_0 = abs(math.cos(math.radians(theta_0)))

    if theta_0 != theta:
        factor1 = mu_0*math.exp(-tau_max/mu_0)/(mu_0 - mu_obs)
        arg_inter = tau_max*(1.0/mu_obs - 1.0/mu_0)
        factor2 = 1.0 - math.exp(-arg_inter)
        return factor1 * factor2

    return tau_max*math.exp(-tau_max/mu_0)/mu_obs


def intensity_at_ground(theta_0, theta, phi, tau_max):
    ''' Give the scattered intensity of
        light at the bottom of the atmosphere
        in the direction theta, phi. '''

    u_0 = photon_unit_vector(theta_0, 0.0)
    u_1 = photon_unit_vector(theta, phi)

    scatter_angle = scattering_angle(u_0, u_1)

    base_intensity = scattered_intensity_transmitted(theta_0, theta, tau_max) * \
                     rayleigh(scatter_angle)

    return base_intensity


def intensity_at_top(theta_0, theta, phi, tau_max):
    ''' Give the scattered intensity of
        light at the top of the atmosphere
        in the direction theta, phi. '''

    u_0 = photon_unit_vector(theta_0, 0.0)
    u_1 = photon_unit_vector(theta, phi)

    scatter_angle = scattering_angle(u_0, u_1)

    top_intensity = scattered_intensity_reflected(theta_0, theta, tau_max) * \
                     rayleigh(scatter_angle)

    return top_intensity


def scattered_intensity_at_tau(theta_0, theta_obs, tau, tau_max):
    ''' Give the scattered intensity (without the phase function) of
        light in the atmosphere. '''

    mu_obs = abs(math.cos(math.radians(theta_obs)))
    mu_0 = abs(math.cos(math.radians(theta_0)))

    factor0 = math.exp(-(tau_max - tau)/mu_0)

    if theta_obs < 90.0:

        factor1 = mu_0/(mu_0 + mu_obs)
        arg_inter = (tau)*(1.0/mu_obs + 1.0/mu_0)
        factor2 = 1.0 - math.exp(-arg_inter)

        return factor0 * factor1 * factor2
    
    if theta_0 != theta_obs:
            factor1 = mu_0*math.exp(-tau/mu_0)/(mu_0 - mu_obs)
            arg_inter = (tau_max - tau)*(1.0/mu_obs - 1.0/mu_0)
            factor2 = 1.0 - math.exp(-arg_inter)
            return factor0 * factor1 * factor2

    return (tau_max - tau) * math.exp(-(tau_max-tau)/mu_obs)/mu_obs 


def intensity_at_tau(theta_0, theta_obs, phi_obs, tau, tau_max):
    ''' Give the scattered intensity of
        light in the atmosphere
        in the direction theta, phi. '''

    u_0 = photon_unit_vector(theta_0, 0.0)
    u_1 = photon_unit_vector(theta_obs, phi_obs)

    if phi_obs == 0.0 and theta_0 == theta_obs:
        (phi, psi, theta) = (0.0, 0.0, 0.0)
    else:
        (phi, psi, theta) = rotation_angles(u_0, u_1)

    intensity = scattered_intensity_at_tau(theta_0, theta_obs, tau, tau_max) * \
                     rayleigh(theta)
    
    # (phi, psi, theta) = rotation_angles(u_0, u_1)
    mueller_scatter = rayleigh_matrix(theta)
    mueller_rotation = rotation_matrix(psi)
    stokes_init = np.array([1.0, 0.0, 0.0, 0.0])
    scattered_stokes = np.dot(mueller_rotation, np.dot(mueller_scatter, stokes_init))

    return intensity * scattered_stokes


def main():
    ''' Function main. '''

    theta_0 = 130.0
    theta_1 = 120.0
    phi_0 = 60.0
    phi_1 = 70.0

    u_0 = photon_unit_vector(theta_0, phi_0)
    u_1 = photon_unit_vector(theta_1, phi_1)

    print("angle between is ", scattering_angle(u_0, u_1))

    print("Test scattered intensity")
    print(scattered_intensity_transmitted(theta_0, theta_1 - 0.1, 0.5))
    print(scattered_intensity_transmitted(theta_0, theta_1, 0.5))
    print(scattered_intensity_transmitted(theta_0, theta_1 + 0.1, 0.5))

    print(intensity_at_ground(130.0, 131.0, 30.0, 0.5))
    print(intensity_at_ground(170.0, 170.0, 180.0, 1.0))
    print(intensity_at_ground(170.0, 157.0, 180.0, 1.0))

    theta_0 = 150.0
    tau = 0.125
    phi = 0.0

    print("\nResults")

    for i in range(0, 70, 5):
        theta = float(i)
        intensity = intensity_at_ground(theta_0, 180.0 - theta, phi, tau)
        print(theta_0, theta, phi, intensity)

    for i in range(71, 80, 1):
        theta = float(i)
        intensity = intensity_at_ground(theta_0, 180.0 - theta, phi, tau)
        print(theta_0, theta, phi, intensity)

    for i in [82, 84, 86, 88, 89]:
        theta = float(i)
        intensity = intensity_at_ground(theta_0, 180.0 - theta, phi, tau)
        print(theta_0, theta, phi, intensity)


if __name__ == '__main__':

    main()
