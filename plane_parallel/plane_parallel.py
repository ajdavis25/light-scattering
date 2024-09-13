#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 15 15:23:36 2022

refactored version: best practices applied

This code simulates the downwelling light field for an observer
standing at the base of a plane-parallel atmosphere. This is a scalar
model using a standard Rayleigh phase function for scattering.
The ground is completely absorbing.

@author: jadams, suit_bingus
"""

import math, numpy as np
from vector_operations import rayleigh_matrix, rotation_angles, rotation_matrix


# constants
PI = math.pi


def rayleigh_phase_function(psi: float) -> float:
    """
    calculates the rayleigh phase function for a given scattering angle

    args:
        psi: the scattering angle in degrees

    returns:
        the value of the rayleigh phase function at the given scattering angle
    """
    psi_rad = math.radians(psi)
    phase_function = (3.0 / (16.0 * PI)) * (1.0 + math.cos(psi_rad) ** 2)
    
    return phase_function


def scattering_angle(vector_1: np.ndarray, vector_2: np.ndarray) -> float:
    """
    calculate the scattering angle between two unit vectors

    args:
        vector_1: first 3D unit vector
        vector_2: second 3D unit vector

    returns:
        the scattering angle in degrees between the two vectors in degrees
    """
    dot_product = np.dot(vector_1, vector_2)

    # clamp the dot product to avoid floating-point inaccuracies outside [-1, 1]
    dot_product = np.clip(dot_product, -1.0, 1.0)
    psi_rad = math.acos(dot_product)

    return math.degrees(psi_rad)


def photon_unit_vector(theta: float, phi: float) -> np.ndarray:
    """
    create the direction vector for a photon heading in the theta, phi direction

    args:
        theta: angle in degrees from the z-axis
        phi: azimuthal angle in degrees in the xy-plane

    returns:
        a 3D unit vector as a NumPy array
    """
    theta_rad = math.radians(theta)
    phi_rad = math.radians(phi)

    u_x = math.sin(theta_rad) * math.cos(phi_rad)
    u_y = math.sin(theta_rad) * math.sin(phi_rad)
    u_z = math.cos(theta_rad)

    return np.array([u_x, u_y, u_z])


def scattered_intensity_reflected(theta_0: float, theta: float, tau_max: float) -> float:
    """
    calculate the scattered intensity (without phase function) at the top of the atmosphere

    args:
        theta_0: angle of incoming light
        theta: scattering angle
        tau_max: maximum optical depth

    returns:
        the scattered intensity at the top of the atmosphere
    """
    mu_0 = abs(math.cos(math.radians(theta_0)))
    mu_obs = abs(math.cos(math.radians(theta)))

    if theta_0 != theta:
        factor1 = mu_0 / (mu_0 + mu_obs)
        arg_inter = tau_max * (1.0 / mu_obs + 1.0 / mu_0)
        factor2 = 1.0 - math.exp(-arg_inter)
        return factor1 * factor2

    return tau_max * math.exp(-tau_max / mu_0) / mu_obs


def scattered_intensity_transmitted(theta_0: float, theta: float, tau_max: float) -> float:
    """
    calculate the scattered intensity (without phase function) at the bottom of the atmosphere

    args:
        theta_0: angle of incoming light
        theta: scattering angle
        tau_max: maximum optical depth

    returns:
        the scattered intensity at the bottom of the atmosphere
    """
    mu_0 = abs(math.cos(math.radians(theta_0)))
    mu_obs = abs(math.cos(math.radians(theta)))

    if theta_0 != theta:
        factor1 = mu_0 * math.exp(-tau_max / mu_0) / (mu_0 - mu_obs)
        arg_inter = tau_max * (1.0 / mu_obs - 1.0 / mu_0)
        factor2 = 1.0 - math.exp(-arg_inter)
        return factor1 * factor2

    return tau_max * math.exp(-tau_max / mu_0) / mu_obs


def intensity_at_ground(theta_0: float, theta: float, phi: float, tau_max: float) -> float:
    """
    calculate the scattered intensity at the bottom of the atmosphere

    args:
        theta_0: angle of the incoming light
        theta: scattering angle
        phi: azimuthal angle
        tau_max: maximum optical depth

    returns:
        the intensity of scattered light at the ground level
    """
    u_0 = photon_unit_vector(theta_0, 0.0)
    u_1 = photon_unit_vector(theta, phi)

    scatter_angle_value = scattering_angle(u_0, u_1)
    base_intensity = scattered_intensity_transmitted(theta_0, theta, tau_max) * \
                     rayleigh_phase_function(scatter_angle_value)

    return base_intensity


def intensity_at_top(theta_0: float, theta: float, phi: float, tau_max: float) -> float:
    """
    calculate the scattered intensity at the top of the atmosphere

    args:
        theta_0: angle of the incoming light
        theta: scattering angle
        phi: azimuthal angle
        tau_max: maximum optical depth

    returns:
        the intensity of scattered light at the top of the atmosphere
    """
    u_0 = photon_unit_vector(theta_0, 0.0)
    u_1 = photon_unit_vector(theta, phi)

    scatter_angle_value = scattering_angle(u_0, u_1)
    top_intensity = scattered_intensity_reflected(theta_0, theta, tau_max) * \
                    rayleigh_phase_function(scatter_angle_value)

    return top_intensity


def scattered_intensity_at_tau(theta_0: float, theta_obs: float, tau: float, tau_max: float) -> float:
    """
    calculate the scattered intensity (without phase function) at a given altitude tau

    args:
        theta_0: angle of the incoming light
        theta_obs: observation angle
        tau: current optical depth
        tau_max: maximum optical depth

    returns:
        the scattered intensity at a given altitude
    """
    mu_0 = abs(math.cos(math.radians(theta_0)))
    mu_obs = abs(math.cos(math.radians(theta_obs)))

    factor0 = math.exp(-(tau_max - tau) / mu_0)

    if theta_obs < 90.0:
        factor1 = mu_0 / (mu_0 + mu_obs)
        arg_inter = tau * (1.0 / mu_obs + 1.0 / mu_0)
        factor2 = 1.0 - math.exp(-arg_inter)
        return factor0 * factor1 * factor2

    if theta_0 != theta_obs:
        factor1 = mu_0 * math.exp(-tau / mu_0) / (mu_0 - mu_obs)
        arg_inter = (tau_max - tau) * (1.0 / mu_obs - 1.0 / mu_0)
        factor2 = 1.0 - math.exp(-arg_inter)
        return factor0 * factor1 * factor2

    return (tau_max - tau) * math.exp(-(tau_max - tau) / mu_obs) / mu_obs


def intensity_at_tau(theta_0: float, theta_obs: float, phi_obs: float, tau: float, tau_max: float) -> np.ndarray:
    """
    calculate the scattered intensity of light at a given altitude tau

    args:
        theta_0: angle of the incoming light
        theta_obs: observation angle
        phi_obs: azimuthal observation angle
        tau: current optical depth
        tau_max: maximum optical depth

    returns:
        the scattered intensity vector of light at a given altitude
    """
    u_0 = photon_unit_vector(theta_0, 0.0)
    u_1 = photon_unit_vector(theta_obs, phi_obs)

    # calculate rotation angles
    if phi_obs == 0.0 and theta_0 == theta_obs:
        phi, psi, theta = 0.0, 0.0, 0.0
    else:
        phi, psi, theta = rotation_angles(u_0, u_1)

    # scattering intensity at specific optical depth
    intensity = scattered_intensity_at_tau(theta_0, theta_obs, tau, tau_max) * \
                rayleigh_phase_function(theta)

    # mueller matrix and rotation matrix calculations for stokes parameters
    mueller_scatter = rayleigh_matrix(theta)
    mueller_rotation = rotation_matrix(psi)
    stokes_init = np.array([1.0, 0.0, 0.0, 0.0])
    scattered_stokes = np.dot(mueller_rotation, np.dot(mueller_scatter, stokes_init))

    # add contributions from both the top and bottom of the atmosphere
    if tau == 0: # top of the atmosphere
        return intensity_at_top(theta_0, theta_obs, phi_obs, tau_max) * scattered_stokes
    if tau == tau_max: # bottom of the atmosphere
        return intensity_at_ground(theta_0, theta_obs, phi_obs, tau_max) * scattered_stokes

    # general case for in-between tau
    return intensity * scattered_stokes
