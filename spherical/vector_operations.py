#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jul 25, 2024

refactored version: best practices applied

This script tests matrix operations including Mueller and rotation matrices 
in the context of light scattering simulations.

@author: jadams, suit_bingus
"""

import math, numpy as np
from typing import Tuple


def direction_vector(theta: float, phi: float) -> np.ndarray:
    """
    return the direction vector corresponding to spherical coordinates

    args:
        theta: the polar angle in degrees
        phi: the azimuthal angle in degrees

    returns:
        a 3D unit direction vector as a NumPy array
    """
    theta_rad = math.radians(theta)
    phi_rad = math.radians(phi)

    x = math.sin(theta_rad) * math.cos(phi_rad)
    y = math.sin(theta_rad) * math.sin(phi_rad)
    z = math.cos(theta_rad)

    return np.array([x, y, z])


def dot_product(a, b):
    ''' Take the dot product of two vectors. '''
    
    # Flatten the vectors to ensure they are compatible
    a = np.ravel(a)
    b = np.ravel(b)
    
    return np.dot(a, b)


def normalize(vec_in: np.ndarray) -> np.ndarray:
    """
    normalize an input vector

    args:
        vec_in: input vector

    returns:
        a unit vector in the same direction as vec_in
    """
    magnitude = np.linalg.norm(vec_in)
    return vec_in / magnitude if magnitude != 0 else vec_in


def rayleigh_matrix(scattering_angle_degrees: float) -> np.ndarray:
    """
    calculate the rayleigh mueller matrix for a given scattering angle

    args:
        scattering_angle_degrees: scattering angle in degrees

    returns:
        the 4x4 rayleigh mueller matrix as a NumPy array
    """
    cos_ang = math.cos(math.radians(scattering_angle_degrees))

    mueller = np.zeros((4, 4))
    mueller[0, 0] = 1 + cos_ang ** 2
    mueller[1, 1] = mueller[0, 0]
    mueller[0, 1] = -(1.0 - cos_ang ** 2)
    mueller[1, 0] = mueller[0, 1]
    mueller[2, 2] = 2.0 * cos_ang
    mueller[3, 3] = mueller[2, 2]

    return mueller / (1.0 + cos_ang ** 2)


def rotation_matrix(rotation_angle_degrees: float) -> np.ndarray:
    """
    calculate the rotation matrix for a given rotation angle

    args:
        rotation_angle_degrees: rotation angle in degrees

    returns:
        the 4x4 rotation matrix as a NumPy array
    """
    cos_ang = math.cos(math.radians(rotation_angle_degrees))
    sin_ang = math.sin(math.radians(rotation_angle_degrees))

    mueller = np.zeros((4, 4))
    mueller[0, 0] = 1.0
    mueller[3, 3] = 1.0
    mueller[1, 1] = cos_ang ** 2 - sin_ang ** 2
    mueller[2, 2] = mueller[1, 1]
    mueller[1, 2] = 2.0 * cos_ang * sin_ang
    mueller[2, 1] = -mueller[1, 2]

    return mueller


def rotation_angles(k_i, k_f):
    ''' Calculate rotation angles and scattering angle for two
        direction vectors. '''
    
    # scattering angle first
    theta = math.degrees(math.acos(np.clip(dot_product(k_i, k_f), -1.0, 1.0)))

    z = np.array([(0.0, 0.0, 1.0)])

    # calculate first rotation phi
    c_i = normalize(np.cross(k_i, z))
    t_i = normalize(np.cross(c_i, k_i))

    d_u = np.cross(np.cross(k_i, k_f), k_i)
    d_i = normalize(d_u)

    dot_dt = dot_product(d_i, t_i)
    phi = 0.0
    if abs(dot_dt) < 1.0:
        phi = math.degrees(math.acos(np.clip(dot_product(d_i, t_i), -1.0, 1.0)))

    test_out = normalize(np.cross(d_u, t_i))
    if dot_product(test_out, k_i) < 0:
        phi = -phi

    # calculate second rotation psi
    c_f = normalize(np.cross(k_f, z))
    t_f = normalize(np.cross(c_f, k_f))

    d_u = -np.cross(np.cross(k_f, k_i), k_f)
    d_f = normalize(d_u)

    dot_dt = dot_product(d_f, t_f)
    psi = 0.0
    if abs(dot_dt) > 1.0:
        if dot_dt > 1.0:
            psi = 180.0
        else:
            psi = -180.0
    else:
        psi = math.degrees(math.acos(np.clip(dot_dt, -1.0, 1.0)))

    test_out = normalize(np.cross(d_f, t_f))
    if dot_product(test_out, k_f) > 0:
        psi = -psi

    return (phi, psi, theta)
