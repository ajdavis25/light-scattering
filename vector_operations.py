#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on  Jul 25 15:23:36 2024

@author: jadams

testing matrix stuff... mueller, rotation, etc
"""

import math
import numpy as np

def direction_vector(theta, phi):
    ''' Return the direction vector corresponding to a direction. '''

    theta_rad = math.radians(theta)
    phi_rad = math.radians(phi)

    x = math.sin(theta_rad)*math.cos(phi_rad)
    y = math.sin(theta_rad)*math.sin(phi_rad)
    z = math.cos(theta_rad)

    return np.array([(x, y, z)])

def dot_product(a, b):
    ''' Take the dot product of two vectors. '''

    return np.sum(a*b)

def normalize(vec_in):
    ''' return a normalized unit vector from an input np vector '''

    magnitude = np.linalg.norm(vec_in)

    vec_out = (1.0/magnitude) * vec_in

    return vec_out

def rayleigh_matrix(scattering_angle_degrees):
    ''' give the Rayleigh mueller matrix for this angle '''

    cos_ang = math.cos(math.radians(scattering_angle_degrees))

    mueller = np.zeros((4, 4))
    mueller[0,0] = 1 + math.pow(cos_ang, 2)
    mueller[1,1] = mueller[0,0]
    mueller[0,1] = -(1.0 - math.pow(cos_ang, 2))
    mueller[1,0] = mueller[0,1]
    mueller[2,2] = 2.0 * cos_ang
    mueller[3,3] = mueller[2,2]

    mueller /= (1.0 + math.pow(cos_ang, 2))

    return mueller

def rotation_matrix(rotation_angle_degrees):
    ''' give the rotation matrix for this angle '''

    cos_ang = math.cos(math.radians(rotation_angle_degrees))
    sin_ang = math.sin(math.radians(rotation_angle_degrees))

    mueller = np.zeros((4, 4))
    mueller[0,0] = 1.0 
    mueller[3,3] = 1.0
    mueller[1,1] = cos_ang * cos_ang - sin_ang * sin_ang
    mueller[2,2] = mueller[1,1]
    mueller[1,2] = 2.0 * cos_ang * sin_ang
    mueller[2,1] = -mueller[1,2]

    return mueller

def rotation_angles(k_i, k_f):
    ''' Calculate rotation angles and scattering angle for two
        direction vectors. '''
    
    # scattering angle first
    theta = math.degrees(math.acos(dot_product(k_i, k_f)))

    z = np.array([(0.0, 0.0, 1.0)])

    # calculate first rotation phi
    c_i = normalize(np.cross(k_i, z))
    t_i = normalize(np.cross(c_i, k_i))

    d_u = np.cross(np.cross(k_i, k_f), k_i)
    d_i = normalize(d_u)

    dot_dt = dot_product(d_i, t_i)
    phi = 0.0
    if abs(dot_dt) < 1.0:
        phi = math.degrees(math.acos(dot_product(d_i, t_i)))

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
        psi = math.degrees(math.acos(dot_dt))

    test_out = normalize(np.cross(d_f, t_f))
    if dot_product(test_out, k_f) > 0:
        psi = -psi

    return (phi, psi, theta)

def my_main():
    ''' do stuff here... '''

    # stokes_vector = np.array([1.0, 0.1, -0.1, 0.5])
    stokes_vector = np.array([1.0, 0.0, -0.0, 0.0])

    rotation_angle = 30.0
    scattering_angle = 30.0

    rotation = rotation_matrix(rotation_angle)
    scattering = rayleigh_matrix(scattering_angle)

    print("\nStokes Vector\n", stokes_vector)
    print("\nRotation Matrix\n", rotation)
    print("\nScattering Matrix\n", scattering)

    rotated_stokes = np.dot(rotation, stokes_vector)
    scattered_stokes = np.dot(scattering, stokes_vector)

    print("\nStokes Vector Rotated\n", rotated_stokes)
    print("\nStokes Vector Scattered\n", scattered_stokes)

    # direction vector takes theta and phi
    theta_1 = 45.0
    phi_1   = 45.0
    theta_2 = 46.0
    phi_2   = 46.0
    theta_3 = 46.0
    phi_3   = 44.0
    u_1 = direction_vector(theta_1, phi_1)
    u_2 = direction_vector(theta_2, phi_2)
    u_3 = direction_vector(theta_3, phi_3)
    print("1 ", theta_1, phi_1)
    print("2 ", theta_2, phi_2)
    print("3 ", theta_3, phi_3)
    print("angle : ", rotation_angles(u_1, u_2))
    print("angle : ", rotation_angles(u_1, u_3))


if __name__ == '__main__':

    my_main()