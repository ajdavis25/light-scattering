#!/usr/bin/python3
''' This code produces a view of the sky at one wavelength. '''

import os, math, numpy as np, matplotlib.pyplot as plt
from typing import List
from matplotlib import cm
from matplotlib.patches import Circle
from plane_parallel import intensity_at_tau
# from PIL import Image, ImageDraw


def calculate_quantity(
    intensity: List[float], optical_quantity: str
) -> float:
    """
    calculate the desired optical quantity based on intensity

    args:
        intensity: list of four stokes parameters [I, Q, U, V]
        optical_quantity: optical quantity to calculate ('I', 'Pol', 'Q', 'U', 'n_lines')

    returns:
        calculated optical quantity
    """
    if optical_quantity == 'I':
        return intensity[0]
    elif optical_quantity == 'Pol':
        if intensity[0] > 0.0:
            return math.sqrt(intensity[1] ** 2 + intensity[2] ** 2) / intensity[0]
    elif optical_quantity == 'Q':
        return intensity[1]
    elif optical_quantity == 'U':
        return intensity[2]
    elif optical_quantity == 'n_lines':
        if optical_quantity[1] > 0:
            return 4.0 if intensity[2] > 0 else 1.0
        else:
            return 3.0 if intensity[2] > 0 else 2.0
    return 0.0


def make_contour_plot(
    optical_quantity: str,
    theta_0: float,
    tau_obs: float,
    tau_atm: float,
    direction: str,
    save_name: str,
    save_path: str = './plots'
) -> None:
    '''
    generate a polar contour polot of the optical quantity

    args:
        optical_quantity: the type of optical quantity to plot ('I', 'Pol', 'Q', 'U', 'n_lines')
        theta_0: solar zenity angle (in degrees)
        tau_obs: optical depty at observation level
        tau_atm: optical depth of the atmosphere
        direction: direction of the light ('downwelling' or 'upwelling')
    '''
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # generate data for the plot using linspace so that the endpoint of 360 is included
    azimuths = np.radians(np.linspace(0, 360, 181)) # 0 to 360 degrees
    zeniths = np.arange(0, 91, 2) # 0 to 90 degrees (zenith)
    r, theta = np.meshgrid(zeniths, azimuths)
    values = np.zeros((azimuths.size, zeniths.size))

    # compute the optical quantity for each (zenith, azimuth) pair
    for azimuth in range(len(azimuths)):
        for zenith in range(len(zeniths)):
            phi = math.degrees(azimuths[azimuth])
            theta_p = 180.0 - zeniths[zenith] if direction == "downwelling" else zeniths[zenith]
            intensity = intensity_at_tau(theta_0, theta_p, phi, tau_obs, tau_atm)
            values[azimuth, zenith] = calculate_quantity(intensity, optical_quantity)

    # create polar contour plot
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))

    # ax.contourf(theta, r, values, 30, cmap=cm.gray)
    ax.contourf(theta, r, values, 70, cmap=cm.bone)
    # ax.contour(theta, r, values, 40)

    ax.set_theta_zero_location("N")
    # or ax.set_theta_offset(pi)
    ax.set(xticks = np.arange(0, 2 * math.pi, math.pi / 6.0),
           yticks = np.arange(0, 91, 30) )
    
    # add the circle of the sun
    sun_position = (0, (180 - theta_0) if direction == "downwelling" else - (180 - theta_0))
    sun = Circle(sun_position, 1.5, transform=ax.transData._b, color='yellow' if direction == "downwelling" else 'red')
    ax.add_artist(sun)

    # save and display the plot
    plt.savefig(os.path.join(save_path, save_name))
    plt.title(f'Single Scatter Plot - {optical_quantity}')
    plt.show()
