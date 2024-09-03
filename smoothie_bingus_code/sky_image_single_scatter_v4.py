#!/usr/bin/python3
''' This code produces a view of the sky at one wavelength. '''

import math
# from PIL import Image, ImageDraw
from plane_parallel_v3 import intensity_at_tau
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
# from matplotlib.patches import Circle

def make_contour_plot(optical_quantity, theta_0, tau_obs, tau_atm, direction):
    ''' use matplotlib contour capability '''

    #-- Generate Data -----------------------------------------
    # Using linspace so that the endpoint of 360 is included...
    azimuths = np.radians(np.linspace(0, 360, 181))
    zeniths = np.arange(0, 91, 2)
    r, theta = np.meshgrid(zeniths, azimuths)
    values = np.random.random((azimuths.size, zeniths.size))


    for i in range(len(azimuths)):
        for j in range(len(zeniths)):
            phi = math.degrees(azimuths[i])
            if direction == "downwelling":
                theta_p = 180.0 - zeniths[j]
            else:
                theta_p = zeniths[j]
            quantity = 0.0
            intensity = intensity_at_tau(theta_0, theta_p, phi, tau_obs, tau_atm)
            if optical_quantity == 'I':
                quantity = intensity[0]
            elif optical_quantity == 'Pol':
                if intensity[0] > 0.0:
                    quantity = math.sqrt(math.pow(intensity[1], 2) + math.pow(intensity[2], 2))/intensity[0]
            elif optical_quantity == 'Q':
                quantity = intensity[1]
            elif optical_quantity == 'U':
                quantity = intensity[2]
            elif optical_quantity == 'n_lines':
                if intensity[1] > 0:
                    if intensity[2] > 0:
                        quantity = 4.0
                    else:
                        quantity = 1.0
                else:
                    if intensity[2] > 0:
                        quantity = 3.0
                    else:
                        quantity = 2.0
            values[i,j] = quantity

    #-- Plot... ------------------------------------------------
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))

    # ax.contourf(theta, r, values, 30, cmap=cm.gray)
    ax.contourf(theta, r, values, 70, cmap=cm.bone)
    # ax.contour(theta, r, values, 40)

    ax.set_theta_zero_location("N")
    # or ax.set_theta_offset(pi)
    ax.set(xticks = np.arange(0, 2*math.pi, math.pi/6.0),
           yticks = np.arange(0, 91, 30) )
    
    # Try to add the circle of the sun
    sun = plt.Circle((0, (180 - theta_0)), 1.5, transform=ax.transData._b, color='yellow')
    if direction == 'upwelling':
        sun = plt.Circle((0, -(180 - theta_0)), 1.5, transform=ax.transData._b, color='red')
    ax.add_artist(sun)

    plt.show()
    # plt.savefig('sky_intensity.png', dpi = 200)

def my_main():
    ''' Main function '''

    tau_atm = 0.5
    tau_obs = 0.00
    solar_zenith_angle = 35.0

    theta_0 = 180.0 - solar_zenith_angle

    # make_contour_plot('I', theta_0, tau_obs, tau_atm, "downwelling")
    # make_contour_plot('I', theta_0, tau_obs, tau_atm, "upwelling")
    # make_contour_plot('Pol', theta_0, tau_obs, tau_atm, "downwelling")
    # make_contour_plot('Pol', theta_0, tau_atm, tau_atm, "upwelling")
    # make_contour_plot('n_lines', theta_0, tau_obs, tau_atm, "downwelling")

    make_contour_plot('I', theta_0, tau_obs, tau_atm, "downwelling")
    make_contour_plot('Pol', theta_0, tau_obs, tau_atm, "downwelling")
    # make_contour_plot('Q', theta_0, tau_obs, tau_atm, "downwelling")
    # make_contour_plot('U', theta_0, tau_obs, tau_atm, "downwelling")
    # make_contour_plot('n_lines', theta_0, tau_obs, tau_atm, "downwelling")


    # print(intensity_at_tau(135.0, 120.0, 30.0, 0.0, 0.125))
    # print(intensity_at_tau(135.0, 120.0, -30.0, 0.0, 0.125))

    return 0

if __name__ == '__main__':

    my_main()