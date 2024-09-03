# this will generate a plot showing the intensity of light vs the theta angle which the observer looks into the sky


import os, numpy as np, matplotlib.pyplot as plt
from plane_parallel import intensity_at_ground


def plot_intensity_vs_theta(tau_atm=0.5, theta_0=145, save_path='./plots', save_name='intensity_vs_theta.png'):

    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    # theta_obs1 intensity1 generation
    phi = 0
    theta_obs1 = list()
    intensity1 = list()
    for theta in range (90, 180):
        theta_obs1.append(-(180-theta))
        intensity1.append(intensity_at_ground(theta_0, float(theta), phi, tau_atm))
    phi = 180
    for theta in range (180, 90, -1):
        theta_obs1.append(180-theta)
        intensity1.append(intensity_at_ground(theta_0, float(theta), phi, tau_atm))
    
    # theta_obs2 intensity2 generation
    phi = 45
    theta_obs2 = list()
    intensity2 = list()
    for theta in range (90,180):
        theta_obs2.append(-(180-theta))
        intensity2.append(intensity_at_ground(theta_0, float(theta), phi, tau_atm))
    phi = 225
    for theta in range (180, 90, -1):
        theta_obs2.append(180-theta)
        intensity2.append(intensity_at_ground(theta_0, float(theta), phi, tau_atm))

    # theta_obs3 intensity3 generation
    phi = 90
    theta_obs3 = list()
    intensity3 = list()
    for theta in range (90, 180):
        theta_obs3.append(-(180-theta))
        intensity3.append(intensity_at_ground(theta_0, float(theta), phi, tau_atm))
    phi = 270
    for theta in range (180, 90, -1):
        theta_obs3.append(180-theta)
        intensity3.append(intensity_at_ground(theta_0, float(theta), phi, tau_atm))    

    plt.plot(theta_obs1,intensity1)
    plt.plot(theta_obs2,intensity2)
    plt.plot(theta_obs3,intensity3)
    plt.legend(['0-180','45-225','90-270'])
    plt.title('Intensity Vs. Theta')
    plt.xlabel('Theta')
    plt.ylabel('Intensity')
    plt.draw()
    plt.savefig(os.path.join(save_path, save_name))
    plt.show()
