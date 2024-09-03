import os, math, numpy as np, matplotlib.pyplot as plt
from plane_parallel import intensity_at_tau


def plot_polarization_vs_theta(tau_atm=0.00, theta_0=145, tau_max=0.1, save_path='./plots', save_name='polarization_vs_theta.png'):

    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # theta_obs1 intensity1 generation
    phi_obs = 0
    theta_obs1 = list()
    intensity1 = list()

    for theta in range(90, 180):
        theta_obs1.append(-(180-theta))
        intensity = intensity_at_tau(theta_0, float(theta), phi_obs, tau_atm, tau_max)
        pol = math.sqrt(math.pow(intensity[1], 2) + math.pow(intensity[2], 2)) / intensity[0]
        intensity1.append(pol)

    phi_obs = 180
    for theta in range(180, 90, -1):
        theta_obs1.append(180-theta)
        intensity = intensity_at_tau(theta_0, float(theta), phi_obs, tau_atm, tau_max)
        pol = math.sqrt(math.pow(intensity[1], 2) + math.pow(intensity[2], 2)) / intensity[0]
        intensity1.append(pol)

    # theta_obs2 intensity2 generation
    phi_obs = 45
    theta_obs2 = list()
    intensity2 = list()
    for theta in range(90, 180):
        theta_obs2.append(-(180-theta))
        intensity = intensity_at_tau(theta_0, float(theta), phi_obs, tau_atm, tau_max)
        pol = math.sqrt(math.pow(intensity[1], 2) + math.pow(intensity[2], 2)) / intensity[0]
        intensity2.append(pol)

    phi_obs = 225
    for theta in range(180, 90, -1):
        theta_obs2.append(180-theta)
        intensity = intensity_at_tau(theta_0, float(theta), phi_obs, tau_atm, tau_max)
        pol = math.sqrt(math.pow(intensity[1], 2) + math.pow(intensity[2], 2)) / intensity[0]
        intensity2.append(pol)

    # theta_obs3 intensity3 generation
    phi_obs = 90
    theta_obs3 = list()
    intensity3 = list()
    for theta in range(90, 180):
        theta_obs3.append(-(180-theta))
        intensity = intensity_at_tau(theta_0, float(theta), phi_obs, tau_atm, tau_max)
        pol = math.sqrt(math.pow(intensity[1], 2) + math.pow(intensity[2], 2)) / intensity[0]
        intensity3.append(pol)

    phi_obs = 270
    for theta in range(180, 90, -1):
        theta_obs3.append(180-theta)
        intensity = intensity_at_tau(theta_0, float(theta), phi_obs, tau_atm, tau_max)
        pol = math.sqrt(math.pow(intensity[1], 2) + math.pow(intensity[2], 2)) / intensity[0]
        intensity3.append(pol)

    plt.plot(theta_obs1, intensity1)
    plt.plot(theta_obs2, intensity2)
    plt.plot(theta_obs3, intensity3)
    plt.legend(['0-180', '45-225', '90-270'])
    plt.title('Polarization Vs. Theta')
    plt.xlabel('Theta')
    plt.ylabel('Polarization')
    plt.savefig(os.path.join(save_path, save_name))
    plt.show()
