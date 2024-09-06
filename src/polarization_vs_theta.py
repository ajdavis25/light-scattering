import os, math, matplotlib.pyplot as plt
from typing import List, Tuple
from plane_parallel import intensity_at_tau


def calculate_polarization(intensity: List[float]) -> float:
    """
    calculates the degree of polarization from a list of stokes parameters

    args:
        intensity: a list of four stokes parameters (I, Q, U, V)

    returns:
        the degree of polarization as a float
    """

    # ensure the intensity list has four elements (I, Q, U, V)
    if len(intensity) != 4:
        raise ValueError("Intensity list must contain 4 elements (I, Q, U, V)")

    # calculate the degree of polarization using the formula
    polarization = math.sqrt(intensity[1] ** 2 + intensity[2] ** 2) / intensity[0]

    # check for division by zero and return 0 if necessary
    if intensity[0] == 0:
        return 0

    return polarization


def generate_polarization_data(
    theta_0: float, phi_obs_angles: Tuple[int, int], tau_atm: float, tau_max: float
) -> Tuple[List[float], List[float]]:
    """
    generates the polarization data for a range of observation angles (theta)
    
    args:
        theta_0: initial angle of light (degrees)
        phi_obs_angles: tuple of two phi observation angles
        tau_atm: atmospheric tau value
        tau_max: maximum tau value

    returns:
        a tuple containing a list of theta values and corresponding polarization values
    """
    theta_obs = []
    polarization = []

    # first half: theta from 90 to 180 degrees
    for theta in range(90, 180):
        theta_obs.append(-(180 - theta))
        intensity = intensity_at_tau(theta_0, float(theta), phi_obs_angles[0], tau_atm, tau_max)
        try:
            pol = calculate_polarization(intensity)
        except ValueError as e:
            print(f"Error calculating polarization at theta {theta}: {e}")
            pol = 0 # handle invalid polarization data gracefully
        polarization.append(pol)

    # second half: theta from 180 to 90 degrees
    for theta in range(180, 90, -1):
        theta_obs.append(180 - theta)
        intensity = intensity_at_tau(theta_0, float(theta), phi_obs_angles[1], tau_atm, tau_max)
        try:
            pol = calculate_polarization(intensity)
        except ValueError as e:
            print(f"Error calculating polarization at theta {theta}: {e}")
            pol = 0 # handle invalid polarization data gracefully
        polarization.append(pol)

    return theta_obs, polarization


def plot_polarization_vs_theta(
    tau_atm: float = 0.00, 
    theta_0: float = 145, 
    tau_max: float = 0.1, 
    save_path: str = './plots', 
    save_name: str = 'polarization_vs_theta.png'
) -> None:
    """
    plots polarization versus theta for different observation angles and saves the plot

    args:
        tau_atm: atmospheric tau value (default: 0.00)
        theta_0: initial angle of light (default: 145)
        tau_max: maximum tau value (default: 0.1)
        save_path: directory to save the plot (default: './plots')
        save_name: filename for the saved plot (default: 'polarization_vs_theta.png')
    """

    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # generate polarization data for different observation angles
    theta_obs1, polarization1 = generate_polarization_data(theta_0, (0, 180), tau_atm, tau_max)
    theta_obs2, polarization2 = generate_polarization_data(theta_0, (45, 225), tau_atm, tau_max)
    theta_obs3, polarization3 = generate_polarization_data(theta_0, (90, 270), tau_atm, tau_max)

    # plot the data
    plt.plot(theta_obs1, polarization1, label='0-180')
    plt.plot(theta_obs2, polarization2, label='45-225')
    plt.plot(theta_obs3, polarization3, label='90-270')

    # add labels and title
    plt.legend()
    plt.title('Polarization Vs. Theta')
    plt.xlabel('Theta')
    plt.ylabel('Polarization')

    # save and display the plot
    plt.savefig(os.path.join(save_path, save_name))
    plt.show()
