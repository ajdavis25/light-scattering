import os, numpy as np, matplotlib.pyplot as plt
from typing import List, Tuple
from plane_parallel import intensity_at_ground


def generate_intensity(
    theta_0: float, phi_angles: Tuple[float, float], tau_atm: float
) -> Tuple[List[float], List[float]]:
    """
    generate the intensity of light at ground level for a range of theta angles and a given phi angle

    args:
        theta_0: initial angle of light
        phi_angles: a tuple containing two phi observation angles
        tau_atm: atmospheric tau value

    returns:
        a tuple of lists containing theta observation angles and corresponding intensity values
    """
    theta_obs = []
    intensities = []

    # first half: theta from 90 to 180 degrees
    for theta in range(90, 180):
        theta_obs.append(-(180 - theta))
        intensities.append(intensity_at_ground(theta_0, float(theta), phi_angles[0], tau_atm))

    # second half: theta from 180 to 90 degrees
    for theta in range(180, 90, -1):
        theta_obs.append(180 - theta)
        intensities.append(intensity_at_ground(theta_0, float(theta), phi_angles[1], tau_atm))

    return theta_obs, intensities


def plot_intensity_vs_theta(
    tau_atm: float = 0.5, 
    theta_0: float = 145, 
    save_path: str = './plots', 
    save_name: str = 'intensity_vs_theta.png'
) -> None:
    """
    plot intensity vs. theta for different observation angles and save the plot

    Args:
        tau_atm: the atmospheric tau value (default: 0.5)
        theta_0: initial angle of light (default: 145)
        save_path: directory to save the plot (default: './plots')
        save_name: file name for the saved plot (default: 'intensity_vs_theta.png')
    """

    # create directory if it doesn't exist
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # generate intensity data for different observation angles
    theta_obs1, intensity1 = generate_intensity(theta_0, (0, 180), tau_atm)
    theta_obs2, intensity2 = generate_intensity(theta_0, (45, 225), tau_atm)
    theta_obs3, intensity3 = generate_intensity(theta_0, (90, 270), tau_atm)

    # plot the data
    plt.plot(theta_obs1, intensity1, label='0-180')
    plt.plot(theta_obs2, intensity2, label='45-225')
    plt.plot(theta_obs3, intensity3, label='90-270')

    plt.legend()
    plt.title('Intensity Vs. Theta')
    plt.xlabel('Theta')
    plt.ylabel('Intensity')

    # save and display the plot
    plt.savefig(os.path.join(save_path, save_name))
    plt.show()


if __name__ == "__main__":
    plot_intensity_vs_theta()
