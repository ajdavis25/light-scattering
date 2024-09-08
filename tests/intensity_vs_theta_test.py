import os, pytest
from src.intensity_vs_theta import (
    generate_intensity,
    plot_intensity_vs_theta
)


def test_generate_intensity():
    """
    test the generate_intensity function to ensure it returns the correct output structure
    """
    theta_0 = 145
    phi_angles = (0, 180)
    tau_atm = 0.5

    theta_obs, intensities = generate_intensity(theta_0, phi_angles, tau_atm)

    # check if the function returns lists
    assert isinstance(theta_obs, list), "Theta observations should be a list."
    assert isinstance(intensities, list), "Intensities should be a list."

    # check if both lists have the same length
    assert len(theta_obs) == len(intensities), "Theta observations and intensities should have the same length."

    # check for expected number of entries
    assert len(theta_obs) == 179, "Expected 179 theta observations."

    # check values in the list
    assert all(isinstance(x, (float, int)) for x in theta_obs), "Theta observations should contain numerical values."
    assert all(isinstance(x, (float, int)) for x in intensities), "Intensities should contain numerical values."

    # optional: check values more rigorously if needed (e.g., ranges of theta values)
    assert theta_obs[0] == -90, "The first theta value should be -90 degrees."
    assert theta_obs[-1] == 90, "The last theta value should be 90 degrees."


@pytest.mark.parametrize("tau_atm, theta_0", [(0.5, 145), (0.25, 130), (0.75, 160)])
def test_plot_intensity_vs_theta(tmpdir, tau_atm, theta_0):
    """
    test the plot_intensity_vs_theta function to check if the plot is correctly saved
    """
    save_path = tmpdir.mkdir("plots")
    save_name = "test_intensity_vs_theta.png"
    
    # call the plot function
    plot_intensity_vs_theta(tau_atm=tau_atm, theta_0=theta_0, save_path=str(save_path), save_name=save_name)
    
    # check if the file is created
    assert os.path.isfile(os.path.join(str(save_path), save_name)), "Plot file should have been created."


@pytest.mark.parametrize("theta_0, phi_angles, tau_atm", [
    (145, (0, 180), 0.5),
    (130, (45, 225), 0.25),
    (160, (90, 270), 0.75)
])
def test_generate_intensity_various_cases(theta_0, phi_angles, tau_atm):
    """
    test the generate_intensity function with various parameter inputs
    """
    theta_obs, intensities = generate_intensity(theta_0, phi_angles, tau_atm)

    # basic checks
    assert isinstance(theta_obs, list), "Theta observations should be a list."
    assert isinstance(intensities, list), "Intensities should be a list."
    assert len(theta_obs) == len(intensities), "Theta observations and intensities should have the same length."
    
    # check the first and last elements of theta
    assert theta_obs[0] == -90, "The first theta value should be -90 degrees."
    assert theta_obs[-1] == 90, "The last theta value should be 90 degrees."
