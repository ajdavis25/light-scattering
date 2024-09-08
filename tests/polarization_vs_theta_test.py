import os, pytest, numpy as np
from src.polarization_vs_theta import (
    calculate_polarization,
    generate_polarization_data,
    plot_polarization_vs_theta
)


# test for calculate_polarization
@pytest.mark.parametrize("intensity, expected", [
    ([1.0, 0.5, 0.5, 0.0], 0.70710678118),
    ([2.0, 1.0, 1.0, 0.0], 0.70710678118),
    ([1.0, 0.0, 0.0, 0.0], 0.0), # no polarization
    ([0.0, 0.0, 0.0, 0.0], 0.0) # division by zero case
])
def test_calculate_polarization(intensity, expected):
    result = calculate_polarization(intensity)
    assert np.isclose(result, expected, atol=1e-6), f"Expected {expected}, got {result}"


# test for calculate_polarization with invalid input
def test_calculate_polarization_invalid_input():
    with pytest.raises(ValueError, match="Intensity list must contain 4 elements"):
        calculate_polarization([1.0, 0.5, 0.5]) # only 3 elements instead of 4


# test for generate_polarization_data
def test_generate_polarization_data():
    theta_0 = 145
    phi_obs_angles = (0, 180)
    tau_atm = 0.5
    tau_max = 1.0

    # generate polarization data
    theta_obs, polarization = generate_polarization_data(theta_0, phi_obs_angles, tau_atm, tau_max)

    # ensure the lists are of equal length
    assert len(theta_obs) == len(polarization), "Theta observation and polarization lists must have the same length."

    # check that the lists are not empty
    assert len(theta_obs) > 0, "Theta observation list should not be empty."
    assert len(polarization) > 0, "Polarization list should not be empty."

    # check that the polarization values are valid numbers
    assert all(isinstance(x, (float, int)) for x in polarization), "Polarization values should be numerical."


# test for plot_polarization_vs_theta (ensuring plot is saved)
def test_plot_polarization_vs_theta(tmpdir):
    save_path = tmpdir.mkdir("plots")
    save_name = "polarization_vs_theta_test.png"

    # call the plot function
    plot_polarization_vs_theta(tau_atm=0.5, theta_0=145, tau_max=1.0, save_path=str(save_path), save_name=save_name)

    # check if the plot file was created
    assert os.path.isfile(os.path.join(str(save_path), save_name)), "Plot file should have been created."
