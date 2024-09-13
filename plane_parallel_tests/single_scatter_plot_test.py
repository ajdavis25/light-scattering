import os, pytest, numpy as np
from plane_parallel.single_scatter_plot import (
    calculate_quantity,
    make_contour_plot
)


# test for calculate_quantity
@pytest.mark.parametrize("intensity, optical_quantity, expected", [
    ([1.0, 0.5, 0.5, 0.0], 'I', 1.0), # test intensity I
    ([1.0, 0.5, 0.5, 0.0], 'Pol', 0.70710678118), # test polarization Pol
    ([2.0, 1.0, 0.5, 0.0], 'Q', 1.0), # test Q component
    ([2.0, 1.0, 0.5, 0.0], 'U', 0.5), # test U component
    ([2.0, 1.0, 0.5, 0.0], 'n_lines', 4.0), # test n_lines
    ([0.0, 0.5, 0.5, 0.0], 'Pol', 0.0), # test polarization when intensity is zero
])
def test_calculate_quantity(intensity, optical_quantity, expected):
    result = calculate_quantity(intensity, optical_quantity)
    assert np.isclose(result, expected, atol=1e-6), f"Expected {expected}, got {result}"


# test for calculate_quantity with invalid optical_quantity
def test_calculate_quantity_invalid_input():
    intensity = [1.0, 0.5, 0.5, 0.0]
    invalid_quantity = "Invalid"

    # test that an invalid optical quantity defaults to 0.0
    result = calculate_quantity(intensity, invalid_quantity)
    assert result == 0.0, f"Expected 0.0 for invalid optical quantity, got {result}"


# test for make_contour_plot (ensuring plot is saved)
def test_make_contour_plot(tmpdir):
    save_path = tmpdir.mkdir("plots")
    save_name = "single_scatter_plot_test.png"
    
    # call the function with example parameters
    make_contour_plot(
        optical_quantity='I',
        theta_0=45.0,
        tau_obs=0.5,
        tau_atm=1.0,
        direction='downwelling',
        save_name=save_name,
        save_path=str(save_path)
    )

    # check if the plot file was created
    assert os.path.isfile(os.path.join(str(save_path), save_name)), "Plot file should have been created."


# test for make_contour_plot with upwelling light direction
def test_make_contour_plot_upwelling(tmpdir):
    save_path = tmpdir.mkdir("plots")
    save_name = "single_scatter_plot_upwelling_test.png"
    
    # call the function with upwelling direction
    make_contour_plot(
        optical_quantity='Pol',
        theta_0=60.0,
        tau_obs=0.3,
        tau_atm=1.0,
        direction='upwelling',
        save_name=save_name,
        save_path=str(save_path)
    )

    # check if the plot file was created
    assert os.path.isfile(os.path.join(str(save_path), save_name)), "Plot file should have been created for upwelling."
