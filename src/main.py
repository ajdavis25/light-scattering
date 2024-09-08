import numpy as np
from intensity_vs_theta import plot_intensity_vs_theta
from polarization_vs_theta import plot_polarization_vs_theta
from plane_parallel import intensity_at_ground
from single_scatter_plot import make_contour_plot
from vector_operations import (
    rotation_matrix,
    rayleigh_matrix,
    direction_vector,
    rotation_angles
)


def main():
    """
    the main function demonstrates the functionalities of various imported modules

    this script performs the following tasks:

    1. calculates and prints scattered intensity at the ground for a range of angles
    2. demonstrates rotation and scattering of a stokes vector using matrices
    3. calculates direction vectors based on theta and phi angles
    4. calculates rotation angles between different direction vectors
    5. calls external functions to plot intensity and polarization vs. theta, and generates contour plots
    """
    # define simulation parameters
    theta_0 = 130.0 # initial scattering angle (degrees)
    tau = 0.5 # optical depth
    phi = 0.0 # azimuthal angle (degrees)

    # calculate and print scattered intensity at ground for various theta angles
    print("Test scattered intensity at ground:")
    for theta in np.linspace(0, 89, 10):
        intensity = intensity_at_ground(theta_0, 180.0 - theta, phi, tau)
        print(f"Theta_0: {theta_0}, Theta: {theta}, Phi: {phi}, Intensity: {intensity}")

    # define a sample Stokes vector and rotation/scattering angles
    stokes_vector = np.array([1.0, 0.0, 0.0, 0.0])
    rotation_angle = 30.0
    scattering_angle = 30.0

    # calculate rotation and scattering matrices
    rotation = rotation_matrix(rotation_angle)
    scattering = rayleigh_matrix(scattering_angle)

    print("\nStokes Vector\n", stokes_vector)
    print("\nRotation Matrix\n", rotation)
    print("\nScattering Matrix\n", scattering)

    rotated_stokes = np.dot(rotation, stokes_vector)
    scattered_stokes = np.dot(scattering, stokes_vector)

    print("\nStokes Vector Rotated\n", rotated_stokes)
    print("\nStokes Vector Scattered\n", scattered_stokes)

    # test direction vectors and rotation angles
    theta_1, phi_1 = 45.0, 45.0
    theta_2, phi_2 = 46.0, 46.0
    theta_3, phi_3 = 46.0, 44.0

    u_1 = direction_vector(theta_1, phi_1)
    u_2 = direction_vector(theta_2, phi_2)
    u_3 = direction_vector(theta_3, phi_3)

    print(f"\nDirection Vector 1: Theta = {theta_1}, Phi = {phi_1}")
    print(f"Direction Vector 2: Theta = {theta_2}, Phi = {phi_2}")
    print(f"Direction Vector 3: Theta = {theta_3}, Phi = {phi_3}")
    
    print("\nRotation Angles between Vector 1 and 2: ", rotation_angles(u_1, u_2))
    print("Rotation Angles between Vector 1 and 3: ", rotation_angles(u_1, u_3))

    plot_intensity_vs_theta()
    plot_polarization_vs_theta()

    # call make_contour_plot with relevant parameters
    tau_obs = 0.0
    tau_atm = 0.5
    solar_zenith_angle = 35.0
    theta_0 = 180.0 - solar_zenith_angle

    # make_contour_plot('I', theta_0, tau_obs, tau_atm, "downwelling")
    # make_contour_plot('I', theta_0, tau_obs, tau_atm, "upwelling")
    # make_contour_plot('Pol', theta_0, tau_obs, tau_atm, "downwelling")
    # make_contour_plot('Pol', theta_0, tau_atm, tau_atm, "upwelling")
    # make_contour_plot('n_lines', theta_0, tau_obs, tau_atm, "downwelling")

    make_contour_plot('I', theta_0, tau_obs, tau_atm, "downwelling", save_name='single_scatter_plot_I')
    make_contour_plot('Pol', theta_0, tau_obs, tau_atm, "downwelling", save_name='single_scatter_plot_Pol')
    # make_contour_plot('Q', theta_0, tau_obs, tau_atm, "downwelling")
    # make_contour_plot('U', theta_0, tau_obs, tau_atm, "downwelling")
    # make_contour_plot('n_lines', theta_0, tau_obs, tau_atm, "downwelling")

    # print(intensity_at_tau(135.0, 120.0, 30.0, 0.0, 0.125))
    # print(intensity_at_tau(135.0, 120.0, -30.0, 0.0, 0.125))


if __name__ == '__main__':
    main()
