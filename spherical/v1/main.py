import numpy as np
from spherical.v1.energy import (
    calculate_energy_loss
)
from spherical_model import (
    calculate_sun_position
)
from spherical.v1.fisheye_projection import (
    generate_unpolarized_fisheye_data,
    generate_polarized_fisheye_data
)
from spherical.v1.plot_utils import (
    plot_intensity_vs_theta_spherical,
    plot_unpolarized_fisheye,
    plot_polarized_fisheye
)


def main():
    """
    the main function that orchestrates the simulation of twilight sky intensity and polarization

    this script performs the following steps:

    1. plots a graph of intensity vs. theta for the spherical earth model
    2. generates unpolarized fisheye data and plots it
    3. generates polarized fisheye data (including I, Q, U components) and plots it

    parameters:
        latitude: the latitude of the observation point in degrees
        longitude: the longitude of the observation point in degrees
        solar_declination: the solar declination in degrees
        tau_atm: the atmospheric optical depth
        num_layers: the number of atmospheric layers to consider in the calculation
        num_zenith_points: the number of points to generate along the zenith angle axis
        num_azimuth_points: the number of points to generate along the azimuth angle axis
        hour_angle: the hour angle of the sun in degrees
    """
    # plot_intensity_vs_theta_spherical()

    latitude = 45.0
    longitude = 0.0
    solar_declination = 23.5 # approximate for summer solstice
    tau_atm = 0.2
    num_layers = 50
    num_zenith_points = 90
    num_azimuth_points = 180
    hour_angle = 90 # ±90° to ±120° for twilight

    # generate fisheye intensity data (unpolarized)
    intensity_data = generate_unpolarized_fisheye_data(
        latitude=latitude,
        longitude=longitude,
        solar_declination=solar_declination,
        tau_atm=tau_atm,
        num_layers=num_layers,
        hour_angle=hour_angle,
        num_zenith_points=num_zenith_points,
        num_azimuth_points=num_azimuth_points
    )

    print("--------------------")
    print(f"Unpolarized Intensity data shape: {intensity_data.shape}")
    print(f"Unpolarized Intensity data min: {np.min(intensity_data)}, max: {np.max(intensity_data)}")
    print("--------------------")

    # calculate energy loss for unpolarized data
    zenith_angles = np.linspace(0, 90, num_zenith_points)
    azimuth_angles = np.linspace(0, 2 * np.pi, num_azimuth_points)
    sun_zenith, _ = calculate_sun_position(latitude, solar_declination, hour_angle)

    E_in, E_out, energy_loss_ratio = calculate_energy_loss(
        intensities=intensity_data,
        sun_zenith=sun_zenith,
        zenith_angles=zenith_angles,
        azimuth_angles=azimuth_angles
    )
    print("--------------------")
    print(f"Energy Loss for Unpolarized Fisheye: E_in={E_in}, E_out={E_out}, Loss Ratio={energy_loss_ratio}")
    print("--------------------")

    # plot the unpolarized fisheye twilight sky
    plot_unpolarized_fisheye(
        intensities=intensity_data,
        latitude=latitude,
        solar_declination=solar_declination,
        hour_angle=hour_angle
    )

    # generate fisheye intensity data with polarization (polarized)
    I_data, Q_data, U_data = generate_polarized_fisheye_data(
        latitude=latitude,
        longitude=longitude,
        solar_declination=solar_declination,
        tau_atm=tau_atm,
        num_layers=num_layers,
        hour_angle=hour_angle,
        num_zenith_points=num_zenith_points,
        num_azimuth_points=num_azimuth_points
    )

    print("--------------------")
    print(f"Polarized Intensity data shape: {I_data.shape}")
    print(f"Polarized Intensity data min: {np.min(I_data)}, max: {np.max(I_data)}")
    print(f"Q_data min: {np.min(Q_data)}, max: {np.max(Q_data)}")
    print(f"U_data min: {np.min(U_data)}, max: {np.max(U_data)}")
    print("--------------------")

    # plot the polarized fisheye twilight sky
    plot_polarized_fisheye(
        longitude=longitude,
        latitude=latitude,
        solar_declination=solar_declination,
        tau_atm=tau_atm,
        num_layers=num_layers,
        hour_angle=hour_angle,
        num_zenith_points=num_zenith_points,
        num_azimuth_points=num_azimuth_points
    )


if __name__ == '__main__':
    main()
