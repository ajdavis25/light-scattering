import numpy as np
from intensity_vs_theta_spherical import (
    plot_intensity_vs_theta_spherical
)
from fisheye_projection import (
    generate_unpolarized_fisheye_data,
    generate_polarized_fisheye_data,
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
    plot_intensity_vs_theta_spherical()

    latitude = 30.0
    longitude = 0.0
    solar_declination = 23.5 # approximate for summer solstice
    tau_atm = 0.5
    num_layers = 100
    num_zenith_points = 90
    num_azimuth_points = 180
    hour_angle = 110 # ±90° to ±120° for twilight

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
    
    print(f"Unpolarized Intensity data shape: {intensity_data.shape}")
    print(f"Unpolarized Intensity data min: {np.min(intensity_data)}, max: {np.max(intensity_data)}")

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

    print(f"Polarized Intensity data shape: {I_data.shape}")
    print(f"Polarized Intensity data min: {np.min(I_data)}, max: {np.max(I_data)}")
    print(f"Q_data min: {np.min(Q_data)}, max: {np.max(Q_data)}")
    print(f"U_data min: {np.min(U_data)}, max: {np.max(U_data)}")

    # plot the polarized fisheye twilight sky
    plot_polarized_fisheye(
        latitude=latitude,
        solar_declination=solar_declination,
        tau_atm=tau_atm,
        hour_angle=hour_angle,
        num_zenith_points=num_zenith_points,
        num_azimuth_points=num_azimuth_points,
        roll_amount=num_azimuth_points // 2
    )


if __name__ == '__main__':
    main()
