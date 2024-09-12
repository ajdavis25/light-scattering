from intensity_vs_theta_spherical import plot_intensity_vs_theta_spherical
from single_scatter_plot_spherical import generate_fisheye_data, plot_fisheye_twilight_sky

def main():
    # plot intensity vs. theta using a spherical earth model
    #plot_intensity_vs_theta_spherical()

    # example usage for fisheye twilight plot
    latitude = 30.0
    longitude = 0.0
    solar_declination = 23.5 # approximate for summer solstice
    tau_atm = 0.5
    num_layers = 100

    # adjust the hour angle to simulate twilight
    hour_angle = 120 # ±90° to ±120° for twilight

    # generate fisheye intensity data
    intensity_data = generate_fisheye_data(latitude, longitude, solar_declination, tau_atm, num_layers)

    # plot the fisheye twilight sky
    plot_fisheye_twilight_sky(intensity_data, solar_declination, hour_angle)


if __name__ == '__main__':
    main()
