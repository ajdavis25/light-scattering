#!/usr/bin/env python3
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from utils import create_twilight_colormap


############################################################
#                     PHYSICS CONSTANTS                    #
############################################################


# earth radius (m)
R_EARTH = 6.371e6  
# top of atmosphere, say 100 km above earth's surface
ATM_HEIGHT = 1.0e5  
R_TOA = R_EARTH + ATM_HEIGHT

# scale height for exponential atmosphere (m) - typical ~8 km for rayleigh
H_SCALE = 8.0e3

# reference scattering coefficient at ground level (arbitrary example)
# this would be something like sigma0 ~ n0 * sigma_Rayleigh
SIGMA_0 = 1e-5  

# simple rayleigh phase function amplitude factor:
#   for unpolarized single scattering, 
#   phase function ~ (3/4)*(1 + cos^2(scatter_angle))
def phase_function_rayleigh(scatter_angle):
    return 0.75 * (1.0 + math.cos(scatter_angle)**2)


############################################################
#               SINGLE SCATTERING INTEGRATOR              #
############################################################


def atmospheric_density(r):
    """
    simple exponential atmospheric density as a function of radial distance r
    r is distance from Earth's center, so altitude = r - R_EARTH
    return dimensionless or relative density, used for scattering
    """
    altitude = r - R_EARTH
    if altitude < 0.0:
        return 0.0
    return math.exp(-altitude / H_SCALE)

def scattering_coefficient(r):
    """
    return the local scattering coefficient at radius r
    we choose sigma(r) = SIGMA_0 * exp(- (r - R_EARTH)/H_SCALE)
    """
    return SIGMA_0 * atmospheric_density(r)

def intensity_spherical_single_scatter(
    zenith_view, azimuth_view,
    zenith_sun, azimuth_sun,
    step_size=1000.0
    ):
    """
    compute single-scattered radiance along a line of sight given by
      (zenith_view, azimuth_view)
    from an observer on the ground (r=R_EARTH) out to r=R_TOA
    
    for each small step:
      1) find the local scattering coefficient sigma(r)
      2) compute the scattering angle between view direction and sun direction
      3) dI = I_sun * sigma(r)*density(r)*phase(scatter_angle)* e^{-tau_sun} * e^{-tau_view} ds
         (in a more detailed derivation you would have geometry factors, 
          but here we keep it very simple)
    
    note: this is a toy example. real codes do more precise geometry 
          (e.g., how the sun ray enters the atmosphere)
    """

    # convert angles to radians
    zv = math.radians(zenith_view)
    av = math.radians(azimuth_view)
    zs = math.radians(zenith_sun)
    as_ = math.radians(azimuth_sun)

    # unit vector in the viewing direction (in spherical coords)
    # we'll define spherical coordinates:
    #   theta = zenith angle measured from +Z
    #   phi   = azimuth angle measured from +X in XY plane
    #   x = sin(theta)*cos(phi), y = sin(theta)*sin(phi), z = cos(theta)
    vx = math.sin(zv)*math.cos(av)
    vy = math.sin(zv)*math.sin(av)
    vz = math.cos(zv)
    
    # unit vector in the sun direction
    sx = math.sin(zs)*math.cos(as_)
    sy = math.sin(zs)*math.sin(as_)
    sz = math.cos(zs)

    # we will integrate from r=R_EARTH to r=R_TOA
    # parametrize the line-of-sight by r, and find position in 3D space
    # the direction is fixed by (vx, vy, vz)
    # we do a simple step integration from r0 -> r1 in increments of step_size
    r0 = R_EARTH
    r1 = R_TOA
    n_steps = int((r1 - r0)//step_size) + 1
    dr = (r1 - r0) / n_steps

    # we'll accumulate scattered intensity along the path
    # for demonstration, assume the direct solar intensity "I_sun" is 1.0 
    # at the top of atmosphere, ignoring the actual solar flux
    # you might want a more physically correct value for actual brightness
    total_intensity = 0.0

    # the optical depth factors for extinction along the line of sight
    # and along the solar path can be included with exp(-tau_view)*exp(-tau_sun)
    # this example has a simplistic approach: 
    # we assume no extinction except at the scattering event 
    # (typical single-scatter approximation with e.g. beer’s law)
    # for better realism, one would integrate or approximate those exponentials 
    # along both the sun’s path and the view path

    # we keep track of a small "transmittance" factor that accounts for
    # e^{-integral_of_sigma ds} for the line of sight behind the scatter point
    # similarly for the sun's path, we need the path from top-of-atmosphere 
    # to that scattering point

    # for speed, we'll do the (very approximate) approach: 
    #   each step: 
    #       1) find local altitude -> sigma_loc
    #       2) find approximate transmittance from ground to "r"
    #       3) find approximate transmittance along sun path
    #       4) scatter contribution
    # real codes do more exact geometry, but let's keep it short

    # precompute the direction cosines for the sun
    #   so we can quickly do geometry for where sun intersects the atmosphere, etc.
    # for single scattering, the scattering angle is the angle between v=(vx,vy,vz) 
    # and s=(sx,sy,sz)

    cos_scatter = vx*sx + vy*sy + vz*sz
    scatter_angle = math.acos(cos_scatter)

    # very naive integration
    # path optical depth in the view direction from ground up to r:
    tau_view = 0.0
    # path optical depth in sun direction from top of atmosphere 
    # (or from local entry to the atmosphere) to the scattering point:
    # we’ll do a simplified approach: compute an approximate airmass 
    # from that altitude to the top of atmosphere
    # this is quite rough, but suffices for demonstration
    
    for i in range(n_steps):
        # current radius
        r = r0 + (i+0.5)*dr
        sig = scattering_coefficient(r)
        # increment optical depth in the viewing path
        dtau_view = sig * dr
        tau_view += dtau_view

        # approximate solar path optical depth from r to top-of-atmosphere
        # this is extremely approximate, but shows the idea
        # if we assume the solar ray is from "above," 
        # we can find a local "airmass" from altitude to top-of-atmosphere:
        # a better approach: do a short line integral in the solar direction
        # for demonstration:
        #   thickness ~ (R_TOA - r)/cos(zs)
        #   so tau_sun ~ sigma_averaged * thickness
        # but let's do something simplistic:

        # the local density is atmospheric_density(r)
        # the total column above that altitude ~ H_SCALE * density(r)
        # we multiply by some geometry factor ~ 1/cos(zs)
        # if the sun is near the horizon, you'd integrate tangentially, 
        # which is more complicated
        # we'll assume the sun is not too close to horizon:
        tau_sun_approx = scattering_coefficient(r)* (R_TOA - r) * (1.0 / max(1e-3, math.cos(zs)))

        # local transmittance in view direction
        T_view = math.exp(-tau_view)
        # local transmittance in sun direction
        T_sun = math.exp(-tau_sun_approx)

        # single-scattered radiance dI
        #   dI ~ I_sun * sigma(r)*phase(scatter_angle) * T_view * T_sun * dr
        # I_sun = 1.0, let's keep it as 1.0 for demonstration
        dI = 1.0 * sig * phase_function_rayleigh(scatter_angle) * T_view * T_sun * dr

        total_intensity += dI

    return total_intensity


############################################################
#         MAKE CONTOUR PLOT FOR SPHERICAL GEOMETRY         #
############################################################


def make_contour_plot_spherical(zenith_sun, azimuth_sun):
    """
    make a polar plot (zenith, azimuth) of the sky brightness 
    using single scattering in a spherical atmosphere
    """
    # using linspace for azimuth from 0->360
    azimuths_deg = np.linspace(0, 360, 181)
    zeniths_deg = np.arange(0, 91, 2)

    # create mesh
    A, Z = np.meshgrid(azimuths_deg, zeniths_deg)

    # prepare array for intensities
    intensities = np.zeros_like(A)

    # fill intensities
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            azi_view = A[i,j]
            zen_view = Z[i,j]
            # compute single-scattered intensity
            val = intensity_spherical_single_scatter(
                zen_view, azi_view,
                zenith_sun, azimuth_sun,
                step_size=2000.0
            )
            intensities[i,j] = val

    # now plot with polar coordinates: 
    # in matplotlib’s polar plot:
    #   "theta" is the angle around the plot (azimuth),
    #   "r" is the radial coordinate
    # we’ll put the zenith angle as r, and the azimuth as theta
    # convert degrees to radians for plotting
    THETA = np.radians(A)
    R = Z

    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(8,6))
    # use a contourf
    # you might want to log-scale intensities or use another colormap
    twilight_colormap = create_twilight_colormap()
    c = ax.contourf(THETA, R, intensities, 60, cmap=twilight_colormap)

    # add colorbar
    plt.colorbar(c, ax=ax, orientation='vertical', label='Single-scatter intensity (arb units)')

    # tweak polar appearance
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1) # so azimuth increases clockwise or counterclockwise
    ax.set_rticks([0, 30, 60, 90])
    ax.set_rlabel_position(135)

    # optionally mark the sun
    # if the sun is at (zenith_sun, azimuth_sun), on the plot we have 
    #    r_sun = zenith_sun 
    #    theta_sun = azimuth_sun (in radians, but flipped sign if you want it around)
    # for simplicity, do something small:
    sun_theta = math.radians(azimuth_sun)
    sun_r = zenith_sun
    ax.plot([sun_theta], [sun_r], 'o', color='red', markersize=10, label='Sun')

    plt.legend(loc='best')
    plt.title(f'Sky Map - Spherical Single Scattering\nSun z={zenith_sun}°, a={azimuth_sun}°')
    plt.show()


############################################################
#                           MAIN                          #
############################################################


def main():
    """
    main driver for testing
    we define a sun position by zenith angle and azimuth angle, 
    then create a contour plot of sky brightness
    """
    # for example, sun at zenith_sun=35°, azimuth_sun=0° (due south or north, etc.)
    zenith_sun = 35.0
    azimuth_sun = 0.0

    make_contour_plot_spherical(zenith_sun, azimuth_sun)

if __name__ == '__main__':
    main()
