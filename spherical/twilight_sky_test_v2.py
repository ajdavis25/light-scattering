#!/usr/bin/env python3
import math
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
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

def atmospheric_density(r):
    """
    simple exponential atmospheric density as a function of radial distance r,
    where r is distance from earth's center (m)
    altitude = r - R_EARTH
    """
    altitude = r - R_EARTH
    if altitude < 0.0:
        return 0.0
    return math.exp(-altitude / H_SCALE)

def scattering_coefficient(r):
    """
    local scattering coefficient at radius r
    here we choose: sigma(r) = SIGMA_0 * exp(-(r - R_EARTH)/H_SCALE)
    """
    return SIGMA_0 * atmospheric_density(r)

def phase_function_rayleigh(scatter_angle):
    """
    rayleigh phase function (unpolarized) ~ (3/4)*(1 + cos^2(scatter_angle))
    """
    return 0.75 * (1.0 + math.cos(scatter_angle)**2)


############################################################
#                SOLAR TRANSMITTANCE FUNCTION             #
############################################################

def solar_transmittance(pos, sun_dir, R_E=R_EARTH, R_TOA=R_TOA, step_sun=500.0):
    """
    return exp(-tau_sun), where tau_sun is the optical depth 
    from the scattering point 'pos' to the top of the atmosphere 
    in direction 'sun_dir'
    """
    r_p = np.linalg.norm(pos) # radial distance of the scattering point

    if r_p >= R_TOA:
        return 1.0

    # solve for intersection with sphere r=R_TOA
    rp_dot_s = pos[0]*sun_dir[0] + pos[1]*sun_dir[1] + pos[2]*sun_dir[2]
    A = 1.0  # sun_dir is assumed unit length
    B = 2.0 * rp_dot_s
    C = r_p*r_p - R_TOA*R_TOA
    disc = B*B - 4.0*A*C
    if disc <= 0.0:
        return 1.0
    sqrt_disc = math.sqrt(disc)

    t1 = -0.5*(B - sqrt_disc)
    t2 = -0.5*(B + sqrt_disc)

    t_candidates = []
    if t1 > 0.0:
        t_candidates.append(t1)
    if t2 > 0.0:
        t_candidates.append(t2)

    if len(t_candidates) == 0:
        return 1.0

    t_max = min(t_candidates)

    n_steps = int(t_max // step_sun) + 1
    dt = t_max / max(n_steps, 1)
    
    tau_sun = 0.0
    for i in range(n_steps):
        t_mid = (i + 0.5) * dt
        xyz = pos + t_mid * sun_dir
        r_loc = np.linalg.norm(xyz)
        sig_loc = scattering_coefficient(r_loc)
        tau_sun += sig_loc * dt

    return math.exp(-tau_sun)


############################################################
#       SINGLE-SCATTER INTENSITY FOR ONE VIEW DIRECTION    #
############################################################

def intensity_spherical_single_scatter(
    zenith_view_deg, azimuth_view_deg,
    zenith_sun_deg, azimuth_sun_deg,
    step_size=2000.0
):
    """
    compute single-scattered radiance along a line of sight given by
    (zenith_view_deg, azimuth_view_deg),
    from an observer at r=R_EARTH out to r=R_TOA,
    including a more precise solar transmittance for each scattering point
    """

    zv = math.radians(zenith_view_deg)
    av = math.radians(azimuth_view_deg)
    zs = math.radians(zenith_sun_deg)
    as_ = math.radians(azimuth_sun_deg)

    vx = math.sin(zv)*math.cos(av)
    vy = math.sin(zv)*math.sin(av)
    vz = math.cos(zv)
    view_dir = np.array([vx, vy, vz], dtype=float)

    sx = math.sin(zs)*math.cos(as_)
    sy = math.sin(zs)*math.sin(as_)
    sz = math.cos(zs)
    sun_dir = np.array([sx, sy, sz], dtype=float)

    r0 = R_EARTH
    r1 = R_TOA
    path_length = r1 - r0
    n_steps = int(path_length // step_size) + 1
    dr = path_length / n_steps

    tau_view = 0.0
    total_intensity = 0.0

    cos_scatter = vx*sx + vy*sy + vz*sz
    cos_scatter = max(-1.0, min(1.0, cos_scatter))
    scatter_angle = math.acos(cos_scatter)

    for i in range(n_steps):
        r_mid = r0 + (i + 0.5)*dr
        sig_loc = scattering_coefficient(r_mid)

        dtau_view = sig_loc * dr
        tau_view += dtau_view
        T_view = math.exp(-tau_view)

        # place the observer at earth's surface on the z-axis
        observer_pos = np.array([0.0, 0.0, R_EARTH], dtype=float)

        # solve for s in |obs + s v| = r_mid
        obs_norm = np.linalg.norm(observer_pos)
        dot_ov = observer_pos.dot(view_dir)

        A = 1.0
        B = 2.0*dot_ov
        C = obs_norm*obs_norm - r_mid*r_mid
        disc = B*B - 4.0*A*C
        if disc < 0.0:
            continue

        s1 = (-B + math.sqrt(disc))/2.0
        s2 = (-B - math.sqrt(disc))/2.0
        s_candidates = []
        for sol in (s1, s2):
            if sol > 0:
                s_candidates.append(sol)
        if len(s_candidates) == 0:
            continue
        s_final = min(s_candidates)
        pos_i = observer_pos + s_final*view_dir

        T_sun = solar_transmittance(pos_i, sun_dir, R_E=R_EARTH, R_TOA=R_TOA, step_sun=500.0)

        dI = (1.0 * sig_loc * phase_function_rayleigh(scatter_angle)
              * T_view * T_sun * dr)
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
    # We can reduce resolution if we need faster testing:
    azimuths_deg = np.linspace(0, 360, 181)  # 0..360
    zeniths_deg = np.arange(0, 91, 2)        # 0..90

    # mesh grid
    A, Z = np.meshgrid(azimuths_deg, zeniths_deg)
    intensities = np.zeros_like(A)

    # loop over sky directions
    for i in tqdm(range(A.shape[0]), desc="Scanning Azimuth Rows"):
        for j in range(A.shape[1]):
            azi_view = A[i, j]
            zen_view = Z[i, j]

            # compute single-scattered intensity
            val = intensity_spherical_single_scatter(
                zen_view, azi_view,
                zenith_sun, azimuth_sun,
                step_size=2000.0
            )
            intensities[i, j] = val

    # convert to polar coords: theta=azimuth, r=zenith
    THETA = np.radians(A)
    R = Z

    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(8, 6))

    twilight_colormap = create_twilight_colormap() # custom colormap
    c = ax.contourf(THETA, R, intensities, 60, cmap=twilight_colormap)
    plt.colorbar(c, ax=ax, orientation='vertical', label='Single-scatter intensity (arb. units)')

    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rticks([0, 30, 60, 90])
    ax.set_rlabel_position(135)

    # mark the sun on the plot
    sun_theta = math.radians(azimuth_sun)
    sun_r = zenith_sun
    ax.plot([sun_theta], [sun_r], 'o', color='red', markersize=10, label='Sun')

    plt.legend(loc='best')
    plt.title(f"Sky Map (Single Scatter, Spherical)\nSun: zenith={zenith_sun}°, azimuth={azimuth_sun}°")
    plt.show()


############################################################
#                           MAIN                          #
############################################################

def main():
    """
    example usage: render a contour plot of the sky for a near-horizon sun
    """
    # for near-horizon scenario, e.g. sun ~ 5° above horizon
    zenith_sun = 85.0
    azimuth_sun = 0.0
    make_contour_plot_spherical(zenith_sun, azimuth_sun)

if __name__ == '__main__':
    main()
