#!/usr/bin/env python3
import math
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from utils import create_twilight_colormap


############################################################
#              PHYSICS CONSTANTS / SIMPLE PARAMS          #
############################################################

R_EARTH = 6.371e6
ATM_HEIGHT = 1.0e5
R_TOA = R_EARTH + ATM_HEIGHT
H_SCALE = 8.0e3
SIGMA_0 = 1e-5

# "fudge" factor to emulate some multiple scattering
# e.g., 1.0 => pure single scattering, 1.3 => mild brightening
MULTIPLE_SCAT_FACTOR = 1.3

# "aerosol fraction" for a simple forward-scattering addition
# e.g., 0.0 => pure rayleigh, 0.3 => 30% extra forward scatter
AEROSOL_FRACTION = 0.3

def atmospheric_density(r):
    """simple exponential density, ignoring negativity below earth surface"""
    altitude = r - R_EARTH
    if altitude < 0.0:
        return 0.0
    return math.exp(-altitude / H_SCALE)

def scattering_coefficient(r):
    """rayleigh base + exponential, we could add aerosol term too"""
    return SIGMA_0 * atmospheric_density(r)


############################################################
#           SCATTERING PHASE FUNCTIONS (SIMPLE)           #
############################################################

def phase_function_rayleigh(scatter_angle):
    """
    rayleigh phase function (unpolarized) ~ (3/4)*(1 + cos^2(scatter_angle))
    we'll blend in a bit of "forward scattering" for aerosols if desired
    """
    # pure rayleigh
    pr = 0.75 * (1.0 + math.cos(scatter_angle)**2)

    # a simple aerosol forward-lobe:
    # e.g. a henyey-greenstein with g=0.7, or just a "peak" near cos=1
    # we'll do a trivial shape: 1 + 5*cos(scatter_angle)
    # real aerosol phase is more complicated, this is just an example
    pa = max(0.0, 1.0 + 5.0*math.cos(scatter_angle))

    # weighted combination: (1 - AEROSOL_FRACTION)*Rayleigh + AEROSOL_FRACTION*(aerosol)
    return (1.0 - AEROSOL_FRACTION)*pr + (AEROSOL_FRACTION)*pa

def rayleigh_polarization_fraction(scatter_angle):
    """
    single-scatter rayleigh fractional polarization:
      Pol = sin^2(theta) / (1 + cos^2(theta))
    ignoring aerosol polarization here for simplicity
    """
    cos_th = math.cos(scatter_angle)
    denominator = 1.0 + cos_th*cos_th
    if abs(denominator) < 1e-12:
        return 0.0
    numerator = 1.0 - cos_th*cos_th
    return numerator / denominator


############################################################
#                SOLAR TRANSMITTANCE FUNCTION             #
############################################################

def solar_transmittance(pos, sun_dir, R_E=R_EARTH, R_TOA=R_TOA, step_sun=500.0):
    """
    returns exp(-tau_sun) from 'pos' to top-of-atmosphere in direction sun_dir
    now includes a check if path dips below earth (r < R_E)
    """
    r_p = np.linalg.norm(pos)
    if r_p >= R_TOA:
        return 1.0

    # solve for intersection with sphere r=R_TOA
    rp_dot_s = np.dot(pos, sun_dir)
    A = 1.0
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
    if not t_candidates:
        return 1.0

    t_max = min(t_candidates)
    n_steps = int(t_max // step_sun) + 1
    dt = t_max / max(n_steps, 1)

    tau_sun = 0.0
    for i in range(n_steps):
        t_mid = (i + 0.5)*dt
        xyz = pos + t_mid*sun_dir
        r_loc = np.linalg.norm(xyz)
        # if we've gone below earth, no solar illumination => break
        if r_loc < R_E:
            tau_sun = 999.0 # effectively zero transmittance
            break

        sig_loc = scattering_coefficient(r_loc)
        tau_sun += sig_loc * dt

    return math.exp(-tau_sun)


############################################################
#       SINGLE-SCATTER I & POL FOR ONE VIEW DIRECTION      #
############################################################

def intensity_and_pol_spherical_single_scatter(
    zenith_view_deg, azimuth_view_deg,
    zenith_sun_deg, azimuth_sun_deg,
    step_size=2000.0
):
    """
    return (I, Pol) for single scattering along a line of sight:
      I : scalar intensity (with a 'multiple-scattering' fudge factor)
      Pol : fractional polarization (rayleigh formula)
    """

    # convert angles to radians
    zv = math.radians(zenith_view_deg)
    av = math.radians(azimuth_view_deg)
    zs = math.radians(zenith_sun_deg)
    as_ = math.radians(azimuth_sun_deg)

    # view direction (unit vector)
    vx = math.sin(zv)*math.cos(av)
    vy = math.sin(zv)*math.sin(av)
    vz = math.cos(zv)
    view_dir = np.array([vx, vy, vz], dtype=float)

    # sun direction (unit vector)
    sx = math.sin(zs)*math.cos(as_)
    sy = math.sin(zs)*math.sin(as_)
    sz = math.cos(zs)
    sun_dir = np.array([sx, sy, sz], dtype=float)

    # setup line-of-sight integration from ground to top-of-atmosphere
    r0 = R_EARTH
    r1 = R_TOA
    path_length = r1 - r0
    n_steps = int(path_length // step_size) + 1
    dr = path_length / n_steps

    tau_view = 0.0
    total_I = 0.0

    # single-scatter angle for entire path (approx)
    cos_scatter = max(-1.0, min(1.0, np.dot(view_dir, sun_dir)))
    scatter_angle = math.acos(cos_scatter)
    pol_fraction = rayleigh_polarization_fraction(scatter_angle)

    for i in range(n_steps):
        r_mid = r0 + (i + 0.5)*dr
        sig_loc = scattering_coefficient(r_mid)

        dtau_view = sig_loc * dr
        tau_view += dtau_view
        T_view = math.exp(-tau_view)

        # 3d position along the view ray:
        # solve for s so that |obs + s*view_dir| = r_mid
        observer_pos = np.array([0.0, 0.0, R_EARTH], dtype=float)
        obs_norm = np.linalg.norm(observer_pos)
        dot_ov = np.dot(observer_pos, view_dir)

        # quadratic: s^2 + 2(dot_ov)*s + (obs_norm^2 - r_mid^2) = 0
        A = 1.0
        B = 2.0 * dot_ov
        C = obs_norm*obs_norm - r_mid*r_mid
        disc = B*B - 4.0*A*C
        if disc < 0.0:
            continue

        s1 = (-B + math.sqrt(disc))/2.0
        s2 = (-B - math.sqrt(disc))/2.0
        s_candidates = [sol for sol in (s1, s2) if sol > 0]
        if not s_candidates:
            continue
        s_final = min(s_candidates)
        pos_i = observer_pos + s_final*view_dir

        # if we've dipped below earth, skip
        if np.linalg.norm(pos_i) < R_EARTH:
            continue

        T_sun = solar_transmittance(pos_i, sun_dir, R_EARTH, R_TOA, step_sun=500.0)

        # single-scatter from this segment
        dI = (1.0
              * sig_loc
              * phase_function_rayleigh(scatter_angle)
              * T_view
              * T_sun
              * dr)

        total_I += dI

    # apply a "multiple-scattering" fudge factor
    total_I *= MULTIPLE_SCAT_FACTOR

    return total_I, pol_fraction


############################################################
#         MAKE CONTOUR PLOT FOR SPHERICAL GEOMETRY         #
############################################################

def make_contour_plot_spherical(zenith_sun, azimuth_sun):
    """
    we'll do 0..90 for zenith so we see from overhead (0 deg) 
    down to the horizon (90 deg)
    """
    azimuths_deg = np.linspace(0, 360, 181) # 0..360
    # restrict zenith to 0..90 (the visible dome from ground)
    zeniths_deg = np.arange(0, 91, 5)

    A, Z = np.meshgrid(azimuths_deg, zeniths_deg)
    intensity_array = np.zeros_like(A)
    pol_array = np.zeros_like(A)

    # loop over directions
    for i in tqdm(range(A.shape[0]), desc="Scanning Azimuth Rows"):
        for j in range(A.shape[1]):
            azi_view = A[i, j]
            zen_view = Z[i, j]
            I_val, pol_val = intensity_and_pol_spherical_single_scatter(
                zen_view, azi_view,
                zenith_sun, azimuth_sun,
                step_size=2000.0
            )
            intensity_array[i, j] = I_val
            pol_array[i, j] = pol_val

    # convert to polar coords: theta=azimuth, r=zenith
    THETA = np.radians(A)
    R = Z

    # --- plot intensity ---
    fig1, ax1 = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(8, 6))
    twilight_colormap = create_twilight_colormap()
    c1 = ax1.contourf(THETA, R, intensity_array, 60, cmap=twilight_colormap)
    plt.colorbar(c1, ax=ax1, orientation='vertical', label='I (Single-scatter + approx)')

    ax1.set_theta_zero_location("N")
    ax1.set_theta_direction(-1)
    ax1.set_rticks([0, 30, 60, 90])
    ax1.set_rlabel_position(135)

    sun_theta = math.radians(azimuth_sun)
    sun_r = zenith_sun
    ax1.plot([sun_theta], [sun_r], 'o', color='red', markersize=10, label='Sun')
    ax1.legend(loc='best')
    ax1.set_title(f"I map: 0–90 zenith\nSun: z={zenith_sun}°, a={azimuth_sun}°")

    # --- plot polarization fraction ---
    fig2, ax2 = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(8, 6))
    c2 = ax2.contourf(THETA, R, pol_array, 60, cmap=twilight_colormap)
    plt.colorbar(c2, ax=ax2, orientation='vertical', label='Pol (fraction 0..1)')

    ax2.set_theta_zero_location("N")
    ax2.set_theta_direction(-1)
    ax2.set_rticks([0, 30, 60, 90])
    ax2.set_rlabel_position(135)

    ax2.plot([sun_theta], [sun_r], 'o', color='red', markersize=10, label='Sun')
    ax2.legend(loc='best')
    ax2.set_title(f"Polarization fraction: 0–90 zenith\nSun: z={zenith_sun}°, a={azimuth_sun}°")

    plt.show()


############################################################
#                           MAIN                          #
############################################################

def main():
    """
    example usage: see sky from 0..90 deg zenith for near-horizon sun
    we also apply a small multiple-scattering factor and an aerosol fraction
    """
    zenith_sun = 85.0
    azimuth_sun = 0.0
    make_contour_plot_spherical(zenith_sun, azimuth_sun)

if __name__ == '__main__':
    main()
