from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Tuple

import numpy as np


def _clamp(value: float, lower: float = -1.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, value))


@dataclass(frozen=True)
class AtmosphereConfig:
    earth_radius_m: float = 6.371e6
    atmosphere_height_m: float = 1.0e5
    rayleigh_scale_height_m: float = 8.0e3
    rayleigh_scattering_ground_m_inv: float = 1.0e-5
    absorption_ground_m_inv: float = 2.5e-6
    absorption_scale_height_m: float = 7.0e3
    integration_step_m: float = 2.0e3

    @property
    def top_of_atmosphere_radius_m(self) -> float:
        return self.earth_radius_m + self.atmosphere_height_m


@dataclass(frozen=True)
class ObserverConfig:
    latitude_deg: float = 45.0
    longitude_deg: float = 0.0
    altitude_m: float = 0.0


@dataclass(frozen=True)
class SolarGeometry:
    zenith_deg: float
    azimuth_deg: float

    @classmethod
    def from_hour_angle(
        cls,
        latitude_deg: float,
        solar_declination_deg: float,
        hour_angle_deg: float,
    ) -> "SolarGeometry":
        zenith_deg, azimuth_deg = solar_angles_from_hour_angle(
            latitude_deg=latitude_deg,
            solar_declination_deg=solar_declination_deg,
            hour_angle_deg=hour_angle_deg,
        )
        return cls(zenith_deg=zenith_deg, azimuth_deg=azimuth_deg)


@dataclass(frozen=True)
class SkySimulationResult:
    zenith_grid_deg: np.ndarray
    azimuth_grid_deg: np.ndarray
    intensity: np.ndarray
    q: np.ndarray
    u: np.ndarray
    sun_zenith_deg: float
    sun_azimuth_deg: float

    @property
    def degree_of_polarization(self) -> np.ndarray:
        polarized = np.sqrt(self.q**2 + self.u**2)
        safe_intensity = np.where(self.intensity > 0.0, self.intensity, np.nan)
        dop = polarized / safe_intensity
        return np.nan_to_num(dop, nan=0.0, posinf=0.0, neginf=0.0)


def solar_angles_from_hour_angle(
    latitude_deg: float,
    solar_declination_deg: float,
    hour_angle_deg: float,
) -> Tuple[float, float]:
    lat = math.radians(latitude_deg)
    dec = math.radians(solar_declination_deg)
    hour = math.radians(hour_angle_deg)

    cos_zenith = (
        math.sin(lat) * math.sin(dec)
        + math.cos(lat) * math.cos(dec) * math.cos(hour)
    )
    zenith = math.degrees(math.acos(_clamp(cos_zenith)))

    sin_zenith = math.sqrt(max(0.0, 1.0 - cos_zenith * cos_zenith))
    if sin_zenith < 1e-12:
        return zenith, 0.0

    sin_azimuth = -math.cos(dec) * math.sin(hour) / sin_zenith
    cos_azimuth = (
        math.sin(dec) - math.sin(lat) * cos_zenith
    ) / max(1e-12, math.cos(lat) * sin_zenith)
    azimuth = math.degrees(math.atan2(sin_azimuth, cos_azimuth)) % 360.0
    return zenith, azimuth


def observer_position(observer: ObserverConfig, atmosphere: AtmosphereConfig) -> np.ndarray:
    lat = math.radians(observer.latitude_deg)
    lon = math.radians(observer.longitude_deg)
    radius = atmosphere.earth_radius_m + observer.altitude_m
    return np.array(
        [
            radius * math.cos(lat) * math.cos(lon),
            radius * math.cos(lat) * math.sin(lon),
            radius * math.sin(lat),
        ],
        dtype=float,
    )


def local_basis(observer: ObserverConfig, atmosphere: AtmosphereConfig) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    up = observer_position(observer, atmosphere)
    up = up / np.linalg.norm(up)

    lon = math.radians(observer.longitude_deg)
    east = np.array([-math.sin(lon), math.cos(lon), 0.0], dtype=float)
    east = east / np.linalg.norm(east)

    north = np.cross(up, east)
    north = north / np.linalg.norm(north)
    return north, east, up


def direction_from_zenith_azimuth(
    zenith_deg: float,
    azimuth_deg: float,
    north: np.ndarray,
    east: np.ndarray,
    up: np.ndarray,
) -> np.ndarray:
    zenith_rad = math.radians(zenith_deg)
    azimuth_rad = math.radians(azimuth_deg)
    direction = (
        math.sin(zenith_rad) * math.cos(azimuth_rad) * north
        + math.sin(zenith_rad) * math.sin(azimuth_rad) * east
        + math.cos(zenith_rad) * up
    )
    return direction / np.linalg.norm(direction)


def sun_direction(
    observer: ObserverConfig,
    atmosphere: AtmosphereConfig,
    sun: SolarGeometry,
) -> np.ndarray:
    north, east, up = local_basis(observer, atmosphere)
    return direction_from_zenith_azimuth(
        zenith_deg=sun.zenith_deg,
        azimuth_deg=sun.azimuth_deg,
        north=north,
        east=east,
        up=up,
    )


def distance_to_sphere_exit(origin: np.ndarray, direction: np.ndarray, sphere_radius_m: float) -> float:
    b = float(np.dot(origin, direction))
    c = float(np.dot(origin, origin) - sphere_radius_m * sphere_radius_m)
    discriminant = b * b - c
    if discriminant < 0.0:
        return 0.0
    root = math.sqrt(discriminant)
    near = -b - root
    far = -b + root
    candidates = [distance for distance in (near, far) if distance > 1e-9]
    return min(candidates) if candidates else 0.0


def altitude_m(position: np.ndarray, atmosphere: AtmosphereConfig) -> float:
    return float(np.linalg.norm(position) - atmosphere.earth_radius_m)


def rayleigh_scattering_coeff(altitude_value_m: float, atmosphere: AtmosphereConfig) -> float:
    if altitude_value_m < 0.0 or altitude_value_m > atmosphere.atmosphere_height_m:
        return 0.0
    return atmosphere.rayleigh_scattering_ground_m_inv * math.exp(
        -altitude_value_m / atmosphere.rayleigh_scale_height_m
    )


def absorption_coeff(altitude_value_m: float, atmosphere: AtmosphereConfig) -> float:
    if altitude_value_m < 0.0 or altitude_value_m > atmosphere.atmosphere_height_m:
        return 0.0
    return atmosphere.absorption_ground_m_inv * math.exp(
        -altitude_value_m / atmosphere.absorption_scale_height_m
    )


def extinction_coeff(altitude_value_m: float, atmosphere: AtmosphereConfig) -> float:
    return rayleigh_scattering_coeff(altitude_value_m, atmosphere) + absorption_coeff(
        altitude_value_m, atmosphere
    )


def rayleigh_phase_function(cos_scatter_angle: float) -> float:
    return (3.0 / (16.0 * math.pi)) * (1.0 + cos_scatter_angle * cos_scatter_angle)


def rayleigh_degree_of_polarization(cos_scatter_angle: float) -> float:
    numerator = max(0.0, 1.0 - cos_scatter_angle * cos_scatter_angle)
    denominator = 1.0 + cos_scatter_angle * cos_scatter_angle
    return numerator / max(1e-12, denominator)


def optical_depth_to_space(
    start_pos: np.ndarray,
    direction: np.ndarray,
    atmosphere: AtmosphereConfig,
) -> float:
    max_distance = distance_to_sphere_exit(
        origin=start_pos,
        direction=direction,
        sphere_radius_m=atmosphere.top_of_atmosphere_radius_m,
    )
    if max_distance <= 0.0:
        return 0.0

    steps = max(1, int(math.ceil(max_distance / atmosphere.integration_step_m)))
    ds = max_distance / steps
    optical_depth = 0.0
    for step_index in range(steps):
        distance = (step_index + 0.5) * ds
        position = start_pos + distance * direction
        radius = float(np.linalg.norm(position))
        if radius <= atmosphere.earth_radius_m:
            return math.inf
        optical_depth += extinction_coeff(altitude_m(position, atmosphere), atmosphere) * ds
    return optical_depth


def polarization_reference_angle(
    view_direction: np.ndarray,
    sun_direction_vector: np.ndarray,
    north: np.ndarray,
    up: np.ndarray,
) -> float:
    scattering_normal = np.cross(sun_direction_vector, view_direction)
    scatter_norm = float(np.linalg.norm(scattering_normal))
    if scatter_norm < 1e-12:
        return 0.0

    scattering_normal /= scatter_norm
    electric_vector = np.cross(scattering_normal, view_direction)
    electric_norm = float(np.linalg.norm(electric_vector))
    if electric_norm < 1e-12:
        return 0.0

    electric_vector /= electric_norm
    reference_vertical = up - float(np.dot(up, view_direction)) * view_direction
    ref_norm = float(np.linalg.norm(reference_vertical))
    if ref_norm < 1e-12:
        reference_vertical = north - float(np.dot(north, view_direction)) * view_direction
        ref_norm = float(np.linalg.norm(reference_vertical))
    if ref_norm < 1e-12:
        return 0.0

    reference_vertical /= ref_norm
    reference_horizontal = np.cross(view_direction, reference_vertical)
    reference_horizontal /= max(1e-12, float(np.linalg.norm(reference_horizontal)))
    return math.atan2(
        float(np.dot(electric_vector, reference_horizontal)),
        float(np.dot(electric_vector, reference_vertical)),
    )


def single_scatter_stokes(
    view_zenith_deg: float,
    view_azimuth_deg: float,
    observer: ObserverConfig,
    sun: SolarGeometry,
    atmosphere: AtmosphereConfig,
) -> Tuple[float, float, float]:
    obs_pos = observer_position(observer, atmosphere)
    north, east, up = local_basis(observer, atmosphere)
    view_dir = direction_from_zenith_azimuth(
        zenith_deg=view_zenith_deg,
        azimuth_deg=view_azimuth_deg,
        north=north,
        east=east,
        up=up,
    )
    sun_dir = sun_direction(observer, atmosphere, sun)

    max_distance = distance_to_sphere_exit(
        origin=obs_pos,
        direction=view_dir,
        sphere_radius_m=atmosphere.top_of_atmosphere_radius_m,
    )
    if max_distance <= 0.0:
        return 0.0, 0.0, 0.0

    cos_scatter_angle = _clamp(float(np.dot(view_dir, sun_dir)))
    phase_value = rayleigh_phase_function(cos_scatter_angle)
    polarization_fraction = rayleigh_degree_of_polarization(cos_scatter_angle)
    rotation_angle = polarization_reference_angle(view_dir, sun_dir, north, up)
    cos_2chi = math.cos(2.0 * rotation_angle)
    sin_2chi = math.sin(2.0 * rotation_angle)

    steps = max(1, int(math.ceil(max_distance / atmosphere.integration_step_m)))
    ds = max_distance / steps
    tau_view = 0.0
    total_i = 0.0
    total_q = 0.0
    total_u = 0.0

    for step_index in range(steps):
        distance = (step_index + 0.5) * ds
        scatter_pos = obs_pos + distance * view_dir
        altitude_value = altitude_m(scatter_pos, atmosphere)
        sigma_s = rayleigh_scattering_coeff(altitude_value, atmosphere)
        sigma_t = extinction_coeff(altitude_value, atmosphere)
        if sigma_s <= 0.0 or sigma_t <= 0.0:
            continue

        tau_mid = tau_view + 0.5 * sigma_t * ds
        tau_sun = optical_depth_to_space(scatter_pos, sun_dir, atmosphere)
        if not math.isfinite(tau_sun):
            tau_view += sigma_t * ds
            continue

        contribution_i = sigma_s * phase_value * math.exp(-(tau_mid + tau_sun)) * ds
        contribution_pol = contribution_i * polarization_fraction
        total_i += contribution_i
        total_q += contribution_pol * cos_2chi
        total_u += contribution_pol * sin_2chi
        tau_view += sigma_t * ds

    return total_i, total_q, total_u


def generate_sky_grid(num_zenith_points: int, num_azimuth_points: int) -> Tuple[np.ndarray, np.ndarray]:
    zenith_grid_deg = np.linspace(0.0, 90.0, num_zenith_points)
    azimuth_grid_deg = np.linspace(0.0, 360.0, num_azimuth_points, endpoint=False)
    return zenith_grid_deg, azimuth_grid_deg


def simulate_twilight_sky(
    observer: ObserverConfig,
    sun: SolarGeometry,
    atmosphere: AtmosphereConfig | None = None,
    num_zenith_points: int = 31,
    num_azimuth_points: int = 72,
) -> SkySimulationResult:
    atmosphere = atmosphere or AtmosphereConfig()
    zenith_grid_deg, azimuth_grid_deg = generate_sky_grid(
        num_zenith_points=num_zenith_points,
        num_azimuth_points=num_azimuth_points,
    )
    intensity = np.zeros((num_zenith_points, num_azimuth_points), dtype=float)
    q = np.zeros_like(intensity)
    u = np.zeros_like(intensity)

    for zenith_index, zenith_deg in enumerate(zenith_grid_deg):
        for azimuth_index, azimuth_deg in enumerate(azimuth_grid_deg):
            i_value, q_value, u_value = single_scatter_stokes(
                view_zenith_deg=float(zenith_deg),
                view_azimuth_deg=float(azimuth_deg),
                observer=observer,
                sun=sun,
                atmosphere=atmosphere,
            )
            intensity[zenith_index, azimuth_index] = i_value
            q[zenith_index, azimuth_index] = q_value
            u[zenith_index, azimuth_index] = u_value

    return SkySimulationResult(
        zenith_grid_deg=zenith_grid_deg,
        azimuth_grid_deg=azimuth_grid_deg,
        intensity=intensity,
        q=q,
        u=u,
        sun_zenith_deg=sun.zenith_deg,
        sun_azimuth_deg=sun.azimuth_deg,
    )


def estimate_hemispheric_flux(result: SkySimulationResult) -> float:
    zenith_rad = np.radians(result.zenith_grid_deg)
    azimuth_rad = np.radians(result.azimuth_grid_deg)
    if result.intensity.size == 0 or zenith_rad.size < 2 or azimuth_rad.size < 2:
        return 0.0

    d_zenith = np.diff(zenith_rad)
    d_azimuth = np.diff(np.append(azimuth_rad, 2.0 * math.pi))
    sin_zenith = np.sin(zenith_rad[:-1])
    cos_zenith = np.cos(zenith_rad[:-1])
    return float(
        np.sum(
            result.intensity[:-1, :]
            * sin_zenith[:, None]
            * cos_zenith[:, None]
            * d_zenith[:, None]
            * d_azimuth[None, :]
        )
    )
