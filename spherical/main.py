from __future__ import annotations

import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from spherical.measurement_plots import save_measurement_case_plots
from spherical.production_io import ProductionSkyResult, load_production_result, write_result_bundle_npz
from spherical.twilight_baseline import (
    AtmosphereConfig,
    ObserverConfig,
    SolarGeometry,
    estimate_hemispheric_flux,
    single_scatter_stokes,
)
from spherical.utils import create_twilight_colormap

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "monte_carlo_cpp" / "config" / "default_clear_sky.cfg"
DEFAULT_MEASUREMENT_CONFIG = (
    ROOT
    / "monte_carlo_cpp"
    / "config"
    / "paper_cases"
    / "frozen_marseille_twilight_20220815_191413z_measurement_interactive.cfg"
)
LEGACY_FULLSKY_MEASUREMENT_CONFIG = (
    ROOT / "monte_carlo_cpp" / "config" / "measurement_gal_lapland_fullsky_450nm_dolp.cfg"
)
DEFAULT_CASE_DIR = ROOT / "monte_carlo_cpp" / "results" / "default_clear_sky"
DEFAULT_PLOT_DIR = ROOT / "plots" / "current"


@dataclass(frozen=True)
class AnalyticReferenceResult:
    intensity: np.ndarray
    q: np.ndarray
    u: np.ndarray

    @property
    def degree_of_polarization(self) -> np.ndarray:
        polarized = np.sqrt(self.q**2 + self.u**2)
        safe_i = np.where(self.intensity > 0.0, self.intensity, np.nan)
        return np.nan_to_num(polarized / safe_i, nan=0.0, posinf=0.0, neginf=0.0)


def load_key_value_config(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def find_solver_executable() -> Path:
    candidates = [
        ROOT / "monte_carlo_cpp" / "build_current" / "MonteCarloCPP.exe",
        ROOT / "monte_carlo_cpp" / "build" / "MonteCarloCPP.exe",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Could not find a built MonteCarloCPP.exe in build_current or build.")


def ensure_production_result(config_path: Path = DEFAULT_CONFIG, case_dir: Path = DEFAULT_CASE_DIR) -> None:
    csv_path = case_dir / "sky_result.csv"
    metadata_path = case_dir / "sky_result_metadata.json"
    if csv_path.exists() and metadata_path.exists():
        return

    executable = find_solver_executable()
    subprocess.run([str(executable), str(config_path)], cwd=ROOT, check=True)


def plot_fisheye(
    zenith_grid_deg: np.ndarray,
    azimuth_grid_deg: np.ndarray,
    values: np.ndarray,
    title: str,
    colorbar_label: str,
    save_path: Path,
    sun_zenith_deg: float,
    sun_azimuth_deg: float,
    cmap: str | None = None,
) -> None:
    save_path.parent.mkdir(parents=True, exist_ok=True)
    theta, radius = np.meshgrid(np.radians(azimuth_grid_deg), zenith_grid_deg)
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(7, 7))
    mesh = ax.pcolormesh(
        theta,
        radius,
        values,
        shading="auto",
        cmap=cmap or create_twilight_colormap(),
    )
    fig.colorbar(mesh, ax=ax, label=colorbar_label)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0.0, 90.0)
    ax.scatter(
        math.radians(sun_azimuth_deg),
        sun_zenith_deg,
        color="yellow",
        s=100,
        edgecolors="black",
        linewidths=0.5,
        label="Sun",
    )
    ax.set_title(title)
    ax.legend(loc="upper right")
    fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)


def build_analytic_reference(
    production: ProductionSkyResult,
    observer: ObserverConfig,
    sun: SolarGeometry,
) -> AnalyticReferenceResult:
    atmosphere = AtmosphereConfig(integration_step_m=5.0e3)
    intensity = np.zeros_like(production.intensity)
    q = np.zeros_like(production.q)
    u = np.zeros_like(production.u)

    for i_zen, zenith_deg in enumerate(production.zenith_grid_deg):
        for i_azi, azimuth_deg in enumerate(production.azimuth_grid_deg):
            i_val, q_val, u_val = single_scatter_stokes(
                view_zenith_deg=float(zenith_deg),
                view_azimuth_deg=float(azimuth_deg),
                observer=observer,
                sun=sun,
                atmosphere=atmosphere,
            )
            intensity[i_zen, i_azi] = i_val
            q[i_zen, i_azi] = q_val
            u[i_zen, i_azi] = u_val

    return AnalyticReferenceResult(intensity=intensity, q=q, u=u)


def save_production_outputs(
    case_dir: Path = DEFAULT_CASE_DIR,
    plot_dir: Path = DEFAULT_PLOT_DIR,
    config_path: Path = DEFAULT_CONFIG,
) -> ProductionSkyResult:
    ensure_production_result(config_path=config_path, case_dir=case_dir)
    production = load_production_result(case_dir)
    config_values = load_key_value_config(config_path)

    observer = ObserverConfig(
        latitude_deg=float(config_values.get("observer_latitude_deg", 45.0)),
        longitude_deg=float(config_values.get("observer_longitude_deg", 0.0)),
        altitude_m=float(config_values.get("observer_altitude_m", 0.0)),
    )
    sun = SolarGeometry(
        zenith_deg=float(production.metadata["sun_zenith_deg"]),
        azimuth_deg=float(production.metadata["sun_azimuth_deg"]),
    )
    analytic = build_analytic_reference(production=production, observer=observer, sun=sun)

    plot_fisheye(
        production.zenith_grid_deg,
        production.azimuth_grid_deg,
        production.intensity,
        "Production Twilight Intensity",
        "Radiance (arb. units)",
        plot_dir / "production_intensity.png",
        float(production.metadata["sun_zenith_deg"]),
        float(production.metadata["sun_azimuth_deg"]),
    )
    plot_fisheye(
        production.zenith_grid_deg,
        production.azimuth_grid_deg,
        production.degree_of_polarization,
        "Production Twilight Degree Of Polarization",
        "DoLP",
        plot_dir / "production_dop.png",
        float(production.metadata["sun_zenith_deg"]),
        float(production.metadata["sun_azimuth_deg"]),
    )
    plot_fisheye(
        production.zenith_grid_deg,
        production.azimuth_grid_deg,
        production.angle_of_polarization_rad,
        "Production Twilight Angle Of Polarization",
        "AoP (rad)",
        plot_dir / "production_aop.png",
        float(production.metadata["sun_zenith_deg"]),
        float(production.metadata["sun_azimuth_deg"]),
        cmap="coolwarm",
    )
    plot_fisheye(
        production.zenith_grid_deg,
        production.azimuth_grid_deg,
        analytic.intensity,
        "Analytic Single-Scatter Reference Intensity",
        "Radiance (arb. units)",
        plot_dir / "analytic_reference_intensity.png",
        float(production.metadata["sun_zenith_deg"]),
        float(production.metadata["sun_azimuth_deg"]),
    )

    production_norm = production.intensity / max(1.0e-12, float(np.max(production.intensity)))
    analytic_norm = analytic.intensity / max(1.0e-12, float(np.max(analytic.intensity)))
    plot_fisheye(
        production.zenith_grid_deg,
        production.azimuth_grid_deg,
        production_norm - analytic_norm,
        "Production Minus Analytic Intensity",
        "Normalized Difference",
        plot_dir / "comparison_intensity_difference.png",
        float(production.metadata["sun_zenith_deg"]),
        float(production.metadata["sun_azimuth_deg"]),
        cmap="coolwarm",
    )
    plot_fisheye(
        production.zenith_grid_deg,
        production.azimuth_grid_deg,
        production.degree_of_polarization - analytic.degree_of_polarization,
        "Production Minus Analytic DoLP",
        "DoLP Difference",
        plot_dir / "comparison_dop_difference.png",
        float(production.metadata["sun_zenith_deg"]),
        float(production.metadata["sun_azimuth_deg"]),
        cmap="coolwarm",
    )
    plot_fisheye(
        production.zenith_grid_deg,
        production.azimuth_grid_deg,
        np.sqrt(np.maximum(production.var_i, 0.0)),
        "Production Intensity Standard Deviation",
        "Std(I)",
        plot_dir / "production_intensity_std.png",
        float(production.metadata["sun_zenith_deg"]),
        float(production.metadata["sun_azimuth_deg"]),
        cmap="viridis",
    )

    bundle_path = write_result_bundle_npz(production, case_dir / "sky_result_bundle.npz")
    comparison_metadata = {
        "analytic_reference": "spherical.twilight_baseline.single_scatter_stokes",
        "case_dir": str(case_dir),
        "bundle_path": str(bundle_path),
        "plot_dir": str(plot_dir),
        "production_peak_intensity": float(np.max(production.intensity)),
        "production_peak_dop": float(np.max(production.degree_of_polarization)),
        "analytic_peak_intensity": float(np.max(analytic.intensity)),
        "analytic_peak_dop": float(np.max(analytic.degree_of_polarization)),
        "analytic_hemispheric_flux_estimate": float(
            estimate_hemispheric_flux(
                type(
                    "ReferenceResult",
                    (),
                    {
                        "zenith_grid_deg": production.zenith_grid_deg,
                        "azimuth_grid_deg": production.azimuth_grid_deg,
                        "intensity": analytic.intensity,
                    },
                )()
            )
        ),
    }
    (case_dir / "comparison_metadata.json").write_text(json.dumps(comparison_metadata, indent=2))
    return production


def main() -> None:
    production = save_production_outputs()
    measurement_summary = save_measurement_case_plots(
        config_path=DEFAULT_MEASUREMENT_CONFIG,
        plot_root=DEFAULT_PLOT_DIR / "measurement_cases",
    )
    print(f"Production case: {production.case_dir}")
    print(f"Peak intensity: {float(np.max(production.intensity)):.6e}")
    print(f"Peak DoLP: {float(np.max(production.degree_of_polarization)):.6f}")
    print(f"Bundle: {production.case_dir / 'sky_result_bundle.npz'}")
    print(f"Measurement case: {measurement_summary['case_id']}")
    print(f"Measurement plots: {measurement_summary['plot_dir']}")


if __name__ == "__main__":
    main()
