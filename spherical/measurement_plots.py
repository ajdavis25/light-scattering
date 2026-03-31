from __future__ import annotations

import csv
import json
import math
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def load_key_value_config(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def find_measurement_runner() -> Path:
    candidates = [
        ROOT / "monte_carlo_cpp" / "build_current" / "MeasurementCaseRunner.exe",
        ROOT / "monte_carlo_cpp" / "build" / "MeasurementCaseRunner.exe",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Could not find a built MeasurementCaseRunner.exe in build_current or build.")


def measurement_case_paths(config_path: Path) -> tuple[str, Path, Path]:
    values = load_key_value_config(config_path)
    case_id = values["case_id"]
    output_dir = (config_path.parent / values.get("output_dir", "../results")).resolve()
    report_dir = output_dir / "measurement_case_reports"
    return case_id, report_dir / f"{case_id}.txt", report_dir / f"{case_id}_comparison.csv"


def ensure_measurement_case_report(config_path: Path) -> tuple[Path, Path]:
    case_id, report_path, comparison_path = measurement_case_paths(config_path)
    if report_path.exists() and comparison_path.exists():
        return report_path, comparison_path

    runner = find_measurement_runner()
    subprocess.run([str(runner), str(config_path)], cwd=ROOT, check=True)
    if not report_path.exists() or not comparison_path.exists():
        raise FileNotFoundError(
            f"Measurement runner completed but did not write expected artifacts for {case_id}: "
            f"{report_path} and {comparison_path}"
        )
    return report_path, comparison_path


def load_measurement_comparison(csv_path: Path) -> dict[str, np.ndarray]:
    rows = list(csv.DictReader(csv_path.open(newline="")))
    if not rows:
        raise ValueError(f"Measurement comparison CSV is empty: {csv_path}")

    columns = {name: [] for name in rows[0].keys()}
    for row in rows:
        for key, value in row.items():
            columns[key].append(float(value))
    return {key: np.asarray(values, dtype=float) for key, values in columns.items()}


def load_measurement_report(report_path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in report_path.read_text().splitlines():
        line = raw_line.strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def has_column(comparison: dict[str, np.ndarray], key: str) -> bool:
    return key in comparison and comparison[key].size > 0


def plot_measurement_scatter(
    *,
    zenith_deg: np.ndarray,
    azimuth_deg: np.ndarray,
    values: np.ndarray,
    title: str,
    colorbar_label: str,
    save_path: Path,
    sun_zenith_deg: float,
    sun_azimuth_deg: float,
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
) -> None:
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(7, 7))
    scatter = ax.scatter(
        np.radians(azimuth_deg),
        zenith_deg,
        c=values,
        s=180.0,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        edgecolors="black",
        linewidths=0.25,
    )
    fig.colorbar(scatter, ax=ax, label=colorbar_label)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0.0, 90.0)
    ax.scatter(
        math.radians(sun_azimuth_deg),
        sun_zenith_deg,
        color="yellow",
        s=120,
        edgecolors="black",
        linewidths=0.5,
        label="Sun",
    )
    ax.set_title(title)
    ax.legend(loc="upper right")
    fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)


def save_measurement_case_plots(config_path: Path, plot_root: Path) -> dict[str, str]:
    report_path, comparison_path = ensure_measurement_case_report(config_path)
    config_values = load_key_value_config(config_path)
    report_values = load_measurement_report(report_path)
    comparison = load_measurement_comparison(comparison_path)

    case_plot_dir = plot_root / config_values["case_id"]
    case_plot_dir.mkdir(parents=True, exist_ok=True)
    sun_zenith_deg = float(config_values["solar_zenith_deg"])
    sun_azimuth_deg = float(config_values["solar_azimuth_deg"])

    if np.any(comparison["reference_dop"] > 0.0):
        shared_max = float(
            max(
                np.max(comparison["reference_dop"]),
                np.max(comparison["model_dop"]),
                1.0e-12,
            )
        )
        plot_measurement_scatter(
            zenith_deg=comparison["zenith_deg"],
            azimuth_deg=comparison["absolute_azimuth_deg"],
            values=comparison["reference_dop"],
            title="Measurement Reference DoLP",
            colorbar_label="DoLP",
            save_path=case_plot_dir / "reference_dop.png",
            sun_zenith_deg=sun_zenith_deg,
            sun_azimuth_deg=sun_azimuth_deg,
            cmap="viridis",
            vmin=0.0,
            vmax=shared_max,
        )
        plot_measurement_scatter(
            zenith_deg=comparison["zenith_deg"],
            azimuth_deg=comparison["absolute_azimuth_deg"],
            values=comparison["model_dop"],
            title="Model DoLP",
            colorbar_label="DoLP",
            save_path=case_plot_dir / "model_dop.png",
            sun_zenith_deg=sun_zenith_deg,
            sun_azimuth_deg=sun_azimuth_deg,
            cmap="viridis",
            vmin=0.0,
            vmax=shared_max,
        )
        plot_measurement_scatter(
            zenith_deg=comparison["zenith_deg"],
            azimuth_deg=comparison["absolute_azimuth_deg"],
            values=comparison["dop_abs_error"],
            title="Measurement Absolute DoLP Error",
            colorbar_label="|DoLP model - reference|",
            save_path=case_plot_dir / "dop_abs_error.png",
            sun_zenith_deg=sun_zenith_deg,
            sun_azimuth_deg=sun_azimuth_deg,
            cmap="magma",
            vmin=0.0,
            vmax=float(np.max(comparison["dop_abs_error"])),
        )
        if has_column(comparison, "signed_dolp_bias"):
            max_abs_bias = float(np.max(np.abs(comparison["signed_dolp_bias"])))
            plot_measurement_scatter(
                zenith_deg=comparison["zenith_deg"],
                azimuth_deg=comparison["absolute_azimuth_deg"],
                values=comparison["signed_dolp_bias"],
                title="Measurement Signed DoLP Bias",
                colorbar_label="DoLP model - reference",
                save_path=case_plot_dir / "dop_signed_bias.png",
                sun_zenith_deg=sun_zenith_deg,
                sun_azimuth_deg=sun_azimuth_deg,
                cmap="coolwarm",
                vmin=-max_abs_bias,
                vmax=max_abs_bias,
            )

    if has_column(comparison, "reference_aop_deg") and has_column(comparison, "model_aop_deg"):
        shared_aop_min = float(
            min(
                np.min(comparison["reference_aop_deg"]),
                np.min(comparison["model_aop_deg"]),
            )
        )
        shared_aop_max = float(
            max(
                np.max(comparison["reference_aop_deg"]),
                np.max(comparison["model_aop_deg"]),
            )
        )
        if shared_aop_max - shared_aop_min < 1.0e-9:
            shared_aop_min = -90.0
            shared_aop_max = 90.0

        plot_measurement_scatter(
            zenith_deg=comparison["zenith_deg"],
            azimuth_deg=comparison["absolute_azimuth_deg"],
            values=comparison["reference_aop_deg"],
            title="Measurement Reference AoP",
            colorbar_label="AoP (deg)",
            save_path=case_plot_dir / "reference_aop_deg.png",
            sun_zenith_deg=sun_zenith_deg,
            sun_azimuth_deg=sun_azimuth_deg,
            cmap="twilight_shifted",
            vmin=shared_aop_min,
            vmax=shared_aop_max,
        )
        plot_measurement_scatter(
            zenith_deg=comparison["zenith_deg"],
            azimuth_deg=comparison["absolute_azimuth_deg"],
            values=comparison["model_aop_deg"],
            title="Model AoP",
            colorbar_label="AoP (deg)",
            save_path=case_plot_dir / "model_aop_deg.png",
            sun_zenith_deg=sun_zenith_deg,
            sun_azimuth_deg=sun_azimuth_deg,
            cmap="twilight_shifted",
            vmin=shared_aop_min,
            vmax=shared_aop_max,
        )
        if has_column(comparison, "aop_abs_error_deg"):
            plot_measurement_scatter(
                zenith_deg=comparison["zenith_deg"],
                azimuth_deg=comparison["absolute_azimuth_deg"],
                values=comparison["aop_abs_error_deg"],
                title="Measurement Absolute AoP Error",
                colorbar_label="|AoP model - reference| (deg)",
                save_path=case_plot_dir / "aop_abs_error_deg.png",
                sun_zenith_deg=sun_zenith_deg,
                sun_azimuth_deg=sun_azimuth_deg,
                cmap="magma",
                vmin=0.0,
                vmax=float(np.max(comparison["aop_abs_error_deg"])),
            )

    if np.any(comparison["reference_intensity"] > 0.0):
        plot_measurement_scatter(
            zenith_deg=comparison["zenith_deg"],
            azimuth_deg=comparison["absolute_azimuth_deg"],
            values=comparison["normalized_reference"],
            title="Measurement Reference Normalized Intensity",
            colorbar_label="Normalized intensity",
            save_path=case_plot_dir / "reference_intensity_norm.png",
            sun_zenith_deg=sun_zenith_deg,
            sun_azimuth_deg=sun_azimuth_deg,
            cmap="viridis",
            vmin=0.0,
            vmax=1.0,
        )
        plot_measurement_scatter(
            zenith_deg=comparison["zenith_deg"],
            azimuth_deg=comparison["absolute_azimuth_deg"],
            values=comparison["normalized_model"],
            title="Model Normalized Intensity",
            colorbar_label="Normalized intensity",
            save_path=case_plot_dir / "model_intensity_norm.png",
            sun_zenith_deg=sun_zenith_deg,
            sun_azimuth_deg=sun_azimuth_deg,
            cmap="viridis",
            vmin=0.0,
            vmax=1.0,
        )
        plot_measurement_scatter(
            zenith_deg=comparison["zenith_deg"],
            azimuth_deg=comparison["absolute_azimuth_deg"],
            values=comparison["normalized_model"] - comparison["normalized_reference"],
            title="Measurement Intensity Difference",
            colorbar_label="Normalized model - reference",
            save_path=case_plot_dir / "intensity_difference.png",
            sun_zenith_deg=sun_zenith_deg,
            sun_azimuth_deg=sun_azimuth_deg,
            cmap="coolwarm",
        )

    order_fraction_specs = [
        ("first_frac", "First-Order Fraction", "first_order_fraction.png"),
        ("second_frac", "Second-Order Fraction", "second_order_fraction.png"),
        ("higher_frac", "Third+ Fraction", "higher_order_fraction.png"),
        ("second_rr_frac", "Second-Order Rayleigh→Rayleigh Fraction", "second_rr_fraction.png"),
        ("second_ar_frac", "Second-Order Aerosol→Rayleigh Fraction", "second_ar_fraction.png"),
        ("second_ra_frac", "Second-Order Rayleigh→Aerosol Fraction", "second_ra_fraction.png"),
        ("second_aa_frac", "Second-Order Aerosol→Aerosol Fraction", "second_aa_fraction.png"),
    ]
    for key, title, filename in order_fraction_specs:
        if not has_column(comparison, key):
            continue
        if not np.any(comparison[key] > 0.0):
            continue
        plot_measurement_scatter(
            zenith_deg=comparison["zenith_deg"],
            azimuth_deg=comparison["absolute_azimuth_deg"],
            values=comparison[key],
            title=title,
            colorbar_label="Fraction",
            save_path=case_plot_dir / filename,
            sun_zenith_deg=sun_zenith_deg,
            sun_azimuth_deg=sun_azimuth_deg,
            cmap="viridis",
            vmin=0.0,
            vmax=1.0,
        )

    summary = {
        "case_id": config_values["case_id"],
        "config_path": str(config_path),
        "report_path": str(report_path),
        "comparison_path": str(comparison_path),
        "plot_dir": str(case_plot_dir),
        "metrics": report_values,
    }
    summary_path = case_plot_dir / "measurement_plot_metadata.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    return summary
