import argparse
import json
import math
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO_ROOT / "monte_carlo_cpp" / "config" / "paper_cases" / "frozen_marseille_twilight_20220815_191413z_measurement_interactive.cfg"


def find_default_runner() -> Path:
    for runner_name in ("MeasurementCaseRunner.exe", "MeasurementCaseRunner"):
        for build_dir in (
            REPO_ROOT / "monte_carlo_cpp" / "build_current",
            REPO_ROOT / "monte_carlo_cpp" / "build",
        ):
            candidate = build_dir / runner_name
            if candidate.exists():
                return candidate
    return REPO_ROOT / "monte_carlo_cpp" / "build_current" / "MeasurementCaseRunner"


DEFAULT_RUNNER = find_default_runner()


def parse_config(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def write_config(path: Path, values: dict[str, str]) -> None:
    lines = [f"{key}={value}" for key, value in values.items()]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_key_value_text(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def as_float(values: dict[str, str], key: str) -> float:
    raw = values.get(key, "nan").strip()
    if raw.lower() == "nan":
        return math.nan
    return float(raw)


def run_case(runner: Path, config_path: Path) -> tuple[dict[str, str], str]:
    completed = subprocess.run(
        [str(runner), str(config_path)],
        cwd=runner.parent,
        capture_output=True,
        text=True,
        check=True,
    )
    stdout = completed.stdout
    report_path = None
    for line in stdout.splitlines():
        if line.startswith("report_txt="):
            report_path = Path(line.split("=", 1)[1].strip())
            break
    if report_path is None:
        raise RuntimeError(f"Measurement runner did not emit report_txt for {config_path}")
    if not report_path.is_absolute():
        report_path = (runner.parent / report_path).resolve()
    return parse_key_value_text(report_path), stdout


def main() -> None:
    parser = argparse.ArgumentParser(description="Run land-vs-ocean Marseille surface sensitivity")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--ocean-wind-speed", type=float, default=5.0)
    args = parser.parse_args()

    config_path = args.config.resolve()
    runner_path = args.runner.resolve()
    if not config_path.exists():
        raise FileNotFoundError(config_path)
    if not runner_path.exists():
        raise FileNotFoundError(runner_path)

    base_values = parse_config(config_path)
    config_dir = config_path.parent
    output_dir = (config_dir / base_values["output_dir"]).resolve()
    report_dir = output_dir / "measurement_case_reports"
    report_dir.mkdir(parents=True, exist_ok=True)

    case_id = base_values.get("case_id", config_path.stem)
    temp_land = config_dir / f"{case_id}_surface_land.tmp.cfg"
    temp_ocean = config_dir / f"{case_id}_surface_ocean.tmp.cfg"

    land_values = dict(base_values)
    land_values["case_id"] = f"{case_id}_surface_land"
    land_values["surface_model"] = "lambertian_land"

    ocean_values = dict(base_values)
    ocean_values["case_id"] = f"{case_id}_surface_ocean"
    ocean_values["surface_model"] = "coxmunk_ocean"
    ocean_values["ocean_wind_speed_m_s"] = f"{args.ocean_wind_speed:.3f}"

    try:
        write_config(temp_land, land_values)
        write_config(temp_ocean, ocean_values)

        land_report, _ = run_case(runner_path, temp_land)
        ocean_report, _ = run_case(runner_path, temp_ocean)
    finally:
        temp_land.unlink(missing_ok=True)
        temp_ocean.unlink(missing_ok=True)

    metric_keys = [
        "normalized_rmse",
        "median_dolp_abs",
        "p95_dolp_abs",
        "median_aop_deg",
        "p95_aop_deg",
        "solar_vertical_signed_dolp_bias",
        "region_bright_horizon_arc_mean_normalized_model",
        "region_bright_horizon_arc_median_dolp_abs",
        "region_bright_horizon_arc_mean_signed_dolp_bias",
        "region_solar_vertical_lowelev_mean_normalized_model",
        "region_solar_vertical_lowelev_median_dolp_abs",
        "region_solar_vertical_lowelev_mean_signed_dolp_bias",
    ]

    comparison: dict[str, dict[str, float]] = {}
    max_abs_delta = 0.0
    for key in metric_keys:
        land_value = as_float(land_report, key)
        ocean_value = as_float(ocean_report, key)
        delta = ocean_value - land_value
        if math.isfinite(delta):
            max_abs_delta = max(max_abs_delta, abs(delta))
        comparison[key] = {
            "land": land_value,
            "ocean": ocean_value,
            "delta_ocean_minus_land": delta,
        }

    negligible = max_abs_delta < 1.0e-6
    conclusion = (
        "surface sensitivity is negligible for the current Marseille twilight path; keep lambertian_land as the default"
        if negligible
        else "surface sensitivity is non-negligible; inspect bright-horizon and low-elevation region deltas before changing the default surface model"
    )

    summary = {
        "base_config": str(config_path),
        "runner": str(runner_path),
        "land_case_id": land_values["case_id"],
        "ocean_case_id": ocean_values["case_id"],
        "ocean_wind_speed_m_s": args.ocean_wind_speed,
        "max_abs_delta": max_abs_delta,
        "negligible": negligible,
        "conclusion": conclusion,
        "comparison": comparison,
    }

    summary_json = report_dir / f"{case_id}_surface_sensitivity_summary.json"
    summary_txt = report_dir / f"{case_id}_surface_sensitivity_summary.txt"
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    with summary_txt.open("w", encoding="utf-8") as output:
        output.write(f"base_config={config_path}\n")
        output.write(f"land_case_id={land_values['case_id']}\n")
        output.write(f"ocean_case_id={ocean_values['case_id']}\n")
        output.write(f"ocean_wind_speed_m_s={args.ocean_wind_speed:.3f}\n")
        output.write(f"max_abs_delta={max_abs_delta:.12g}\n")
        output.write(f"negligible={str(negligible).lower()}\n")
        output.write(f"conclusion={conclusion}\n")
        for key in metric_keys:
            values = comparison[key]
            output.write(f"{key}_land={values['land']}\n")
            output.write(f"{key}_ocean={values['ocean']}\n")
            output.write(f"{key}_delta_ocean_minus_land={values['delta_ocean_minus_land']}\n")

    print(f"summary_json={summary_json}")
    print(f"summary_txt={summary_txt}")
    print(f"max_abs_delta={max_abs_delta}")
    print(f"negligible={str(negligible).lower()}")
    print(f"conclusion={conclusion}")


if __name__ == "__main__":
    main()
