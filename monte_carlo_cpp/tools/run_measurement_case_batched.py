from __future__ import annotations

import argparse
from collections import deque
import csv
import json
import math
import os
import selectors
import subprocess
import sys
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
COMPARISON_FIELDS = [
    "index",
    "zenith_deg",
    "relative_azimuth_deg",
    "absolute_azimuth_deg",
    "reference_intensity",
    "model_intensity",
    "normalized_reference",
    "normalized_model",
    "reference_dop",
    "model_dop",
    "dop_abs_error",
    "reference_aop_deg",
    "model_aop_deg",
    "aop_abs_error_deg",
    "signed_dolp_bias",
    "first_frac",
    "second_frac",
    "higher_frac",
    "second_rr_frac",
    "second_ar_frac",
    "second_ra_frac",
    "second_aa_frac",
    "model_q",
    "model_u",
    "model_v",
]
REGION_FIELDS = [
    "region_name",
    "count",
    "dop_count",
    "aop_count",
    "median_dolp_abs",
    "p95_dolp_abs",
    "median_aop_deg",
    "p95_aop_deg",
    "mean_signed_dolp_bias",
    "mean_first_frac",
    "mean_second_frac",
    "mean_higher_frac",
    "mean_second_rr_frac",
    "mean_second_ar_frac",
    "mean_second_ra_frac",
    "mean_second_aa_frac",
    "mean_normalized_reference",
    "mean_normalized_model",
    "max_dolp_abs",
    "max_dolp_zenith_deg",
    "max_dolp_relative_azimuth_deg",
]


def candidate_build_dirs() -> list[Path]:
    build_dirs: list[Path] = []
    env_build_dir = os.environ.get("MONTE_CARLO_BUILD_DIR", "").strip()
    if env_build_dir:
        build_dirs.append(Path(env_build_dir).expanduser().resolve())
    build_dirs.extend([
        REPO_ROOT / "monte_carlo_cpp" / "build_current",
        REPO_ROOT / "monte_carlo_cpp" / "build",
    ])

    unique: list[Path] = []
    seen: set[str] = set()
    for build_dir in build_dirs:
        key = str(build_dir)
        if key in seen:
            continue
        seen.add(key)
        unique.append(build_dir)
    return unique


def find_runner(name: str) -> Path:
    names = [name]
    if name.endswith(".exe"):
        names.append(name[:-4])
    else:
        names.append(f"{name}.exe")

    for build_dir in candidate_build_dirs():
        for runner_name in names:
            candidate = build_dir / runner_name
            if candidate.exists():
                return candidate
    raise FileNotFoundError(
        f"Could not find {name} or platform-specific variant in {', '.join(str(path) for path in candidate_build_dirs())}."
    )


def parse_config(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def rewrite_config(base_text: str, overrides: dict[str, str]) -> str:
    replaced: set[str] = set()
    output_lines: list[str] = []
    for raw_line in base_text.splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#") or "=" not in raw_line:
            output_lines.append(raw_line)
            continue
        key, _ = raw_line.split("=", 1)
        key = key.strip()
        if key in overrides:
            output_lines.append(f"{key}={overrides[key]}")
            replaced.add(key)
        else:
            output_lines.append(raw_line)
    for key, value in overrides.items():
        if key not in replaced:
            output_lines.append(f"{key}={value}")
    return "\n".join(output_lines) + "\n"


def load_reference_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise RuntimeError(f"Reference CSV is missing a header: {path}")
        fieldnames = list(reader.fieldnames)
        rows = [dict(row) for row in reader]
    if "original_index" not in fieldnames:
        fieldnames.append("original_index")
    cleaned_rows: list[dict[str, str]] = []
    for index, row in enumerate(rows):
        cleaned = {key: (value if value is not None else "") for key, value in row.items()}
        original_index_raw = cleaned.get("original_index", "").strip()
        if original_index_raw:
            try:
                original_index = int(float(original_index_raw))
            except ValueError as exc:
                raise RuntimeError(
                    f"Invalid original_index value in reference CSV {path}: {original_index_raw!r}"
                ) from exc
        else:
            original_index = index
        cleaned["original_index"] = str(original_index)
        cleaned_rows.append(cleaned)
    return fieldnames, cleaned_rows


def resolve_config_path(config_path: Path, value: str) -> Path:
    candidate = Path(value)
    if candidate.is_absolute():
        return candidate
    return (config_path.parent / candidate).resolve()


def read_comparison_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return [dict(row) for row in reader]


def load_completed_indices(path: Path) -> set[int]:
    if not path.exists() or path.stat().st_size == 0:
        return set()
    completed: set[int] = set()
    for row in read_comparison_rows(path):
        completed.add(int(float(row["index"])))
    return completed


def append_rows(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    write_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COMPARISON_FIELDS)
        if write_header:
            writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in COMPARISON_FIELDS})


def float_value(row: dict[str, str], key: str, default: float = 0.0) -> float:
    value = row.get(key, "")
    if value in ("", None):
        return default
    return float(value)


def wrap_half_turn_deg(value: float) -> float:
    while value <= -90.0:
        value += 180.0
    while value > 90.0:
        value -= 180.0
    return value


def aop_difference_deg(lhs: float, rhs: float) -> float:
    return abs(wrap_half_turn_deg(lhs - rhs))


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    index = int(max(0, min(len(ordered) - 1, fraction * (len(ordered) - 1))))
    return ordered[index]


def angular_separation_deg(zenith_a: float, azimuth_a: float, zenith_b: float, azimuth_b: float) -> float:
    def vector(zenith_deg: float, azimuth_deg: float) -> tuple[float, float, float]:
        zen = math.radians(zenith_deg)
        az = math.radians(azimuth_deg)
        return (
            math.sin(zen) * math.cos(az),
            math.sin(zen) * math.sin(az),
            math.cos(zen),
        )

    ax, ay, az = vector(zenith_a, azimuth_a)
    bx, by, bz = vector(zenith_b, azimuth_b)
    cosine = max(-1.0, min(1.0, ax * bx + ay * by + az * bz))
    return math.degrees(math.acos(cosine))


def in_solar_vertical_midzen(row: dict[str, float]) -> bool:
    return 240.0 <= row["relative_azimuth_deg"] <= 300.0 and 30.0 <= row["zenith_deg"] <= 70.0


def in_solar_vertical_lowelev(row: dict[str, float]) -> bool:
    return 240.0 <= row["relative_azimuth_deg"] <= 300.0 and 70.0 <= row["zenith_deg"] <= 88.0


def in_antisolar_midzen(row: dict[str, float]) -> bool:
    return (row["relative_azimuth_deg"] <= 30.0 or row["relative_azimuth_deg"] >= 330.0) and 30.0 <= row["zenith_deg"] <= 70.0


def in_bright_horizon_arc(row: dict[str, float]) -> bool:
    return 320.0 <= row["relative_azimuth_deg"] <= 350.0 and 65.0 <= row["zenith_deg"] <= 88.0


def in_near_zenith(row: dict[str, float]) -> bool:
    return row["zenith_deg"] <= 25.0


def in_aop_flip_45(row: dict[str, float]) -> bool:
    return 30.0 <= row["relative_azimuth_deg"] <= 60.0 and 5.0 <= row["zenith_deg"] <= 85.0


def in_aop_flip_215(row: dict[str, float]) -> bool:
    return 200.0 <= row["relative_azimuth_deg"] <= 230.0 and 5.0 <= row["zenith_deg"] <= 85.0


def summarize_region(name: str, rows: list[dict[str, float]], predicate) -> dict[str, float | str]:
    selected = [row for row in rows if predicate(row)]
    dop_rows = [row for row in selected if row["has_reference_dop"] > 0.5]
    aop_rows = [row for row in selected if row["has_reference_aop"] > 0.5 and row["reference_dop"] >= 0.15]
    worst = max(dop_rows, key=lambda row: row["dop_abs_error"], default=None)
    return {
        "region_name": name,
        "count": len(selected),
        "dop_count": len(dop_rows),
        "aop_count": len(aop_rows),
        "median_dolp_abs": percentile([row["dop_abs_error"] for row in dop_rows], 0.5),
        "p95_dolp_abs": percentile([row["dop_abs_error"] for row in dop_rows], 0.95),
        "median_aop_deg": percentile([row["aop_abs_error_deg"] for row in aop_rows], 0.5),
        "p95_aop_deg": percentile([row["aop_abs_error_deg"] for row in aop_rows], 0.95),
        "mean_signed_dolp_bias": (
            sum(row["signed_dolp_bias"] for row in dop_rows) / len(dop_rows) if dop_rows else math.nan
        ),
        "mean_first_frac": sum(row["first_frac"] for row in selected) / len(selected) if selected else math.nan,
        "mean_second_frac": sum(row["second_frac"] for row in selected) / len(selected) if selected else math.nan,
        "mean_higher_frac": sum(row["higher_frac"] for row in selected) / len(selected) if selected else math.nan,
        "mean_second_rr_frac": sum(row["second_rr_frac"] for row in selected) / len(selected) if selected else math.nan,
        "mean_second_ar_frac": sum(row["second_ar_frac"] for row in selected) / len(selected) if selected else math.nan,
        "mean_second_ra_frac": sum(row["second_ra_frac"] for row in selected) / len(selected) if selected else math.nan,
        "mean_second_aa_frac": sum(row["second_aa_frac"] for row in selected) / len(selected) if selected else math.nan,
        "mean_normalized_reference": (
            sum(row["normalized_reference"] for row in selected) / len(selected) if selected else math.nan
        ),
        "mean_normalized_model": sum(row["normalized_model"] for row in selected) / len(selected) if selected else math.nan,
        "max_dolp_abs": worst["dop_abs_error"] if worst else math.nan,
        "max_dolp_zenith_deg": worst["zenith_deg"] if worst else math.nan,
        "max_dolp_relative_azimuth_deg": worst["relative_azimuth_deg"] if worst else math.nan,
    }


def write_progress(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def read_checkpoint_progress(path: Path) -> tuple[int, int, int]:
    if not path.exists():
        return (0, 0, 0)
    has_first_order = 0
    has_second_order = 0
    completed = 0
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        if key == "has_first_order":
            has_first_order = 1 if int(float(value.strip())) != 0 else 0
        elif key == "has_second_order":
            has_second_order = 1 if int(float(value.strip())) != 0 else 0
        elif key == "higher_completed_samples":
            completed = int(float(value.strip()))
    return (has_first_order, has_second_order, completed)


def single_direction_case_id(parent_case_id: str, original_index: int) -> str:
    return f"{parent_case_id}__batch_{original_index:04d}_{original_index:04d}"


def make_direction_work_item(
    parent_case_id: str,
    row: dict[str, str],
    report_dir: Path,
    work_dir: Path,
) -> dict[str, object]:
    original_index = int(row["original_index"])
    case_id = single_direction_case_id(parent_case_id, original_index)
    return {
        "row": row,
        "original_index": original_index,
        "case_id": case_id,
        "reference_path": work_dir / f"{case_id}_reference.csv",
        "config_path": work_dir / f"{case_id}.cfg",
        "comparison_csv": report_dir / f"{case_id}_comparison.csv",
        "report_txt": report_dir / f"{case_id}.txt",
        "region_csv": report_dir / f"{case_id}_region_summary.csv",
        "checkpoint_path": work_dir / f"{case_id}_checkpoint.txt",
        "prepared": False,
    }


def materialize_direction_work_item(
    item: dict[str, object],
    fieldnames: list[str],
    config_path: Path,
    config_text: str,
    config_values: dict[str, str],
    output_dir: Path,
    *,
    reset_outputs: bool,
) -> None:
    reference_path = item["reference_path"]
    config_path_out = item["config_path"]
    comparison_csv = item["comparison_csv"]
    report_txt = item["report_txt"]
    region_csv = item["region_csv"]
    checkpoint_path = item["checkpoint_path"]
    row = item["row"]
    case_id = item["case_id"]

    if reset_outputs:
        for stale_path in (comparison_csv, report_txt, region_csv, checkpoint_path):
            if stale_path.exists():
                stale_path.unlink()

    with Path(reference_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(row)

    overrides = {
        "case_id": str(case_id),
        "measurement_reference_csv": str(reference_path),
        "output_dir": str(output_dir.resolve()),
    }
    for key, value in config_values.items():
        if not value:
            continue
        if key in overrides:
            continue
        if key.endswith("_csv") or key.endswith("_json") or key.endswith("_cfg") or key.endswith("_config"):
            overrides[key] = str(resolve_config_path(config_path, value))

    Path(config_path_out).write_text(
        rewrite_config(config_text, overrides),
        encoding="utf-8",
    )
    item["prepared"] = True


def launch_direction_worker(
    measurement_runner: Path,
    item: dict[str, object],
    higher_order_block_size: int,
) -> subprocess.Popen[str]:
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    process = subprocess.Popen(
        [
            str(measurement_runner),
            str(item["config_path"]),
            "--checkpoint-state",
            str(item["checkpoint_path"]),
            "--higher-order-block-size",
            str(higher_order_block_size),
        ],
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        env=env,
    )
    assert process.stdout is not None
    return process


def finalize_outputs(
    rows: list[dict[str, str]],
    config_values: dict[str, str],
    report_dir: Path,
    case_id: str,
    partial_rows_csv: Path,
    started: float,
) -> None:
    parsed_rows: list[dict[str, float]] = []
    for row in rows:
        has_reference_dop = row.get("reference_dop", "").strip() != ""
        has_reference_aop = row.get("reference_aop_deg", "").strip() != ""
        parsed_rows.append({
            "index": float_value(row, "index"),
            "zenith_deg": float_value(row, "zenith_deg"),
            "relative_azimuth_deg": float_value(row, "relative_azimuth_deg"),
            "absolute_azimuth_deg": float_value(row, "absolute_azimuth_deg"),
            "reference_intensity": float_value(row, "reference_intensity"),
            "model_intensity": float_value(row, "model_intensity"),
            "reference_dop": float_value(row, "reference_dop"),
            "model_dop": float_value(row, "model_dop"),
            "reference_aop_deg": float_value(row, "reference_aop_deg"),
            "model_aop_deg": float_value(row, "model_aop_deg"),
            "signed_dolp_bias": float_value(row, "signed_dolp_bias"),
            "first_frac": float_value(row, "first_frac"),
            "second_frac": float_value(row, "second_frac"),
            "higher_frac": float_value(row, "higher_frac"),
            "second_rr_frac": float_value(row, "second_rr_frac"),
            "second_ar_frac": float_value(row, "second_ar_frac"),
            "second_ra_frac": float_value(row, "second_ra_frac"),
            "second_aa_frac": float_value(row, "second_aa_frac"),
            "model_q": float_value(row, "model_q"),
            "model_u": float_value(row, "model_u"),
            "model_v": float_value(row, "model_v"),
            "has_reference_dop": 1.0 if has_reference_dop else 0.0,
            "has_reference_aop": 1.0 if has_reference_aop else 0.0,
        })

    parsed_rows.sort(key=lambda row: int(row["index"]))
    reference_peak = max((row["reference_intensity"] for row in parsed_rows), default=0.0)
    model_peak = max((row["model_intensity"] for row in parsed_rows), default=0.0)

    measurement_mask_fraction = float(config_values.get("measurement_mask_fraction_of_peak", "0.05"))
    dop_errors: list[float] = []
    aop_errors: list[float] = []
    solar_vertical_bias: list[float] = []
    sum_squared = 0.0
    count = 0.0
    brightest_reference = None
    brightest_model = None

    comparison_csv = report_dir / f"{case_id}_comparison.csv"
    with comparison_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COMPARISON_FIELDS)
        writer.writeheader()
        for row in parsed_rows:
            normalized_reference = row["reference_intensity"] / max(1.0e-12, reference_peak) if reference_peak > 0.0 else 0.0
            normalized_model = row["model_intensity"] / max(1.0e-12, model_peak) if model_peak > 0.0 else 0.0
            dop_abs_error = abs(row["model_dop"] - row["reference_dop"])
            aop_abs_error = aop_difference_deg(row["model_aop_deg"], row["reference_aop_deg"])
            row["normalized_reference"] = normalized_reference
            row["normalized_model"] = normalized_model
            row["dop_abs_error"] = dop_abs_error
            row["aop_abs_error_deg"] = aop_abs_error
            if brightest_reference is None or normalized_reference > brightest_reference["normalized_reference"]:
                brightest_reference = row
            if brightest_model is None or normalized_model > brightest_model["normalized_model"]:
                brightest_model = row
            if normalized_reference >= measurement_mask_fraction:
                diff = normalized_model - normalized_reference
                sum_squared += diff * diff
                count += 1.0
                dop_errors.append(dop_abs_error)
                if 30.0 <= row["zenith_deg"] <= 70.0 and 240.0 <= row["relative_azimuth_deg"] <= 300.0:
                    solar_vertical_bias.append(row["signed_dolp_bias"])
                if row["reference_dop"] >= 0.15:
                    aop_errors.append(aop_abs_error)
            writer.writerow({
                "index": int(row["index"]),
                "zenith_deg": row["zenith_deg"],
                "relative_azimuth_deg": row["relative_azimuth_deg"],
                "absolute_azimuth_deg": row["absolute_azimuth_deg"],
                "reference_intensity": row["reference_intensity"],
                "model_intensity": row["model_intensity"],
                "normalized_reference": normalized_reference,
                "normalized_model": normalized_model,
                "reference_dop": row["reference_dop"],
                "model_dop": row["model_dop"],
                "dop_abs_error": dop_abs_error,
                "reference_aop_deg": row["reference_aop_deg"],
                "model_aop_deg": row["model_aop_deg"],
                "aop_abs_error_deg": aop_abs_error,
                "signed_dolp_bias": row["signed_dolp_bias"],
                "first_frac": row["first_frac"],
                "second_frac": row["second_frac"],
                "higher_frac": row["higher_frac"],
                "second_rr_frac": row["second_rr_frac"],
                "second_ar_frac": row["second_ar_frac"],
                "second_ra_frac": row["second_ra_frac"],
                "second_aa_frac": row["second_aa_frac"],
                "model_q": row["model_q"],
                "model_u": row["model_u"],
                "model_v": row["model_v"],
            })

    region_rows = [
        summarize_region("solar_vertical_midzen", parsed_rows, in_solar_vertical_midzen),
        summarize_region("solar_vertical_lowelev", parsed_rows, in_solar_vertical_lowelev),
        summarize_region("antisolar_midzen", parsed_rows, in_antisolar_midzen),
        summarize_region("bright_horizon_arc", parsed_rows, in_bright_horizon_arc),
        summarize_region("near_zenith", parsed_rows, in_near_zenith),
        summarize_region("aop_flip_relaz_45_sector", parsed_rows, in_aop_flip_45),
        summarize_region("aop_flip_relaz_215_sector", parsed_rows, in_aop_flip_215),
    ]
    region_summary_csv = report_dir / f"{case_id}_region_summary.csv"
    with region_summary_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REGION_FIELDS)
        writer.writeheader()
        writer.writerows(region_rows)

    report_path = report_dir / f"{case_id}.txt"
    lines = [
        f"case_id={case_id}",
        f"reference_csv={config_values.get('measurement_reference_csv', '')}",
        f"reference_points={len(parsed_rows)}",
        f"partial_rows_csv={partial_rows_csv}",
        f"comparison_csv={comparison_csv}",
        f"region_summary_csv={region_summary_csv}",
    ]
    if count > 0.0 and brightest_reference and brightest_model:
        lines.append(f"normalized_rmse={math.sqrt(sum_squared / count)}")
        lines.append(
            "brightest_location_deg="
            + str(
                angular_separation_deg(
                    brightest_reference["zenith_deg"],
                    brightest_reference["absolute_azimuth_deg"],
                    brightest_model["zenith_deg"],
                    brightest_model["absolute_azimuth_deg"],
                )
            )
        )
    lines.append(f"median_dolp_abs={percentile(dop_errors, 0.5)}")
    lines.append(f"p95_dolp_abs={percentile(dop_errors, 0.95)}")
    if parsed_rows:
        worst = max(parsed_rows, key=lambda row: row["dop_abs_error"])
        lines.append(f"max_dolp_abs={worst['dop_abs_error']}")
        lines.append(f"max_dolp_zenith_deg={worst['zenith_deg']}")
        lines.append(f"max_dolp_relative_azimuth_deg={worst['relative_azimuth_deg']}")
        lines.append(f"max_dolp_absolute_azimuth_deg={worst['absolute_azimuth_deg']}")
        lines.append(f"max_dolp_reference={worst['reference_dop']}")
        lines.append(f"max_dolp_model={worst['model_dop']}")
    lines.append(f"median_aop_deg={percentile(aop_errors, 0.5)}")
    lines.append(f"p95_aop_deg={percentile(aop_errors, 0.95)}")
    lines.append(
        "solar_vertical_signed_dolp_bias="
        + str(sum(solar_vertical_bias) / len(solar_vertical_bias) if solar_vertical_bias else math.nan)
    )
    for region in region_rows:
        prefix = f"region_{region['region_name']}_"
        for key, value in region.items():
            if key == "region_name":
                continue
            lines.append(f"{prefix}{key}={value}")
    lines.append(f"timing_total_runtime_seconds={time.time() - started}")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a measurement case in resumable exact-direction batches.")
    parser.add_argument("config", type=Path, help="Measurement-case config to run.")
    parser.add_argument(
        "--batch-size",
        type=int,
        default=8,
        help="Exact directions per batch, or max concurrent single-direction workers when checkpointing is enabled.",
    )
    parser.add_argument(
        "--higher-order-block-size",
        type=int,
        default=32,
        help="Higher-order samples per checkpointed single-direction worker invocation.",
    )
    parser.add_argument("--resume", action="store_true", help="Resume from an existing partial rows CSV if present.")
    parser.add_argument("--max-batches", type=int, default=0, help="Optional limit on batches to execute this run.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    started = time.time()
    measurement_runner = find_runner("MeasurementCaseRunner.exe")
    config_path = args.config.resolve()
    config_text = config_path.read_text(encoding="utf-8")
    config_values = parse_config(config_path)
    case_id = config_values["case_id"]
    output_dir = resolve_config_path(config_path, config_values["output_dir"])
    reference_csv = resolve_config_path(config_path, config_values["measurement_reference_csv"])
    report_dir = output_dir / "measurement_case_reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    work_dir = report_dir / "_batched_work"
    work_dir.mkdir(parents=True, exist_ok=True)
    partial_rows_csv = report_dir / f"{case_id}_batched_partial_rows.csv"
    progress_path = report_dir / f"{case_id}_batched_progress.json"
    log_path = report_dir / f"{case_id}_batched.log"

    if not args.resume:
        for path in (partial_rows_csv, progress_path, log_path):
            if path.exists():
                path.unlink()

    fieldnames, reference_rows = load_reference_rows(reference_csv)
    completed = load_completed_indices(partial_rows_csv) if args.resume else set()
    pending_rows = [row for row in reference_rows if int(row["original_index"]) not in completed]
    checkpoint_parallel_mode = args.higher_order_block_size > 0

    progress = {
        "case_id": case_id,
        "config": str(config_path),
        "measurement_runner": str(measurement_runner),
        "batch_size": args.batch_size,
        "mode": "checkpoint_parallel" if checkpoint_parallel_mode else "grouped_batch",
        "resume": args.resume,
        "completed_points": len(completed),
        "total_points": len(reference_rows),
        "remaining_points": len(pending_rows),
        "partial_rows_csv": str(partial_rows_csv),
        "log": str(log_path),
        "status": "in_progress" if pending_rows else "complete",
    }
    write_progress(progress_path, progress)

    with log_path.open("a", encoding="utf-8", buffering=1) as log_file:
        if checkpoint_parallel_mode:
            concurrency = max(1, args.batch_size)
            selected_rows = pending_rows
            if args.max_batches > 0:
                selected_rows = selected_rows[:args.max_batches * concurrency]

            queue: deque[dict[str, object]] = deque(
                make_direction_work_item(case_id, row, report_dir, work_dir) for row in selected_rows
            )
            selector = selectors.DefaultSelector()
            active: dict[str, dict[str, object]] = {}

            recovered: deque[dict[str, object]] = deque()
            while queue:
                item = queue.popleft()
                comparison_csv = Path(item["comparison_csv"])
                if args.resume and comparison_csv.exists():
                    recovered_rows = read_comparison_rows(comparison_csv)
                    append_rows(partial_rows_csv, recovered_rows)
                    completed.update(int(float(row["index"])) for row in recovered_rows)
                    progress.update({
                        "completed_points": len(completed),
                        "remaining_points": len(reference_rows) - len(completed),
                        "last_completed_batch": str(item["case_id"]),
                        "status": "in_progress" if len(completed) < len(reference_rows) else "complete",
                    })
                    write_progress(progress_path, progress)
                    recovered_line = (
                        f"[batched-measurement] recovered_completed_direction case_id={item['case_id']} "
                        f"completed_points={len(completed)}/{len(reference_rows)}"
                    )
                    print(recovered_line, flush=True)
                    log_file.write(recovered_line + "\n")
                    log_file.flush()
                    continue
                recovered.append(item)
            queue = recovered

            while active or queue:
                while len(active) < concurrency and queue:
                    item = queue.popleft()
                    if not bool(item["prepared"]):
                        materialize_direction_work_item(
                            item,
                            fieldnames,
                            config_path,
                            config_text,
                            config_values,
                            output_dir,
                            reset_outputs=not args.resume,
                        )

                    checkpoint_path = Path(item["checkpoint_path"])
                    previous_checkpoint_progress = read_checkpoint_progress(checkpoint_path)
                    previous_completed_samples = previous_checkpoint_progress[2]
                    process = launch_direction_worker(
                        measurement_runner,
                        item,
                        args.higher_order_block_size,
                    )
                    item["process"] = process
                    item["stdout_closed"] = False
                    item["last_completed_samples"] = previous_completed_samples
                    item["last_checkpoint_progress"] = previous_checkpoint_progress
                    active[str(item["case_id"])] = item
                    selector.register(process.stdout, selectors.EVENT_READ, data=str(item["case_id"]))

                    launch_line = (
                        f"=== direction worker started case_id={item['case_id']} "
                        f"direction_index={item['original_index']} "
                        f"completed_samples={previous_completed_samples} "
                        f"block_size={args.higher_order_block_size} ==="
                    )
                    print(launch_line, flush=True)
                    log_file.write(launch_line + "\n")
                    log_file.flush()

                if not active:
                    break

                for key, _ in selector.select(timeout=1.0):
                    active_case_id = str(key.data)
                    item = active.get(active_case_id)
                    if item is None:
                        continue
                    stream = key.fileobj
                    line = stream.readline()
                    if line:
                        print(line, end="", flush=True)
                        log_file.write(line)
                        log_file.flush()
                    else:
                        selector.unregister(stream)
                        stream.close()
                        item["stdout_closed"] = True

                finished_case_ids: list[str] = []
                for active_case_id, item in list(active.items()):
                    process = item["process"]
                    if process.poll() is None or not bool(item["stdout_closed"]):
                        continue

                    exit_line = (
                        f"=== direction worker exit_code={process.returncode} case_id={item['case_id']} "
                        f"direction_index={item['original_index']} ==="
                    )
                    print(exit_line, flush=True)
                    log_file.write(exit_line + "\n")
                    log_file.flush()
                    if process.returncode != 0:
                        return int(process.returncode)

                    comparison_csv = Path(item["comparison_csv"])
                    if comparison_csv.exists():
                        batch_comparison_rows = read_comparison_rows(comparison_csv)
                        append_rows(partial_rows_csv, batch_comparison_rows)
                        completed.update(int(float(row["index"])) for row in batch_comparison_rows)
                        progress.update({
                            "completed_points": len(completed),
                            "remaining_points": len(reference_rows) - len(completed),
                            "last_completed_batch": str(item["case_id"]),
                            "status": "in_progress" if len(completed) < len(reference_rows) else "complete",
                        })
                        write_progress(progress_path, progress)
                        completed_line = (
                            f"[batched-measurement] completed_points={len(completed)}/{len(reference_rows)} "
                            f"last_batch={item['case_id']} partial_rows_csv={partial_rows_csv}"
                        )
                        print(completed_line, flush=True)
                        log_file.write(completed_line + "\n")
                        log_file.flush()
                    else:
                        checkpoint_path = Path(item["checkpoint_path"])
                        checkpoint_progress = read_checkpoint_progress(checkpoint_path)
                        completed_samples = checkpoint_progress[2]
                        checkpoint_line = (
                            f"[batched-measurement] checkpoint_first_order={checkpoint_progress[0]} "
                            f"checkpoint_second_order={checkpoint_progress[1]} "
                            f"checkpoint_samples={completed_samples} "
                            f"block_size={args.higher_order_block_size} batch_case_id={item['case_id']}"
                        )
                        print(checkpoint_line, flush=True)
                        log_file.write(checkpoint_line + "\n")
                        log_file.flush()
                        if checkpoint_progress == tuple(item["last_checkpoint_progress"]):
                            print(
                                f"Checkpoint did not advance for {item['case_id']}: "
                                f"{checkpoint_progress} "
                                f"(previous {item['last_checkpoint_progress']}).",
                                file=sys.stderr,
                            )
                            return 1
                        item["last_completed_samples"] = completed_samples
                        progress.update({
                            "last_checkpoint_batch": str(item["case_id"]),
                            "last_checkpoint_completed_samples": completed_samples,
                            "status": "in_progress",
                        })
                        write_progress(progress_path, progress)
                        queue.append(item)

                    finished_case_ids.append(active_case_id)

                for active_case_id in finished_case_ids:
                    active.pop(active_case_id, None)
        else:
            executed_batches = 0
            for batch_index, start in enumerate(range(0, len(pending_rows), args.batch_size), start=1):
                if args.max_batches > 0 and executed_batches >= args.max_batches:
                    break
                batch_rows = pending_rows[start:start + args.batch_size]
                batch_indices = [int(row["original_index"]) for row in batch_rows]
                batch_case_id = f"batch_{batch_indices[0]:04d}_{batch_indices[-1]:04d}"
                batch_reference_path = work_dir / f"{batch_case_id}_reference.csv"
                batch_config_path = work_dir / f"{batch_case_id}.cfg"

                with batch_reference_path.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=fieldnames)
                    writer.writeheader()
                    writer.writerows(batch_rows)

                batch_overrides = {
                    "case_id": batch_case_id,
                    "measurement_reference_csv": str(batch_reference_path),
                    "output_dir": str(output_dir.resolve()),
                }
                for key, value in config_values.items():
                    if not value:
                        continue
                    if key in batch_overrides:
                        continue
                    if key.endswith("_csv") or key.endswith("_json") or key.endswith("_cfg") or key.endswith("_config"):
                        batch_overrides[key] = str(resolve_config_path(config_path, value))

                batch_config_path.write_text(
                    rewrite_config(config_text, batch_overrides),
                    encoding="utf-8",
                )

                log_file.write(f"=== batch {batch_index} started indices={batch_indices[0]}..{batch_indices[-1]} ===\n")
                log_file.flush()
                batch_comparison_csv = report_dir / f"{batch_case_id}_comparison.csv"
                batch_report_txt = report_dir / f"{batch_case_id}.txt"
                batch_region_csv = report_dir / f"{batch_case_id}_region_summary.csv"
                checkpoint_path = work_dir / f"{batch_case_id}_checkpoint.txt"
                if not args.resume:
                    for stale_path in (batch_comparison_csv, batch_report_txt, batch_region_csv, checkpoint_path):
                        if stale_path.exists():
                            stale_path.unlink()

                command = [str(measurement_runner), str(batch_config_path)]
                process = subprocess.Popen(
                    command,
                    cwd=REPO_ROOT,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                )
                assert process.stdout is not None
                for line in process.stdout:
                    print(line, end="", flush=True)
                    log_file.write(line)
                    log_file.flush()
                return_code = process.wait()
                log_file.write(f"=== batch {batch_index} exit_code={return_code} ===\n")
                log_file.flush()
                if return_code != 0:
                    return return_code

                batch_comparison_rows = read_comparison_rows(batch_comparison_csv)
                append_rows(partial_rows_csv, batch_comparison_rows)
                completed.update(int(float(row["index"])) for row in batch_comparison_rows)
                executed_batches += 1
                progress.update({
                    "completed_points": len(completed),
                    "remaining_points": len(reference_rows) - len(completed),
                    "last_completed_batch": batch_case_id,
                    "status": "in_progress" if len(completed) < len(reference_rows) else "complete",
                })
                write_progress(progress_path, progress)
                completed_line = (
                    f"[batched-measurement] completed_points={len(completed)}/{len(reference_rows)} "
                    f"last_batch={batch_case_id} partial_rows_csv={partial_rows_csv}"
                )
                print(completed_line, flush=True)
                log_file.write(completed_line + "\n")
                log_file.flush()

    all_rows = read_comparison_rows(partial_rows_csv)
    deduped = {int(float(row["index"])): row for row in all_rows}
    if len(deduped) != len(reference_rows):
        progress.update({
            "completed_points": len(deduped),
            "remaining_points": len(reference_rows) - len(deduped),
            "status": "incomplete",
        })
        write_progress(progress_path, progress)
        print(f"Incomplete batched run: {len(deduped)}/{len(reference_rows)} points available.", flush=True)
        return 0

    finalize_outputs(
        [deduped[index] for index in sorted(deduped)],
        config_values,
        report_dir,
        case_id,
        partial_rows_csv,
        started,
    )
    progress.update({
        "completed_points": len(reference_rows),
        "remaining_points": 0,
        "status": "complete",
    })
    write_progress(progress_path, progress)
    print(f"Final report written for case_id={case_id}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
