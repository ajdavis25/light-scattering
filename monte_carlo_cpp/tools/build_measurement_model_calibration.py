from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def wrap_half_turn_deg(value: float) -> float:
    while value <= -90.0:
        value += 180.0
    while value > 90.0:
        value -= 180.0
    return value


def dolp_from_stokes(intensity: float, q_value: float, u_value: float) -> float:
    if intensity <= 0.0:
        return 0.0
    return math.hypot(q_value, u_value) / max(intensity, 1.0e-12)


def aop_from_stokes(q_value: float, u_value: float) -> float:
    return wrap_half_turn_deg(0.5 * math.degrees(math.atan2(u_value, q_value)))


def float_value(row: dict[str, str], key: str, default: float = 0.0) -> float:
    value = row.get(key, "")
    if value in ("", None):
        return default
    return float(value)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a frozen row-wise model-quality calibration from a completed measurement comparison CSV. "
            "This is a calibration artifact, not an independent validation reference."
        )
    )
    parser.add_argument("comparison_csv", type=Path, help="Completed raw comparison CSV to calibrate against.")
    parser.add_argument("output_csv", type=Path, help="Calibration CSV to write.")
    parser.add_argument(
        "--max-dolp-scale",
        type=float,
        default=25.0,
        help="Safety cap for DoLP scaling when the raw model DoLP is tiny.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    with args.comparison_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise RuntimeError(f"Comparison CSV is missing a header: {args.comparison_csv}")
        rows = [dict(row) for row in reader]

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "index",
        "zenith_deg",
        "relative_azimuth_deg",
        "intensity_gain",
        "dolp_scale",
        "aop_offset_deg",
        "source_reference_intensity",
        "source_model_intensity",
        "source_reference_dop",
        "source_model_dop",
        "source_reference_aop_deg",
        "source_model_aop_deg",
    ]
    with args.output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            index = int(float(row["index"]))
            reference_intensity = float_value(row, "reference_intensity")
            model_intensity = float_value(row, "raw_model_intensity", float_value(row, "model_intensity"))
            model_q = float_value(row, "raw_model_q", float_value(row, "model_q"))
            model_u = float_value(row, "raw_model_u", float_value(row, "model_u"))
            model_dop = dolp_from_stokes(model_intensity, model_q, model_u)
            model_aop = aop_from_stokes(model_q, model_u)
            reference_dop = float_value(row, "reference_dop")
            reference_aop = float_value(row, "reference_aop_deg")
            intensity_gain = reference_intensity / max(model_intensity, 1.0e-300)
            if model_dop <= 1.0e-12:
                dolp_scale = args.max_dolp_scale if reference_dop > 0.0 else 1.0
            else:
                dolp_scale = reference_dop / model_dop
            dolp_scale = max(0.0, min(args.max_dolp_scale, dolp_scale))
            writer.writerow({
                "index": index,
                "zenith_deg": row["zenith_deg"],
                "relative_azimuth_deg": row["relative_azimuth_deg"],
                "intensity_gain": f"{intensity_gain:.17e}",
                "dolp_scale": f"{dolp_scale:.17e}",
                "aop_offset_deg": f"{wrap_half_turn_deg(reference_aop - model_aop):.17e}",
                "source_reference_intensity": f"{reference_intensity:.17e}",
                "source_model_intensity": f"{model_intensity:.17e}",
                "source_reference_dop": f"{reference_dop:.17e}",
                "source_model_dop": f"{model_dop:.17e}",
                "source_reference_aop_deg": f"{reference_aop:.17e}",
                "source_model_aop_deg": f"{model_aop:.17e}",
            })

    print(f"wrote_calibration_csv={args.output_csv}")
    print(f"rows={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
