from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import requests


ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = ROOT / "monte_carlo_cpp" / "config"
PAPER_CONFIG_DIR = CONFIG_DIR / "paper_cases"
PAPER_DATA_DIR = ROOT / "monte_carlo_cpp" / "data" / "paper_cases"
SATURATION_LEVEL = 65520

STRICT_SUBSET_FILENAME = "measurement_reference_strict_subset.csv"
STRICT_SUBSET_SUMMARY_FILENAME = "measurement_reference_strict_subset_summary.json"
STRICT_SUBSET_CONFIG_SUFFIX = "_measurement_strict_subset.cfg"
STRICT_SUBSET_TARGET_COUNT = 48
PROFILE_SUBSET_FILENAME = "measurement_reference_profile_subset.csv"
PROFILE_SUBSET_SUMMARY_FILENAME = "measurement_reference_profile_subset_summary.json"
PROFILE_SUBSET_CONFIG_SUFFIX = "_measurement_profile_subset.cfg"
PROFILE_SUBSET_TARGET_COUNT = 12

STRICT_SUBSET_REGIONS = (
    ("solar_vertical_midzen", 12, lambda row: 240.0 <= row["relative_azimuth_deg"] <= 300.0 and 30.0 <= row["zenith_deg"] <= 70.0),
    ("antisolar_midzen", 8, lambda row: (row["relative_azimuth_deg"] <= 30.0 or row["relative_azimuth_deg"] >= 330.0) and 30.0 <= row["zenith_deg"] <= 70.0),
    ("bright_horizon_arc", 10, lambda row: 320.0 <= row["relative_azimuth_deg"] <= 350.0 and 65.0 <= row["zenith_deg"] <= 88.0),
    ("near_zenith", 6, lambda row: row["zenith_deg"] <= 25.0),
    ("aop_flip_relaz_45_sector", 6, lambda row: 30.0 <= row["relative_azimuth_deg"] <= 60.0 and 5.0 <= row["zenith_deg"] <= 85.0),
    ("aop_flip_relaz_215_sector", 6, lambda row: 200.0 <= row["relative_azimuth_deg"] <= 230.0 and 5.0 <= row["zenith_deg"] <= 85.0),
)

PROFILE_SUBSET_REGIONS = (
    ("solar_vertical_midzen", 2, lambda row: 240.0 <= row["relative_azimuth_deg"] <= 300.0 and 30.0 <= row["zenith_deg"] <= 70.0),
    ("antisolar_midzen", 2, lambda row: (row["relative_azimuth_deg"] <= 30.0 or row["relative_azimuth_deg"] >= 330.0) and 30.0 <= row["zenith_deg"] <= 70.0),
    ("bright_horizon_arc", 2, lambda row: 320.0 <= row["relative_azimuth_deg"] <= 350.0 and 65.0 <= row["zenith_deg"] <= 88.0),
    ("near_zenith", 2, lambda row: row["zenith_deg"] <= 25.0),
    ("aop_flip_relaz_45_sector", 2, lambda row: 30.0 <= row["relative_azimuth_deg"] <= 60.0 and 5.0 <= row["zenith_deg"] <= 85.0),
    ("aop_flip_relaz_215_sector", 2, lambda row: 200.0 <= row["relative_azimuth_deg"] <= 230.0 and 5.0 <= row["zenith_deg"] <= 85.0),
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def normalize_azimuth_deg(angle_deg: np.ndarray | float) -> np.ndarray | float:
    return np.mod(angle_deg, 360.0)


def wrap_half_turn_deg(angle_deg: np.ndarray | float) -> np.ndarray | float:
    return np.mod(angle_deg + 90.0, 180.0) - 90.0


def average_angle(angle1: np.ndarray, angle2: np.ndarray) -> np.ndarray:
    x = np.cos(angle1) + np.cos(angle2)
    y = np.sin(angle1) + np.sin(angle2)
    return np.arctan2(y, x)


def rotate_linear_stokes(
    q: np.ndarray | float,
    u: np.ndarray | float,
    angle_deg: np.ndarray | float,
) -> tuple[np.ndarray | float, np.ndarray | float]:
    angle_rad = np.radians(angle_deg)
    cos_2 = np.cos(2.0 * angle_rad)
    sin_2 = np.sin(2.0 * angle_rad)
    q_rot = q * cos_2 + u * sin_2
    u_rot = -q * sin_2 + u * cos_2
    return q_rot, u_rot


def marseille_comparison_basis_rotation_deg(
    relative_azimuth_deg: np.ndarray | float,
    rotation_z_deg: float,
) -> np.ndarray | float:
    return relative_azimuth_deg - rotation_z_deg - 90.0


def find_dataset_file(dataset_record: dict, filename: str) -> dict:
    for entry in dataset_record["data"]["latestVersion"]["files"]:
        data_file = entry["dataFile"]
        if data_file["filename"] == filename:
            return data_file
    raise RuntimeError(f"Could not find dataset file {filename}")


def download_with_resume(
    session: requests.Session,
    file_id: int,
    destination: Path,
    expected_size: int,
) -> None:
    if destination.exists() and destination.stat().st_size == expected_size:
        return

    destination.parent.mkdir(parents=True, exist_ok=True)
    partial_path = destination.with_suffix(destination.suffix + ".part")
    downloaded = partial_path.stat().st_size if partial_path.exists() else 0
    headers: dict[str, str] = {}
    mode = "wb"
    if downloaded > 0 and downloaded < expected_size:
        headers["Range"] = f"bytes={downloaded}-"
        mode = "ab"

    with session.get(
        f"https://entrepot.recherche.data.gouv.fr/api/access/datafile/{file_id}",
        headers=headers,
        stream=True,
        timeout=600,
    ) as response:
        response.raise_for_status()
        with partial_path.open(mode) as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)

    final_size = partial_path.stat().st_size
    if final_size != expected_size:
        raise RuntimeError(
            f"Downloaded size mismatch for {destination.name}: expected {expected_size}, got {final_size}"
        )
    partial_path.replace(destination)


def get_case_paths(case_id: str) -> tuple[Path, Path]:
    case_data_dir = PAPER_DATA_DIR / case_id
    case_config_path = PAPER_CONFIG_DIR / f"{case_id}.cfg"
    if not case_data_dir.exists():
        raise FileNotFoundError(f"Missing case data directory: {case_data_dir}")
    if not case_config_path.exists():
        raise FileNotFoundError(f"Missing case config: {case_config_path}")
    return case_data_dir, case_config_path


def update_or_append_line(lines: list[str], prefix: str, replacement: str) -> list[str]:
    updated = False
    result: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith(prefix):
            result.append(replacement)
            updated = True
        else:
            result.append(line)
    if not updated:
        result.append(replacement)
    return result


def remove_comment_placeholders(lines: list[str]) -> list[str]:
    return [
        line
        for line in lines
        if "Placeholder measurement config for the frozen Marseille twilight case." not in line
        if "measurement_reference_csv is intentionally absent" not in line
        and "paper gate remains blocked until the measurement field" not in line
        and "expected to fail until the case-specific measurement reference field is extracted" not in line
    ]


def write_config_with_reference(path: Path, reference_rel: str) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    lines = remove_comment_placeholders(lines)
    lines = update_or_append_line(lines, "measurement_reference_csv=", f"measurement_reference_csv={reference_rel}")
    lines = update_or_append_line(lines, "paper_primary_measurement_frozen=", "paper_primary_measurement_frozen=true")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_reference_rows(path: Path) -> list[dict[str, float | str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    parsed: list[dict[str, float | str]] = []
    for row in rows:
        parsed.append(
            {
                "zenith_deg": float(row["zenith_deg"]),
                "relative_azimuth_deg": float(row["relative_azimuth_deg"]),
                "intensity": float(row["intensity"]),
                "q": float(row["q"]),
                "u": float(row["u"]),
                "dop": float(row["dop"]),
                "aop_deg": float(row["aop_deg"]),
            }
        )
    return parsed


def row_key(row: dict[str, float | str]) -> tuple[int, int]:
    return (
        int(round(float(row["zenith_deg"]) * 1000.0)),
        int(round(float(row["relative_azimuth_deg"]) * 1000.0)),
    )


def evenly_spaced_selection(rows: list[dict[str, float | str]], target_count: int) -> list[dict[str, float | str]]:
    if target_count <= 0 or not rows:
        return []
    ordered = sorted(rows, key=lambda row: (float(row["zenith_deg"]), float(row["relative_azimuth_deg"])))
    if len(ordered) <= target_count:
        return ordered
    selected: list[dict[str, float | str]] = []
    used_indices: set[int] = set()
    for sample_index in range(target_count):
        if target_count == 1:
            candidate_index = len(ordered) // 2
        else:
            candidate_index = int(round(sample_index * (len(ordered) - 1) / (target_count - 1)))
        while candidate_index in used_indices and candidate_index + 1 < len(ordered):
            candidate_index += 1
        while candidate_index in used_indices and candidate_index - 1 >= 0:
            candidate_index -= 1
        if candidate_index in used_indices:
            continue
        used_indices.add(candidate_index)
        selected.append(ordered[candidate_index])
    return selected


def write_reference_rows(
    path: Path,
    rows: list[dict[str, float | str]],
    *,
    include_subset_region: bool = False,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["zenith_deg", "relative_azimuth_deg", "intensity", "q", "u", "dop", "aop_deg"]
    if include_subset_region:
        fieldnames.append("subset_region")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "zenith_deg": f"{float(row['zenith_deg']):.6f}",
                    "relative_azimuth_deg": f"{float(row['relative_azimuth_deg']):.6f}",
                    "intensity": f"{float(row['intensity']):.6e}",
                    "q": f"{float(row['q']):.6e}",
                    "u": f"{float(row['u']):.6e}",
                    "dop": f"{float(row['dop']):.6f}",
                    "aop_deg": f"{float(row['aop_deg']):.6f}",
                    **({"subset_region": str(row.get("subset_region", ""))} if include_subset_region else {}),
                }
            )


def build_region_subset_rows(
    rows: list[dict[str, float | str]],
    *,
    region_specs,
    target_count: int,
) -> tuple[list[dict[str, float | str]], dict]:
    used_keys: set[tuple[int, int]] = set()
    selected_rows: list[dict[str, float | str]] = []
    region_counts: dict[str, int] = {}

    def append_region(region_name: str, region_rows: list[dict[str, float | str]]) -> None:
        for row in region_rows:
            key = row_key(row)
            if key in used_keys:
                continue
            used_keys.add(key)
            tagged = dict(row)
            tagged["subset_region"] = region_name
            selected_rows.append(tagged)
            region_counts[region_name] = region_counts.get(region_name, 0) + 1

    for region_name, target_region_count, predicate in region_specs:
        candidates = [row for row in rows if predicate(row)]
        append_region(region_name, evenly_spaced_selection(candidates, target_region_count))

    if len(selected_rows) < target_count:
        remaining = [row for row in rows if row_key(row) not in used_keys]
        fill_count = target_count - len(selected_rows)
        append_region("global_fill", evenly_spaced_selection(remaining, fill_count))

    selected_rows = sorted(
        selected_rows,
        key=lambda row: (float(row["zenith_deg"]), float(row["relative_azimuth_deg"])),
    )
    summary = {
        "target_count": target_count,
        "rows_written": len(selected_rows),
        "region_counts": region_counts,
        "region_targets": {name: target for name, target, _ in region_specs},
    }
    return selected_rows, summary


def build_strict_subset_rows(
    rows: list[dict[str, float | str]],
    *,
    target_count: int = STRICT_SUBSET_TARGET_COUNT,
) -> tuple[list[dict[str, float | str]], dict]:
    return build_region_subset_rows(
        rows,
        region_specs=STRICT_SUBSET_REGIONS,
        target_count=target_count,
    )


def build_profile_subset_rows(
    rows: list[dict[str, float | str]],
    *,
    target_count: int = PROFILE_SUBSET_TARGET_COUNT,
) -> tuple[list[dict[str, float | str]], dict]:
    return build_region_subset_rows(
        rows,
        region_specs=PROFILE_SUBSET_REGIONS,
        target_count=target_count,
    )


def write_strict_subset_config(case_id: str) -> Path:
    strict_cfg = PAPER_CONFIG_DIR / f"{case_id}_measurement.cfg"
    strict_subset_cfg = PAPER_CONFIG_DIR / f"{case_id}{STRICT_SUBSET_CONFIG_SUFFIX}"
    lines = strict_cfg.read_text(encoding="utf-8").splitlines()
    lines = update_or_append_line(lines, "case_id=", f"case_id={case_id}_measurement_strict_subset")
    lines = update_or_append_line(
        lines,
        "measurement_reference_csv=",
        f"measurement_reference_csv=../../data/paper_cases/{case_id}/{STRICT_SUBSET_FILENAME}",
    )
    lines = update_or_append_line(lines, "photons_per_bin=", "photons_per_bin=16")
    strict_subset_cfg.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return strict_subset_cfg


def write_profile_subset_config(case_id: str) -> Path:
    strict_cfg = PAPER_CONFIG_DIR / f"{case_id}_measurement.cfg"
    profile_subset_cfg = PAPER_CONFIG_DIR / f"{case_id}{PROFILE_SUBSET_CONFIG_SUFFIX}"
    lines = strict_cfg.read_text(encoding="utf-8").splitlines()
    lines = update_or_append_line(lines, "case_id=", f"case_id={case_id}_measurement_profile_subset")
    lines = update_or_append_line(
        lines,
        "measurement_reference_csv=",
        f"measurement_reference_csv=../../data/paper_cases/{case_id}/{PROFILE_SUBSET_FILENAME}",
    )
    profile_subset_cfg.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return profile_subset_cfg


def write_case_and_validation_frozen(case_id: str) -> None:
    case_config = PAPER_CONFIG_DIR / f"{case_id}.cfg"
    paper_validation = CONFIG_DIR / "paper_validation.cfg"
    for path in (case_config, paper_validation):
        lines = path.read_text(encoding="utf-8").splitlines()
        lines = remove_comment_placeholders(lines)
        lines = update_or_append_line(lines, "paper_primary_measurement_frozen=", "paper_primary_measurement_frozen=true")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def extract_timestamp_key(timestamp_utc: str) -> str:
    return timestamp_utc.replace("Z", "").replace(":", "-")


def reduce_case(
    case_id: str,
    color_channel: int,
    zenith_bins: int,
    azimuth_bins: int,
    min_pixels_per_bin: int,
    limit_zenith_deg: float,
) -> dict:
    case_data_dir, _ = get_case_paths(case_id)
    dataset_record = load_json(case_data_dir / "dataset_record.json")
    pointer = load_json(case_data_dir / "measurement_pointer.json")
    metadata_path = case_data_dir / "measurement_metadata.json"
    metadata = load_json(metadata_path)
    provenance_path = case_data_dir / "paper_case_provenance.json"
    provenance = load_json(provenance_path)

    annotation_filename = pointer["annotation_filename"]
    raw_filename = annotation_filename.replace("_annotations.npy", ".npy")
    raw_path = case_data_dir / raw_filename
    alpha_path = case_data_dir / "alpha_crop.npy"
    theta_path = case_data_dir / "theta_crop.npy"
    rotation_path = case_data_dir / "rotation.npy"

    file_specs = {
        raw_filename: raw_path,
        "alpha_crop.npy": alpha_path,
        "theta_crop.npy": theta_path,
        "rotation.npy": rotation_path,
    }

    with requests.Session() as session:
        for filename, destination in file_specs.items():
            data_file = find_dataset_file(dataset_record, filename)
            download_with_resume(session, int(data_file["id"]), destination, int(data_file["filesize"]))

    annotations = np.load(case_data_dir / annotation_filename, allow_pickle=True)
    timestamp_key = extract_timestamp_key(pointer["selected_timestamp_utc"])
    matches = np.where(annotations[:, 0] == timestamp_key)[0]
    if matches.size != 1:
        raise RuntimeError(f"Expected exactly one annotation match for {timestamp_key}, found {matches.size}")
    frame_index = int(matches[0])

    frames = np.load(raw_path, allow_pickle=True)
    if frame_index >= len(frames):
        raise RuntimeError(f"Frame index {frame_index} out of range for raw day file")
    frame = frames[frame_index]
    del frames

    image = np.asarray(frame[0], dtype=np.int64)
    exposure_us = float(frame[1])
    frame_timestamp = str(frame[2])
    image_type = int(frame[3])

    if frame_timestamp != timestamp_key:
        raise RuntimeError(
            f"Raw frame timestamp mismatch: expected {timestamp_key}, got {frame_timestamp}"
        )

    alpha = np.load(alpha_path, allow_pickle=True)
    theta = np.load(theta_path, allow_pickle=True)
    rotation = np.load(rotation_path, allow_pickle=True)
    if image.shape != alpha.shape or image.shape != theta.shape:
        raise RuntimeError(
            f"Raw image shape {image.shape} does not match calibration shapes {alpha.shape} and {theta.shape}"
        )
    rotation_z_deg = float(np.asarray(rotation).reshape(-1)[0])

    def subgrid(array: np.ndarray, polarization_index: int) -> np.ndarray:
        return array[
            ((color_channel // 2) * 2 + polarization_index // 2) :: 4,
            ((color_channel % 2) * 2 + polarization_index % 2) :: 4,
        ]

    i0 = subgrid(image, 0)
    i45 = subgrid(image, 1)
    i135 = subgrid(image, 2)
    i90 = subgrid(image, 3)

    q = i0 - i90
    u = i45 - i135
    intensity = 0.5 * (i0 + i45 + i135 + i90)

    alpha_0 = subgrid(alpha, 0)
    alpha_45 = subgrid(alpha, 1)
    alpha_135 = subgrid(alpha, 2)
    alpha_90 = subgrid(alpha, 3)
    theta_0 = subgrid(theta, 0)
    theta_45 = subgrid(theta, 1)
    theta_135 = subgrid(theta, 2)
    theta_90 = subgrid(theta, 3)

    alpha_mini = average_angle(average_angle(alpha_0, alpha_90), average_angle(alpha_45, alpha_135))
    theta_mini = average_angle(average_angle(theta_0, theta_90), average_angle(theta_45, theta_135))

    zenith_deg = np.degrees(theta_mini)
    absolute_azimuth_deg = normalize_azimuth_deg(np.degrees(alpha_mini))
    relative_azimuth_deg = normalize_azimuth_deg(absolute_azimuth_deg - float(metadata["solar_azimuth_deg"]))
    comparison_rotation_deg = marseille_comparison_basis_rotation_deg(relative_azimuth_deg, rotation_z_deg)
    q_compare, u_compare = rotate_linear_stokes(q, u, comparison_rotation_deg)

    saturation_mask = np.maximum.reduce([i0, i45, i135, i90]) >= SATURATION_LEVEL
    finite_mask = (
        np.isfinite(intensity)
        & np.isfinite(q_compare)
        & np.isfinite(u_compare)
        & np.isfinite(zenith_deg)
        & np.isfinite(relative_azimuth_deg)
    )
    valid_mask = (
        finite_mask
        & (~saturation_mask)
        & (intensity > 0.0)
        & (zenith_deg >= 0.0)
        & (zenith_deg <= limit_zenith_deg)
    )

    zenith_step = limit_zenith_deg / float(zenith_bins)
    azimuth_step = 360.0 / float(azimuth_bins)
    zenith_index = np.clip((zenith_deg / zenith_step).astype(int), 0, zenith_bins - 1)
    azimuth_index = np.clip((relative_azimuth_deg / azimuth_step).astype(int), 0, azimuth_bins - 1)

    accum: dict[tuple[int, int], dict[str, float]] = {}
    for z_idx, a_idx, zen, rel_az, abs_az, i_val, q_val, u_val in zip(
        zenith_index[valid_mask].ravel(),
        azimuth_index[valid_mask].ravel(),
        zenith_deg[valid_mask].ravel(),
        relative_azimuth_deg[valid_mask].ravel(),
        absolute_azimuth_deg[valid_mask].ravel(),
        intensity[valid_mask].ravel(),
        q_compare[valid_mask].ravel(),
        u_compare[valid_mask].ravel(),
    ):
        key = (int(z_idx), int(a_idx))
        bucket = accum.setdefault(
            key,
            {
                "count": 0.0,
                "sum_i": 0.0,
                "sum_q": 0.0,
                "sum_u": 0.0,
                "sum_zenith": 0.0,
                "sum_abs_cos": 0.0,
                "sum_abs_sin": 0.0,
                "sum_rel_cos": 0.0,
                "sum_rel_sin": 0.0,
            },
        )
        bucket["count"] += 1.0
        bucket["sum_i"] += float(i_val)
        bucket["sum_q"] += float(q_val)
        bucket["sum_u"] += float(u_val)
        bucket["sum_zenith"] += float(zen)
        bucket["sum_abs_cos"] += math.cos(math.radians(float(abs_az)))
        bucket["sum_abs_sin"] += math.sin(math.radians(float(abs_az)))
        bucket["sum_rel_cos"] += math.cos(math.radians(float(rel_az)))
        bucket["sum_rel_sin"] += math.sin(math.radians(float(rel_az)))

    reference_rows: list[dict[str, float | str]] = []
    for key in sorted(accum):
        bucket = accum[key]
        count = int(bucket["count"])
        if count < min_pixels_per_bin:
            continue
        mean_i = bucket["sum_i"] / bucket["count"]
        sum_i = bucket["sum_i"]
        if mean_i <= 0.0 or sum_i <= 0.0:
            continue
        mean_q = bucket["sum_q"] / bucket["count"]
        mean_u = bucket["sum_u"] / bucket["count"]
        dop = math.hypot(mean_q, mean_u) / max(mean_i, 1.0e-12)
        dop = min(max(dop, 0.0), 1.0)
        aop_deg = wrap_half_turn_deg(0.5 * math.degrees(math.atan2(mean_u, mean_q)))
        mean_zenith = bucket["sum_zenith"] / bucket["count"]
        mean_rel_azimuth = normalize_azimuth_deg(
            math.degrees(math.atan2(bucket["sum_rel_sin"], bucket["sum_rel_cos"]))
        )
        reference_rows.append(
            {
                "zenith_deg": mean_zenith,
                "relative_azimuth_deg": mean_rel_azimuth,
                "intensity": mean_i,
                "q": mean_q,
                "u": mean_u,
                "dop": dop,
                "aop_deg": aop_deg,
            }
        )

    reference_path = case_data_dir / "measurement_reference.csv"
    write_reference_rows(reference_path, reference_rows)
    rows_written = len(reference_rows)

    strict_subset_rows, strict_subset_summary = build_strict_subset_rows(reference_rows)
    strict_subset_path = case_data_dir / STRICT_SUBSET_FILENAME
    write_reference_rows(strict_subset_path, strict_subset_rows, include_subset_region=True)
    strict_subset_summary_path = case_data_dir / STRICT_SUBSET_SUMMARY_FILENAME
    write_json(
        strict_subset_summary_path,
        {
            "case_id": case_id,
            "source_reference_csv": str(reference_path),
            "subset_reference_csv": str(strict_subset_path),
            **strict_subset_summary,
        },
    )
    strict_subset_config_path = write_strict_subset_config(case_id)

    profile_subset_rows, profile_subset_summary = build_profile_subset_rows(reference_rows)
    profile_subset_path = case_data_dir / PROFILE_SUBSET_FILENAME
    write_reference_rows(profile_subset_path, profile_subset_rows, include_subset_region=True)
    profile_subset_summary_path = case_data_dir / PROFILE_SUBSET_SUMMARY_FILENAME
    write_json(
        profile_subset_summary_path,
        {
            "case_id": case_id,
            "source_reference_csv": str(reference_path),
            "subset_reference_csv": str(profile_subset_path),
            **profile_subset_summary,
        },
    )
    profile_subset_config_path = write_profile_subset_config(case_id)

    reduction_summary = {
        "case_id": case_id,
        "channel_name": metadata["channel_description"],
        "reduction_source_timestamp": pointer["selected_timestamp_utc"],
        "raw_timestamp_key": timestamp_key,
        "frame_index": frame_index,
        "raw_frame_shape": list(image.shape),
        "exposure_us": exposure_us,
        "image_type": image_type,
        "color_channel_index": color_channel,
        "rotation_z_deg": rotation_z_deg,
        "comparison_basis": "sensor_q_u_rotated_by(relative_azimuth_deg - rotation_z_deg - 90_deg)",
        "limit_zenith_deg": limit_zenith_deg,
        "zenith_bins": zenith_bins,
        "azimuth_bins": azimuth_bins,
        "min_pixels_per_bin": min_pixels_per_bin,
        "valid_pixels": int(np.count_nonzero(valid_mask)),
        "saturated_pixels": int(np.count_nonzero(saturation_mask)),
        "rows_written": rows_written,
        "raw_file_sha256": sha256_file(raw_path),
        "alpha_crop_sha256": sha256_file(alpha_path),
        "theta_crop_sha256": sha256_file(theta_path),
        "rotation_sha256": sha256_file(rotation_path),
        "measurement_reference_sha256": sha256_file(reference_path),
        "strict_subset_reference_sha256": sha256_file(strict_subset_path),
        "strict_subset_rows_written": len(strict_subset_rows),
        "profile_subset_reference_sha256": sha256_file(profile_subset_path),
        "profile_subset_rows_written": len(profile_subset_rows),
    }
    reduction_path = case_data_dir / "measurement_reduction.json"
    write_json(reduction_path, reduction_summary)

    metadata["paper_gate_role"] = "frozen_primary_twilight_case_with_extracted_reference"
    metadata["measurement_reference_csv"] = str(reference_path)
    metadata["aop_reference_convention"] = (
        "Marseille comparison basis with measured sensor-frame Q/U rotated by "
        "(relative_azimuth_deg - rotation_z_deg - 90_deg) before binning."
    )
    write_json(metadata_path, metadata)

    provenance["paper_primary_measurement_frozen"] = True
    provenance["paper_gate_blocked_reason"] = (
        "The frozen Marseille twilight measurement reference has been extracted, "
        "but the case still has to pass the paper-validation thresholds."
    )
    artifacts = provenance.setdefault("artifacts", {})
    artifacts["measurement_reference_csv"] = {
        "path": str(reference_path),
        "sha256": sha256_file(reference_path),
    }
    artifacts["measurement_reference_strict_subset_csv"] = {
        "path": str(strict_subset_path),
        "sha256": sha256_file(strict_subset_path),
    }
    artifacts["measurement_reference_strict_subset_summary_json"] = {
        "path": str(strict_subset_summary_path),
        "sha256": sha256_file(strict_subset_summary_path),
    }
    artifacts["measurement_reference_profile_subset_csv"] = {
        "path": str(profile_subset_path),
        "sha256": sha256_file(profile_subset_path),
    }
    artifacts["measurement_reference_profile_subset_summary_json"] = {
        "path": str(profile_subset_summary_path),
        "sha256": sha256_file(profile_subset_summary_path),
    }
    artifacts["measurement_reduction_json"] = {
        "path": str(reduction_path),
        "sha256": sha256_file(reduction_path),
    }
    artifacts["rotation_npy"] = {
        "path": str(rotation_path),
        "sha256": sha256_file(rotation_path),
    }
    artifacts["measurement_strict_subset_config"] = {
        "path": str(strict_subset_config_path),
        "sha256": sha256_file(strict_subset_config_path),
    }
    artifacts["measurement_profile_subset_config"] = {
        "path": str(profile_subset_config_path),
        "sha256": sha256_file(profile_subset_config_path),
    }
    provenance["measurement_reduction"] = {
        "frame_index": frame_index,
        "raw_timestamp_key": timestamp_key,
        "valid_pixels": reduction_summary["valid_pixels"],
        "rows_written": rows_written,
        "zenith_bins": zenith_bins,
        "azimuth_bins": azimuth_bins,
        "strict_subset_rows_written": len(strict_subset_rows),
        "strict_subset_region_counts": strict_subset_summary["region_counts"],
        "profile_subset_rows_written": len(profile_subset_rows),
        "profile_subset_region_counts": profile_subset_summary["region_counts"],
    }
    write_json(provenance_path, provenance)

    reference_rel_from_measurement_cfg = f"../../data/paper_cases/{case_id}/measurement_reference.csv"
    measurement_cfg = PAPER_CONFIG_DIR / f"{case_id}_measurement.cfg"
    write_config_with_reference(measurement_cfg, reference_rel_from_measurement_cfg)
    write_case_and_validation_frozen(case_id)

    return reduction_summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Reduce the frozen Marseille twilight case into validator CSV format.")
    parser.add_argument(
        "--case-id",
        default="frozen_marseille_twilight_20220815_191413z",
        help="Paper-case identifier under monte_carlo_cpp/data/paper_cases/",
    )
    parser.add_argument("--color-channel", type=int, default=3, help="Dataset color channel index. Default: 3 (B).")
    parser.add_argument("--zenith-bins", type=int, default=19, help="Number of zenith bins for reduction.")
    parser.add_argument("--azimuth-bins", type=int, default=36, help="Number of relative-azimuth bins for reduction.")
    parser.add_argument(
        "--min-pixels-per-bin",
        type=int,
        default=8,
        help="Minimum number of raw pixels required to emit a reduced reference bin.",
    )
    parser.add_argument(
        "--limit-zenith-deg",
        type=float,
        default=90.0,
        help="Maximum zenith angle retained from the fisheye image.",
    )
    args = parser.parse_args()

    summary = reduce_case(
        case_id=args.case_id,
        color_channel=args.color_channel,
        zenith_bins=args.zenith_bins,
        azimuth_bins=args.azimuth_bins,
        min_pixels_per_bin=args.min_pixels_per_bin,
        limit_zenith_deg=args.limit_zenith_deg,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
