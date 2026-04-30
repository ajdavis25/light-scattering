from __future__ import annotations

import argparse
import csv
import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

from run_measurement_case_batched import (
    load_reference_rows,
    parse_config,
    resolve_config_path,
    rewrite_config,
    single_direction_case_id,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
PATH_SUFFIXES = ("_csv", "_json", "_cfg", "_config")


def split_rows(reference_rows: list[dict[str, str]], shard_count: int, layout: str) -> list[list[dict[str, str]]]:
    if shard_count <= 1:
        return [list(reference_rows)]

    shards: list[list[dict[str, str]]] = [[] for _ in range(shard_count)]
    if layout == "contiguous":
        chunk_size = math.ceil(len(reference_rows) / shard_count)
        for shard_index in range(shard_count):
            start = shard_index * chunk_size
            stop = min(len(reference_rows), start + chunk_size)
            shards[shard_index].extend(reference_rows[start:stop])
    else:
        for row_index, row in enumerate(reference_rows):
            shards[row_index % shard_count].append(row)

    return [rows for rows in shards if rows]


def make_path_overrides(config_path: Path, config_values: dict[str, str]) -> dict[str, str]:
    overrides: dict[str, str] = {}
    for key, value in config_values.items():
        if not value:
            continue
        if key.endswith(PATH_SUFFIXES):
            overrides[key] = str(resolve_config_path(config_path, value))
    return overrides


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create sharded configs for a measurement-case batch run.")
    parser.add_argument("config", type=Path, help="Base measurement config to shard.")
    parser.add_argument("--shards", type=int, default=8, help="Requested shard count.")
    parser.add_argument(
        "--layout",
        choices=("strided", "contiguous"),
        default="strided",
        help="Direction-to-shard assignment layout.",
    )
    parser.add_argument(
        "--case-id-suffix",
        default="arc_sharded",
        help="Suffix appended to the parent case_id to avoid clobbering existing outputs.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="Optional manifest output path. Defaults under results/measurement_case_reports/_sharded.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config_path = args.config.resolve()
    config_text = config_path.read_text(encoding="utf-8")
    config_values = parse_config(config_path)
    base_case_id = config_values["case_id"]
    parent_case_id = f"{base_case_id}__{args.case_id_suffix}" if args.case_id_suffix else base_case_id
    output_dir = resolve_config_path(config_path, config_values["output_dir"]).resolve()
    reference_csv = resolve_config_path(config_path, config_values["measurement_reference_csv"]).resolve()
    report_dir = output_dir / "measurement_case_reports"
    report_dir.mkdir(parents=True, exist_ok=True)

    fieldnames, reference_rows = load_reference_rows(reference_csv)
    shard_count = max(1, min(args.shards, len(reference_rows)))
    shard_rows = split_rows(reference_rows, shard_count, args.layout)

    shard_root = report_dir / "_sharded" / parent_case_id
    reference_root = shard_root / "reference_csv"
    config_root = shard_root / "configs"
    reference_root.mkdir(parents=True, exist_ok=True)
    config_root.mkdir(parents=True, exist_ok=True)

    if args.manifest is None:
        manifest_path = shard_root / "manifest.json"
    else:
        manifest_path = args.manifest.resolve()
        manifest_path.parent.mkdir(parents=True, exist_ok=True)

    path_overrides = make_path_overrides(config_path, config_values)
    manifest = {
        "prepared_utc": datetime.now(timezone.utc).isoformat(),
        "prepared_epoch_seconds": time.time(),
        "repo_root": str(REPO_ROOT),
        "base_config": str(config_path),
        "base_case_id": base_case_id,
        "parent_case_id": parent_case_id,
        "output_dir": str(output_dir),
        "report_dir": str(report_dir),
        "reference_csv": str(reference_csv),
        "layout": args.layout,
        "requested_shards": args.shards,
        "shard_count": len(shard_rows),
        "parent_partial_rows_csv": str(report_dir / f"{parent_case_id}_batched_partial_rows.csv"),
        "parent_progress_json": str(report_dir / f"{parent_case_id}_batched_progress.json"),
        "parent_comparison_csv": str(report_dir / f"{parent_case_id}_comparison.csv"),
        "parent_report_txt": str(report_dir / f"{parent_case_id}.txt"),
        "parent_region_summary_csv": str(report_dir / f"{parent_case_id}_region_summary.csv"),
        "shards": [],
    }

    for shard_index, rows in enumerate(shard_rows):
        shard_case_id = f"{parent_case_id}__shard_{shard_index:03d}"
        shard_reference_path = reference_root / f"{shard_case_id}.csv"
        shard_config_path = config_root / f"{shard_case_id}.cfg"
        with shard_reference_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow({key: row.get(key, "") for key in fieldnames})

        shard_overrides = dict(path_overrides)
        shard_overrides["case_id"] = shard_case_id
        shard_overrides["measurement_reference_csv"] = str(shard_reference_path.resolve())
        shard_overrides["output_dir"] = str(output_dir)
        shard_config_path.write_text(rewrite_config(config_text, shard_overrides), encoding="utf-8")

        original_indices = sorted(int(float(row["original_index"])) for row in rows)
        checkpoint_paths = [
            str((report_dir / "_batched_work" / f"{single_direction_case_id(shard_case_id, original_index)}_checkpoint.txt").resolve())
            for original_index in original_indices
        ]
        manifest["shards"].append({
            "index": shard_index,
            "case_id": shard_case_id,
            "config": str(shard_config_path),
            "reference_csv": str(shard_reference_path),
            "direction_count": len(rows),
            "min_original_index": original_indices[0],
            "max_original_index": original_indices[-1],
            "original_indices": original_indices,
            "partial_rows_csv": str(report_dir / f"{shard_case_id}_batched_partial_rows.csv"),
            "progress_json": str(report_dir / f"{shard_case_id}_batched_progress.json"),
            "comparison_csv": str(report_dir / f"{shard_case_id}_comparison.csv"),
            "checkpoint_paths": checkpoint_paths,
            "report_txt": str(report_dir / f"{shard_case_id}.txt"),
            "region_summary_csv": str(report_dir / f"{shard_case_id}_region_summary.csv"),
        })

    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"manifest={manifest_path}")
    print(f"parent_case_id={parent_case_id}")
    print(f"shard_count={len(shard_rows)}")
    for shard in manifest["shards"]:
        print(
            "shard="
            f"{shard['index']:03d}"
            f" directions={shard['direction_count']}"
            f" case_id={shard['case_id']}"
            f" original_index_range={shard['min_original_index']}..{shard['max_original_index']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
