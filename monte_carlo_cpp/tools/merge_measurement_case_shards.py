from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from run_measurement_case_batched import (
    append_rows,
    finalize_outputs,
    load_reference_rows,
    parse_config,
    read_checkpoint_progress,
    read_comparison_rows,
    resolve_config_path,
    single_direction_case_id,
    write_progress,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge completed measurement shard outputs into one parent report.")
    parser.add_argument("manifest", type=Path, help="Shard manifest produced by prepare_measurement_case_shards.py.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    config_path = Path(manifest["base_config"]).resolve()
    config_values = parse_config(config_path)
    report_dir = Path(manifest["report_dir"]).resolve()
    parent_case_id = manifest["parent_case_id"]
    partial_rows_csv = Path(manifest["parent_partial_rows_csv"]).resolve()
    progress_path = Path(manifest["parent_progress_json"]).resolve()

    reference_csv = resolve_config_path(config_path, config_values["measurement_reference_csv"]).resolve()
    _, reference_rows = load_reference_rows(reference_csv)

    deduped: dict[int, dict[str, str]] = {}
    missing_shards: list[str] = []
    resumable_shards: list[str] = []
    resumable_checkpoint_progress: dict[str, list[tuple[int, int, int, int, int]]] = {}
    for shard in manifest["shards"]:
        shard_partial_rows = Path(shard["partial_rows_csv"]).resolve()
        shard_comparison_csv = Path(shard["comparison_csv"]).resolve()
        if shard_partial_rows.exists() and shard_partial_rows.stat().st_size > 0:
            rows = read_comparison_rows(shard_partial_rows)
        elif shard_comparison_csv.exists() and shard_comparison_csv.stat().st_size > 0:
            rows = read_comparison_rows(shard_comparison_csv)
        else:
            checkpoint_paths = shard.get("checkpoint_paths")
            if checkpoint_paths is None:
                checkpoint_paths = [
                    str(
                        (
                            report_dir
                            / "_batched_work"
                            / f"{single_direction_case_id(shard['case_id'], int(original_index))}_checkpoint.txt"
                        ).resolve()
                    )
                    for original_index in shard.get("original_indices", [])
                ]
            checkpoint_progresses = [
                read_checkpoint_progress(Path(checkpoint_path).resolve())
                for checkpoint_path in checkpoint_paths
            ]
            has_resumable_checkpoint = any(progress != (0, 0, 0, 0, 0) for progress in checkpoint_progresses)
            if has_resumable_checkpoint:
                resumable_shards.append(shard["case_id"])
                resumable_checkpoint_progress[shard["case_id"]] = checkpoint_progresses
            else:
                missing_shards.append(shard["case_id"])
            continue
        for row in rows:
            deduped[int(float(row["index"]))] = row

    if missing_shards:
        status = "incomplete"
    elif len(deduped) == len(reference_rows):
        status = "complete"
    elif resumable_shards:
        status = "in_progress"
    else:
        status = "incomplete"

    progress = {
        "case_id": parent_case_id,
        "config": str(config_path),
        "measurement_runner": "",
        "batch_size": 0,
        "mode": "sharded_merge",
        "resume": True,
        "completed_points": len(deduped),
        "total_points": len(reference_rows),
        "remaining_points": len(reference_rows) - len(deduped),
        "partial_rows_csv": str(partial_rows_csv),
        "log": "",
        "status": status,
        "resumable_shards": resumable_shards,
        "missing_shards": missing_shards,
    }

    if partial_rows_csv.exists():
        partial_rows_csv.unlink()
    append_rows(partial_rows_csv, [deduped[index] for index in sorted(deduped)])
    write_progress(progress_path, progress)

    if missing_shards:
        print("missing_shards=" + ",".join(missing_shards))
    if resumable_shards:
        print("resumable_shards=" + ",".join(resumable_shards))
        for shard_case_id in resumable_shards:
            summaries = [
                f"({progress[0]},{progress[1]},{progress[2]},{progress[3]},{progress[4]})"
                for progress in resumable_checkpoint_progress.get(shard_case_id, [])
            ]
            print(f"resumable_checkpoint_progress[{shard_case_id}]=" + ",".join(summaries))
    if len(deduped) != len(reference_rows) and missing_shards:
        print(f"incomplete_merge={len(deduped)}/{len(reference_rows)}")
        return 1
    if len(deduped) != len(reference_rows):
        print(f"resumable_merge={len(deduped)}/{len(reference_rows)}")
        return 0

    started = float(manifest.get("prepared_epoch_seconds", time.time()))
    finalize_outputs(
        [deduped[index] for index in sorted(deduped)],
        config_values,
        config_path,
        report_dir,
        parent_case_id,
        partial_rows_csv,
        started,
    )
    progress["completed_points"] = len(reference_rows)
    progress["remaining_points"] = 0
    progress["status"] = "complete"
    write_progress(progress_path, progress)
    print(f"merged_case_id={parent_case_id}")
    print(f"comparison_csv={manifest['parent_comparison_csv']}")
    print(f"report_txt={manifest['parent_report_txt']}")
    print(f"region_summary_csv={manifest['parent_region_summary_csv']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
