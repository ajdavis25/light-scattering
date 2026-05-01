#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Set


REPO_ROOT = Path(__file__).resolve().parents[1]
PAPER_DIR = REPO_ROOT / "paper"
MANIFEST_PATH = PAPER_DIR / "MANIFEST.csv"
SUMMARY_PATH = REPO_ROOT / "notebooks" / "marseille_calibrated_validation_summary_2026-04-29.json"
MEASUREMENT_CONFIG = (
    REPO_ROOT
    / "monte_carlo_cpp"
    / "config"
    / "paper_cases"
    / "frozen_marseille_twilight_20220815_191413z_measurement.cfg"
)
REPORT_CASE_ID = "frozen_marseille_twilight_20220815_191413z_measurement__full_branchcap_robust_r2"

PUBLICATION_DOCS = [
    PAPER_DIR / "README.md",
    PAPER_DIR / "CLAIMS.md",
    PAPER_DIR / "provenance.md",
    PAPER_DIR / "text_snippets.md",
    REPO_ROOT / "notebooks" / "MARSEILLE_CALIBRATED_PAPER_PACKAGE_2026-04-29.md",
]

REQUIRED_PHRASES = [
    "calibrated pipeline validation",
    "frozen row-wise measurement-model calibration",
]

DANGEROUS_PHRASES = [
    "independent raw first-principles closure",
    "raw first-principles closure",
    "independent raw predictive validation",
    "raw predictive physics closure",
    "independently predicts the marseille",
]

NEGATION_MARKERS = [
    "not",
    "do not",
    "does not",
    "no ",
    "forbidden",
    "rather than",
    "must not",
    "should not",
    "not_claimed",
]


def repo_relative(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def resolve_repo_path(path_text: str) -> Path:
    path = Path(path_text)
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest() -> List[Dict[str, str]]:
    with MANIFEST_PATH.open(newline="") as handle:
        return list(csv.DictReader(handle))


def validate_json(path: Path) -> None:
    with path.open() as handle:
        json.load(handle)


def validate_csv(path: Path) -> None:
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
    if not header:
        raise ValueError(f"CSV has no header: {repo_relative(path)}")


def verify_manifest() -> None:
    rows = load_manifest()
    if not rows:
        raise ValueError("MANIFEST.csv is empty")
    seen: Set[str] = set()
    for row in rows:
        artifact_id = row["artifact_id"]
        if artifact_id in seen:
            raise ValueError(f"Duplicate artifact id in manifest: {artifact_id}")
        seen.add(artifact_id)

        source_path = resolve_repo_path(row["source_path"])
        if row["required"].lower() == "true" and not source_path.exists():
            raise FileNotFoundError(f"Required artifact missing: {row['source_path']}")
        if not source_path.exists():
            continue

        expected_sha = row["sha256"].strip()
        actual_sha = sha256_file(source_path)
        if expected_sha and actual_sha != expected_sha:
            raise ValueError(
                f"Checksum mismatch for {row['source_path']}: expected {expected_sha} got {actual_sha}"
            )

        suffix = source_path.suffix.lower()
        if suffix == ".json":
            validate_json(source_path)
        elif suffix == ".csv":
            validate_csv(source_path)


def verify_summary_scope() -> None:
    with SUMMARY_PATH.open() as handle:
        summary = json.load(handle)
    metrics = summary.get("metrics", {})
    if summary.get("claim_scope") != "calibrated_pipeline_validation":
        raise ValueError("Summary claim_scope must be calibrated_pipeline_validation")
    if summary.get("not_claimed_scope") != "independent_raw_first_principles_closure":
        raise ValueError("Summary not_claimed_scope must be independent_raw_first_principles_closure")
    if summary.get("gate_status") != "pass_after_calibration":
        raise ValueError("Summary gate_status must be pass_after_calibration")
    if metrics.get("measurement_model_calibration_applied") is not True:
        raise ValueError("Summary must record measurement_model_calibration_applied=true")


def sentences(text: str) -> List[str]:
    return [item.strip() for item in re.split(r"(?<=[.!?])\s+", text) if item.strip()]


def has_negation(sentence: str) -> bool:
    lower = sentence.lower()
    return any(marker in lower for marker in NEGATION_MARKERS)


def verify_claim_language() -> None:
    combined = "\n".join(path.read_text() for path in PUBLICATION_DOCS)
    lower_combined = combined.lower()
    for phrase in REQUIRED_PHRASES:
        if phrase not in lower_combined:
            raise ValueError(f"Required claim phrase missing from publication docs: {phrase}")

    for path in PUBLICATION_DOCS:
        text = path.read_text()
        for sentence in sentences(text):
            lower = sentence.lower()
            for phrase in DANGEROUS_PHRASES:
                if phrase in lower and not has_negation(lower):
                    raise ValueError(
                        f"Potentially unsafe unqualified claim in {repo_relative(path)}: {sentence}"
                    )


def compact_validation_table_markdown(summary: Dict[str, object]) -> str:
    metrics = summary["metrics"]
    rows = [
        ("Reference sky directions", metrics["reference_points"], ""),
        ("Measurement-model calibration applied", str(metrics["measurement_model_calibration_applied"]).lower(), ""),
        ("Normalized intensity RMSE", metrics["normalized_rmse"], f"<= {metrics['normalized_rmse_gate']}"),
        ("Brightest-location error", f"{metrics['brightest_location_deg']} deg", f"<= {metrics['brightest_location_deg_gate']} deg"),
        ("Median DoLP absolute error", metrics["median_dolp_abs"], f"<= {metrics['median_dolp_abs_gate']}"),
        ("p95 DoLP absolute error", metrics["p95_dolp_abs"], f"<= {metrics['p95_dolp_abs_gate']}"),
        ("Median AoP error", f"{metrics['median_aop_deg']} deg", f"<= {metrics['median_aop_deg_gate']} deg"),
        ("p95 AoP error", f"{metrics['p95_aop_deg']} deg", f"<= {metrics['p95_aop_deg_gate']} deg"),
        (
            "Solar-vertical signed DoLP bias",
            metrics["solar_vertical_signed_dolp_bias"],
            f"<= {metrics['solar_vertical_signed_dolp_bias_gate']}",
        ),
        ("Interpretation", "calibrated row-wise closure, not independent raw physics", "caveat required"),
    ]

    lines = [
        "# Marseille Calibrated Validation Table",
        "",
        "| Quantity | Calibrated value | Gate | Status |",
        "| --- | ---: | ---: | --- |",
    ]
    for quantity, value, gate in rows:
        lines.append(f"| {quantity} | {value} | {gate} | pass |")
    lines.append("")
    lines.append(
        "This table is generated from `notebooks/marseille_calibrated_validation_summary_2026-04-29.json`."
    )
    lines.append(
        "The interpretation row is a required claim-boundary caveat, not an additional numerical gate."
    )
    lines.append("")
    return "\n".join(lines)


def export_manifest_artifacts(skip_artifact_ids=None) -> None:
    skip_artifact_ids = set(skip_artifact_ids or [])
    for row in load_manifest():
        if row["artifact_id"] in skip_artifact_ids:
            continue
        export_path_text = row["export_path"].strip()
        if not export_path_text:
            continue
        source_path = resolve_repo_path(row["source_path"])
        export_path = resolve_repo_path(export_path_text)
        if not source_path.exists():
            raise FileNotFoundError(f"Cannot export missing artifact: {row['source_path']}")
        export_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, export_path)
        print(f"exported {repo_relative(export_path)}")


def export_compact_validation_table() -> None:
    with SUMMARY_PATH.open() as handle:
        summary = json.load(handle)
    output_path = PAPER_DIR / "tables" / "marseille_calibrated_validation_table.md"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(compact_validation_table_markdown(summary))
    print(f"generated {repo_relative(output_path)}")


def verify_exports() -> None:
    missing: List[str] = []
    for row in load_manifest():
        export_path_text = row["export_path"].strip()
        if not export_path_text:
            continue
        export_path = resolve_repo_path(export_path_text)
        if not export_path.exists():
            missing.append(export_path_text)
    table_path = PAPER_DIR / "tables" / "marseille_calibrated_validation_table.md"
    if not table_path.exists():
        missing.append(repo_relative(table_path))
    if missing:
        raise FileNotFoundError("Missing exported artifacts: " + ", ".join(missing))


def quicklook_artifact_ids() -> Set[str]:
    ids = set()
    for row in load_manifest():
        if row["source_path"].startswith("plots/current/measurement_cases/"):
            ids.add(row["artifact_id"])
    return ids


def copy_regenerated_plots(regenerated_case_dir: Path) -> None:
    for row in load_manifest():
        if row["artifact_id"] not in quicklook_artifact_ids():
            continue
        export_path = resolve_repo_path(row["export_path"])
        regenerated_path = regenerated_case_dir / Path(row["source_path"]).name
        if not regenerated_path.exists():
            raise FileNotFoundError(f"Regenerated plot artifact missing: {regenerated_path}")
        export_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(regenerated_path, export_path)
        print(f"exported regenerated {repo_relative(export_path)}")


def regenerate_plots_from_frozen_artifacts(plot_root: Path) -> Path:
    sys.path.insert(0, str(REPO_ROOT))
    from spherical.measurement_plots import save_measurement_case_plots

    summary = save_measurement_case_plots(
        MEASUREMENT_CONFIG,
        plot_root,
        report_case_id=REPORT_CASE_ID,
    )
    print(f"regenerated plots in {summary['plot_dir']}")
    return Path(summary["plot_dir"])


def find_newer_python() -> str:
    for executable in ("python3.11", "python3.10", "python3.9", "python3.8"):
        path = shutil.which(executable)
        if path:
            return path
    raise RuntimeError(
        "--regenerate-plots requires a Python interpreter new enough to import the spherical plotting package"
    )


def run_plot_regeneration(plot_root: Path) -> Path:
    if sys.version_info >= (3, 7):
        return regenerate_plots_from_frozen_artifacts(plot_root)

    python_executable = find_newer_python()
    try:
        subprocess.run(
            [
                python_executable,
                str(Path(__file__).resolve()),
                "--regenerate-plots-worker",
                str(plot_root),
            ],
            cwd=str(REPO_ROOT),
            check=True,
        )
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            "Plot regeneration failed under {}. Install numpy and matplotlib for that interpreter, "
            "or use --export to mirror the already-frozen plots.".format(python_executable)
        ) from error
    return plot_root / REPORT_CASE_ID


def run_verify() -> None:
    verify_manifest()
    verify_summary_scope()
    verify_claim_language()
    print("verification passed")


def run_export() -> None:
    run_verify()
    export_manifest_artifacts()
    export_compact_validation_table()
    verify_exports()
    print("export passed")


def run_regenerate_plots() -> None:
    run_verify()
    with tempfile.TemporaryDirectory(prefix="marseille_plots_") as tmp_dir:
        regenerated_case_dir = run_plot_regeneration(Path(tmp_dir))
        export_manifest_artifacts(skip_artifact_ids=quicklook_artifact_ids())
        copy_regenerated_plots(regenerated_case_dir)
        export_compact_validation_table()
        verify_exports()
    print("regenerate-plots export passed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify and export the Marseille calibrated pipeline validation paper package."
    )
    parser.add_argument("--verify-only", action="store_true", help="Verify frozen artifacts and claim discipline.")
    parser.add_argument("--export", action="store_true", help="Export figures and tables from frozen artifacts.")
    parser.add_argument(
        "--regenerate-plots",
        action="store_true",
        help="Regenerate Marseille plots from frozen report/comparison artifacts before export.",
    )
    parser.add_argument("--regenerate-plots-worker", help=argparse.SUPPRESS)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.regenerate_plots_worker:
        regenerate_plots_from_frozen_artifacts(Path(args.regenerate_plots_worker))
        return 0

    if not (args.verify_only or args.export or args.regenerate_plots):
        args.verify_only = True

    if args.regenerate_plots:
        run_regenerate_plots()
    elif args.export:
        run_export()
    elif args.verify_only:
        run_verify()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
