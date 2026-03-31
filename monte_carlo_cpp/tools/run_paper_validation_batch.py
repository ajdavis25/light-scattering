from __future__ import annotations

import argparse
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO_ROOT / "monte_carlo_cpp" / "config" / "paper_validation.cfg"
DEFAULT_LOG = REPO_ROOT / "monte_carlo_cpp" / "results" / "validation" / "paper_validation_batch.log"


def find_validation_runner() -> Path:
    candidates = [
        REPO_ROOT / "monte_carlo_cpp" / "build_current" / "ValidationRunner.exe",
        REPO_ROOT / "monte_carlo_cpp" / "build" / "ValidationRunner.exe",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Could not find ValidationRunner.exe in build_current or build.")


def write_log_header(log_path: Path, runner: Path, config: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        f"started_utc={datetime.now(timezone.utc).isoformat()}",
        f"runner={runner}",
        f"config={config}",
        f"cwd={REPO_ROOT}",
        "",
    ]
    log_path.write_text("\n".join(header))


def run_streaming(runner: Path, config: Path, log_path: Path) -> int:
    with log_path.open("a", encoding="utf-8") as log_file:
        process = subprocess.Popen(
            [str(runner), str(config)],
            cwd=REPO_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="")
            log_file.write(line)
        return process.wait()


def run_detached(runner: Path, config: Path, log_path: Path) -> int:
    creationflags = 0
    if sys.platform == "win32":
        creationflags = subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.DETACHED_PROCESS

    with log_path.open("a", encoding="utf-8") as log_file:
        process = subprocess.Popen(
            [str(runner), str(config)],
            cwd=REPO_ROOT,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            text=True,
            creationflags=creationflags,
        )

    print(f"Started detached validation run with PID {process.pid}")
    print(f"Log: {log_path}")
    print("Use the log file and monte_carlo_cpp/results/validation/validation_report.json to inspect completion.")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the strict paper validation suite as a logged batch job."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Validation config to run. Defaults to the strict paper-validation config.",
    )
    parser.add_argument(
        "--log",
        type=Path,
        default=DEFAULT_LOG,
        help="Path to the batch log file.",
    )
    parser.add_argument(
        "--detached",
        action="store_true",
        help="Launch ValidationRunner in the background and return immediately.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the resolved runner/config/log paths without launching validation.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    runner = find_validation_runner()
    config = args.config.resolve()
    log_path = args.log.resolve()

    if args.dry_run:
        print(f"runner={runner}")
        print(f"config={config}")
        print(f"log={log_path}")
        return 0

    write_log_header(log_path, runner, config)
    if args.detached:
        return run_detached(runner, config, log_path)
    return run_streaming(runner, config, log_path)


if __name__ == "__main__":
    raise SystemExit(main())
