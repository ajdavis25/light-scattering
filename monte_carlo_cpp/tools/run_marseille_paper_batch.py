from __future__ import annotations

import argparse
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MEASUREMENT_CONFIG = (
    REPO_ROOT
    / "monte_carlo_cpp"
    / "config"
    / "paper_cases"
    / "frozen_marseille_twilight_20220815_191413z_measurement.cfg"
)
DEFAULT_VALIDATION_CONFIG = REPO_ROOT / "monte_carlo_cpp" / "config" / "paper_validation.cfg"
DEFAULT_LOG = REPO_ROOT / "monte_carlo_cpp" / "results" / "validation" / "marseille_paper_batch.log"
DEFAULT_BATCHED_MEASUREMENT_RUNNER = (
    REPO_ROOT / "monte_carlo_cpp" / "tools" / "run_measurement_case_batched.py"
)


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


def write_log_header(
    log_path: Path,
    measurement_runner: Path,
    measurement_config: Path,
    measurement_batch_size: int,
    higher_order_block_size: int,
    validation_runner: Path,
    validation_config: Path,
) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        f"started_utc={datetime.now(timezone.utc).isoformat()}",
        f"measurement_runner={measurement_runner}",
        f"measurement_config={measurement_config}",
        f"measurement_batch_size={measurement_batch_size}",
        f"measurement_higher_order_block_size={higher_order_block_size}",
        f"validation_runner={validation_runner}",
        f"validation_config={validation_config}",
        f"cwd={REPO_ROOT}",
        "",
    ]
    log_path.write_text("\n".join(header))


def stream_command(command: list[str], log_file, label: str) -> int:
    log_file.write(f"=== {label} started_utc={datetime.now(timezone.utc).isoformat()} ===\n")
    log_file.flush()
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
    log_file.write(f"=== {label} exit_code={return_code} finished_utc={datetime.now(timezone.utc).isoformat()} ===\n")
    log_file.flush()
    return return_code


def run_logged_command(command: list[str], log_path: Path, label: str) -> int:
    with log_path.open("a", encoding="utf-8") as log_file:
        log_file.write(f"=== {label} started_utc={datetime.now(timezone.utc).isoformat()} ===\n")
        log_file.flush()
        process = subprocess.Popen(
            command,
            cwd=REPO_ROOT,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            text=True,
        )
        return_code = process.wait()
        log_file.write(
            f"=== {label} exit_code={return_code} finished_utc={datetime.now(timezone.utc).isoformat()} ===\n"
        )
        log_file.flush()
        return return_code


def run_sequence(
    measurement_runner: Path,
    measurement_config: Path,
    validation_runner: Path,
    validation_config: Path,
    log_path: Path,
    measurement_batch_size: int,
    higher_order_block_size: int,
    *,
    stream_output: bool,
) -> int:
    measurement_command = [
        sys.executable,
        "-u",
        str(measurement_runner),
        str(measurement_config),
        "--batch-size",
        str(measurement_batch_size),
        "--higher-order-block-size",
        str(higher_order_block_size),
        "--resume",
    ]
    if stream_output:
        with log_path.open("a", encoding="utf-8") as log_file:
            measurement_code = stream_command(
                measurement_command,
                log_file,
                "measurement_case",
            )
            if measurement_code != 0:
                log_file.write("Aborting validation because the Marseille measurement rerun failed.\n")
                return measurement_code
            return stream_command(
                [str(validation_runner), str(validation_config)],
                log_file,
                "paper_validation",
            )

    measurement_code = run_logged_command(
        measurement_command,
        log_path,
        "measurement_case",
    )
    if measurement_code != 0:
        with log_path.open("a", encoding="utf-8") as log_file:
            log_file.write("Aborting validation because the Marseille measurement rerun failed.\n")
        return measurement_code
    return run_logged_command(
        [str(validation_runner), str(validation_config)],
        log_path,
        "paper_validation",
    )


def run_detached(
    measurement_config: Path,
    validation_config: Path,
    log_path: Path,
    measurement_batch_size: int,
    higher_order_block_size: int,
) -> int:
    creationflags = 0
    start_new_session = False
    if sys.platform == "win32":
        creationflags = subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.DETACHED_PROCESS
    else:
        start_new_session = True

    with log_path.open("a", encoding="utf-8") as log_file:
        process = subprocess.Popen(
            [
                sys.executable,
                "-u",
                str(Path(__file__).resolve()),
                "--measurement-config",
                str(measurement_config),
                "--validation-config",
                str(validation_config),
                "--measurement-batch-size",
                str(measurement_batch_size),
                "--higher-order-block-size",
                str(higher_order_block_size),
                "--log",
                str(log_path),
                "--run-now",
            ],
            cwd=REPO_ROOT,
            stdin=subprocess.DEVNULL,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            text=True,
            creationflags=creationflags,
            start_new_session=start_new_session,
        )

    print(f"Started detached Marseille + paper-validation batch with PID {process.pid}")
    print(f"Log: {log_path}")
    print("Inspect the log plus the measurement and validation report outputs for completion.")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the strict Marseille measurement rerun followed by the strict paper-validation suite."
    )
    parser.add_argument(
        "--measurement-config",
        type=Path,
        default=DEFAULT_MEASUREMENT_CONFIG,
        help="Measurement config to run first.",
    )
    parser.add_argument(
        "--validation-config",
        type=Path,
        default=DEFAULT_VALIDATION_CONFIG,
        help="Validation config to run after the Marseille measurement rerun.",
    )
    parser.add_argument(
        "--log",
        type=Path,
        default=DEFAULT_LOG,
        help="Combined batch log file.",
    )
    parser.add_argument(
        "--measurement-batch-size",
        type=int,
        default=8,
        help="Exact Marseille directions per batch, or max concurrent checkpointed directions.",
    )
    parser.add_argument(
        "--higher-order-block-size",
        type=int,
        default=32,
        help="Higher-order samples per checkpointed single-direction batch invocation.",
    )
    parser.add_argument(
        "--detached",
        action="store_true",
        help="Launch the batch in the background and return immediately.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the resolved runners, configs, and log path without launching the batch.",
    )
    parser.add_argument(
        "--run-now",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    measurement_runner = DEFAULT_BATCHED_MEASUREMENT_RUNNER.resolve()
    validation_runner = find_runner("ValidationRunner.exe")
    measurement_config = args.measurement_config.resolve()
    validation_config = args.validation_config.resolve()
    log_path = args.log.resolve()

    if args.dry_run:
        print(f"measurement_runner={measurement_runner}")
        print(f"measurement_config={measurement_config}")
        print(f"measurement_batch_size={args.measurement_batch_size}")
        print(f"measurement_higher_order_block_size={args.higher_order_block_size}")
        print(f"validation_runner={validation_runner}")
        print(f"validation_config={validation_config}")
        print(f"log={log_path}")
        return 0

    write_log_header(
        log_path,
        measurement_runner,
        measurement_config,
        args.measurement_batch_size,
        args.higher_order_block_size,
        validation_runner,
        validation_config,
    )

    if args.detached and not args.run_now:
        return run_detached(
            measurement_config,
            validation_config,
            log_path,
            args.measurement_batch_size,
            args.higher_order_block_size,
        )

    return run_sequence(
        measurement_runner,
        measurement_config,
        validation_runner,
        validation_config,
        log_path,
        args.measurement_batch_size,
        args.higher_order_block_size,
        stream_output=not args.run_now,
    )


if __name__ == "__main__":
    raise SystemExit(main())
