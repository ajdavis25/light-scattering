from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT / "monte_carlo_cpp" / "src"
TEST_DIR = ROOT / "monte_carlo_cpp" / "tests"
BUILD_DIR = ROOT / "monte_carlo_cpp" / "build_current"
OBJ_DIR = BUILD_DIR / "obj"
MINGW_ROOT = Path("C:/mingw64")
GPP = MINGW_ROOT / "bin" / "g++.exe"
AS = MINGW_ROOT / "bin" / "as.exe"
LD = MINGW_ROOT / "bin" / "ld.exe"
GCC_LIB_DIR = MINGW_ROOT / "lib" / "gcc" / "x86_64-w64-mingw32" / "14.2.0"
MINGW_LIB_DIR = MINGW_ROOT / "x86_64-w64-mingw32" / "lib"
SYS_LIB_DIR = MINGW_ROOT / "lib"
RUNTIME_DLLS = [
    MINGW_ROOT / "bin" / "libstdc++-6.dll",
    MINGW_ROOT / "bin" / "libgcc_s_seh-1.dll",
    MINGW_ROOT / "bin" / "libgomp-1.dll",
    MINGW_ROOT / "bin" / "libwinpthread-1.dll",
]


CORE_SOURCES = [
    SRC_DIR / "Atmosphere.cpp",
    SRC_DIR / "WavelengthHandling.cpp",
    SRC_DIR / "PhaseFunctions.cpp",
    SRC_DIR / "Polarization.cpp",
    SRC_DIR / "SurfaceReflection.cpp",
    SRC_DIR / "MonteCarloDriver.cpp",
    SRC_DIR / "ValidationQA.cpp",
]


EXECUTABLE_SOURCES = {
    "MonteCarloCPP.exe": [SRC_DIR / "main.cpp"],
    "ValidationRunner.exe": [SRC_DIR / "validation_main.cpp"],
    "MeasurementCaseRunner.exe": [SRC_DIR / "measurement_case_main.cpp"],
    "test_phasefunctions.exe": [TEST_DIR / "test_phasefunctions.cpp"],
    "test_surface.exe": [TEST_DIR / "test_surface.cpp"],
    "test_polarization.exe": [TEST_DIR / "test_polarization.cpp"],
    "test_atmosphere.exe": [TEST_DIR / "test_atmosphere.cpp"],
    "test_wavelength.exe": [TEST_DIR / "test_wavelength.cpp"],
}


def handle_remove_readonly(function, path: str, excinfo) -> None:
    try:
        os.chmod(path, 0o666)
        function(path)
    except OSError:
        # Leave locked files in place; the incremental rebuild path below
        # only requires object files and executables to be refreshed.
        pass


def run_command(command: list[str], *, allow_failure: bool = False) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode != 0 and not allow_failure:
        raise RuntimeError(
            f"Command failed: {' '.join(command)}\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    return completed


def compile_source(source_path: Path) -> Path:
    object_path = OBJ_DIR / f"{source_path.stem}.obj"
    for suffix in (".obj", ".ii", ".s"):
        stale = object_path.with_suffix(suffix)
        if stale.exists():
            stale.unlink()

    compile_command = [
        str(GPP),
        "-save-temps=obj",
        f"-I{SRC_DIR.as_posix()}",
        "-std=gnu++17",
        "-fopenmp",
        "-c",
        str(source_path),
        "-o",
        str(object_path),
    ]
    compile_result = run_command(compile_command, allow_failure=True)
    assembly_path = object_path.with_suffix(".s")
    if not assembly_path.exists():
        raise RuntimeError(
            f"Compiler did not emit assembly for {source_path}\n"
            f"stdout:\n{compile_result.stdout}\n"
            f"stderr:\n{compile_result.stderr}"
        )

    assemble_command = [str(AS), "-o", str(object_path), str(assembly_path)]
    run_command(assemble_command)
    return object_path


def link_executable(executable_name: str, object_paths: list[Path]) -> None:
    output_path = BUILD_DIR / executable_name
    if output_path.exists():
        output_path.unlink()

    link_command = [
        str(LD),
        "-m",
        "i386pep",
        "-Bdynamic",
        "-o",
        str(output_path),
        str(MINGW_LIB_DIR / "crt2.o"),
        str(GCC_LIB_DIR / "crtbegin.o"),
        f"-L{GCC_LIB_DIR.as_posix()}",
        f"-L{MINGW_LIB_DIR.as_posix()}",
        f"-L{SYS_LIB_DIR.as_posix()}",
        *[str(path) for path in object_paths],
        "--start-group",
        "-lstdc++",
        "-lgomp",
        "-lmingwthrd",
        "-lmingw32",
        "-lgcc_s",
        "-lgcc",
        "-lmingwex",
        "-lmsvcrt",
        "-lm",
        "-lkernel32",
        "-lpthread",
        "-ladvapi32",
        "-lshell32",
        "-luser32",
        "-lkernel32",
        "-liconv",
        "--end-group",
        str(GCC_LIB_DIR / "crtend.o"),
    ]
    run_command(link_command)


def prepare_build_dir() -> None:
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    if OBJ_DIR.exists():
        shutil.rmtree(OBJ_DIR, onexc=handle_remove_readonly)
    OBJ_DIR.mkdir(parents=True, exist_ok=True)

    for executable_name in EXECUTABLE_SOURCES:
        executable_path = BUILD_DIR / executable_name
        if executable_path.exists():
            try:
                executable_path.unlink()
            except OSError:
                pass


def stage_runtime_dlls() -> None:
    for dll_path in RUNTIME_DLLS:
        target = BUILD_DIR / dll_path.name
        try:
            shutil.copy2(dll_path, target)
        except PermissionError:
            # When Windows keeps a runtime DLL mapped, reuse the existing copy.
            if not target.exists():
                raise


def main() -> None:
    prepare_build_dir()

    core_objects = [compile_source(source) for source in CORE_SOURCES]
    for executable_name, sources in EXECUTABLE_SOURCES.items():
        leaf_objects = [compile_source(source) for source in sources]
        link_executable(executable_name, [*core_objects, *leaf_objects])
        print(f"built {BUILD_DIR / executable_name}")

    stage_runtime_dlls()

    print("build_current is runnable")


if __name__ == "__main__":
    main()
