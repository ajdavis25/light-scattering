from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OPTICS_DIR = ROOT / "monte_carlo_cpp" / "data" / "optics"


def write_table(filename: str, header: tuple[str, str], rows: list[tuple[float, float]]) -> None:
    path = OPTICS_DIR / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(header)
        for wavelength_nm, value in rows:
            writer.writerow([f"{wavelength_nm:.1f}", f"{value:.6e}"])


def o2_rows() -> list[tuple[float, float]]:
    return [
        (350.0, 0.0),
        (500.0, 0.0),
        (680.0, 5.0e-31),
        (720.0, 5.0e-30),
        (740.0, 2.0e-29),
        (760.0, 1.2e-28),
        (770.0, 3.0e-29),
        (780.0, 5.0e-30),
        (800.0, 1.0e-30),
    ]


def o4_rows() -> list[tuple[float, float]]:
    # Reference O4 values are commonly tabulated in cm^5 molecule^-2.
    # The runtime expects m^5, so convert by 1e-10.
    return [
        (350.0, 4.0e-57),
        (400.0, 8.0e-57),
        (450.0, 1.4e-56),
        (477.0, 1.7e-56),
        (530.0, 7.0e-57),
        (577.0, 8.0e-57),
        (630.0, 9.0e-57),
        (700.0, 3.0e-57),
        (800.0, 1.0e-57),
    ]


def h2o_rows() -> list[tuple[float, float]]:
    return [
        (350.0, 0.0),
        (600.0, 0.0),
        (680.0, 3.0e-30),
        (700.0, 1.0e-28),
        (720.0, 4.0e-28),
        (740.0, 1.5e-28),
        (760.0, 8.0e-29),
        (780.0, 2.0e-28),
        (800.0, 4.0e-28),
    ]


def no2_rows() -> list[tuple[float, float]]:
    return [
        (350.0, 7.0e-27),
        (400.0, 5.0e-27),
        (450.0, 3.0e-27),
        (500.0, 1.4e-27),
        (550.0, 7.0e-28),
        (600.0, 2.5e-28),
        (650.0, 8.0e-29),
        (700.0, 1.5e-29),
        (800.0, 0.0),
    ]


def main() -> None:
    write_table("o2_cross_section_reference.csv", ("wavelength_nm", "o2_cross_section_m2"), o2_rows())
    write_table("o4_cross_section_reference.csv", ("wavelength_nm", "o4_cross_section_m5"), o4_rows())
    write_table("h2o_cross_section_reference.csv", ("wavelength_nm", "h2o_cross_section_m2"), h2o_rows())
    write_table("no2_cross_section_reference.csv", ("wavelength_nm", "no2_cross_section_m2"), no2_rows())
    print("gas cross-section reference tables written")


if __name__ == "__main__":
    main()
