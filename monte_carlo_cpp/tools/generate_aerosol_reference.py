from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ATMOSPHERE_PATH = ROOT / "monte_carlo_cpp" / "data" / "atmosphere" / "clear_sky_midlatitude.csv"
PHASE_PATH = ROOT / "monte_carlo_cpp" / "data" / "optics" / "aerosol_phase_matrix_reference.csv"
METADATA_PATH = ROOT / "monte_carlo_cpp" / "data" / "optics" / "aerosol_reference_metadata.json"


PROFILE_FIELDNAMES = [
    "altitude_m",
    "pressure_pa",
    "temperature_k",
    "molecular_number_density_m3",
    "ozone_number_density_m3",
    "aerosol_extinction_550_m_inv",
    "aerosol_single_scattering_albedo",
    "aerosol_asymmetry",
    "aerosol_scattering_angstrom_exponent",
    "aerosol_absorption_angstrom_exponent",
]


PHASE_FIELDNAMES = [
    "wavelength_nm",
    "angle_deg",
    "f11",
    "f12",
    "f22",
    "f33",
    "f34",
    "f44",
]


def read_base_profile() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with ATMOSPHERE_PATH.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                {
                    "altitude_m": float(row["altitude_m"]),
                    "pressure_pa": float(row["pressure_pa"]),
                    "temperature_k": float(row["temperature_k"]),
                    "molecular_number_density_m3": float(row["molecular_number_density_m3"]),
                    "ozone_number_density_m3": float(row["ozone_number_density_m3"]),
                }
            )
    return rows


def hg_phase(mu: float, g: float) -> float:
    denom = max(1.0e-12, 1.0 + g * g - 2.0 * g * mu)
    return (1.0 - g * g) / (4.0 * math.pi * denom ** 1.5)


def aerosol_profile_components(altitude_m: float) -> dict[str, float]:
    z = max(0.0, altitude_m)
    boundary_fine_scattering = 5.2e-5 * math.exp(-z / 1800.0)
    boundary_coarse_scattering = 5.2e-5 * math.exp(-z / 900.0)
    elevated_fine_scattering = 6.1e-6 * math.exp(-z / 6000.0)
    background_sulfate_scattering = 3.5e-7 * math.exp(-max(0.0, z - 18000.0) / 8000.0)

    black_carbon_absorption = 3.7e-6 * math.exp(-z / 1800.0)
    residual_absorption = 4.0e-7 * math.exp(-z / 6500.0)

    scattering_components = {
        "boundary_fine": boundary_fine_scattering,
        "boundary_coarse": boundary_coarse_scattering,
        "elevated_fine": elevated_fine_scattering,
        "background_sulfate": background_sulfate_scattering,
    }
    absorption_components = {
        "black_carbon": black_carbon_absorption,
        "residual_absorber": residual_absorption,
    }

    scattering = sum(scattering_components.values())
    absorption = sum(absorption_components.values())
    extinction = scattering + absorption

    def weighted_average(values: dict[str, float], weights: dict[str, float]) -> float:
        total = sum(weights.values())
        if total <= 0.0:
            return 0.0
        return sum(values[key] * weights[key] for key in weights) / total

    asymmetry = weighted_average(
        {
            "boundary_fine": 0.68,
            "boundary_coarse": 0.84,
            "elevated_fine": 0.63,
            "background_sulfate": 0.58,
        },
        scattering_components,
    )
    scattering_angstrom = weighted_average(
        {
            "boundary_fine": 1.75,
            "boundary_coarse": 0.20,
            "elevated_fine": 1.30,
            "background_sulfate": 2.05,
        },
        scattering_components,
    )
    absorption_angstrom = weighted_average(
        {
            "black_carbon": 1.00,
            "residual_absorber": 1.35,
        },
        absorption_components,
    )

    return {
        "aerosol_extinction_550_m_inv": extinction,
        "aerosol_single_scattering_albedo": scattering / max(1.0e-16, extinction),
        "aerosol_asymmetry": asymmetry,
        "aerosol_scattering_angstrom_exponent": scattering_angstrom,
        "aerosol_absorption_angstrom_exponent": absorption_angstrom,
    }


def write_profile() -> float:
    rows = read_base_profile()
    output_rows: list[dict[str, str]] = []
    aod_550 = 0.0
    previous_altitude = None
    previous_extinction = None

    for base in rows:
        aerosol = aerosol_profile_components(base["altitude_m"])
        combined = {**base, **aerosol}
        output_rows.append(
            {
                "altitude_m": f"{combined['altitude_m']:.1f}",
                "pressure_pa": f"{combined['pressure_pa']:.6e}",
                "temperature_k": f"{combined['temperature_k']:.3f}",
                "molecular_number_density_m3": f"{combined['molecular_number_density_m3']:.6e}",
                "ozone_number_density_m3": f"{combined['ozone_number_density_m3']:.6e}",
                "aerosol_extinction_550_m_inv": f"{combined['aerosol_extinction_550_m_inv']:.6e}",
                "aerosol_single_scattering_albedo": f"{combined['aerosol_single_scattering_albedo']:.6f}",
                "aerosol_asymmetry": f"{combined['aerosol_asymmetry']:.6f}",
                "aerosol_scattering_angstrom_exponent": f"{combined['aerosol_scattering_angstrom_exponent']:.6f}",
                "aerosol_absorption_angstrom_exponent": f"{combined['aerosol_absorption_angstrom_exponent']:.6f}",
            }
        )

        if previous_altitude is not None and previous_extinction is not None:
            dz = combined["altitude_m"] - previous_altitude
            aod_550 += 0.5 * (previous_extinction + combined["aerosol_extinction_550_m_inv"]) * dz
        previous_altitude = combined["altitude_m"]
        previous_extinction = combined["aerosol_extinction_550_m_inv"]

    with ATMOSPHERE_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=PROFILE_FIELDNAMES)
        writer.writeheader()
        writer.writerows(output_rows)

    return aod_550


def aerosol_phase_parameters(wavelength_nm: float) -> dict[str, float]:
    wavelength_t = (wavelength_nm - 350.0) / (800.0 - 350.0)
    return {
        "coarse_weight": 0.40 + 0.08 * wavelength_t,
        "backscatter_weight": 0.04,
        "g_coarse": 0.83 - 0.01 * wavelength_t,
        "g_fine": 0.72 - 0.10 * wavelength_t,
        "g_backscatter": -0.25,
        "p_max": 0.085 - 0.020 * wavelength_t,
    }


def aerosol_phase_row(wavelength_nm: float, angle_deg: float) -> dict[str, str]:
    params = aerosol_phase_parameters(wavelength_nm)
    mu = math.cos(math.radians(angle_deg))
    fine_weight = 1.0 - params["coarse_weight"] - params["backscatter_weight"]

    f11 = (
        params["coarse_weight"] * hg_phase(mu, params["g_coarse"])
        + fine_weight * hg_phase(mu, params["g_fine"])
        + params["backscatter_weight"] * hg_phase(mu, params["g_backscatter"])
    )

    forward_damp = 1.0 - math.exp(-((angle_deg / 18.0) ** 2))
    backward_damp = 1.0 - math.exp(-((((180.0 - angle_deg) / 24.0) ** 2)))
    polarization_shape = (1.0 - mu * mu) / (1.0 + 0.25 * mu * mu)
    polarization_fraction = params["p_max"] * polarization_shape * forward_damp * backward_damp
    polarization_fraction = min(max(polarization_fraction, 0.0), 0.20)

    f12 = -polarization_fraction * f11
    orthogonal_term = f11 * math.sqrt(max(0.0, 1.0 - polarization_fraction * polarization_fraction))

    return {
        "wavelength_nm": f"{wavelength_nm:.0f}",
        "angle_deg": f"{angle_deg:.0f}",
        "f11": f"{f11:.9e}",
        "f12": f"{f12:.9e}",
        "f22": f"{f11:.9e}",
        "f33": f"{orthogonal_term:.9e}",
        "f34": f"{0.0:.9e}",
        "f44": f"{orthogonal_term:.9e}",
    }


def write_phase_table() -> None:
    wavelengths = list(range(350, 801, 50))
    angles = list(range(0, 181, 5))

    with PHASE_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=PHASE_FIELDNAMES)
        writer.writeheader()
        for wavelength_nm in wavelengths:
            for angle_deg in angles:
                writer.writerow(aerosol_phase_row(float(wavelength_nm), float(angle_deg)))


def write_metadata(aod_550: float) -> None:
    metadata = {
        "description": (
            "Repository-owned clear-sky continental aerosol reference used for the bundled midlatitude profile "
            "and aerosol phase/Mueller lookup table."
        ),
        "intended_use": "Development-stage clear-sky reference input, not a study-specific paper dataset.",
        "column_aod_550": aod_550,
        "profile_model": {
            "scattering_components": {
                "boundary_fine": {
                    "extinction_550_ground_m_inv": 5.2e-5,
                    "scale_height_m": 1800.0,
                    "asymmetry": 0.68,
                    "angstrom_exponent": 1.75,
                },
                "boundary_coarse": {
                    "extinction_550_ground_m_inv": 5.2e-5,
                    "scale_height_m": 900.0,
                    "asymmetry": 0.84,
                    "angstrom_exponent": 0.20,
                },
                "elevated_fine": {
                    "extinction_550_ground_m_inv": 6.1e-6,
                    "scale_height_m": 6000.0,
                    "asymmetry": 0.63,
                    "angstrom_exponent": 1.30,
                },
                "background_sulfate": {
                    "extinction_550_ground_m_inv": 3.5e-7,
                    "onset_altitude_m": 18000.0,
                    "scale_height_m": 8000.0,
                    "asymmetry": 0.58,
                    "angstrom_exponent": 2.05,
                },
            },
            "absorption_components": {
                "black_carbon": {
                    "extinction_550_ground_m_inv": 3.7e-6,
                    "scale_height_m": 1800.0,
                    "angstrom_exponent": 1.00,
                },
                "residual_absorber": {
                    "extinction_550_ground_m_inv": 4.0e-7,
                    "scale_height_m": 6500.0,
                    "angstrom_exponent": 1.35,
                },
            },
        },
        "phase_matrix_model": {
            "f11_shape": "wavelength-dependent coarse/fine/backscatter Henyey-Greenstein mixture",
            "polarization_fraction": "bounded broadside envelope with suppressed forward and backward scattering",
            "mueller_constraints": "f22=f11, f33=f44=f11*sqrt(1-p^2), f34=0",
        },
    }
    METADATA_PATH.write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def main() -> None:
    aod_550 = write_profile()
    write_phase_table()
    write_metadata(aod_550)
    print(f"wrote {ATMOSPHERE_PATH}")
    print(f"wrote {PHASE_PATH}")
    print(f"wrote {METADATA_PATH}")
    print(f"column_aod_550={aod_550:.6f}")


if __name__ == "__main__":
    main()
