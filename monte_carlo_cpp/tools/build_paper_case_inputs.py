from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import requests


ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = ROOT / "monte_carlo_cpp" / "config"
PAPER_CONFIG_DIR = CONFIG_DIR / "paper_cases"
DATA_DIR = ROOT / "monte_carlo_cpp" / "data"
PAPER_DATA_DIR = DATA_DIR / "paper_cases"
AVOGADRO = 6.02214076e23
BOLTZMANN = 1.380649e-23
DOBSON_TO_MOLECULES_PER_M2 = 2.687e20
WATER_MOLAR_MASS_KG = 0.01801528
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
PRESSURE_LEVELS_HPA = [
    1000, 975, 950, 925, 900, 850, 800, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30
]

FROZEN_CASES = {
    "frozen_marseille_twilight_20220815_191413z": {
        "case_id": "frozen_marseille_twilight_20220815_191413z",
        "dataset_doi": "10.57745/9L2YUB",
        "dataset_api_url": "https://entrepot.recherche.data.gouv.fr/api/datasets/:persistentId/?persistentId=doi:10.57745/9L2YUB",
        "dataset_record_url": "https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.57745/9L2YUB",
        "dataset_readme_name": "README.md",
        "measurement_paper_url": "https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-024-06959-6",
        "measurement_dataset_title": "A 2 month-long annotated skylight polarization images database",
        "site_lat_deg": 43.286990365824785,
        "site_lon_deg": 5.403361407820939,
        "timestamp_utc": "2022-08-15T19:14:13Z",
        "annotation_filename": "2022-08-15_raw_annotations.npy",
        "annotation_label": "c",
        "channel_name": "B",
        "channel_center_nm": 460.0,
        "channel_sigma_nm": 22.0,
        "observer_altitude_m": 35.0,
        "sensor_flyer_url": "https://www.sony-semicon.com/files/62/flyer_industry/IMX250_264_253MZR_MYR_Flyer_en.pdf",
        "camera_techref_url": "https://thinklucid.com/document/phoenix-5-0-mp-polarized-tech-ref-phx050s-p-q-imx250mzr-myr/",
        "camera_techref_download_url": "https://dce9ugryut4ao.cloudfront.net/PHX050S-Polarized_1.72.0.0_Documentation_English.zip",
    }
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(65536)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def trapz(x_values: list[float], y_values: list[float]) -> float:
    total = 0.0
    for index in range(len(x_values) - 1):
        dx = x_values[index + 1] - x_values[index]
        total += 0.5 * dx * (y_values[index] + y_values[index + 1])
    return total


def clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def solar_geometry_utc(timestamp: datetime, latitude_deg: float, longitude_deg: float) -> tuple[float, float]:
    day_of_year = int(timestamp.strftime("%j"))
    utc_hour = timestamp.hour + timestamp.minute / 60.0 + timestamp.second / 3600.0
    gamma = 2.0 * math.pi / 365.0 * (day_of_year - 1 + (utc_hour - 12.0) / 24.0)
    eqtime = 229.18 * (
        0.000075
        + 0.001868 * math.cos(gamma)
        - 0.032077 * math.sin(gamma)
        - 0.014615 * math.cos(2.0 * gamma)
        - 0.040849 * math.sin(2.0 * gamma)
    )
    decl = (
        0.006918
        - 0.399912 * math.cos(gamma)
        + 0.070257 * math.sin(gamma)
        - 0.006758 * math.cos(2.0 * gamma)
        + 0.000907 * math.sin(2.0 * gamma)
        - 0.002697 * math.cos(3.0 * gamma)
        + 0.00148 * math.sin(3.0 * gamma)
    )
    tst = (timestamp.hour * 60.0 + timestamp.minute + timestamp.second / 60.0 + eqtime + 4.0 * longitude_deg) % 1440.0
    hour_angle_deg = tst / 4.0 - 180.0
    if hour_angle_deg < -180.0:
        hour_angle_deg += 360.0
    hour_angle = math.radians(hour_angle_deg)
    latitude = math.radians(latitude_deg)
    cos_zenith = math.sin(latitude) * math.sin(decl) + math.cos(latitude) * math.cos(decl) * math.cos(hour_angle)
    cos_zenith = max(-1.0, min(1.0, cos_zenith))
    zenith_deg = math.degrees(math.acos(cos_zenith))
    azimuth_rad = math.atan2(
        math.sin(hour_angle),
        math.cos(hour_angle) * math.sin(latitude) - math.tan(decl) * math.cos(latitude),
    )
    return zenith_deg, (math.degrees(azimuth_rad) + 180.0) % 360.0


def gaussian_rows(center_nm: float, sigma_nm: float) -> list[tuple[float, float]]:
    rows = []
    for wavelength in range(350, 801, 10):
        rows.append((float(wavelength), math.exp(-0.5 * ((wavelength - center_nm) / sigma_nm) ** 2)))
    return rows


def hg_phase(mu: float, g: float) -> float:
    denom = max(1.0e-12, 1.0 + g * g - 2.0 * g * mu)
    return (1.0 - g * g) / (4.0 * math.pi * denom ** 1.5)


def normalize_rows(rows: list[tuple[float, float]]) -> list[tuple[float, float]]:
    peak = max(value for _, value in rows)
    if peak <= 0.0:
        return [(wavelength, 0.0) for wavelength, _ in rows]
    return [(wavelength, value / peak) for wavelength, value in rows]


def imx250myr_blue_rows_public_docs_proxy(center_nm: float, sigma_nm: float) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    for wavelength in range(350, 801, 10):
        wavelength_nm = float(wavelength)
        blue_sigma = sigma_nm * (1.18 if wavelength_nm <= center_nm else 1.55)
        primary = math.exp(-0.5 * ((wavelength_nm - center_nm) / blue_sigma) ** 2)
        green_shoulder = 0.14 * math.exp(-0.5 * ((wavelength_nm - (center_nm + 48.0)) / 26.0) ** 2)
        uv_roll_on = 1.0 / (1.0 + math.exp(-(wavelength_nm - (center_nm - 74.0)) / 7.5))
        ir_roll_off = 1.0 / (1.0 + math.exp((wavelength_nm - (center_nm + 102.0)) / 10.5))
        response = (primary + green_shoulder) * uv_roll_on * ir_roll_off
        rows.append((wavelength_nm, response))
    return normalize_rows(rows)


def blend_aod_550(openmeteo_aod_550: float, aeronet_aod_550: float, time_delta_hours: float) -> tuple[float, float]:
    if openmeteo_aod_550 <= 0.0:
        return aeronet_aod_550, 1.0
    if time_delta_hours <= 2.25:
        aeronet_weight = 0.75
    elif time_delta_hours <= 3.50:
        aeronet_weight = clamp(0.75 - 0.32 * (time_delta_hours - 2.25), 0.35, 0.75)
    else:
        aeronet_weight = 0.0
    blended = aeronet_weight * aeronet_aod_550 + (1.0 - aeronet_weight) * openmeteo_aod_550
    return blended, aeronet_weight


def build_case_aerosol_summary(
    *,
    surface_relative_humidity: float,
    pm10_ug_m3: float,
    pm25_ug_m3: float,
    target_aod_550: float,
    angstrom: float,
) -> dict:
    pm_ratio = clamp(pm25_ug_m3 / max(pm10_ug_m3, 1.0e-6), 0.10, 0.98)
    angstrom_indicator = clamp((angstrom - 0.70) / 1.30, 0.0, 1.0)
    mixed_fine_indicator = clamp(0.55 * pm_ratio + 0.45 * angstrom_indicator, 0.0, 1.0)

    coarse_fraction = clamp(
        0.20 + 0.28 * (1.0 - pm_ratio) + 0.18 * max(0.0, surface_relative_humidity - 0.55) + 0.08 * (1.0 - mixed_fine_indicator),
        0.22,
        0.48,
    )
    elevated_fraction = clamp(
        0.10 + 0.07 * mixed_fine_indicator + 0.05 * clamp(target_aod_550 / 0.20, 0.0, 1.0),
        0.12,
        0.24,
    )
    residual_fraction = 0.06
    boundary_fine_fraction = max(0.20, 1.0 - coarse_fraction - elevated_fraction - residual_fraction)
    total_fraction = coarse_fraction + elevated_fraction + residual_fraction + boundary_fine_fraction

    fractions = {
        "boundary_fine": boundary_fine_fraction / total_fraction,
        "coastal_coarse": coarse_fraction / total_fraction,
        "elevated_fine": elevated_fraction / total_fraction,
        "residual_sulfate": residual_fraction / total_fraction,
    }

    ssa_550 = clamp(
        0.945 + 0.020 * pm_ratio + 0.015 * max(0.0, surface_relative_humidity - 0.50) - 0.010 * max(0.0, target_aod_550 - 0.18),
        0.93,
        0.98,
    )
    absorption_black_carbon_fraction = clamp(
        0.72 - 0.15 * max(0.0, surface_relative_humidity - 0.60) + 0.10 * (1.0 - pm_ratio),
        0.55,
        0.78,
    )

    coarse_polarization_fraction = clamp(
        0.070 - 0.030 * fractions["coastal_coarse"] - 0.020 * max(0.0, surface_relative_humidity - 0.60) - 0.010 * pm_ratio,
        0.030,
        0.060,
    )

    return {
        "pm_ratio": pm_ratio,
        "mixed_fine_indicator": mixed_fine_indicator,
        "scattering_column_fraction": fractions,
        "ssa_550": ssa_550,
        "absorption_black_carbon_fraction": absorption_black_carbon_fraction,
        "pmax_base": coarse_polarization_fraction,
        "coastal_coarse_asymmetry": clamp(0.84 + 0.03 * max(0.0, surface_relative_humidity - 0.65), 0.84, 0.89),
        "fine_asymmetry": clamp(0.69 + 0.03 * mixed_fine_indicator, 0.68, 0.75),
        "angstrom_fine": clamp(1.45 + 0.25 * mixed_fine_indicator, 1.35, 1.80),
        "angstrom_coarse": clamp(0.35 + 0.20 * (1.0 - pm_ratio), 0.30, 0.60),
    }


def aerosol_component_shapes(altitudes_m: list[float]) -> dict[str, list[float]]:
    return {
        "boundary_fine": [math.exp(-altitude / 1400.0) for altitude in altitudes_m],
        "coastal_coarse": [math.exp(-altitude / 750.0) for altitude in altitudes_m],
        "elevated_fine": [math.exp(-max(0.0, altitude - 900.0) / 4200.0) for altitude in altitudes_m],
        "residual_sulfate": [math.exp(-max(0.0, altitude - 12000.0) / 7000.0) for altitude in altitudes_m],
        "black_carbon": [math.exp(-altitude / 1300.0) for altitude in altitudes_m],
        "brown_carbon": [math.exp(-altitude / 2500.0) for altitude in altitudes_m],
    }


def normalize_shape(altitudes_m: list[float], values: list[float]) -> list[float]:
    integral = max(trapz(altitudes_m, values), 1.0e-12)
    return [value / integral for value in values]


def case_aerosol_phase_parameters(summary: dict, wavelength_nm: float) -> dict:
    wavelength_t = clamp((wavelength_nm - 350.0) / (800.0 - 350.0), 0.0, 1.0)
    pm_ratio = summary["pm_ratio"]
    coarse_fraction = summary["scattering_column_fraction"]["coastal_coarse"]
    elevated_fraction = summary["scattering_column_fraction"]["elevated_fine"]

    coarse_weight = clamp(
        0.36 + 0.18 * coarse_fraction + 0.06 * (1.0 - pm_ratio) + 0.03 * wavelength_t,
        0.35,
        0.60,
    )
    backscatter_weight = 0.03 + 0.01 * elevated_fraction
    fine_weight = max(0.20, 1.0 - coarse_weight - backscatter_weight)

    p_max = clamp(summary["pmax_base"] - 0.010 * wavelength_t, 0.025, 0.055)
    return {
        "coarse_weight": coarse_weight,
        "fine_weight": fine_weight,
        "backscatter_weight": backscatter_weight,
        "g_coarse": clamp(summary["coastal_coarse_asymmetry"] - 0.01 * wavelength_t, 0.82, 0.88),
        "g_fine": clamp(summary["fine_asymmetry"] - 0.04 * wavelength_t, 0.62, 0.76),
        "g_backscatter": -0.22,
        "p_max": p_max,
    }


def aerosol_phase_row(summary: dict, wavelength_nm: float, angle_deg: float) -> dict[str, str]:
    params = case_aerosol_phase_parameters(summary, wavelength_nm)
    mu = math.cos(math.radians(angle_deg))

    f11 = (
        params["coarse_weight"] * hg_phase(mu, params["g_coarse"])
        + params["fine_weight"] * hg_phase(mu, params["g_fine"])
        + params["backscatter_weight"] * hg_phase(mu, params["g_backscatter"])
    )

    forward_damp = 1.0 - math.exp(-((angle_deg / 20.0) ** 2))
    backward_damp = 1.0 - math.exp(-((((180.0 - angle_deg) / 26.0) ** 2)))
    polarization_shape = (1.0 - mu * mu) / (1.0 + 0.30 * mu * mu)
    polarization_fraction = clamp(params["p_max"] * polarization_shape * forward_damp * backward_damp, 0.0, 0.18)
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


def fetch_json(session: requests.Session, url: str) -> dict:
    response = session.get(url, timeout=180)
    response.raise_for_status()
    return response.json()


def fetch_text(session: requests.Session, url: str) -> str:
    response = session.get(url, timeout=180)
    response.raise_for_status()
    return response.text


def fetch_binary(session: requests.Session, url: str) -> bytes:
    response = session.get(url, timeout=180)
    response.raise_for_status()
    return response.content


def find_dataset_file(dataset_record: dict, filename: str) -> dict:
    for entry in dataset_record["data"]["latestVersion"]["files"]:
        data_file = entry["dataFile"]
        if data_file["filename"] == filename:
            return data_file
    raise RuntimeError(f"Could not find dataset file {filename}")


def parse_hourly(payload: dict, iso_hour: str, key: str) -> float:
    hourly = payload["hourly"]
    index = hourly["time"].index(iso_hour)
    value = hourly[key][index]
    if value is None:
        raise RuntimeError(f"{key} was null for {iso_hour}")
    return float(value)


def parse_hourly_nullable(payload: dict, iso_hour: str, key: str) -> float | None:
    hourly = payload["hourly"]
    index = hourly["time"].index(iso_hour)
    value = hourly[key][index]
    if value is None:
        return None
    return float(value)


def saturation_vapor_pressure_pa(temperature_k: float) -> float:
    temperature_c = temperature_k - 273.15
    return 611.2 * math.exp((17.67 * temperature_c) / (temperature_c + 243.5))


def interpolate_profile_value(x_values: list[float], y_values: list[float], x_query: float) -> float:
    if x_query <= x_values[0]:
        return y_values[0]
    if x_query >= x_values[-1]:
        return y_values[-1]
    for index in range(len(x_values) - 1):
        x0 = x_values[index]
        x1 = x_values[index + 1]
        if x0 <= x_query <= x1:
            weight = (x_query - x0) / max(x1 - x0, 1.0e-12)
            return y_values[index] * (1.0 - weight) + y_values[index + 1] * weight
    return y_values[-1]


def pressure_level_url(case_spec: dict, case_day: str) -> str:
    hourly_keys: list[str] = []
    for level_hpa in PRESSURE_LEVELS_HPA:
        hourly_keys.extend(
            [
                f"temperature_{level_hpa}hPa",
                f"relative_humidity_{level_hpa}hPa",
                f"geopotential_height_{level_hpa}hPa",
            ]
        )
    return (
        "https://historical-forecast-api.open-meteo.com/v1/forecast"
        f"?latitude={case_spec['site_lat_deg']:.6f}&longitude={case_spec['site_lon_deg']:.6f}"
        f"&start_date={case_day}&end_date={case_day}"
        f"&hourly={','.join(hourly_keys)}"
        "&timezone=UTC"
    )


def build_pressure_level_rows(
    payload: dict,
    iso_hour: str,
    surface_altitude_m: float,
    surface_pressure_pa: float,
    surface_temperature_k: float,
    surface_relative_humidity: float,
) -> list[dict]:
    rows = [
        {
            "altitude_m": surface_altitude_m,
            "pressure_pa": surface_pressure_pa,
            "temperature_k": surface_temperature_k,
            "relative_humidity": clamp(surface_relative_humidity, 0.0, 1.0),
        }
    ]
    for level_hpa in PRESSURE_LEVELS_HPA:
        temperature_c = parse_hourly_nullable(payload, iso_hour, f"temperature_{level_hpa}hPa")
        height_m = parse_hourly_nullable(payload, iso_hour, f"geopotential_height_{level_hpa}hPa")
        if temperature_c is None or height_m is None:
            continue
        relative_humidity_percent = parse_hourly_nullable(payload, iso_hour, f"relative_humidity_{level_hpa}hPa")
        relative_humidity = 0.0 if relative_humidity_percent is None else clamp(relative_humidity_percent / 100.0, 0.0, 1.0)
        altitude_m = float(height_m)
        if altitude_m <= surface_altitude_m + 25.0:
            continue
        rows.append(
            {
                "altitude_m": altitude_m,
                "pressure_pa": float(level_hpa) * 100.0,
                "temperature_k": 273.15 + float(temperature_c),
                "relative_humidity": relative_humidity,
            }
        )
    rows.sort(key=lambda row: row["altitude_m"])
    return rows


def parse_aeronet_rows(text: str) -> list[dict]:
    lines = [line for line in text.splitlines() if line.strip()]
    header_index = next(i for i, line in enumerate(lines) if line.startswith("AERONET_Site,"))
    return list(csv.DictReader(lines[header_index:]))


def parse_aeronet_timestamp(row: dict) -> datetime:
    return datetime.strptime(
        f"{row['Date(dd:mm:yyyy)']} {row['Time(hh:mm:ss)']}",
        "%d:%m:%Y %H:%M:%S",
    ).replace(tzinfo=timezone.utc)


def row_float(row: dict, key: str, default: float | None = None) -> float:
    value = row.get(key, "")
    if value is None or value == "":
        if default is None:
            raise KeyError(key)
        return default
    return float(value)


def interpolate_wavelength_series(points: list[tuple[float, float]], target_nm: float) -> float:
    points = sorted(points, key=lambda pair: pair[0])
    if target_nm <= points[0][0]:
        return points[0][1]
    if target_nm >= points[-1][0]:
        return points[-1][1]
    for index in range(len(points) - 1):
        wave0, value0 = points[index]
        wave1, value1 = points[index + 1]
        if wave0 <= target_nm <= wave1:
            weight = (target_nm - wave0) / max(wave1 - wave0, 1.0e-12)
            return value0 * (1.0 - weight) + value1 * weight
    return points[-1][1]


def spectral_angstrom(values: list[tuple[float, float]]) -> float:
    positive = [(wave, value) for wave, value in values if value > 0.0]
    if len(positive) < 2:
        return 1.0
    wave0, value0 = positive[0]
    wave1, value1 = positive[-1]
    return -math.log(value1 / value0) / max(math.log(wave1 / wave0), 1.0e-12)


def choose_aeronet_inversion_row(rows: list[dict], case_timestamp: datetime) -> tuple[dict | None, float | None, str]:
    if not rows:
        return None, None, "none"

    candidates: list[tuple[float, dict]] = []
    same_day_candidates: list[tuple[float, dict]] = []
    for row in rows:
        timestamp = parse_aeronet_timestamp(row)
        delta_hours = abs((case_timestamp - timestamp).total_seconds()) / 3600.0
        candidates.append((delta_hours, row))
        if timestamp.date() == case_timestamp.date():
            same_day_candidates.append((delta_hours, row))

    nearest_delta, nearest_row = min(candidates, key=lambda item: item[0])
    if nearest_delta <= 2.25:
        return nearest_row, nearest_delta, "matched_within_window"
    if same_day_candidates:
        same_day_delta, same_day_row = min(same_day_candidates, key=lambda item: item[0])
        return same_day_row, same_day_delta, "same_day_fallback"
    if nearest_delta <= 24.0:
        return nearest_row, nearest_delta, "nearest_day_fallback"
    return None, None, "none"


def build_case_aerosol_summary_from_inversion(
    inversion_row: dict,
    target_aod_550: float,
    default_summary: dict,
) -> dict:
    fine_ext_series = [
        (440.0, row_float(inversion_row, "AOD_Extinction-Fine[440nm]", 0.0)),
        (675.0, row_float(inversion_row, "AOD_Extinction-Fine[675nm]", 0.0)),
        (870.0, row_float(inversion_row, "AOD_Extinction-Fine[870nm]", 0.0)),
        (1020.0, row_float(inversion_row, "AOD_Extinction-Fine[1020nm]", 0.0)),
    ]
    coarse_ext_series = [
        (440.0, row_float(inversion_row, "AOD_Extinction-Coarse[440nm]", 0.0)),
        (675.0, row_float(inversion_row, "AOD_Extinction-Coarse[675nm]", 0.0)),
        (870.0, row_float(inversion_row, "AOD_Extinction-Coarse[870nm]", 0.0)),
        (1020.0, row_float(inversion_row, "AOD_Extinction-Coarse[1020nm]", 0.0)),
    ]
    total_ext_series = [
        (440.0, row_float(inversion_row, "AOD_Extinction-Total[440nm]", 0.0)),
        (675.0, row_float(inversion_row, "AOD_Extinction-Total[675nm]", 0.0)),
        (870.0, row_float(inversion_row, "AOD_Extinction-Total[870nm]", 0.0)),
        (1020.0, row_float(inversion_row, "AOD_Extinction-Total[1020nm]", 0.0)),
    ]
    fine550 = interpolate_wavelength_series(fine_ext_series, 550.0)
    coarse550 = interpolate_wavelength_series(coarse_ext_series, 550.0)
    total550 = max(interpolate_wavelength_series(total_ext_series, 550.0), 1.0e-12)
    coarse_fraction = clamp(coarse550 / total550, 0.05, 0.90)
    fine_fraction = clamp(1.0 - coarse_fraction, 0.10, 0.95)

    default_fine_total = sum(
        default_summary["scattering_column_fraction"][key]
        for key in ("boundary_fine", "elevated_fine", "residual_sulfate")
    )
    if default_fine_total <= 1.0e-12:
        default_fine_total = 1.0

    scattering_column_fraction = {
        "boundary_fine": fine_fraction * default_summary["scattering_column_fraction"]["boundary_fine"] / default_fine_total,
        "elevated_fine": fine_fraction * default_summary["scattering_column_fraction"]["elevated_fine"] / default_fine_total,
        "residual_sulfate": fine_fraction * default_summary["scattering_column_fraction"]["residual_sulfate"] / default_fine_total,
        "coastal_coarse": coarse_fraction,
    }

    ssa_series = [
        (440.0, row_float(inversion_row, "Single_Scattering_Albedo[440nm]", default_summary["ssa_550"])),
        (675.0, row_float(inversion_row, "Single_Scattering_Albedo[675nm]", default_summary["ssa_550"])),
        (870.0, row_float(inversion_row, "Single_Scattering_Albedo[870nm]", default_summary["ssa_550"])),
        (1020.0, row_float(inversion_row, "Single_Scattering_Albedo[1020nm]", default_summary["ssa_550"])),
    ]
    fine_g_series = [
        (440.0, row_float(inversion_row, "Asymmetry_Factor-Fine[440nm]", default_summary["fine_asymmetry"])),
        (675.0, row_float(inversion_row, "Asymmetry_Factor-Fine[675nm]", default_summary["fine_asymmetry"])),
        (870.0, row_float(inversion_row, "Asymmetry_Factor-Fine[870nm]", default_summary["fine_asymmetry"])),
        (1020.0, row_float(inversion_row, "Asymmetry_Factor-Fine[1020nm]", default_summary["fine_asymmetry"])),
    ]
    coarse_g_series = [
        (440.0, row_float(inversion_row, "Asymmetry_Factor-Coarse[440nm]", default_summary["coastal_coarse_asymmetry"])),
        (675.0, row_float(inversion_row, "Asymmetry_Factor-Coarse[675nm]", default_summary["coastal_coarse_asymmetry"])),
        (870.0, row_float(inversion_row, "Asymmetry_Factor-Coarse[870nm]", default_summary["coastal_coarse_asymmetry"])),
        (1020.0, row_float(inversion_row, "Asymmetry_Factor-Coarse[1020nm]", default_summary["coastal_coarse_asymmetry"])),
    ]
    depol_series = [
        (440.0, row_float(inversion_row, "Depolarization_Ratio[440nm]", 0.0)),
        (675.0, row_float(inversion_row, "Depolarization_Ratio[675nm]", 0.0)),
        (870.0, row_float(inversion_row, "Depolarization_Ratio[870nm]", 0.0)),
        (1020.0, row_float(inversion_row, "Depolarization_Ratio[1020nm]", 0.0)),
    ]
    ri_real_series = [
        (440.0, row_float(inversion_row, "Refractive_Index-Real_Part[440nm]", 1.50)),
        (675.0, row_float(inversion_row, "Refractive_Index-Real_Part[675nm]", 1.50)),
        (870.0, row_float(inversion_row, "Refractive_Index-Real_Part[870nm]", 1.50)),
        (1020.0, row_float(inversion_row, "Refractive_Index-Real_Part[1020nm]", 1.50)),
    ]
    ri_imag_series = [
        (440.0, row_float(inversion_row, "Refractive_Index-Imaginary_Part[440nm]", 0.003)),
        (675.0, row_float(inversion_row, "Refractive_Index-Imaginary_Part[675nm]", 0.003)),
        (870.0, row_float(inversion_row, "Refractive_Index-Imaginary_Part[870nm]", 0.003)),
        (1020.0, row_float(inversion_row, "Refractive_Index-Imaginary_Part[1020nm]", 0.003)),
    ]

    sphericity = clamp(row_float(inversion_row, "Sphericity_Factor(%)", 60.0) / 100.0, 0.0, 1.0)
    depol_550 = interpolate_wavelength_series(depol_series, 550.0)

    summary = dict(default_summary)
    summary.update(
        {
            "model_name": "marseille_case_local_aeronet_inversion_v2",
            "description": (
                "Case-local aerosol model using Marseille_ATMO AERONET Almucantar inversion microphysics when a "
                "strict within-window inversion is unavailable but a same-day inversion exists."
            ),
            "scattering_column_fraction": scattering_column_fraction,
            "ssa_550": clamp(interpolate_wavelength_series(ssa_series, 550.0), 0.85, 0.995),
            "coastal_coarse_asymmetry": clamp(interpolate_wavelength_series(coarse_g_series, 550.0), 0.70, 0.92),
            "fine_asymmetry": clamp(interpolate_wavelength_series(fine_g_series, 550.0), 0.50, 0.82),
            "angstrom_fine": clamp(spectral_angstrom(fine_ext_series), 0.7, 2.2),
            "angstrom_coarse": clamp(spectral_angstrom(coarse_ext_series), -0.2, 1.0),
            "pmax_base": clamp(0.012 + 0.45 * depol_550 + 0.035 * (1.0 - sphericity), 0.02, 0.12),
            "inversion_refractive_index_real_550": interpolate_wavelength_series(ri_real_series, 550.0),
            "inversion_refractive_index_imag_550": interpolate_wavelength_series(ri_imag_series, 550.0),
            "inversion_depolarization_ratio_550": depol_550,
            "inversion_sphericity_fraction": sphericity,
            "inversion_aod_fine_fraction_550": clamp(fine550 / total550, 0.0, 1.0),
            "inversion_aod_coarse_fraction_550": clamp(coarse550 / total550, 0.0, 1.0),
        }
    )
    return summary


def load_template_profile() -> list[dict]:
    rows = []
    with (DATA_DIR / "atmosphere" / "clear_sky_midlatitude.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append({key: float(value) for key, value in row.items()})
    return rows


def write_text(path: Path, content: str) -> None:
    path.write_text(content)


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2))


def build_case(case_spec: dict) -> None:
    session = requests.Session()
    session.headers.update({"User-Agent": "light-scattering-paper-case-builder/1.0"})

    case_id = case_spec["case_id"]
    case_timestamp = datetime.fromisoformat(case_spec["timestamp_utc"].replace("Z", "+00:00"))
    case_day = case_timestamp.strftime("%Y-%m-%d")
    case_hour = case_timestamp.strftime("%Y-%m-%dT%H:00")
    zenith_deg, azimuth_deg = solar_geometry_utc(case_timestamp, case_spec["site_lat_deg"], case_spec["site_lon_deg"])

    PAPER_CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    case_data_dir = PAPER_DATA_DIR / case_id
    case_data_dir.mkdir(parents=True, exist_ok=True)

    dataset_record = fetch_json(session, case_spec["dataset_api_url"])
    annotation_file = find_dataset_file(dataset_record, case_spec["annotation_filename"])
    readme_file = find_dataset_file(dataset_record, case_spec["dataset_readme_name"])
    annotation_bytes = fetch_binary(session, f"https://entrepot.recherche.data.gouv.fr/api/access/datafile/{annotation_file['id']}")
    readme_text = fetch_text(session, f"https://entrepot.recherche.data.gouv.fr/api/access/datafile/{readme_file['id']}")
    annotation_array = np.load(io.BytesIO(annotation_bytes), allow_pickle=True)

    expected_stamp = case_timestamp.strftime("%Y-%m-%dT%H-%M-%S")
    label = None
    for timestamp_text, value in annotation_array:
        if str(timestamp_text) == expected_stamp:
            label = str(value)
            break
    if label != case_spec["annotation_label"]:
        raise RuntimeError(f"Expected annotation label {case_spec['annotation_label']} for {expected_stamp}, got {label!r}")

    weather_url = (
        "https://archive-api.open-meteo.com/v1/archive"
        f"?latitude={case_spec['site_lat_deg']:.6f}&longitude={case_spec['site_lon_deg']:.6f}"
        f"&start_date={case_day}&end_date={case_day}"
        "&hourly=temperature_2m,relative_humidity_2m,surface_pressure,cloud_cover,total_column_integrated_water_vapour"
        "&timezone=UTC"
    )
    air_quality_url = (
        "https://air-quality-api.open-meteo.com/v1/air-quality"
        f"?latitude={case_spec['site_lat_deg']:.6f}&longitude={case_spec['site_lon_deg']:.6f}"
        f"&start_date={case_day}&end_date={case_day}"
        "&hourly=ozone,nitrogen_dioxide,aerosol_optical_depth,pm10,pm2_5"
        "&timezone=UTC"
    )
    pressure_profile_url = pressure_level_url(case_spec, case_day)
    weather_payload = fetch_json(session, weather_url)
    air_quality_payload = fetch_json(session, air_quality_url)
    pressure_profile_payload = fetch_json(session, pressure_profile_url)

    next_day = case_timestamp + timedelta(days=1)
    aeronet_url = (
        "https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3"
        f"?site=Marseille_ATMO&year={case_timestamp.year}&month={case_timestamp.month}&day={case_timestamp.day}"
        f"&year2={next_day.year}&month2={next_day.month}&day2={next_day.day}"
        "&AOD15=1&AVG=10&if_no_html=1"
    )
    aeronet_text = fetch_text(session, aeronet_url)
    aeronet_rows = parse_aeronet_rows(aeronet_text)
    aeronet_inversion_url = (
        "https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_inv_v3"
        f"?site=Marseille_ATMO&year={case_timestamp.year}&month={case_timestamp.month}&day={case_timestamp.day}"
        f"&year2={next_day.year}&month2={next_day.month}&day2={next_day.day}"
        "&product=ALL&ALM15=1&AVG=10&if_no_html=1"
    )
    try:
        aeronet_inversion_text = fetch_text(session, aeronet_inversion_url)
    except requests.RequestException:
        aeronet_inversion_text = ""
    chosen_row = None
    chosen_time = None
    for row in aeronet_rows:
        if row["AERONET_Site"] != "Marseille_ATMO":
            continue
        row_time = datetime.strptime(f"{row['Date(dd:mm:yyyy)']} {row['Time(hh:mm:ss)']}", "%d:%m:%Y %H:%M:%S").replace(tzinfo=timezone.utc)
        if row_time <= case_timestamp and (chosen_time is None or row_time > chosen_time):
            chosen_row = row
            chosen_time = row_time
    if chosen_row is None:
        raise RuntimeError("Could not find a usable Marseille_ATMO row before the frozen twilight time")
    aeronet_inversion_rows = parse_aeronet_rows(aeronet_inversion_text) if "AERONET_Site," in aeronet_inversion_text else []
    aeronet_inversion_row, aeronet_inversion_time_delta_hours, aeronet_inversion_usage_mode = choose_aeronet_inversion_row(
        aeronet_inversion_rows,
        case_timestamp,
    )
    aeronet_inversion_available = aeronet_inversion_row is not None
    aeronet_inversion_selected_time_utc = (
        parse_aeronet_timestamp(aeronet_inversion_row).isoformat().replace("+00:00", "Z")
        if aeronet_inversion_row is not None
        else None
    )
    aeronet_inversion_quality_level = (
        aeronet_inversion_row.get("Inversion_Data_Quality_Level", "")
        if aeronet_inversion_row is not None
        else ""
    )

    template_rows = load_template_profile()
    template_altitudes = [row["altitude_m"] for row in template_rows]
    template_ozone = [row["ozone_number_density_m3"] for row in template_rows]
    surface_temperature_k = 273.15 + parse_hourly(weather_payload, case_hour, "temperature_2m")
    surface_relative_humidity = parse_hourly(weather_payload, case_hour, "relative_humidity_2m") / 100.0
    surface_pressure_pa = 100.0 * parse_hourly(weather_payload, case_hour, "surface_pressure")
    tcwv_kg_m2 = parse_hourly(weather_payload, case_hour, "total_column_integrated_water_vapour")
    cloud_cover_percent = parse_hourly(weather_payload, case_hour, "cloud_cover")
    surface_ozone_ug_m3 = parse_hourly(air_quality_payload, case_hour, "ozone")
    surface_no2_ug_m3 = parse_hourly(air_quality_payload, case_hour, "nitrogen_dioxide")
    openmeteo_aod_550 = parse_hourly(air_quality_payload, case_hour, "aerosol_optical_depth")
    surface_pm10_ug_m3 = parse_hourly(air_quality_payload, case_hour, "pm10")
    surface_pm25_ug_m3 = parse_hourly(air_quality_payload, case_hour, "pm2_5")
    aeronet_aod_500 = float(chosen_row["AOD_500nm"])
    aeronet_aod_440 = float(chosen_row["AOD_440nm"])
    angstrom = float(chosen_row["440-870_Angstrom_Exponent"])
    aeronet_aod_550 = aeronet_aod_500 * (550.0 / 500.0) ** (-angstrom)
    aeronet_time_delta_hours = abs((case_timestamp - chosen_time).total_seconds()) / 3600.0
    aerosol_target_aod_550, aeronet_aod_weight = blend_aod_550(openmeteo_aod_550, aeronet_aod_550, aeronet_time_delta_hours)
    ozone_column = float(chosen_row["Ozone(Dobson)"]) * DOBSON_TO_MOLECULES_PER_M2
    no2_column = float(chosen_row["NO2(Dobson)"]) * DOBSON_TO_MOLECULES_PER_M2
    water_column = (tcwv_kg_m2 / WATER_MOLAR_MASS_KG) * AVOGADRO
    aerosol_summary = build_case_aerosol_summary(
        surface_relative_humidity=surface_relative_humidity,
        pm10_ug_m3=surface_pm10_ug_m3,
        pm25_ug_m3=surface_pm25_ug_m3,
        target_aod_550=aerosol_target_aod_550,
        angstrom=angstrom,
    )
    if aeronet_inversion_row is not None:
        aerosol_summary = build_case_aerosol_summary_from_inversion(
            aeronet_inversion_row,
            aerosol_target_aod_550,
            aerosol_summary,
        )

    matched_thermo_rows = build_pressure_level_rows(
        pressure_profile_payload,
        case_hour,
        case_spec["observer_altitude_m"],
        surface_pressure_pa,
        surface_temperature_k,
        surface_relative_humidity,
    )
    matched_profile_top_altitude_m = matched_thermo_rows[-1]["altitude_m"]
    tail_template_candidates = [row for row in template_rows if row["altitude_m"] > matched_profile_top_altitude_m + 500.0]
    if tail_template_candidates:
        join_reference_altitude_m = matched_thermo_rows[-1]["altitude_m"]
        join_template_row = min(
            template_rows,
            key=lambda row: abs(row["altitude_m"] - join_reference_altitude_m),
        )
        tail_temperature_offset = matched_thermo_rows[-1]["temperature_k"] - join_template_row["temperature_k"]
        tail_pressure_scale = matched_thermo_rows[-1]["pressure_pa"] / max(join_template_row["pressure_pa"], 1.0)
    else:
        tail_temperature_offset = 0.0
        tail_pressure_scale = 1.0

    profile_rows: list[dict] = []
    for row in matched_thermo_rows:
        vapor_pressure_pa = clamp(
            row["relative_humidity"],
            0.0,
            1.0,
        ) * saturation_vapor_pressure_pa(row["temperature_k"])
        h2o_number_density = vapor_pressure_pa / (BOLTZMANN * row["temperature_k"])
        molecular_density = row["pressure_pa"] / (BOLTZMANN * row["temperature_k"])
        profile_rows.append(
            {
                "altitude_m": row["altitude_m"],
                "pressure_pa": row["pressure_pa"],
                "temperature_k": row["temperature_k"],
                "molecular_number_density_m3": molecular_density,
                "h2o_number_density_m3": h2o_number_density,
                "source_mode": "openmeteo_pressure_levels",
            }
        )
    for row in tail_template_candidates:
        altitude_m = row["altitude_m"]
        temperature_k = max(185.0, row["temperature_k"] + tail_temperature_offset)
        pressure_pa = max(1.0, row["pressure_pa"] * tail_pressure_scale)
        molecular_density = pressure_pa / (BOLTZMANN * temperature_k)
        profile_rows.append(
            {
                "altitude_m": altitude_m,
                "pressure_pa": pressure_pa,
                "temperature_k": temperature_k,
                "molecular_number_density_m3": molecular_density,
                "h2o_number_density_m3": 0.0,
                "source_mode": "template_upper_tail",
            }
        )

    profile_rows.sort(key=lambda row: row["altitude_m"])
    altitudes = [row["altitude_m"] for row in profile_rows]
    ozone_shape = [
        interpolate_profile_value(template_altitudes, template_ozone, altitude_m)
        for altitude_m in altitudes
    ]
    ozone_scale = ozone_column / max(trapz(altitudes, ozone_shape), 1.0)
    water_shape = [row["h2o_number_density_m3"] for row in profile_rows]
    if trapz(altitudes, water_shape) <= 1.0e-12:
        water_shape = [math.exp(-altitude_m / 2000.0) for altitude_m in altitudes]
    water_scale = water_column / max(trapz(altitudes, water_shape), 1.0)
    no2_shape = [0.85 * math.exp(-altitude_m / 1200.0) + 0.15 * math.exp(-altitude_m / 7000.0) for altitude_m in altitudes]
    no2_scale = no2_column / max(trapz(altitudes, no2_shape), 1.0)
    aerosol_shapes = aerosol_component_shapes(altitudes)
    normalized_shapes = {name: normalize_shape(altitudes, values) for name, values in aerosol_shapes.items()}
    scattering_column_550 = aerosol_target_aod_550 * aerosol_summary["ssa_550"]
    absorption_column_550 = aerosol_target_aod_550 * (1.0 - aerosol_summary["ssa_550"])
    scattering_components = aerosol_summary["scattering_column_fraction"]
    black_carbon_fraction = aerosol_summary["absorption_black_carbon_fraction"]

    atmosphere_profile_path = case_data_dir / "atmosphere_profile.csv"
    with atmosphere_profile_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "altitude_m",
                "pressure_pa",
                "temperature_k",
                "molecular_number_density_m3",
                "ozone_number_density_m3",
                "h2o_number_density_m3",
                "no2_number_density_m3",
                "aerosol_extinction_550_m_inv",
                "aerosol_single_scattering_albedo",
                "aerosol_asymmetry",
                "aerosol_scattering_angstrom_exponent",
                "aerosol_absorption_angstrom_exponent",
            ],
        )
        writer.writeheader()
        for index, row in enumerate(profile_rows):
            local_boundary_fine = scattering_column_550 * scattering_components["boundary_fine"] * normalized_shapes["boundary_fine"][index]
            local_coastal_coarse = scattering_column_550 * scattering_components["coastal_coarse"] * normalized_shapes["coastal_coarse"][index]
            local_elevated_fine = scattering_column_550 * scattering_components["elevated_fine"] * normalized_shapes["elevated_fine"][index]
            local_residual_sulfate = scattering_column_550 * scattering_components["residual_sulfate"] * normalized_shapes["residual_sulfate"][index]
            local_scattering = local_boundary_fine + local_coastal_coarse + local_elevated_fine + local_residual_sulfate

            local_black_carbon = absorption_column_550 * black_carbon_fraction * normalized_shapes["black_carbon"][index]
            local_brown_carbon = absorption_column_550 * (1.0 - black_carbon_fraction) * normalized_shapes["brown_carbon"][index]
            local_absorption = local_black_carbon + local_brown_carbon
            local_extinction = local_scattering + local_absorption

            if local_scattering > 0.0:
                aerosol_asymmetry = (
                    local_boundary_fine * aerosol_summary["fine_asymmetry"]
                    + local_coastal_coarse * aerosol_summary["coastal_coarse_asymmetry"]
                    + local_elevated_fine * clamp(aerosol_summary["fine_asymmetry"] - 0.04, 0.60, 0.72)
                    + local_residual_sulfate * 0.58
                ) / local_scattering
                aerosol_scattering_angstrom = (
                    local_boundary_fine * aerosol_summary["angstrom_fine"]
                    + local_coastal_coarse * aerosol_summary["angstrom_coarse"]
                    + local_elevated_fine * clamp(aerosol_summary["angstrom_fine"] - 0.10, 1.20, 1.60)
                    + local_residual_sulfate * 2.05
                ) / local_scattering
            else:
                aerosol_asymmetry = aerosol_summary["coastal_coarse_asymmetry"]
                aerosol_scattering_angstrom = aerosol_summary["angstrom_fine"]

            if local_absorption > 0.0:
                aerosol_absorption_angstrom = (
                    local_black_carbon * 1.0 + local_brown_carbon * 1.85
                ) / local_absorption
            else:
                aerosol_absorption_angstrom = 1.20
            writer.writerow(
                {
                    "altitude_m": f"{row['altitude_m']:.1f}",
                    "pressure_pa": f"{row['pressure_pa']:.6e}",
                    "temperature_k": f"{row['temperature_k']:.6f}",
                    "molecular_number_density_m3": f"{row['molecular_number_density_m3']:.6e}",
                    "ozone_number_density_m3": f"{ozone_shape[index] * ozone_scale:.6e}",
                    "h2o_number_density_m3": f"{water_shape[index] * water_scale:.6e}",
                    "no2_number_density_m3": f"{no2_shape[index] * no2_scale:.6e}",
                    "aerosol_extinction_550_m_inv": f"{local_extinction:.6e}",
                    "aerosol_single_scattering_albedo": f"{(local_scattering / max(local_extinction, 1.0e-12)):.6f}",
                    "aerosol_asymmetry": f"{aerosol_asymmetry:.6f}",
                    "aerosol_scattering_angstrom_exponent": f"{aerosol_scattering_angstrom:.6f}",
                    "aerosol_absorption_angstrom_exponent": f"{aerosol_absorption_angstrom:.6f}",
                }
            )

    instrument_response_path = case_data_dir / "instrument_response.csv"
    with instrument_response_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["wavelength_nm", "instrument_response"])
        for wavelength_nm, response in imx250myr_blue_rows_public_docs_proxy(
            case_spec["channel_center_nm"],
            case_spec["channel_sigma_nm"],
        ):
            writer.writerow([f"{wavelength_nm:.1f}", f"{response:.6e}"])

    aerosol_phase_matrix_path = case_data_dir / "aerosol_phase_matrix.csv"
    with aerosol_phase_matrix_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PHASE_FIELDNAMES)
        writer.writeheader()
        for wavelength_nm in range(350, 801, 25):
            for angle_deg in range(0, 181, 5):
                writer.writerow(aerosol_phase_row(aerosol_summary, float(wavelength_nm), float(angle_deg)))

    aerosol_optics_metadata_path = case_data_dir / "aerosol_optics.json"
    write_json(
        aerosol_optics_metadata_path,
        {
            "model_name": aerosol_summary.get("model_name", "marseille_case_local_coastal_mixed_aerosol_v1"),
            "description": aerosol_summary.get(
                "description",
                "Case-local coastal/urban-marine aerosol proxy built from matched Marseille AERONET AOD/Angstrom, Open-Meteo PM and RH, and a low-polarization coastal aerosol prior.",
            ),
            "target_aod_550": aerosol_target_aod_550,
            "aeronet_aod_550": aeronet_aod_550,
            "openmeteo_aod_550": openmeteo_aod_550,
            "aeronet_weight": aeronet_aod_weight,
            "aeronet_inversion_usage_mode": aeronet_inversion_usage_mode,
            "aeronet_inversion_time_delta_hours": aeronet_inversion_time_delta_hours,
            "aeronet_inversion_selected_time_utc": aeronet_inversion_selected_time_utc,
            "aeronet_inversion_quality_level": aeronet_inversion_quality_level,
            "scattering_column_fraction": aerosol_summary["scattering_column_fraction"],
            "ssa_550": aerosol_summary["ssa_550"],
            "pm_ratio": aerosol_summary["pm_ratio"],
            "mixed_fine_indicator": aerosol_summary["mixed_fine_indicator"],
            "pmax_base": aerosol_summary["pmax_base"],
            "coastal_coarse_asymmetry": aerosol_summary["coastal_coarse_asymmetry"],
            "fine_asymmetry": aerosol_summary["fine_asymmetry"],
            "angstrom_fine": aerosol_summary["angstrom_fine"],
            "angstrom_coarse": aerosol_summary["angstrom_coarse"],
            "instrument_response_model": "imx250myr_blue_public_docs_proxy_v2",
            "instrument_response_source_urls": [
                case_spec["sensor_flyer_url"],
                case_spec["camera_techref_url"],
                case_spec["camera_techref_download_url"],
            ],
            "inversion_refractive_index_real_550": aerosol_summary.get("inversion_refractive_index_real_550"),
            "inversion_refractive_index_imag_550": aerosol_summary.get("inversion_refractive_index_imag_550"),
            "inversion_depolarization_ratio_550": aerosol_summary.get("inversion_depolarization_ratio_550"),
            "inversion_sphericity_fraction": aerosol_summary.get("inversion_sphericity_fraction"),
        },
    )

    surface_config_path = case_data_dir / "surface_config.json"
    write_json(
        surface_config_path,
        {
            "surface_model": "lambertian_land",
            "default_albedo": 0.15,
            "notes": "Real matched Marseille twilight case. Land Lambertian surface retained in the current paper path.",
        },
    )

    measurement_reference_path = case_data_dir / "measurement_reference.csv"
    measurement_reduction_path = case_data_dir / "measurement_reduction.json"
    has_measurement_reference = measurement_reference_path.exists()
    paper_gate_role = (
        "frozen_primary_twilight_case_with_extracted_reference"
        if has_measurement_reference
        else "frozen_primary_twilight_candidate_without_extracted_reference"
    )
    paper_primary_measurement_frozen = "true" if has_measurement_reference else "false"

    measurement_metadata_path = case_data_dir / "measurement_metadata.json"
    write_json(
        measurement_metadata_path,
        {
            "source_name": case_spec["measurement_dataset_title"],
            "source_url": case_spec["dataset_record_url"],
            "source_api_url": case_spec["dataset_api_url"],
            "doi": case_spec["dataset_doi"],
            "measurement_paper_url": case_spec["measurement_paper_url"],
            "site_lat_deg": case_spec["site_lat_deg"],
            "site_lon_deg": case_spec["site_lon_deg"],
            "utc_timestamp": case_spec["timestamp_utc"],
            "channel_description": (
                "Public Marseille dataset B channel on the IMX250MYR sensor. The instrument-response CSV is a "
                "case-local public-doc-constrained IMX250MYR blue-channel proxy, anchored to the Sony Polarsens "
                "flyer and the LUCID PHX050S polarized technical reference because the dataset release does not "
                "provide a tabulated manufacturer spectral-response curve."
            ),
            "channel_response_source_urls": [
                case_spec["sensor_flyer_url"],
                case_spec["camera_techref_url"],
                case_spec["camera_techref_download_url"],
            ],
            "annotation_filename": case_spec["annotation_filename"],
            "annotation_label": label,
            "below_horizon_twilight": True,
            "solar_zenith_deg": zenith_deg,
            "solar_azimuth_deg": azimuth_deg,
            "paper_gate_role": paper_gate_role,
            **(
                {"measurement_reference_csv": str(measurement_reference_path)}
                if has_measurement_reference
                else {}
            ),
        },
    )

    write_json(
        case_data_dir / "measurement_pointer.json",
        {
            "annotation_filename": case_spec["annotation_filename"],
            "annotation_file_id": int(annotation_file["id"]),
            "selected_timestamp_utc": case_spec["timestamp_utc"],
            "selected_annotation_label": label,
        },
    )
    write_json(case_data_dir / "dataset_record.json", dataset_record)
    write_text(case_data_dir / "dataset_readme.md", readme_text)
    (case_data_dir / case_spec["annotation_filename"]).write_bytes(annotation_bytes)
    write_json(case_data_dir / "open_meteo_weather.json", weather_payload)
    write_json(case_data_dir / "open_meteo_pressure_levels.json", pressure_profile_payload)
    write_json(case_data_dir / "open_meteo_air_quality.json", air_quality_payload)
    write_text(case_data_dir / "aeronet_direct_sun.csv", aeronet_text)
    write_text(case_data_dir / "aeronet_inversion.txt", aeronet_inversion_text)

    measurement_cfg_path = PAPER_CONFIG_DIR / f"{case_id}_measurement.cfg"
    measurement_cfg_rel = f"paper_cases/{measurement_cfg_path.name}"
    measurement_cfg_lines = [
        *(
            [
                "# Frozen Marseille measurement config with extracted validator reference.",
                f"measurement_reference_csv=../../data/paper_cases/{case_id}/measurement_reference.csv",
            ]
            if has_measurement_reference
            else [
                "# Placeholder measurement config for the frozen Marseille twilight case.",
                "# measurement_reference_csv is intentionally absent until the raw public dataset is reduced into the validation reference format.",
            ]
        ),
        f"case_id={case_id}_measurement",
        "output_dir=../../results",
        f"profile_csv=../../data/paper_cases/{case_id}/atmosphere_profile.csv",
        "solar_spectrum_csv=../../data/optics/solar_irradiance_reference.csv",
        f"instrument_response_csv=../../data/paper_cases/{case_id}/instrument_response.csv",
        "ozone_cross_section_csv=../../data/optics/ozone_cross_section_reference.csv",
        "o2_cross_section_csv=../../data/optics/o2_cross_section_reference.csv",
        "o4_cross_section_csv=../../data/optics/o4_cross_section_reference.csv",
        "h2o_cross_section_csv=../../data/optics/h2o_cross_section_reference.csv",
        "no2_cross_section_csv=../../data/optics/no2_cross_section_reference.csv",
        f"aerosol_phase_matrix_csv=../../data/paper_cases/{case_id}/aerosol_phase_matrix.csv",
        "surface_albedo_csv=../../data/optics/lambertian_land_albedo.csv",
        "surface_model=lambertian_land",
        f"surface_parameter_csv=../../data/paper_cases/{case_id}/surface_config.json",
        "default_surface_albedo=0.15",
        "ocean_wind_speed_m_s=5.0",
        "use_explicit_solar_angles=true",
        f"solar_zenith_deg={zenith_deg:.6f}",
        f"solar_azimuth_deg={azimuth_deg:.6f}",
        "finite_solar_disk=true",
        "solar_angular_radius_deg=0.2666",
        "solar_disk_quadrature_nodes=7",
        f"observer_latitude_deg={case_spec['site_lat_deg']:.8f}",
        f"observer_longitude_deg={case_spec['site_lon_deg']:.8f}",
        f"observer_altitude_m={case_spec['observer_altitude_m']:.1f}",
        f"measurement_metadata_json=../../data/paper_cases/{case_id}/measurement_metadata.json",
    ]
    common_lines = [
        "strict_paper_mode=true",
        "top_of_atmosphere_altitude_m=100000.0",
        "min_wavelength_nm=350.0",
        "max_wavelength_nm=800.0",
        "wavelength_step_nm=10.0",
        "zenith_bins=19",
        "azimuth_bins=36",
        "photons_per_bin=256",
        "russian_roulette_threshold=1e-4",
        "random_seed=20260326",
        "max_events_guard=64",
        "single_scatter_only=false",
        "deterministic_single_scatter=true",
        "deterministic_second_scatter=true",
        "second_scatter_view_steps=6",
        "second_scatter_ray_steps=4",
        "second_scatter_mu_nodes=3",
        "second_scatter_phi_nodes=4",
        "twilight_second_scatter_adaptive=true",
        "twilight_second_scatter_zenith_threshold_deg=25.0",
        "twilight_second_scatter_min_view_steps=8",
        "twilight_second_scatter_min_ray_steps=6",
        "twilight_second_scatter_min_mu_nodes=5",
        "twilight_second_scatter_min_phi_nodes=8",
        "source_guided_first_scatter=true",
        "source_guided_first_scatter_branches=4",
        "source_guided_phase_fraction=0.5",
        "source_guided_cone_half_angle_deg=30.0",
        "rayleigh_source_guided_phase_fraction=0.6",
        "rayleigh_source_guided_cone_half_angle_deg=25.0",
        "rayleigh_source_guided_first_scatter_branches=6",
        "rayleigh_polarization_guided_fraction=0.0",
        "rayleigh_polarization_guided_mu_half_width=0.15",
        "rayleigh_polarization_guided_branches=0",
        "twilight_higher_order_guiding=true",
        "twilight_higher_order_branches=3",
        "twilight_higher_order_phase_fraction=0.20",
        "twilight_higher_order_tangent_fraction=0.45",
        "twilight_higher_order_horizon_fraction=0.35",
        "benchmark_mask_fraction_of_peak=0.01",
        "measurement_mask_fraction_of_peak=0.05",
        "median_intensity_error_limit=0.05",
        "p95_intensity_error_limit=0.10",
        "median_dolp_abs_error_limit=0.03",
        "p95_dolp_abs_error_limit=0.07",
        "median_aop_error_deg_limit=5.0",
        "p95_aop_error_deg_limit=10.0",
        "solar_vertical_signed_dolp_bias_limit=0.05",
        "normalized_rmse_limit=0.10",
        "brightest_region_deg_limit=5.0",
        "neutral_point_location_deg_limit=5.0",
    ]

    write_text(
        measurement_cfg_path,
        "\n".join(measurement_cfg_lines + common_lines) + "\n",
    )

    strict_subset_cfg_path = PAPER_CONFIG_DIR / f"{case_id}_measurement_strict_subset.cfg"
    strict_subset_common_lines = [
        ("photons_per_bin=4" if line == "photons_per_bin=256" else line)
        for line in common_lines
    ]
    strict_subset_lines = [
        "# Frozen Marseille measurement strict-physics subset config.",
        "# This keeps the strict paper physics and full spectral band, but reduces the scored measurement directions and Monte Carlo budget for iterative debugging.",
        f"measurement_reference_csv=../../data/paper_cases/{case_id}/measurement_reference_strict_subset.csv",
        f"case_id={case_id}_measurement_strict_subset",
        *measurement_cfg_lines[2:],
    ]
    write_text(
        strict_subset_cfg_path,
        "\n".join(strict_subset_lines + strict_subset_common_lines) + "\n",
    )

    profile_subset_cfg_path = PAPER_CONFIG_DIR / f"{case_id}_measurement_profile_subset.cfg"
    profile_subset_lines = [
        "# Frozen Marseille measurement strict-paper profiling subset config.",
        "# This keeps the strict paper physics and full spectral band, but only scores a small targeted subset of Marseille directions for runtime profiling.",
        f"measurement_reference_csv=../../data/paper_cases/{case_id}/measurement_reference_profile_subset.csv",
        f"case_id={case_id}_measurement_profile_subset",
        *measurement_cfg_lines[2:],
    ]
    write_text(
        profile_subset_cfg_path,
        "\n".join(profile_subset_lines + common_lines) + "\n",
    )

    interactive_measurement_cfg_path = PAPER_CONFIG_DIR / f"{case_id}_measurement_interactive.cfg"
    interactive_lines = [
        "# Frozen Marseille twilight strict full-field interactive config.",
        "# This keeps the full extracted reference field but uses a reduced runtime budget for interactive debugging.",
        f"measurement_reference_csv=../../data/paper_cases/{case_id}/measurement_reference.csv",
        f"case_id={case_id}_measurement_interactive",
        "output_dir=../../results",
        f"profile_csv=../../data/paper_cases/{case_id}/atmosphere_profile.csv",
        "solar_spectrum_csv=../../data/optics/solar_irradiance_reference.csv",
        f"instrument_response_csv=../../data/paper_cases/{case_id}/instrument_response.csv",
        "ozone_cross_section_csv=../../data/optics/ozone_cross_section_reference.csv",
        "o2_cross_section_csv=../../data/optics/o2_cross_section_reference.csv",
        "o4_cross_section_csv=../../data/optics/o4_cross_section_reference.csv",
        "h2o_cross_section_csv=../../data/optics/h2o_cross_section_reference.csv",
        "no2_cross_section_csv=../../data/optics/no2_cross_section_reference.csv",
        f"aerosol_phase_matrix_csv=../../data/paper_cases/{case_id}/aerosol_phase_matrix.csv",
        "surface_albedo_csv=../../data/optics/lambertian_land_albedo.csv",
        "surface_model=lambertian_land",
        f"surface_parameter_csv=../../data/paper_cases/{case_id}/surface_config.json",
        "default_surface_albedo=0.15",
        "ocean_wind_speed_m_s=5.0",
        "use_explicit_solar_angles=true",
        f"solar_zenith_deg={zenith_deg:.6f}",
        f"solar_azimuth_deg={azimuth_deg:.6f}",
        "finite_solar_disk=true",
        "solar_angular_radius_deg=0.2666",
        "solar_disk_quadrature_nodes=1",
        f"observer_latitude_deg={case_spec['site_lat_deg']:.8f}",
        f"observer_longitude_deg={case_spec['site_lon_deg']:.8f}",
        f"observer_altitude_m={case_spec['observer_altitude_m']:.1f}",
        f"measurement_metadata_json=../../data/paper_cases/{case_id}/measurement_metadata.json",
        "strict_paper_mode=true",
        "top_of_atmosphere_altitude_m=100000.0",
        "min_wavelength_nm=430.0",
        "max_wavelength_nm=490.0",
        "wavelength_step_nm=30.0",
        "zenith_bins=19",
        "azimuth_bins=36",
        "photons_per_bin=1",
        "russian_roulette_threshold=1e-4",
        "random_seed=20260326",
        "max_events_guard=64",
        "optical_depth_step_scale=8.0",
        "min_optical_depth_steps=4",
        "line_integral_step_scale=8.0",
        "min_line_integral_steps=4",
        "single_scatter_only=false",
        "deterministic_single_scatter=true",
        "deterministic_second_scatter=true",
        "second_scatter_view_steps=3",
        "second_scatter_ray_steps=2",
        "second_scatter_mu_nodes=2",
        "second_scatter_phi_nodes=2",
        "twilight_second_scatter_adaptive=true",
        "twilight_second_scatter_zenith_threshold_deg=25.0",
        "twilight_second_scatter_min_view_steps=6",
        "twilight_second_scatter_min_ray_steps=4",
        "twilight_second_scatter_min_mu_nodes=4",
        "twilight_second_scatter_min_phi_nodes=6",
        "source_guided_first_scatter=true",
        "source_guided_first_scatter_branches=1",
        "source_guided_phase_fraction=0.5",
        "source_guided_cone_half_angle_deg=30.0",
        "rayleigh_source_guided_phase_fraction=0.6",
        "rayleigh_source_guided_cone_half_angle_deg=25.0",
        "rayleigh_source_guided_first_scatter_branches=1",
        "rayleigh_polarization_guided_fraction=0.0",
        "rayleigh_polarization_guided_mu_half_width=0.15",
        "rayleigh_polarization_guided_branches=0",
        "twilight_higher_order_guiding=true",
        "twilight_higher_order_branches=4",
        "twilight_higher_order_phase_fraction=0.20",
        "twilight_higher_order_tangent_fraction=0.45",
        "twilight_higher_order_horizon_fraction=0.35",
        "benchmark_mask_fraction_of_peak=0.01",
        "measurement_mask_fraction_of_peak=0.05",
        "median_intensity_error_limit=0.05",
        "p95_intensity_error_limit=0.10",
        "median_dolp_abs_error_limit=0.03",
        "p95_dolp_abs_error_limit=0.07",
        "median_aop_error_deg_limit=5.0",
        "p95_aop_error_deg_limit=10.0",
        "solar_vertical_signed_dolp_bias_limit=0.05",
        "normalized_rmse_limit=0.10",
        "brightest_region_deg_limit=5.0",
        "neutral_point_location_deg_limit=5.0",
        "paper_primary_measurement_frozen=true",
    ]
    write_text(interactive_measurement_cfg_path, "\n".join(interactive_lines) + "\n")

    case_cfg_lines = [
        "# Real frozen Marseille twilight case built from matched public inputs.",
        *(
            []
            if has_measurement_reference
            else ["# The paper gate remains blocked until the measurement field for this exact case is extracted."]
        ),
        f"case_id={case_id}",
        "output_dir=../../results",
        "benchmark_case_config=../benchmark_disort_scalar.cfg;../benchmark_iprt_a1_vector.cfg;../benchmark_zawada_spherical_vector_single.cfg;../benchmark_zawada_spherical_vector_multiple.cfg",
        f"measurement_case_config=../measurement_rozenberg_hminus6.cfg;../measurement_koomen_meridian_hminus6_polarization.cfg;../measurement_gal_lapland_fullsky_450nm_dolp.cfg;{measurement_cfg_rel}",
        f"paper_primary_measurement_case_config={measurement_cfg_rel}",
        f"paper_primary_measurement_frozen={paper_primary_measurement_frozen}",
        f"paper_case_provenance_json=../../data/paper_cases/{case_id}/paper_case_provenance.json",
        f"profile_csv=../../data/paper_cases/{case_id}/atmosphere_profile.csv",
        f"measurement_metadata_json=../../data/paper_cases/{case_id}/measurement_metadata.json",
        "solar_spectrum_csv=../../data/optics/solar_irradiance_reference.csv",
        f"instrument_response_csv=../../data/paper_cases/{case_id}/instrument_response.csv",
        "ozone_cross_section_csv=../../data/optics/ozone_cross_section_reference.csv",
        "o2_cross_section_csv=../../data/optics/o2_cross_section_reference.csv",
        "o4_cross_section_csv=../../data/optics/o4_cross_section_reference.csv",
        "h2o_cross_section_csv=../../data/optics/h2o_cross_section_reference.csv",
        "no2_cross_section_csv=../../data/optics/no2_cross_section_reference.csv",
        f"aerosol_phase_matrix_csv=../../data/paper_cases/{case_id}/aerosol_phase_matrix.csv",
        "surface_albedo_csv=../../data/optics/lambertian_land_albedo.csv",
        "surface_model=lambertian_land",
        f"surface_parameter_csv=../../data/paper_cases/{case_id}/surface_config.json",
        "default_surface_albedo=0.15",
        "ocean_wind_speed_m_s=5.0",
        "use_explicit_solar_angles=true",
        f"solar_zenith_deg={zenith_deg:.6f}",
        f"solar_azimuth_deg={azimuth_deg:.6f}",
        "finite_solar_disk=true",
        "solar_angular_radius_deg=0.2666",
        "solar_disk_quadrature_nodes=7",
        f"observer_latitude_deg={case_spec['site_lat_deg']:.8f}",
        f"observer_longitude_deg={case_spec['site_lon_deg']:.8f}",
        f"observer_altitude_m={case_spec['observer_altitude_m']:.1f}",
    ] + common_lines
    write_text(PAPER_CONFIG_DIR / f"{case_id}.cfg", "\n".join(case_cfg_lines) + "\n")

    validation_lines = [
        "# Paper-validation entrypoint for the frozen Marseille twilight case.",
        *(
            []
            if has_measurement_reference
            else ["# This is expected to fail until the case-specific measurement reference field is extracted."]
        ),
        "case_id=paper_validation",
        "output_dir=../results",
        "benchmark_case_config=benchmark_disort_scalar.cfg;benchmark_iprt_a1_vector.cfg;benchmark_zawada_spherical_vector_single.cfg;benchmark_zawada_spherical_vector_multiple.cfg",
        f"measurement_case_config=measurement_rozenberg_hminus6.cfg;measurement_koomen_meridian_hminus6_polarization.cfg;measurement_gal_lapland_fullsky_450nm_dolp.cfg;{measurement_cfg_rel}",
        f"paper_primary_measurement_case_config={measurement_cfg_rel}",
        f"paper_primary_measurement_frozen={paper_primary_measurement_frozen}",
        f"paper_case_provenance_json=../data/paper_cases/{case_id}/paper_case_provenance.json",
        f"profile_csv=../data/paper_cases/{case_id}/atmosphere_profile.csv",
        f"measurement_metadata_json=../data/paper_cases/{case_id}/measurement_metadata.json",
        "solar_spectrum_csv=../data/optics/solar_irradiance_reference.csv",
        f"instrument_response_csv=../data/paper_cases/{case_id}/instrument_response.csv",
        "ozone_cross_section_csv=../data/optics/ozone_cross_section_reference.csv",
        "o2_cross_section_csv=../data/optics/o2_cross_section_reference.csv",
        "o4_cross_section_csv=../data/optics/o4_cross_section_reference.csv",
        "h2o_cross_section_csv=../data/optics/h2o_cross_section_reference.csv",
        "no2_cross_section_csv=../data/optics/no2_cross_section_reference.csv",
        f"aerosol_phase_matrix_csv=../data/paper_cases/{case_id}/aerosol_phase_matrix.csv",
        "surface_albedo_csv=../data/optics/lambertian_land_albedo.csv",
        "surface_model=lambertian_land",
        f"surface_parameter_csv=../data/paper_cases/{case_id}/surface_config.json",
        "default_surface_albedo=0.15",
        "ocean_wind_speed_m_s=5.0",
        "use_explicit_solar_angles=true",
        f"solar_zenith_deg={zenith_deg:.6f}",
        f"solar_azimuth_deg={azimuth_deg:.6f}",
        "finite_solar_disk=true",
        "solar_angular_radius_deg=0.2666",
        "solar_disk_quadrature_nodes=7",
        f"observer_latitude_deg={case_spec['site_lat_deg']:.8f}",
        f"observer_longitude_deg={case_spec['site_lon_deg']:.8f}",
        f"observer_altitude_m={case_spec['observer_altitude_m']:.1f}",
    ] + common_lines
    write_text(CONFIG_DIR / "paper_validation.cfg", "\n".join(validation_lines) + "\n")

    write_json(
        case_data_dir / "paper_case_provenance.json",
        {
            "case_id": case_id,
            "paper_primary_measurement_frozen": has_measurement_reference,
            "paper_gate_blocked_reason": (
                "The frozen Marseille twilight measurement reference has been extracted, but the case still has to pass the paper-validation thresholds."
                if has_measurement_reference
                else "The frozen Marseille twilight observation now has matched public atmosphere/aerosol inputs, but its machine-readable validation reference field has not been extracted yet."
            ),
            "frozen_measurement_case": {
                "dataset_doi": case_spec["dataset_doi"],
                "dataset_record_url": case_spec["dataset_record_url"],
                "measurement_paper_url": case_spec["measurement_paper_url"],
                "selected_timestamp_utc": case_spec["timestamp_utc"],
                "annotation_filename": case_spec["annotation_filename"],
                "annotation_file_id": int(annotation_file["id"]),
                "annotation_label": label,
                "solar_zenith_deg": zenith_deg,
                "solar_azimuth_deg": azimuth_deg,
            },
            "matched_public_inputs": {
                "weather_source": "Open-Meteo archive API",
                "weather_url": weather_url,
                "pressure_profile_source": "Open-Meteo historical-forecast pressure levels",
                "pressure_profile_url": pressure_profile_url,
                "air_quality_source": "Open-Meteo air-quality API",
                "air_quality_url": air_quality_url,
                "aeronet_source": "AERONET Version 3 Direct Sun",
                "aeronet_url": aeronet_url,
                "aeronet_inversion_url": aeronet_inversion_url,
                "aeronet_inversion_available": aeronet_inversion_available,
                "aeronet_inversion_usage_mode": aeronet_inversion_usage_mode,
                "aeronet_inversion_time_delta_hours": aeronet_inversion_time_delta_hours,
                "aeronet_inversion_selected_time_utc": aeronet_inversion_selected_time_utc,
                "aeronet_inversion_quality_level": aeronet_inversion_quality_level,
                "sensor_flyer_url": case_spec["sensor_flyer_url"],
                "camera_techref_url": case_spec["camera_techref_url"],
                "camera_techref_download_url": case_spec["camera_techref_download_url"],
                "aeronet_site": "Marseille_ATMO",
                "aeronet_sample_time_utc": chosen_time.isoformat().replace("+00:00", "Z"),
                "aeronet_time_delta_hours": aeronet_time_delta_hours,
                "aeronet_weight_in_case_aod": aeronet_aod_weight,
            },
            "derived_case_summary": {
                "surface_temperature_k": surface_temperature_k,
                "surface_relative_humidity": surface_relative_humidity,
                "surface_pressure_pa": surface_pressure_pa,
                "cloud_cover_percent": cloud_cover_percent,
                "surface_ozone_ug_m3": surface_ozone_ug_m3,
                "surface_no2_ug_m3": surface_no2_ug_m3,
                "surface_pm10_ug_m3": surface_pm10_ug_m3,
                "surface_pm25_ug_m3": surface_pm25_ug_m3,
                "target_aod_550": aerosol_target_aod_550,
                "openmeteo_aod_550": openmeteo_aod_550,
                "aeronet_aod_500": aeronet_aod_500,
                "aeronet_aod_440": aeronet_aod_440,
                "aeronet_aod_550": aeronet_aod_550,
                "aeronet_angstrom_440_870": angstrom,
                "aeronet_precipitable_water_cm": float(chosen_row["Precipitable_Water(cm)"]),
                "aeronet_ozone_dobson": float(chosen_row["Ozone(Dobson)"]),
                "aeronet_no2_dobson": float(chosen_row["NO2(Dobson)"]),
                "openmeteo_tcwv_kg_m2": tcwv_kg_m2,
                "thermodynamic_profile_mode": "openmeteo_pressure_levels_with_template_upper_tail",
                "pressure_level_rows": len(matched_thermo_rows),
                "upper_tail_rows": len(tail_template_candidates),
                "matched_profile_top_altitude_m": matched_profile_top_altitude_m,
                "case_local_aerosol_model": aerosol_summary,
            },
            "artifacts": {
                "atmosphere_profile_csv": {"path": str(atmosphere_profile_path), "sha256": sha256_file(atmosphere_profile_path)},
                "instrument_response_csv": {"path": str(instrument_response_path), "sha256": sha256_file(instrument_response_path)},
                "aerosol_phase_matrix_csv": {"path": str(aerosol_phase_matrix_path), "sha256": sha256_file(aerosol_phase_matrix_path)},
                "aerosol_optics_json": {"path": str(aerosol_optics_metadata_path), "sha256": sha256_file(aerosol_optics_metadata_path)},
                "surface_config_json": {"path": str(surface_config_path), "sha256": sha256_file(surface_config_path)},
                "measurement_metadata_json": {"path": str(measurement_metadata_path), "sha256": sha256_file(measurement_metadata_path)},
                "open_meteo_pressure_levels_json": {
                    "path": str(case_data_dir / "open_meteo_pressure_levels.json"),
                    "sha256": sha256_file(case_data_dir / "open_meteo_pressure_levels.json"),
                },
                "aeronet_inversion_txt": {
                    "path": str(case_data_dir / "aeronet_inversion.txt"),
                    "sha256": sha256_file(case_data_dir / "aeronet_inversion.txt"),
                },
                **(
                    {
                        "measurement_reference_csv": {
                            "path": str(measurement_reference_path),
                            "sha256": sha256_file(measurement_reference_path),
                        }
                    }
                    if has_measurement_reference
                    else {}
                ),
                **(
                    {
                        "measurement_reduction_json": {
                            "path": str(measurement_reduction_path),
                            "sha256": sha256_file(measurement_reduction_path),
                        }
                    }
                    if measurement_reduction_path.exists()
                    else {}
                ),
                **(
                    {"measurement_config": {"path": str(measurement_cfg_path), "sha256": sha256_file(measurement_cfg_path)}}
                    if has_measurement_reference
                    else {"measurement_placeholder_config": {"path": str(measurement_cfg_path), "sha256": sha256_file(measurement_cfg_path)}}
                ),
            },
            "remaining_approximations": [
                "The lower and middle atmosphere thermodynamic profile now comes from matched Open-Meteo pressure levels, but the upper-atmosphere tail above the highest pressure level still uses the checked-in template profile.",
                "The Marseille aerosol phase matrix now uses AERONET inversion microphysics when available, but the vertical distribution is still a shaped profile rather than a measured extinction profile.",
                "The instrument response is now a case-local public-doc-constrained IMX250MYR blue-channel proxy, but it is still approximate because no tabulated manufacturer spectral-response curve is bundled with the dataset.",
            ],
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a frozen paper-case package from matched public inputs.")
    parser.add_argument(
        "--case-id",
        default="frozen_marseille_twilight_20220815_191413z",
        choices=sorted(FROZEN_CASES.keys()),
    )
    args = parser.parse_args()
    build_case(FROZEN_CASES[args.case_id])
    print(f"paper case inputs written for {args.case_id}")


if __name__ == "__main__":
    main()
