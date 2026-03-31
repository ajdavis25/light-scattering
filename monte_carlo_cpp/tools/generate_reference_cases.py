from __future__ import annotations

import csv
import json
from pathlib import Path
import ssl
import urllib.request

import numpy as np
from PythonicDISORT import pydisort
from PythonicDISORT.subroutines import interpolate


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_DIR = ROOT / "monte_carlo_cpp" / "data" / "reference_cases"
ATMOSPHERE_DIR = ROOT / "monte_carlo_cpp" / "data" / "atmosphere"
OPTICS_DIR = ROOT / "monte_carlo_cpp" / "data" / "optics"


def fetch_url_bytes(url: str, *, allow_insecure_ssl: bool = False) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return response.read()
    except Exception:
        if not allow_insecure_ssl:
            raise
        context = ssl._create_unverified_context()
        with urllib.request.urlopen(request, timeout=60, context=context) as response:
            return response.read()


def write_benchmark_case() -> None:
    tau_arr = np.array([0.05], dtype=float)
    omega_arr = np.array([0.999999], dtype=float)
    legendre = np.zeros((1, 32), dtype=float)
    legendre[0, 0] = 1.0
    # PythonicDISORT expects unweighted Legendre coefficients beta_l.
    # For Rayleigh, p(mu) = (1 / 4pi) * [1 + 0.5 P2(mu)], so beta_2 = 0.1
    # because the code applies the (2l + 1) factor internally.
    legendre[0, 2] = 0.1
    mu0 = float(np.cos(np.deg2rad(60.0)))
    _, _, _, _, u = pydisort(
        tau_arr,
        omega_arr,
        32,
        legendre,
        mu0,
        1.0,
        0.0,
        NLeg=32,
        NFourier=16,
    )
    u_interp = interpolate(u)

    zeniths = [5, 15, 25, 35, 45, 55, 65, 75]
    rows = []
    for zenith in zeniths:
        mu = -float(np.cos(np.deg2rad(zenith)))
        intensity = float(u_interp(mu, tau_arr[-1], 0.0))
        rows.append(
            {
                "zenith_deg": zenith,
                "relative_azimuth_deg": 0.0,
                "intensity": intensity,
            }
        )

    csv_path = REFERENCE_DIR / "benchmark_disort_scalar_principal_plane.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["zenith_deg", "relative_azimuth_deg", "intensity"])
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "source_name": "PythonicDISORT 1.6",
        "source_url": "https://pythonic-disort.readthedocs.io/en/latest/Pythonic-DISORT.html",
        "description": "Plane-parallel scalar Rayleigh principal-plane reference generated with PythonicDISORT.",
        "geometry": {
            "solar_zenith_deg": 60.0,
            "relative_azimuth_deg": 0.0,
            "observer_boundary": "surface",
        },
        "optics": {
            "single_layer_optical_depth": 0.05,
            "single_scattering_albedo": 0.999999,
            "rayleigh_legendre_coefficients": [1.0, 0.0, 0.1],
        },
    }
    (REFERENCE_DIR / "benchmark_disort_scalar_metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )


def write_measurement_case() -> None:
    rows = []
    table = {
        0.0: [0.900, 0.700, 0.180, 0.060, 0.033, 0.022],
        22.5: [0.680, 0.670, 0.170, 0.060, 0.033, 0.022],
        45.0: [0.400, 0.350, 0.100, 0.046, 0.031, 0.022],
        90.0: [0.170, 0.155, 0.060, 0.038, 0.028, 0.022],
        180.0: [0.085, 0.085, 0.056, 0.037, 0.028, 0.022],
    }
    altitudes = [0.0, 10.0, 30.0, 50.0, 70.0, 90.0]
    for relative_azimuth, values in table.items():
        for altitude_deg, intensity in zip(altitudes, values):
            rows.append(
                {
                    "altitude_deg": altitude_deg,
                    "relative_azimuth_deg": relative_azimuth,
                    "intensity": intensity,
                }
            )

    csv_path = REFERENCE_DIR / "measurement_rozenberg_table1_hminus6.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["altitude_deg", "relative_azimuth_deg", "intensity"])
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "source_name": "Rozenberg (1952) Table I, Sacramento Peak",
        "source_url": "https://ntrs.nasa.gov/api/citations/19930092270/downloads/19930092270.pdf",
        "description": "Measured twilight sky brightness pattern for solar altitude H=-6 degrees from Table I.",
        "geometry": {
            "solar_altitude_deg": -6.0,
            "solar_zenith_deg": 96.0,
            "relative_azimuth_definition": "bearing from the direction of the Sun",
        },
        "notes": [
            "Brightness values are taken directly from the published table.",
            "Validation uses normalized intensity pattern rather than absolute photometric units.",
        ],
    }
    (REFERENCE_DIR / "measurement_rozenberg_table1_metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )


def write_iprt_vector_benchmark_case() -> None:
    url = (
        "https://www.meteo.physik.uni-muenchen.de/~iprt/lib/exe/fetch.php?"
        "media=intercomparisons:phase_a:a1:iprt_case_a1_mystic.dat"
    )
    raw_text = urllib.request.urlopen(url).read().decode("utf-8", errors="ignore")

    target_vza = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0}
    rows = []
    for line in raw_text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) < 10:
            continue
        depol, altitude, sza, saa, vza, vaa = map(float, fields[:6])
        if depol != 0.0 or altitude != 0.0 or sza != 0.0 or saa != 65.0 or vaa != 0.0 or vza not in target_vza:
            continue
        intensity, q, u, _v = map(float, fields[6:10])
        rows.append(
            {
                "zenith_deg": vza,
                "azimuth_deg": vaa,
                "intensity": intensity,
                "q": q,
                "u": u,
            }
        )

    rows.sort(key=lambda row: row["zenith_deg"])
    csv_path = REFERENCE_DIR / "benchmark_iprt_a1_vector_mystic.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["zenith_deg", "azimuth_deg", "intensity", "q", "u"])
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "source_name": "IPRT Phase A1 benchmark, MYSTIC reference subset",
        "source_url": "https://www.meteo.physik.uni-muenchen.de/~iprt/doku.php?id=intercomparisons:a1_rayleigh",
        "reference_file_url": url,
        "description": (
            "Vector Rayleigh benchmark subset from the IPRT A1 intercomparison. "
            "Subset uses depol=0, altitude=0 km, sza=0 deg, saa=65 deg, vaa=0 deg, vza=10..80 deg."
        ),
        "notes": [
            "This benchmark includes I, Q, and U from an external vector radiative transfer model.",
            "The direct-sun viewing direction vza=0 is excluded because it is not a stable diffuse-sky comparison point.",
            "The subset is intentionally small to keep regression validation runtime bounded.",
        ],
    }
    (REFERENCE_DIR / "benchmark_iprt_a1_vector_metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )


def write_koomen_polarization_measurement_case() -> None:
    west_intensity = {
        10.0: 0.700,
        30.0: 0.180,
        50.0: 0.060,
        70.0: 0.033,
        90.0: 0.022,
    }
    # Digitized from Fig. 3 in Koomen et al. (1952) for H=-6 on the meridian through the Sun.
    # The paper defines p = B_parallel / B_perpendicular. Convert to DoLP via (1 - p) / (1 + p).
    polarization_factor = {
        10.0: 0.93,
        30.0: 0.82,
        50.0: 0.66,
        70.0: 0.24,
        90.0: 0.20,
    }

    rows = []
    for altitude_deg, intensity in west_intensity.items():
        p_factor = polarization_factor[altitude_deg]
        dop = (1.0 - p_factor) / (1.0 + p_factor)
        rows.append(
            {
                "altitude_deg": altitude_deg,
                "relative_azimuth_deg": 0.0,
                "intensity": intensity,
                "dop": dop,
            }
        )

    csv_path = REFERENCE_DIR / "measurement_koomen_meridian_hminus6_polarization.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["altitude_deg", "relative_azimuth_deg", "intensity", "dop"])
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "source_name": "Koomen et al. (1952) Fig. 3 plus Table I, Sacramento Peak",
        "source_url": "https://static1.squarespace.com/static/54694fa6e4b0eaec4530f99d/t/5e876886cdfbaa441f538327/1585932427823/measurements%2Bof%2Bthe%2Bbrightness%2Bof%2Bthe%2Btwilight%2Bsky%2B1952.pdf",
        "description": (
            "Solar-meridian twilight brightness and polarization reference for H=-6 deg at Sacramento Peak. "
            "Intensity values come from Table I (Z=0); polarization factors are digitized from Fig. 3 and converted to DoLP."
        ),
        "notes": [
            "This is a meridian-only polarization validation subset, not a full-sky measurement field.",
            "The DoLP values are approximate because the publication provides the H=-6 polarization data as a plotted curve, not a table.",
        ],
    }
    (REFERENCE_DIR / "measurement_koomen_meridian_hminus6_metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )


def write_gal_fullsky_polarization_measurement_case() -> None:
    import fitz
    from PIL import Image

    pdf_bytes = fetch_url_bytes(
        "https://arago.elte.hu/sites/default/files/SkyPolLapland_PRSLA.pdf",
        allow_insecure_ssl=True,
    )
    document = fitz.open(stream=pdf_bytes, filetype="pdf")
    page = document.load_page(5)
    pixmap = page.get_pixmap(matrix=fitz.Matrix(3.0, 3.0), alpha=False)
    image = Image.frombytes("RGB", (pixmap.width, pixmap.height), pixmap.samples)

    # Figure 2a, first-row degree-of-polarization panel at 2 h.
    panel = image.crop((360, 234, 500, 361)).convert("L")
    gray = np.array(panel, dtype=float)

    center_x = 71.25
    center_y = 60.5
    radius_px = 58.5
    sun_azimuth_deg = 22.18
    patch_half_width = 2

    rows = []
    for zenith_deg in (15.0, 30.0, 45.0, 60.0, 75.0):
        radial_distance = radius_px * zenith_deg / 90.0
        for relative_azimuth_deg in range(0, 360, 30):
            absolute_azimuth_deg = (sun_azimuth_deg + relative_azimuth_deg) % 360.0
            azimuth_rad = np.deg2rad(absolute_azimuth_deg)
            x = center_x - radial_distance * np.sin(azimuth_rad)
            y = center_y - radial_distance * np.cos(azimuth_rad)

            x_index = int(round(x))
            y_index = int(round(y))
            patch = gray[
                max(0, y_index - patch_half_width):min(gray.shape[0], y_index + patch_half_width + 1),
                max(0, x_index - patch_half_width):min(gray.shape[1], x_index + patch_half_width + 1),
            ]
            if patch.size == 0:
                continue

            dop = 1.0 - float(np.median(patch)) / 255.0
            rows.append(
                {
                    "zenith_deg": zenith_deg,
                    "relative_azimuth_deg": float(relative_azimuth_deg),
                    "dop": max(0.0, min(1.0, dop)),
                }
            )

    csv_path = REFERENCE_DIR / "measurement_gal_lapland_fullsky_450nm_dolp.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["zenith_deg", "relative_azimuth_deg", "dop"])
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "source_name": "Gal et al. (2001) Figure 2a, Sodankyla full-sky imaging polarimetry",
        "source_url": "https://arago.elte.hu/sites/default/files/SkyPolLapland_PRSLA.pdf",
        "doi": "10.1098/rspa.2000.0726",
        "description": (
            "Digitized full-sky degree-of-linear-polarization reference from the 2 h panel "
            "of Figure 2a in Gal et al. (2001), measured at 450 nm with calibrated full-sky imaging polarimetry."
        ),
        "geometry": {
            "solar_zenith_deg": 83.1,
            "relative_azimuth_definition": "bearing from the solar direction in the published fisheye map",
            "sampling_grid": {
                "zenith_deg": [15, 30, 45, 60, 75],
                "relative_azimuth_deg": list(range(0, 360, 30)),
            },
        },
        "notes": [
            "This is a digitized coarse full-sky DoLP field derived from a published calibrated fisheye map, not raw instrument output.",
            "The published circular map was sampled assuming approximately linear zenith radius in the rendered figure.",
            "The reference is low-sun daytime sky at solar zenith 83.1 deg, not below-horizon twilight.",
            "Only polarization metrics are enforced for this case because the source figure does not provide a matching calibrated intensity table.",
        ],
    }
    (REFERENCE_DIR / "measurement_gal_lapland_fullsky_450nm_metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )


def write_zawada_spherical_vector_benchmark_case() -> None:
    from netCDF4 import Dataset

    dataset_url = "https://zenodo.org/api/records/4292303/files/zawada_AMT_rtm_comparison_data_v1.nc/content"
    dataset_path = REFERENCE_DIR / "zawada_amt_rtm_comparison_data_v1.nc"
    if not dataset_path.exists():
        dataset_path.write_bytes(fetch_url_bytes(dataset_url))

    with Dataset(dataset_path) as dataset:
        model_data = dataset.groups["model_data"]
        ancillary_data = dataset.groups["ancillary_data"]
        geometry_data = dataset.groups["geometry_data"]

        wavelengths_nm = np.array(model_data.variables["wavelength"][:], dtype=float)
        altitude_km = np.array(model_data.variables["altitude"][:], dtype=float)
        test_cases = np.array(model_data.variables["test_case"][:], dtype=int)
        compositions = np.array(model_data.variables["composition"][:], dtype=str)

        benchmark_wavelength_nm = 351.0
        benchmark_test_case = 0
        benchmark_solar_index = 0
        benchmark_composition = "rayleigh+ozone+aerosol"
        benchmark_albedo_index = 0
        altitude_indices = list(range(0, altitude_km.size, 5))
        if altitude_indices[-1] != altitude_km.size - 1:
            altitude_indices.append(altitude_km.size - 1)

        wavelength_index = int(np.where(np.isclose(wavelengths_nm, benchmark_wavelength_nm))[0][0])
        test_case_index = int(np.where(test_cases == benchmark_test_case)[0][0])
        composition_index = int(np.where(compositions == benchmark_composition)[0][0])

        toa_mu = np.array(geometry_data.variables["toa_mu"][benchmark_solar_index, :], dtype=float)
        toa_sza = np.array(geometry_data.variables["toa_sza"][benchmark_solar_index, :], dtype=float)
        toa_saa = np.array(geometry_data.variables["toa_saa"][benchmark_solar_index, :], dtype=float)
        benchmark_radiance = np.array(
            model_data.variables["mmm"][
                test_case_index,
                benchmark_solar_index,
                composition_index,
                benchmark_albedo_index,
                wavelength_index,
                :,
                :3,
            ],
            dtype=float,
        )

        rows = []
        for altitude_index in altitude_indices:
            viewing_nadir_deg = float(np.degrees(np.arccos(np.clip(toa_mu[altitude_index], -1.0, 1.0))))
            rows.append(
                {
                    "tangent_altitude_km": float(altitude_km[altitude_index]),
                    "zenith_deg": 180.0 - viewing_nadir_deg,
                    "azimuth_deg": 0.0,
                    "solar_zenith_deg": float(toa_sza[altitude_index]),
                    "solar_azimuth_deg": float(toa_saa[altitude_index]),
                    "intensity": float(benchmark_radiance[altitude_index, 0]),
                    "q": float(benchmark_radiance[altitude_index, 1]),
                    "u": float(benchmark_radiance[altitude_index, 2]),
                }
            )

        benchmark_csv = REFERENCE_DIR / "benchmark_zawada_spherical_vector_single.csv"
        with benchmark_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "tangent_altitude_km",
                    "zenith_deg",
                    "azimuth_deg",
                    "solar_zenith_deg",
                    "solar_azimuth_deg",
                    "intensity",
                    "q",
                    "u",
                ],
            )
            writer.writeheader()
            writer.writerows(rows)

        ancillary_altitude_m = np.array(ancillary_data.variables["altitude"][:], dtype=float)
        pressure_pa = np.array(ancillary_data.variables["pressure"][:], dtype=float)
        temperature_k = np.array(ancillary_data.variables["temperature"][:], dtype=float)
        air_numden_cm3 = np.array(ancillary_data.variables["air_numden"][:], dtype=float)
        ozone_numden_cm3 = np.array(ancillary_data.variables["ozone_numden"][:], dtype=float)
        aerosol_numden_cm3 = np.array(ancillary_data.variables["aerosol_numden"][:], dtype=float)
        aerosol_scat_cross_section_cm2 = float(
            ancillary_data.variables["aerosol_scattering_cross_section"][wavelength_index]
        )
        aerosol_abs_cross_section_cm2 = float(
            ancillary_data.variables["aerosol_absorption_cross_section"][wavelength_index]
        )
        aerosol_extinction_m_inv = (
            aerosol_numden_cm3 * (aerosol_scat_cross_section_cm2 + aerosol_abs_cross_section_cm2) * 1.0e2
        )
        aerosol_ssa = np.divide(
            aerosol_numden_cm3 * aerosol_scat_cross_section_cm2,
            aerosol_numden_cm3 * (aerosol_scat_cross_section_cm2 + aerosol_abs_cross_section_cm2),
            out=np.zeros_like(aerosol_numden_cm3),
            where=(aerosol_scat_cross_section_cm2 + aerosol_abs_cross_section_cm2) > 0.0,
        )

        phase_mu = np.array(ancillary_data.variables["cos_scatter_angle"][:], dtype=float)
        phase_matrix = np.array(
            ancillary_data.variables["aerosol_phase_matrix"][wavelength_index, :, :, :],
            dtype=float,
        )
        phase_f11 = phase_matrix[:, 0, 0]
        phase_angles_deg = np.degrees(np.arccos(np.clip(phase_mu, -1.0, 1.0)))
        normalization = np.trapezoid(phase_f11, phase_mu)
        aerosol_asymmetry = float(np.trapezoid(phase_f11 * phase_mu, phase_mu) / max(1.0e-12, normalization))

        atmosphere_csv = ATMOSPHERE_DIR / "benchmark_zawada_limb_351nm.csv"
        with atmosphere_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
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
                ],
            )
            writer.writeheader()
            for index, altitude_m in enumerate(ancillary_altitude_m):
                writer.writerow(
                    {
                        "altitude_m": float(altitude_m),
                        "pressure_pa": float(pressure_pa[index]),
                        "temperature_k": float(temperature_k[index]),
                        "molecular_number_density_m3": float(air_numden_cm3[index] * 1.0e6),
                        "ozone_number_density_m3": float(ozone_numden_cm3[index] * 1.0e6),
                        "aerosol_extinction_550_m_inv": float(aerosol_extinction_m_inv[index]),
                        "aerosol_single_scattering_albedo": float(aerosol_ssa[index]),
                        "aerosol_asymmetry": aerosol_asymmetry,
                        "aerosol_scattering_angstrom_exponent": 0.0,
                        "aerosol_absorption_angstrom_exponent": 0.0,
                    }
                )

        phase_csv = OPTICS_DIR / "aerosol_phase_matrix_zawada_351nm.csv"
        with phase_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["wavelength_nm", "angle_deg", "f11", "f12", "f22", "f33", "f34", "f44"],
            )
            writer.writeheader()
            for angle_deg, coeffs in zip(phase_angles_deg, phase_matrix, strict=True):
                writer.writerow(
                    {
                        "wavelength_nm": benchmark_wavelength_nm,
                        "angle_deg": float(angle_deg),
                        "f11": float(coeffs[0, 0]),
                        "f12": float(coeffs[0, 1]),
                        "f22": float(coeffs[1, 1]),
                        "f33": float(coeffs[2, 2]),
                        "f34": float(coeffs[2, 3]),
                        "f44": float(coeffs[3, 3]),
                    }
                )

        solar_csv = OPTICS_DIR / "solar_unit_351nm.csv"
        with solar_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["wavelength_nm", "irradiance_w_m2_nm"])
            writer.writeheader()
            writer.writerow({"wavelength_nm": benchmark_wavelength_nm, "irradiance_w_m2_nm": 1.0})

        ozone_csv = OPTICS_DIR / "ozone_cross_section_zawada_351nm.csv"
        with ozone_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["wavelength_nm", "cross_section_m2"])
            writer.writeheader()
            writer.writerow(
                {
                    "wavelength_nm": benchmark_wavelength_nm,
                    "cross_section_m2": float(ancillary_data.variables["ozone_absorption_cross_section"][wavelength_index]) * 1.0e-4,
                }
            )

        rayleigh_csv = OPTICS_DIR / "rayleigh_cross_section_zawada_351nm.csv"
        with rayleigh_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["wavelength_nm", "cross_section_m2"])
            writer.writeheader()
            writer.writerow(
                {
                    "wavelength_nm": benchmark_wavelength_nm,
                    "cross_section_m2": float(ancillary_data.variables["rayleigh_scattering_cross_section"][wavelength_index]) * 1.0e-4,
                }
            )

        metadata = {
            "source_name": "Zawada et al. (2020) spherical-vector limb benchmark subset",
            "source_url": "https://doi.org/10.5281/zenodo.4292303",
            "paper_url": "https://doi.org/10.5194/amt-14-3953-2021",
            "reference_file_url": dataset_url,
            "description": (
                "Published spherical vector benchmark subset extracted from the Zawada et al. limb-scattering "
                "intercomparison dataset. The subset uses the multi-model mean for the single-scatter, no-refraction "
                "test case with Rayleigh scattering, ozone absorption, and stratospheric aerosol scattering at 351 nm."
            ),
            "subset_definition": {
                "test_case": "single scatter only, no refraction",
                "solar_condition_index": benchmark_solar_index,
                "tangent_sza_deg": float(geometry_data.variables["tangent_sza"][benchmark_solar_index]),
                "tangent_saa_deg": float(geometry_data.variables["tangent_saa"][benchmark_solar_index]),
                "composition": benchmark_composition,
                "surface_albedo": float(model_data.variables["albedo"][benchmark_albedo_index]),
                "wavelength_nm": benchmark_wavelength_nm,
                "tangent_altitude_sampling_km": [float(altitude_km[index]) for index in altitude_indices],
            },
            "notes": [
                "Reference radiances are the published multi-model mean from the Zenodo benchmark dataset, not an internally generated surrogate.",
                "Each sampled line of sight carries its own observer-local solar zenith and azimuth from the published geometry_data group.",
                "The benchmark atmosphere and aerosol phase table are derived directly from the ancillary_data group at 351 nm to avoid cross-study optics drift.",
                "The benchmark-specific Rayleigh scattering cross section is also taken directly from the ancillary_data group at 351 nm.",
            ],
        }
        (REFERENCE_DIR / "benchmark_zawada_spherical_vector_single_metadata.json").write_text(
            json.dumps(metadata, indent=2),
            encoding="utf-8",
        )

        benchmark_test_case = 1
        test_case_index = int(np.where(test_cases == benchmark_test_case)[0][0])
        benchmark_radiance = np.array(
            model_data.variables["mmm"][
                test_case_index,
                benchmark_solar_index,
                composition_index,
                benchmark_albedo_index,
                wavelength_index,
                :,
                :3,
            ],
            dtype=float,
        )

        rows = []
        for altitude_index in altitude_indices:
            viewing_nadir_deg = float(np.degrees(np.arccos(np.clip(toa_mu[altitude_index], -1.0, 1.0))))
            rows.append(
                {
                    "tangent_altitude_km": float(altitude_km[altitude_index]),
                    "zenith_deg": 180.0 - viewing_nadir_deg,
                    "azimuth_deg": 0.0,
                    "solar_zenith_deg": float(toa_sza[altitude_index]),
                    "solar_azimuth_deg": float(toa_saa[altitude_index]),
                    "intensity": float(benchmark_radiance[altitude_index, 0]),
                    "q": float(benchmark_radiance[altitude_index, 1]),
                    "u": float(benchmark_radiance[altitude_index, 2]),
                }
            )

        benchmark_csv = REFERENCE_DIR / "benchmark_zawada_spherical_vector_multiple.csv"
        with benchmark_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "tangent_altitude_km",
                    "zenith_deg",
                    "azimuth_deg",
                    "solar_zenith_deg",
                    "solar_azimuth_deg",
                    "intensity",
                    "q",
                    "u",
                ],
            )
            writer.writeheader()
            writer.writerows(rows)

        metadata = {
            "source_name": "Zawada et al. (2020) spherical-vector limb benchmark subset",
            "source_url": "https://doi.org/10.5281/zenodo.4292303",
            "paper_url": "https://doi.org/10.5194/amt-14-3953-2021",
            "reference_file_url": dataset_url,
            "description": (
                "Published spherical vector benchmark subset extracted from the Zawada et al. limb-scattering "
                "intercomparison dataset. The subset uses the multi-model mean for the all-orders, no-refraction "
                "test case with Rayleigh scattering, ozone absorption, and stratospheric aerosol scattering at 351 nm."
            ),
            "subset_definition": {
                "test_case": "all orders of scatter, no refraction",
                "solar_condition_index": benchmark_solar_index,
                "tangent_sza_deg": float(geometry_data.variables["tangent_sza"][benchmark_solar_index]),
                "tangent_saa_deg": float(geometry_data.variables["tangent_saa"][benchmark_solar_index]),
                "composition": benchmark_composition,
                "surface_albedo": float(model_data.variables["albedo"][benchmark_albedo_index]),
                "wavelength_nm": benchmark_wavelength_nm,
                "tangent_altitude_sampling_km": [float(altitude_km[index]) for index in altitude_indices],
            },
            "notes": [
                "Reference radiances are the published multi-model mean from the Zenodo benchmark dataset, not an internally generated surrogate.",
                "Each sampled line of sight carries its own observer-local solar zenith and azimuth from the published geometry_data group.",
                "This case is stronger than the single-scatter spherical benchmark because it includes higher-order multiple scattering while still remaining within the current no-refraction solver assumptions.",
                "The benchmark-specific Rayleigh scattering cross section is also taken directly from the ancillary_data group at 351 nm.",
            ],
        }
        (REFERENCE_DIR / "benchmark_zawada_spherical_vector_multiple_metadata.json").write_text(
            json.dumps(metadata, indent=2),
            encoding="utf-8",
        )


def write_reference_manifest() -> None:
    manifest = {
        "benchmark_references": [
            {
                "status": "bundled",
                "case_config": "benchmark_disort_scalar.cfg",
                "reference_csv": "benchmark_disort_scalar_principal_plane.csv",
                "metadata_json": "benchmark_disort_scalar_metadata.json",
                "notes": "External plane-parallel scalar Rayleigh benchmark generated with PythonicDISORT 1.6.",
            },
            {
                "status": "bundled",
                "case_config": "benchmark_iprt_a1_vector.cfg",
                "reference_csv": "benchmark_iprt_a1_vector_mystic.csv",
                "metadata_json": "benchmark_iprt_a1_vector_metadata.json",
                "notes": "External vector Rayleigh benchmark subset from the IPRT A1 intercomparison using MYSTIC results.",
            },
            {
                "status": "bundled",
                "case_config": "benchmark_zawada_spherical_vector_single.cfg",
                "reference_csv": "benchmark_zawada_spherical_vector_single.csv",
                "metadata_json": "benchmark_zawada_spherical_vector_single_metadata.json",
                "notes": "Published spherical-vector limb-scattering benchmark subset from Zawada et al. (2020) / AMT 2021 using the Zenodo multi-model mean at 351 nm.",
            },
            {
                "status": "bundled",
                "case_config": "benchmark_zawada_spherical_vector_multiple.cfg",
                "reference_csv": "benchmark_zawada_spherical_vector_multiple.csv",
                "metadata_json": "benchmark_zawada_spherical_vector_multiple_metadata.json",
                "notes": "Published spherical-vector multiple-scattering limb benchmark subset from Zawada et al. (2020) / AMT 2021 using the Zenodo multi-model mean at 351 nm.",
            },
        ],
        "measurement_references": [
            {
                "status": "bundled",
                "case_config": "measurement_rozenberg_hminus6.cfg",
                "reference_csv": "measurement_rozenberg_table1_hminus6.csv",
                "metadata_json": "measurement_rozenberg_table1_metadata.json",
                "notes": "Published twilight brightness pattern from Rozenberg (1952) Table I, Sacramento Peak, H=-6 degrees.",
            },
            {
                "status": "bundled",
                "case_config": "measurement_koomen_meridian_hminus6_polarization.cfg",
                "reference_csv": "measurement_koomen_meridian_hminus6_polarization.csv",
                "metadata_json": "measurement_koomen_meridian_hminus6_metadata.json",
                "notes": "Published twilight meridian polarization subset from Koomen et al. (1952) Fig. 3 plus Table I, H=-6 degrees.",
            },
            {
                "status": "bundled",
                "case_config": "measurement_gal_lapland_fullsky_450nm_dolp.cfg",
                "reference_csv": "measurement_gal_lapland_fullsky_450nm_dolp.csv",
                "metadata_json": "measurement_gal_lapland_fullsky_450nm_metadata.json",
                "notes": "Digitized full-sky calibrated DoLP field from Gal et al. (2001) Figure 2a, 2 h panel, 450 nm, solar zenith 83.1 degrees.",
            },
        ],
    }
    (REFERENCE_DIR / "reference_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )


def main() -> None:
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    ATMOSPHERE_DIR.mkdir(parents=True, exist_ok=True)
    OPTICS_DIR.mkdir(parents=True, exist_ok=True)
    write_benchmark_case()
    write_iprt_vector_benchmark_case()
    write_zawada_spherical_vector_benchmark_case()
    write_measurement_case()
    write_koomen_polarization_measurement_case()
    write_gal_fullsky_polarization_measurement_case()
    write_reference_manifest()


if __name__ == "__main__":
    main()
