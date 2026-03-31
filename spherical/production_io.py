from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class ProductionSkyResult:
    case_dir: Path
    zenith_grid_deg: np.ndarray
    azimuth_grid_deg: np.ndarray
    intensity: np.ndarray
    q: np.ndarray
    u: np.ndarray
    v: np.ndarray
    var_i: np.ndarray
    var_q: np.ndarray
    var_u: np.ndarray
    var_v: np.ndarray
    metadata: dict

    @property
    def degree_of_polarization(self) -> np.ndarray:
        polarized = np.sqrt(self.q**2 + self.u**2)
        safe_i = np.where(self.intensity > 0.0, self.intensity, np.nan)
        return np.nan_to_num(polarized / safe_i, nan=0.0, posinf=0.0, neginf=0.0)

    @property
    def angle_of_polarization_rad(self) -> np.ndarray:
        return 0.5 * np.arctan2(self.u, self.q)


def load_production_result(case_dir: str | Path) -> ProductionSkyResult:
    case_path = Path(case_dir)
    csv_path = case_path / "sky_result.csv"
    metadata_path = case_path / "sky_result_metadata.json"
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing production sky result CSV: {csv_path}")
    if not metadata_path.exists():
        raise FileNotFoundError(f"Missing production metadata JSON: {metadata_path}")

    metadata = json.loads(metadata_path.read_text())
    rows = list(csv.DictReader(csv_path.open(newline="")))
    zenith_values = sorted({float(row["zenith_deg"]) for row in rows})
    azimuth_values = sorted({float(row["azimuth_deg"]) for row in rows})
    zenith_index = {value: index for index, value in enumerate(zenith_values)}
    azimuth_index = {value: index for index, value in enumerate(azimuth_values)}

    shape = (len(zenith_values), len(azimuth_values))
    intensity = np.zeros(shape, dtype=float)
    q = np.zeros(shape, dtype=float)
    u = np.zeros(shape, dtype=float)
    v = np.zeros(shape, dtype=float)
    var_i = np.zeros(shape, dtype=float)
    var_q = np.zeros(shape, dtype=float)
    var_u = np.zeros(shape, dtype=float)
    var_v = np.zeros(shape, dtype=float)

    for row in rows:
        i_zen = zenith_index[float(row["zenith_deg"])]
        i_azi = azimuth_index[float(row["azimuth_deg"])]
        intensity[i_zen, i_azi] = float(row["I"])
        q[i_zen, i_azi] = float(row["Q"])
        u[i_zen, i_azi] = float(row["U"])
        v[i_zen, i_azi] = float(row["V"])
        var_i[i_zen, i_azi] = float(row["var_I"])
        var_q[i_zen, i_azi] = float(row["var_Q"])
        var_u[i_zen, i_azi] = float(row["var_U"])
        var_v[i_zen, i_azi] = float(row["var_V"])

    return ProductionSkyResult(
        case_dir=case_path,
        zenith_grid_deg=np.asarray(zenith_values, dtype=float),
        azimuth_grid_deg=np.asarray(azimuth_values, dtype=float),
        intensity=intensity,
        q=q,
        u=u,
        v=v,
        var_i=var_i,
        var_q=var_q,
        var_u=var_u,
        var_v=var_v,
        metadata=metadata,
    )


def write_result_bundle_npz(result: ProductionSkyResult, output_path: str | Path) -> Path:
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        target,
        zenith_grid_deg=result.zenith_grid_deg,
        azimuth_grid_deg=result.azimuth_grid_deg,
        intensity=result.intensity,
        q=result.q,
        u=result.u,
        v=result.v,
        var_i=result.var_i,
        var_q=result.var_q,
        var_u=result.var_u,
        var_v=result.var_v,
        dop=result.degree_of_polarization,
        aop_rad=result.angle_of_polarization_rad,
        sun_zenith_deg=float(result.metadata["sun_zenith_deg"]),
        sun_azimuth_deg=float(result.metadata["sun_azimuth_deg"]),
        photons_per_bin=int(result.metadata["photons_per_bin"]),
    )
    return target
