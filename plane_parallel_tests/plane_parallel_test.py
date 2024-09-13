import pytest, numpy as np
from plane_parallel.plane_parallel import (
    rayleigh_phase_function,
    scattering_angle,
    photon_unit_vector,
    scattered_intensity_reflected,
    scattered_intensity_transmitted,
    intensity_at_ground,
    intensity_at_top,
    intensity_at_tau,
)


# test data
@pytest.mark.parametrize("psi, expected", [
    (0, 3.0 / (8.0 * np.pi)),
    (90, 3.0 / (16.0 * np.pi)),
    (180, 3.0 / (8.0 * np.pi)),
])
def test_rayleigh_phase_function(psi, expected):
    result = rayleigh_phase_function(psi)
    assert np.isclose(result, expected, atol=1e-6), f"Expected {expected}, got {result}"


def test_scattering_angle():
    vector_1 = np.array([1, 0, 0])
    vector_2 = np.array([0, 1, 0])
    result = scattering_angle(vector_1, vector_2)
    expected = 90.0
    assert np.isclose(result, expected, atol=1e-6), f"Expected {expected}, got {result}"

    vector_3 = np.array([1, 0, 0])
    vector_4 = np.array([1, 0, 0])
    result = scattering_angle(vector_3, vector_4)
    expected = 0.0
    assert np.isclose(result, expected, atol=1e-6), f"Expected {expected}, got {result}"


@pytest.mark.parametrize("theta, phi, expected", [
    (0, 0, np.array([0, 0, 1])),
    (90, 0, np.array([1, 0, 0])),
    (90, 90, np.array([0, 1, 0])),
])
def test_photon_unit_vector(theta, phi, expected):
    result = photon_unit_vector(theta, phi)
    assert np.allclose(result, expected, atol=1e-6), f"Expected {expected}, got {result}"


def test_scattered_intensity_reflected():
    result = scattered_intensity_reflected(45, 45, 1.0)
    assert result > 0, f"Expected a positive scattered intensity, got {result}"

    result = scattered_intensity_reflected(0, 90, 1.0)
    assert result > 0, f"Expected a positive scattered intensity, got {result}"


def test_scattered_intensity_transmitted():
    result = scattered_intensity_transmitted(45, 45, 1.0)
    assert result > 0, f"Expected a positive scattered intensity, got {result}"

    result = scattered_intensity_transmitted(0, 90, 1.0)
    assert result > 0, f"Expected a positive scattered intensity, got {result}"


def test_intensity_at_ground():
    result = intensity_at_ground(45, 30, 0, 1.0)
    assert result > 0, f"Expected a positive intensity, got {result}"

    result = intensity_at_ground(0, 90, 180, 1.0)
    assert result > 0, f"Expected a positive intensity, got {result}"


def test_intensity_at_top():
    result = intensity_at_top(45, 30, 0, 1.0)
    assert result > 0, f"Expected a positive intensity, got {result}"

    result = intensity_at_top(0, 90, 180, 1.0)
    assert result > 0, f"Expected a positive intensity, got {result}"


def test_intensity_at_tau():
    result = intensity_at_tau(45, 30, 0, 1.0, 1.0)
    assert isinstance(result, np.ndarray), "Expected result to be a NumPy array."
    assert result[0] > 0, f"Expected the first component of the intensity to be positive, got {result[0]}"

    result = intensity_at_tau(0, 90, 180, 1.0, 1.0)
    assert isinstance(result, np.ndarray), "Expected result to be a NumPy array."
    assert result[0] > 0, f"Expected the first component of the intensity to be positive, got {result[0]}"
