import pytest, numpy as np
from src.vector_operations import (
    direction_vector,
    dot_product,
    normalize,
    rayleigh_matrix,
    rotation_matrix,
    rotation_angles
)


# test for direction_vector
@pytest.mark.parametrize("theta, phi, expected", [
    (0, 0, [0, 0, 1]), # along z-axis
    (90, 0, [1, 0, 0]), # along x-axis
    (90, 90, [0, 1, 0]), # along y-axis
    (45, 45, [0.5, 0.5, np.sqrt(2)/2]), # diagonal direction
])
def test_direction_vector(theta, phi, expected):
    result = direction_vector(theta, phi)
    assert np.allclose(result, expected, atol=1e-6), f"Expected {expected}, got {result}"


# test for dot_product
def test_dot_product():
    a = np.array([1, 0, 0])
    b = np.array([0, 1, 0])
    result = dot_product(a, b)
    assert result == 0, f"Expected dot product of orthogonal vectors to be 0, got {result}"

    a = np.array([1, 1, 1])
    b = np.array([1, 1, 1])
    result = dot_product(a, b)
    assert result == 3, f"Expected dot product of [1, 1, 1] and [1, 1, 1] to be 3, got {result}"


# test for normalize
def test_normalize():
    vec = np.array([3, 4, 0])
    result = normalize(vec)
    expected = np.array([0.6, 0.8, 0.0])
    assert np.allclose(result, expected, atol=1e-6), f"Expected {expected}, got {result}"

    vec = np.array([0, 0, 0])
    result = normalize(vec)
    expected = np.array([0, 0, 0]) # should return zero vector
    assert np.allclose(result, expected, atol=1e-6), f"Expected {expected}, got {result}"


# test for rayleigh_matrix
def test_rayleigh_matrix():
    result = rayleigh_matrix(90) # test at 90 degrees
    expected = np.array([
        [1, -1, 0, 0],
        [-1, 1, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0]
    ])
    assert np.allclose(result, expected, atol=1e-6), f"Expected {expected}, got {result}"

    result = rayleigh_matrix(0) # test at 0 degrees
    expected = np.eye(4) # identity matrix at 0 degrees
    assert np.allclose(result, expected, atol=1e-6), f"Expected identity matrix, got {result}"


# test for rotation_matrix
def test_rotation_matrix():
    result = rotation_matrix(90) # test at 90 degrees
    expected = np.array([
        [1, 0, 0, 0],
        [0, 0, 1, 0],
        [0, -1, 0, 0],
        [0, 0, 0, 1]
    ])
    assert np.allclose(result, expected, atol=1e-6), f"Expected {expected}, got {result}"

    result = rotation_matrix(0) # test at 0 degrees
    expected = np.eye(4) # identity matrix at 0 degrees
    assert np.allclose(result, expected, atol=1e-6), f"Expected identity matrix, got {result}"


# test for rotation_angles
def test_rotation_angles():
    k_i = np.array([0, 0, 1]) # initial vector along z-axis
    k_f = np.array([1, 0, 0]) # final vector along x-axis
    phi, psi, theta = rotation_angles(k_i, k_f)

    assert np.isclose(theta, 90), f"Expected scattering angle theta to be 90 degrees, got {theta}"
    assert np.isclose(phi, 90), f"Expected rotation angle phi to be 90 degrees, got {phi}"
    assert np.isclose(psi, -90), f"Expected rotation angle psi to be -90 degrees, got {psi}"

    k_f = np.array([0, 1, 0])  # Final vector along y-axis
    phi, psi, theta = rotation_angles(k_i, k_f)

    assert np.isclose(theta, 90), f"Expected scattering angle theta to be 90 degrees, got {theta}"
    assert np.isclose(phi, 90), f"Expected rotation angle phi to be 90 degrees, got {phi}"
    assert np.isclose(psi, -90), f"Expected rotation angle psi to be -90 degrees, got {psi}"
