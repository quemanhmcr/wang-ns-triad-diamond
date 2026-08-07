import math
import numpy as np
import pytest

from src.ancestor_reservoir_sync import (
    CLEAN_GENERATION_PROGRESS,
    CLEAN_RESERVOIR_GROWTH,
    affine_material_step,
    amortized_service_capacity_upper,
    arb_reservoir_certificate,
    critical_mass_service_ratio_upper,
    kelvin_covector,
    low_band_service_from_physical_energy,
    physical_energy_service_ratio_upper,
    total_amortized_service_upper,
)


def test_exact_service_half_life_constants():
    assert critical_mass_service_ratio_upper() < 1 / 2
    assert physical_energy_service_ratio_upper() < 1 / 2
    assert float(CLEAN_RESERVOIR_GROWTH / CLEAN_GENERATION_PROGRESS) == 21 / 32


def test_kelvin_covector_exact_constant_A():
    A = np.array([[0.2, 0.1, 0.0], [-0.3, -0.1, 0.2], [0.0, 0.1, -0.1]])
    A -= np.trace(A) / 3 * np.eye(3)
    L0 = np.array([[1.2, 0.1, 0.0], [0.0, 0.9, 0.2], [0.1, 0.0, 1.1]])
    k0 = np.array([1.0, -0.4, 0.7])
    L1, k1 = affine_material_step(A, 0.2, L0, k0)
    assert np.linalg.norm(kelvin_covector(L1, k1) - kelvin_covector(L0, k0)) < 1e-12


def test_amortized_capacity_below_clean_half_life():
    base = low_band_service_from_physical_energy(1.0, 4.0, 3.0)
    for q in range(12):
        assert amortized_service_capacity_upper(q, 1.0, 4.0, 3.0) <= base * 0.5**q + 1e-14
    assert total_amortized_service_upper(1.0, 4.0, 3.0) < 2 * base


def test_arb_optional():
    pytest.importorskip("flint")
    cert = arb_reservoir_certificate()
    assert cert["clean_generation_progress"] == ">8/5"
    assert cert["critical_mass_service_clean"] == "<1/2"
