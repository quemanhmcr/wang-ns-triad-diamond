import math
import numpy as np
from src.curvature_sideband_irrep import (
    clean_combined_observability, divergence_trace, envelope_reconstruction,
    hook_M, hook_component, hook_from_M, hook_strain_sideband, random_divfree_curvature,
)
from src.hermite_helicity_ledger import sym3


def test_exact_irrep_split():
    rng=np.random.default_rng(10); B=random_divfree_curvature(rng)
    E=envelope_reconstruction(B); H=hook_component(B)
    assert np.linalg.norm(sym3(E)-sym3(B))<1e-11
    assert np.linalg.norm(sym3(H))<1e-11
    assert np.linalg.norm(divergence_trace(H))<1e-11
    assert abs(np.einsum('abc,abc',E,H))<1e-10


def test_hook_reconstruction_and_norms():
    rng=np.random.default_rng(11); B=random_divfree_curvature(rng); H=hook_component(B); M=hook_M(H)
    assert np.linalg.norm(hook_from_M(M)-H)<1e-11
    assert math.isclose(np.sum(H*H),6*np.sum(M*M),rel_tol=1e-11,abs_tol=1e-11)
    C=hook_strain_sideband(H)
    assert math.isclose(np.sum(C*C),.25*np.sum(H*H),rel_tol=1e-11,abs_tol=1e-11)


def test_combined_observability():
    rng=np.random.default_rng(12)
    for _ in range(100): assert clean_combined_observability(random_divfree_curvature(rng))>-1e-10
