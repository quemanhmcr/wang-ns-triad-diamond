from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Sequence

STATUS = (
    "DRAFT_NATIVE_CURL_CENTERED_DEFECT_LAW__"
    "BELTRAMI_VARIANCE_AND_TRIAD_CURVATURE_CURRENT"
)


@dataclass(frozen=True)
class CurlCenteredState:
    """Spectral invariants of r=(curl-lambda)u for a finite helical state.

    ``signed_frequencies`` are curl eigenvalues a=s|k| and ``modal_energies``
    are the corresponding nonnegative |u_a|^2 masses.  The zero field is
    excluded only to avoid the removable H/E quotient.
    """

    energy: float
    helicity: float
    enstrophy: float
    mean_signed_frequency: float
    defect_energy: float
    defect_variance: float
    curl_defect_energy: float

    def __post_init__(self) -> None:
        vals = (
            self.energy,
            self.helicity,
            self.enstrophy,
            self.mean_signed_frequency,
            self.defect_energy,
            self.defect_variance,
            self.curl_defect_energy,
        )
        if not all(math.isfinite(x) for x in vals):
            raise ValueError("finite curl-centered state required")
        if self.energy <= 0.0:
            raise ValueError("positive energy required")
        scale = max(1.0, self.enstrophy, abs(self.helicity), self.energy)
        if self.defect_energy < -1.0e-12 * scale:
            raise AssertionError("Beltrami defect lost nonnegativity")
        if self.curl_defect_energy < -1.0e-12 * scale:
            raise AssertionError("curl-defect norm lost nonnegativity")


def curl_centered_state(
    signed_frequencies: Sequence[float], modal_energies: Sequence[float]
) -> CurlCenteredState:
    if len(signed_frequencies) != len(modal_energies) or not signed_frequencies:
        raise ValueError("matching nonempty signed-frequency/energy data required")
    a = tuple(float(x) for x in signed_frequencies)
    e = tuple(float(x) for x in modal_energies)
    if not all(math.isfinite(x) for x in a + e):
        raise ValueError("finite spectral data required")
    if any(x < 0.0 for x in e):
        raise ValueError("modal energies must be nonnegative")
    E = math.fsum(e)
    if E <= 0.0:
        raise ValueError("positive total energy required")
    H = math.fsum(ai * ei for ai, ei in zip(a, e))
    Z = math.fsum(ai * ai * ei for ai, ei in zip(a, e))
    lam = H / E
    B_direct = math.fsum((ai - lam) ** 2 * ei for ai, ei in zip(a, e))
    B_invariant = Z - H * H / E
    scale = max(1.0, Z, B_direct, abs(B_invariant))
    if abs(B_direct - B_invariant) > 2.0e-12 * scale:
        raise AssertionError("B=Z-H^2/E failed")
    Cr2 = math.fsum(ai * ai * (ai - lam) ** 2 * ei for ai, ei in zip(a, e))
    return CurlCenteredState(
        energy=E,
        helicity=H,
        enstrophy=Z,
        mean_signed_frequency=lam,
        defect_energy=max(0.0, B_invariant),
        defect_variance=max(0.0, B_direct / E),
        curl_defect_energy=max(0.0, Cr2),
    )


def energy_dissipation_split(state: CurlCenteredState, viscosity: float) -> dict[str, float]:
    """Exact split (1/2)E'=-nu Z=-nu lambda^2 E-nu B."""
    nu = float(viscosity)
    if nu < 0.0 or not math.isfinite(nu):
        raise ValueError("finite nonnegative viscosity required")
    baseline = nu * state.mean_signed_frequency**2 * state.energy
    defect = nu * state.defect_energy
    total = nu * state.enstrophy
    scale = max(1.0, total, baseline + defect)
    if abs(total - baseline - defect) > 2.0e-12 * scale:
        raise AssertionError("curl-centered energy dissipation split failed")
    return {
        "minus_half_energy_rate": total,
        "beltrami_baseline": baseline,
        "defect_dissipation": defect,
    }


def viscous_defect_rate(state: CurlCenteredState, viscosity: float) -> float:
    """Pure-viscous contribution B'=-2 nu ||curl r||_2^2."""
    nu = float(viscosity)
    if nu < 0.0 or not math.isfinite(nu):
        raise ValueError("finite nonnegative viscosity required")
    return -2.0 * nu * state.curl_defect_energy


def defect_balance_rate(
    stretching_inner_product: float, curl_defect_energy: float, viscosity: float
) -> float:
    """Full native identity B'=2<curl r,u x r>-2 nu||curl r||^2.

    The first argument is the physical scalar <curl r,u x r>; this helper does
    not manufacture it from an analysis proxy.
    """
    q = float(stretching_inner_product)
    cr2 = float(curl_defect_energy)
    nu = float(viscosity)
    if not all(math.isfinite(x) for x in (q, cr2, nu)) or cr2 < 0.0 or nu < 0.0:
        raise ValueError("finite stretching data and nonnegative norm/viscosity required")
    return 2.0 * q - 2.0 * nu * cr2


def triad_curvature_current(
    signed_frequencies: Sequence[float], modal_works: Sequence[float]
) -> float:
    """Q_Delta=sum a_i^2 T_i, after enforcing the two Euler null laws."""
    if len(signed_frequencies) != 3 or len(modal_works) != 3:
        raise ValueError("exactly three triad frequencies and works required")
    a = tuple(float(x) for x in signed_frequencies)
    T = tuple(float(x) for x in modal_works)
    if not all(math.isfinite(x) for x in a + T):
        raise ValueError("finite triad data required")
    energy_residual = math.fsum(T)
    helicity_residual = math.fsum(ai * ti for ai, ti in zip(a, T))
    scale = max(1.0, *(abs(x) for x in T), *(abs(ai * ti) for ai, ti in zip(a, T)))
    if abs(energy_residual) > 2.0e-12 * scale:
        raise ValueError("triad works do not satisfy energy conservation")
    if abs(helicity_residual) > 2.0e-12 * scale:
        raise ValueError("triad works do not satisfy signed-helicity conservation")
    return math.fsum(ai * ai * ti for ai, ti in zip(a, T))


def second_divided_difference(
    nodes: Sequence[float], observable: Callable[[float], float]
) -> float:
    """Second divided difference Phi[a0,a1,a2] for pairwise distinct nodes."""
    if len(nodes) != 3:
        raise ValueError("exactly three nodes required")
    a0, a1, a2 = (float(x) for x in nodes)
    if not all(math.isfinite(x) for x in (a0, a1, a2)):
        raise ValueError("finite nodes required")
    scale = max(1.0, abs(a0), abs(a1), abs(a2))
    if min(abs(a0 - a1), abs(a0 - a2), abs(a1 - a2)) <= 1.0e-14 * scale:
        raise ValueError("pairwise distinct nodes required; degenerate cases use the confluent limit")
    p0, p1, p2 = (float(observable(x)) for x in (a0, a1, a2))
    if not all(math.isfinite(x) for x in (p0, p1, p2)):
        raise ValueError("observable must be finite at all nodes")
    return (
        p0 / ((a0 - a1) * (a0 - a2))
        + p1 / ((a1 - a0) * (a1 - a2))
        + p2 / ((a2 - a0) * (a2 - a1))
    )


def triad_observable_response(
    signed_frequencies: Sequence[float],
    modal_works: Sequence[float],
    observable: Callable[[float], float],
) -> dict[str, float]:
    """Exact response sum Phi(a_i)T_i=Q_Delta Phi[a0,a1,a2]."""
    a = tuple(float(x) for x in signed_frequencies)
    T = tuple(float(x) for x in modal_works)
    Q = triad_curvature_current(a, T)
    dd = second_divided_difference(a, observable)
    direct = math.fsum(float(observable(ai)) * ti for ai, ti in zip(a, T))
    represented = Q * dd
    scale = max(1.0, abs(direct), abs(represented))
    residual = direct - represented
    if abs(residual) > 5.0e-12 * scale:
        raise AssertionError("triad second-divided-difference representation failed")
    return {
        "direct_response": direct,
        "curvature_current": Q,
        "second_divided_difference": dd,
        "represented_response": represented,
        "residual": residual,
    }


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "rotational_form": "for divergence-free u, C=curl and C^2=-Delta, so u_t=P(u x Cu)-nu C^2u",
        "curl_gauge": "u x Cu=u x (C-lambda)u for every spatially constant scalar lambda",
        "canonical_center": "lambda=<u,Cu>/||u||_2^2 minimizes ||(C-alpha)u||_2^2",
        "defect": "r=(C-lambda)u and B=||r||_2^2=Z-H^2/E=E Var_mu(a) for signed curl frequency a=s|k|",
        "energy": "(1/2)E'=-nu Z=-nu lambda^2 E-nu B; hence nu int B dt is actual viscous energy dissipation, not an imposed reset currency",
        "nonlinear_force": "P(u x Cu)=P(u x r); B=0 is an exact Beltrami state and turns off the Euler nonlinearity",
        "defect_evolution": "B'=2<curl r,u x r>-2nu||curl r||_2^2",
        "triad_affine_null": "closed-triad work annihilates Phi(a)=c0+c1 a because sum T_i=sum a_i T_i=0",
        "triad_curvature": "for distinct a_i, sum Phi(a_i)T_i=Q_Delta Phi[a0,a1,a2], Q_Delta=sum a_i^2 T_i; degenerate nodes follow by confluent continuity",
        "quadratic_case": "Phi(a)=a^2 has second divided difference 1, so Q_Delta is the nonlinear enstrophy/vortex-stretching current",
        "interpretation": "vortex stretching is signed-curl spectral variance production; viscosity reads the same defect through C^2",
        "temporal_matching_used": False,
        "owner_bookkeeping_used": False,
        "global_regularity_claimed": False,
        "open_question": "identify the sharp intrinsic conversion law from Beltrami defect into canonical positive true-upward spectral transfer, and test whether the r_* geometry is its genuine extremizer",
    }
