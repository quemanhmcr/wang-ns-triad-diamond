from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import combinations
from typing import Sequence


STATUS = (
    "DRAFT_NATIVE_CURL_KRYLOV_CURRENT_LAW__"
    "PAIRWISE_CURL_SHEAR__CRITICAL_GRAM_DETERMINANT__"
    "SECOND_KRYLOV_RESIDUAL__SHARP_THREE_POINT_CURRENT"
)


@dataclass(frozen=True)
class CurlKrylovState:
    """First two intrinsic curl-Krylov levels of a finite spectral state.

    The spectral nodes are signed curl eigenvalues ``a=s|k|`` and the weights
    are physical modal energies.  ``B`` is the centered curl defect,
    ``Delta=E*B=E*Z-H^2`` is the critical two-vector Gram determinant, and
    ``h2`` is the squared norm of the component of ``curl r`` orthogonal to
    ``span{u,r}``.
    """

    energy: float
    helicity: float
    enstrophy: float
    center: float
    defect_energy: float
    critical_determinant: float
    defect_curl_moment: float
    curl_defect_energy: float
    fourth_curl_moment: float
    second_residual_energy: float
    third_hankel_determinant: float
    beta1: float
    alpha1: float
    beta2: float


@dataclass(frozen=True)
class ThreePointCurrentLaw:
    """Sharp three-signed-frequency amplitude law before any recursion."""

    total_energy: float
    defect_energy: float
    critical_determinant: float
    third_hankel_determinant: float
    barycentric_efficiency: float
    sharp_efficiency_bound: float
    waleffe_magnitude: float
    phase_cosine_abs: float
    curvature_current_magnitude: float
    gross_energy_current_magnitude: float
    gross_current_upper: float
    global_gross_current_upper: float
    equality_weights: tuple[float, float, float]
    median_index: int


def _spectral_data(
    signed_frequencies: Sequence[float], modal_energies: Sequence[float]
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    if len(signed_frequencies) != len(modal_energies) or not signed_frequencies:
        raise ValueError("matching nonempty signed-frequency/energy data required")
    a = tuple(float(x) for x in signed_frequencies)
    e = tuple(float(x) for x in modal_energies)
    if not all(math.isfinite(x) for x in a + e):
        raise ValueError("finite spectral data required")
    if any(x < 0.0 for x in e):
        raise ValueError("modal energies must be nonnegative")
    if math.fsum(e) <= 0.0:
        raise ValueError("positive total energy required")
    return a, e


def hankel_vandermonde_determinant(
    signed_frequencies: Sequence[float], modal_energies: Sequence[float], order: int
) -> float:
    """Exact Hankel/Gram determinant by Cauchy--Binet/Vandermonde expansion.

    For ``order=n`` this is the Gram determinant of
    ``{u,Cu,...,C^n u}``.  For a discrete spectral law it equals the sum over
    every ``n+1``-node subset of ``prod energy * Vandermonde^2``.
    """

    a, e = _spectral_data(signed_frequencies, modal_energies)
    n = int(order)
    if n < 0:
        raise ValueError("nonnegative determinant order required")
    m = n + 1
    if m > len(a):
        return 0.0
    terms: list[float] = []
    for ids in combinations(range(len(a)), m):
        weight = math.prod(e[i] for i in ids)
        vand = 1.0
        for p, q in combinations(ids, 2):
            vand *= a[p] - a[q]
        terms.append(weight * vand * vand)
    return max(0.0, math.fsum(terms))


def pairwise_curl_shear_capacity(
    signed_frequencies: Sequence[float], modal_energies: Sequence[float]
) -> float:
    """Return sum_{i<j} e_i e_j (a_i-a_j)^2 = E Z-H^2."""

    return hankel_vandermonde_determinant(signed_frequencies, modal_energies, 1)


def curl_krylov_state(
    signed_frequencies: Sequence[float], modal_energies: Sequence[float]
) -> CurlKrylovState:
    a, e = _spectral_data(signed_frequencies, modal_energies)
    E = math.fsum(e)
    H = math.fsum(ai * ei for ai, ei in zip(a, e))
    Z = math.fsum(ai * ai * ei for ai, ei in zip(a, e))
    lam = H / E
    B = math.fsum((ai - lam) ** 2 * ei for ai, ei in zip(a, e))
    D1 = E * B
    pairwise = pairwise_curl_shear_capacity(a, e)
    scale = max(1.0, abs(D1), abs(pairwise), E * Z, H * H)
    if abs(D1 - (E * Z - H * H)) > 3.0e-12 * scale:
        raise AssertionError("critical determinant EZ-H^2 failed")
    if abs(D1 - pairwise) > 3.0e-12 * scale:
        raise AssertionError("pairwise curl-shear identity failed")

    if B <= 2.0e-14 * max(1.0, Z, abs(H)):
        # Exact Beltrami boundary.  alpha1 is not geometrically needed here;
        # the zero convention keeps the finite certificate total and explicit.
        return CurlKrylovState(
            energy=E,
            helicity=H,
            enstrophy=Z,
            center=lam,
            defect_energy=max(0.0, B),
            critical_determinant=max(0.0, D1),
            defect_curl_moment=0.0,
            curl_defect_energy=0.0,
            fourth_curl_moment=math.fsum(ai**4 * ei for ai, ei in zip(a, e)),
            second_residual_energy=0.0,
            third_hankel_determinant=0.0,
            beta1=0.0,
            alpha1=0.0,
            beta2=0.0,
        )

    J = math.fsum(ai * (ai - lam) ** 2 * ei for ai, ei in zip(a, e))
    Cr2 = math.fsum(ai * ai * (ai - lam) ** 2 * ei for ai, ei in zip(a, e))
    Y = math.fsum(ai**4 * ei for ai, ei in zip(a, e))

    # Near a two-node spectral collapse, the direct difference
    # Cr2-B^2/E-J^2/B suffers catastrophic cancellation.  The Cauchy--Binet
    # Gram volume is the same exact quantity and stays manifestly positive, so
    # use D2/D1 as the canonical residual norm and reconstruct Cr2 afterwards.
    D2 = hankel_vandermonde_determinant(a, e, 2)
    h2 = D2 / D1
    reconstructed_Cr2 = B * B / E + J * J / B + h2
    hscale = max(1.0, Cr2, reconstructed_Cr2)
    if abs(Cr2 - reconstructed_Cr2) > 2.0e-10 * hscale:
        raise AssertionError("curl-defect orthogonal decomposition failed")

    return CurlKrylovState(
        energy=E,
        helicity=H,
        enstrophy=Z,
        center=lam,
        defect_energy=B,
        critical_determinant=D1,
        defect_curl_moment=J,
        curl_defect_energy=Cr2,
        fourth_curl_moment=Y,
        second_residual_energy=h2,
        third_hankel_determinant=D2,
        beta1=math.sqrt(B / E),
        alpha1=J / B,
        beta2=math.sqrt(h2 / B),
    )


def critical_determinant_log_rate(
    state: CurlKrylovState,
    defect_strain_inner: float,
    viscosity: float,
) -> dict[str, float]:
    """Exact log law for Delta=E*B when B>0.

    ``defect_strain_inner`` is the physical scalar ``int r.S r``.  Nothing in
    this helper infers it from capacity, a norm, or an owner label.
    """

    s = float(defect_strain_inner)
    nu = float(viscosity)
    if not math.isfinite(s) or not math.isfinite(nu) or nu < 0.0:
        raise ValueError("finite strain and nonnegative viscosity required")
    if state.defect_energy <= 0.0:
        raise ValueError("positive curl defect required for log Delta")
    nonlinear = 2.0 * s / state.defect_energy
    viscous_rayleigh = state.enstrophy / state.energy + state.curl_defect_energy / state.defect_energy
    viscous = -2.0 * nu * viscous_rayleigh
    lanczos_rayleigh = (
        state.center**2
        + 2.0 * state.beta1**2
        + state.alpha1**2
        + state.beta2**2
    )
    scale = max(1.0, abs(viscous_rayleigh), abs(lanczos_rayleigh))
    if abs(viscous_rayleigh - lanczos_rayleigh) > 5.0e-12 * scale:
        raise AssertionError("Rayleigh/Lanczos viscous rate mismatch")
    return {
        "log_delta_rate": nonlinear + viscous,
        "nonlinear_log_rate": nonlinear,
        "viscous_log_rate": viscous,
        "viscous_rayleigh": viscous_rayleigh,
        "lanczos_rayleigh": lanczos_rayleigh,
    }







def isolated_three_point_euler_law(
    signed_frequencies: Sequence[float],
    amplitudes: Sequence[complex],
    coupling: complex,
) -> dict[str, object]:
    """Exact isolated three-wave Euler source in the repository cyclic gauge.

    With ``a_i=s_i|k_i|`` and the common cyclic Waleffe coefficient ``g``,
    ``dot A_i=2(a_j-a_k)g conjugate(A_j A_k)`` reproduces the certified three
    rooted energy works.  This helper records the intrinsic radial/phase speed
    split and the conserved cubic quadrature; it is not a claim that a selected
    triad is an invariant subspace of the full PDE.
    """

    if len(signed_frequencies) != 3 or len(amplitudes) != 3:
        raise ValueError("exactly three signed frequencies and amplitudes required")
    a = tuple(float(x) for x in signed_frequencies)
    A = tuple(complex(x) for x in amplitudes)
    g = complex(coupling)
    if not all(math.isfinite(x) for x in a):
        raise ValueError("finite signed frequencies required")
    vals = tuple(z.real for z in A) + tuple(z.imag for z in A) + (g.real, g.imag)
    if not all(math.isfinite(x) for x in vals):
        raise ValueError("finite amplitudes/coupling required")
    d = (a[1] - a[2], a[2] - a[0], a[0] - a[1])
    rhs = (
        2.0 * d[0] * g * (A[1] * A[2]).conjugate(),
        2.0 * d[1] * g * (A[2] * A[0]).conjugate(),
        2.0 * d[2] * g * (A[0] * A[1]).conjugate(),
    )
    e = tuple(abs(z) ** 2 for z in A)
    T = tuple(2.0 * (A[i].conjugate() * rhs[i]).real for i in range(3))
    delta = pairwise_curl_shear_capacity(a, e)
    total_speed = math.fsum(abs(z) ** 2 for z in rhs)
    represented_total = 4.0 * abs(g) ** 2 * delta
    scale = max(1.0, total_speed, represented_total)
    if abs(total_speed - represented_total) > 8.0e-12 * scale:
        raise AssertionError("isolated triad total-speed/Delta identity failed")
    er = math.fsum(T)
    hr = math.fsum(ai * ti for ai, ti in zip(a, T))
    if abs(er) > 5.0e-12 * max(1.0, *(abs(x) for x in T)):
        raise AssertionError("isolated triad energy null law failed")
    if abs(hr) > 5.0e-12 * max(1.0, *(abs(ai * ti) for ai, ti in zip(a, T))):
        raise AssertionError("isolated triad helicity null law failed")

    radial_speed = 0.0
    phase_speed = 0.0
    for Ai, Fi in zip(A, rhs):
        if abs(Ai) == 0.0:
            phase_speed += abs(Fi) ** 2
            continue
        z = Ai.conjugate() * Fi / abs(Ai)
        radial_speed += z.real * z.real
        phase_speed += z.imag * z.imag
    if abs(radial_speed + phase_speed - total_speed) > 8.0e-12 * scale:
        raise AssertionError("isolated triad radial/phase Pythagorean law failed")

    P = A[0] * A[1] * A[2] * g.conjugate()
    Pdot = (
        rhs[0] * A[1] * A[2]
        + A[0] * rhs[1] * A[2]
        + A[0] * A[1] * rhs[2]
    ) * g.conjugate()
    if abs(Pdot.imag) > 8.0e-12 * max(1.0, abs(Pdot)):
        raise AssertionError("isolated triad cubic quadrature is not invariant")
    phase_cos = 0.0
    denom = abs(g) * math.sqrt(math.prod(e))
    if denom > 0.0:
        phase_cos = max(-1.0, min(1.0, P.real / denom))
    radial_rep = represented_total * phase_cos * phase_cos
    phase_rep = represented_total - radial_rep
    if abs(radial_speed - radial_rep) > 1.0e-10 * max(1.0, radial_speed, radial_rep):
        raise AssertionError("isolated triad radial speed lost phase-cosine factor")
    if abs(phase_speed - phase_rep) > 1.0e-10 * max(1.0, phase_speed, phase_rep):
        raise AssertionError("isolated triad phase speed lost complementary factor")
    return {
        "rhs": rhs,
        "modal_energy_rates": T,
        "critical_determinant": delta,
        "total_source_speed_squared": total_speed,
        "radial_energy_speed_squared": radial_speed,
        "phase_shape_speed_squared": phase_speed,
        "phase_cosine": phase_cos,
        "cubic_product": P,
        "cubic_product_rate": Pdot,
        "cubic_quadrature_invariant": P.imag,
        "energy_rate_residual": er,
        "helicity_rate_residual": hr,
    }

def three_point_martingale_spread(
    ordered_signed_frequencies: Sequence[float], donor_work: float
) -> dict[str, float]:
    """Unique same-time barycentric spread from the median curl node.

    ``a<m<b`` are ordered signed curl eigenvalues and ``donor_work=q>0`` is
    actual physical energy work leaving the median node.  The returned weights
    are work fractions, not a temporal ancestry probability.
    """

    if len(ordered_signed_frequencies) != 3:
        raise ValueError("exactly three ordered signed frequencies required")
    a, m, b = (float(x) for x in ordered_signed_frequencies)
    q = float(donor_work)
    if not all(math.isfinite(x) for x in (a, m, b, q)) or not (a < m < b) or q < 0.0:
        raise ValueError("finite a<m<b and nonnegative donor work required")
    left = (b - m) / (b - a)
    right = (m - a) / (b - a)
    mean = left * a + right * b
    variance = left * (a - m) ** 2 + right * (b - m) ** 2
    intrinsic_variance = (m - a) * (b - m)
    scale = max(1.0, abs(a), abs(m), abs(b), variance, intrinsic_variance)
    if abs(mean - m) > 3.0e-13 * scale:
        raise AssertionError("martingale barycenter failed")
    if abs(variance - intrinsic_variance) > 4.0e-13 * scale:
        raise AssertionError("martingale quadratic variation failed")
    return {
        "left_recipient_fraction": left,
        "right_recipient_fraction": right,
        "left_recipient_work": q * left,
        "right_recipient_work": q * right,
        "donor_work": q,
        "conditional_mean": mean,
        "conditional_variance": variance,
        "quadratic_variation_rate": q * variance,
    }


def martingale_observable_increment(
    ordered_signed_frequencies: Sequence[float],
    donor_work: float,
    observable_values: Sequence[float],
) -> dict[str, float]:
    """Observable production of the unique three-node martingale spread."""

    if len(observable_values) != 3:
        raise ValueError("three observable values required")
    a, m, b = (float(x) for x in ordered_signed_frequencies)
    phi_a, phi_m, phi_b = (float(x) for x in observable_values)
    spread = three_point_martingale_spread((a, m, b), donor_work)
    direct_per_mass = (
        spread["left_recipient_fraction"] * phi_a
        + spread["right_recipient_fraction"] * phi_b
        - phi_m
    )
    # Second divided difference on ordered nodes.
    dd = (
        phi_a / ((a - m) * (a - b))
        + phi_m / ((m - a) * (m - b))
        + phi_b / ((b - a) * (b - m))
    )
    represented_per_mass = spread["conditional_variance"] * dd
    scale = max(1.0, abs(direct_per_mass), abs(represented_per_mass))
    if abs(direct_per_mass - represented_per_mass) > 5.0e-12 * scale:
        raise AssertionError("martingale divided-difference generator failed")
    return {
        **spread,
        "observable_increment": spread["donor_work"] * direct_per_mass,
        "observable_increment_per_donor_work": direct_per_mass,
        "second_divided_difference": dd,
    }


def heterochiral_frontier_progress_side_bound(
    donor_radius_ratio: float, side_radius_ratio: float
) -> dict[str, float]:
    """Continuous no-free-progress law for a heterochiral strict-UV spread.

    Child radius is normalized to one.  ``D`` is the same-helicity median donor
    radius and ``S`` the opposite-helicity side-recipient radius.  The physical
    triangle requires ``D,S<1`` and ``D+S>1``.  Per common positive current R,
    donor/child/side works are ``1+S``, ``D+S`` and ``1-D``.
    """

    D = float(donor_radius_ratio)
    S = float(side_radius_ratio)
    if not all(math.isfinite(x) for x in (D, S)) or not (0.0 < D < 1.0 and 0.0 < S < 1.0 and D + S > 1.0):
        raise ValueError("strict UV physical triangle requires 0<D,S<1 and D+S>1")
    top = max(D, S)
    donor = 1.0 + S
    child = D + S
    side = 1.0 - D
    progress = child * math.log(1.0 / top)
    # If top=D, child<=2D and D log(1/D)<=1-D.  If top=S,
    # child<=2S and S log(1/S)<=1-S<=1-D.
    upper = 2.0 * side
    if progress > upper + 5.0e-13 * max(1.0, progress, upper):
        raise AssertionError("physical log-progress exceeded twice the compulsory side work")
    return {
        "donor_work_per_common_current": donor,
        "child_work_per_common_current": child,
        "side_work_per_common_current": side,
        "child_log_progress_per_common_current": progress,
        "progress_upper_from_side": upper,
        "high_child_retained_fraction": child / donor,
        "side_fraction": side / donor,
    }

def observable_tangent_gram(
    signed_frequencies: Sequence[float],
    modal_energies: Sequence[float],
    observable_values: Sequence[float],
) -> dict[str, float]:
    """Gram volume of ``u,Cu,phi(C)u`` and the tangent residual of ``phi``.

    The observable is supplied by its finite values on the occupied signed-curl
    nodes.  ``residual_energy`` is the squared norm of the projection of
    ``phi(C)u`` onto the common energy--helicity tangent space.
    """

    a, e = _spectral_data(signed_frequencies, modal_energies)
    if len(observable_values) != len(a):
        raise ValueError("one finite observable value per signed frequency required")
    phi = tuple(float(x) for x in observable_values)
    if not all(math.isfinite(x) for x in phi):
        raise ValueError("finite observable values required")
    E = math.fsum(e)
    H = math.fsum(ai * ei for ai, ei in zip(a, e))
    Z = math.fsum(ai * ai * ei for ai, ei in zip(a, e))
    P = math.fsum(pi * ei for pi, ei in zip(phi, e))
    CP = math.fsum(ai * pi * ei for ai, pi, ei in zip(a, phi, e))
    PP = math.fsum(pi * pi * ei for pi, ei in zip(phi, e))
    delta = E * Z - H * H
    # Cauchy--Binet is the native positive representation of this Gram volume.
    # It avoids catastrophic cancellation when phi is almost affine on a very
    # uneven spectral state.
    gram_terms: list[float] = []
    for i, j, k in combinations(range(len(a)), 3):
        det = (
            phi[i] * (a[j] - a[k])
            + phi[j] * (a[k] - a[i])
            + phi[k] * (a[i] - a[j])
        )
        gram_terms.append(e[i] * e[j] * e[k] * det * det)
    gram = math.fsum(gram_terms)
    residual = 0.0 if delta <= 0.0 else gram / delta
    return {
        "gram_determinant": gram,
        "critical_determinant": max(0.0, delta),
        "residual_energy": residual,
        "observable_stock": P,
        "curl_observable_cross": CP,
        "observable_square_stock": PP,
    }


def three_point_observable_volume_law(
    signed_frequencies: Sequence[float],
    modal_energies: Sequence[float],
    observable_values: Sequence[float],
    *,
    waleffe_magnitude: float,
    phase_cosine_abs: float = 1.0,
) -> dict[str, float]:
    """Universal three-node Waleffe/Gram production law for one observable.

    For a physical closed triad the cyclic work is ``T_i=(a_j-a_k)R``.
    Therefore ``sum phi_i T_i`` is the oriented determinant of the three
    columns ``1,a,phi`` times the same physical scalar ``R``.  Cauchy--Binet
    identifies its amplitude/Vandermonde factor with the Gram volume returned
    above.
    """

    a0, e0 = _spectral_data(signed_frequencies, modal_energies)
    if len(a0) != 3 or len(observable_values) != 3:
        raise ValueError("exactly three signed frequencies/energies/observable values required")
    if any(x <= 0.0 for x in e0):
        raise ValueError("strictly positive three-point energies required")
    phi = tuple(float(x) for x in observable_values)
    if not all(math.isfinite(x) for x in phi):
        raise ValueError("finite observable values required")
    g = float(waleffe_magnitude)
    pc = abs(float(phase_cosine_abs))
    if not math.isfinite(g) or g < 0.0 or g > 0.5 + 1.0e-12:
        raise ValueError("physical Waleffe magnitude must lie in [0,1/2]")
    if not math.isfinite(pc) or pc > 1.0 + 1.0e-12:
        raise ValueError("absolute phase cosine must lie in [0,1]")
    pc = min(1.0, pc)
    a = tuple(a0)
    e = tuple(e0)
    gram = observable_tangent_gram(a, e, phi)
    det = (
        phi[0] * (a[1] - a[2])
        + phi[1] * (a[2] - a[0])
        + phi[2] * (a[0] - a[1])
    )
    common_abs = 4.0 * g * pc * math.sqrt(math.prod(e))
    direct = common_abs * abs(det)
    represented = 4.0 * g * pc * math.sqrt(gram["gram_determinant"])
    scale = max(1.0, direct, represented)
    if abs(direct - represented) > 8.0e-11 * scale:
        raise AssertionError("universal three-point observable volume law failed")
    return {
        **gram,
        "observable_current_magnitude": direct,
        "represented_current_magnitude": represented,
        "waleffe_magnitude": g,
        "phase_cosine_abs": pc,
        "oriented_observable_determinant_abs": abs(det),
    }


def spectral_source_action(
    signed_frequencies: Sequence[float],
    modal_energies: Sequence[float],
    modal_energy_rates: Sequence[float],
) -> dict[str, float]:
    """Physical Fisher--Rao speed of a signed-curl energy source.

    For one curl eigenspace with energy ``E_a`` and nonlinear modal-energy rate
    ``S_a=2<u_a,F_a>``, the component of ``F_a`` parallel to the existing state
    has squared norm ``S_a^2/(4E_a)``.  Summing gives the norm of the Euler
    velocity projected onto the state's curl-cyclic subspace.

    The source is required to satisfy the Euler affine null laws.  A zero-energy
    node must have zero first-order energy rate; newly born amplitude is a
    higher-order state effect and is not forged into this action.
    """

    a, e = _spectral_data(signed_frequencies, modal_energies)
    if len(modal_energy_rates) != len(a):
        raise ValueError("one modal-energy rate per signed curl node required")
    s = tuple(float(x) for x in modal_energy_rates)
    if not all(math.isfinite(x) for x in s):
        raise ValueError("finite modal-energy rates required")
    scale = max(1.0, *(abs(x) for x in s))
    if abs(math.fsum(s)) > 5.0e-11 * scale:
        raise ValueError("spectral source violates Euler energy conservation")
    hscale = max(1.0, *(abs(ai * si) for ai, si in zip(a, s)))
    if abs(math.fsum(ai * si for ai, si in zip(a, s))) > 5.0e-11 * hscale:
        raise ValueError("spectral source violates Euler helicity conservation")
    action_terms: list[float] = []
    for ei, si in zip(e, s):
        if ei <= 0.0:
            if abs(si) > 5.0e-12 * scale:
                raise ValueError("zero-energy node cannot carry first-order physical energy work")
            continue
        action_terms.append(si * si / (4.0 * ei))
    return {
        "spectral_velocity_norm_squared": math.fsum(action_terms),
        "energy_rate_residual": math.fsum(s),
        "helicity_rate_residual": math.fsum(ai * si for ai, si in zip(a, s)),
    }


def observable_source_speed_bound(
    signed_frequencies: Sequence[float],
    modal_energies: Sequence[float],
    modal_energy_rates: Sequence[float],
    observable_values: Sequence[float],
) -> dict[str, float]:
    """All-observable Cauchy law for the one physical spectral Euler velocity."""

    action = spectral_source_action(signed_frequencies, modal_energies, modal_energy_rates)
    gram = observable_tangent_gram(signed_frequencies, modal_energies, observable_values)
    phi = tuple(float(x) for x in observable_values)
    rates = tuple(float(x) for x in modal_energy_rates)
    direct = math.fsum(pi * si for pi, si in zip(phi, rates))
    rhs = 2.0 * math.sqrt(
        gram["residual_energy"] * action["spectral_velocity_norm_squared"]
    )
    if abs(direct) > rhs + 8.0e-11 * max(1.0, abs(direct), rhs):
        raise AssertionError("spectral observable speed exceeded the physical tangent action")
    return {
        **action,
        **gram,
        "observable_rate": direct,
        "observable_rate_upper": rhs,
    }


def three_point_spectral_speed_law(
    signed_frequencies: Sequence[float],
    modal_energies: Sequence[float],
    *,
    waleffe_magnitude: float,
    phase_cosine: float = 1.0,
) -> dict[str, float]:
    """Exact Fisher--Rao speed of one physical three-node cyclic work law."""

    a0, e0 = _spectral_data(signed_frequencies, modal_energies)
    if len(a0) != 3 or any(x <= 0.0 for x in e0):
        raise ValueError("three strictly positive signed-curl energy nodes required")
    g = float(waleffe_magnitude)
    c = float(phase_cosine)
    if not math.isfinite(g) or g < 0.0 or g > 0.5 + 1.0e-12:
        raise ValueError("physical Waleffe magnitude must lie in [0,1/2]")
    if not math.isfinite(c) or abs(c) > 1.0 + 1.0e-12:
        raise ValueError("phase cosine must lie in [-1,1]")
    c = max(-1.0, min(1.0, c))
    a = tuple(a0)
    e = tuple(e0)
    common = 4.0 * g * c * math.sqrt(math.prod(e))
    T = (
        (a[1] - a[2]) * common,
        (a[2] - a[0]) * common,
        (a[0] - a[1]) * common,
    )
    action = spectral_source_action(a, e, T)
    delta = pairwise_curl_shear_capacity(a, e)
    represented = 4.0 * g * g * c * c * delta
    scale = max(1.0, represented, action["spectral_velocity_norm_squared"])
    if abs(action["spectral_velocity_norm_squared"] - represented) > 8.0e-12 * scale:
        raise AssertionError("three-point spectral speed lost Waleffe/Delta factorization")
    return {
        **action,
        "critical_determinant": delta,
        "represented_spectral_speed": represented,
        "speed_efficiency": 2.0 * g * abs(c),
        "modal_energy_rates": T,
    }


def intrinsic_spectral_reynolds_barrier(
    state: CurlKrylovState, spectral_velocity_norm_squared: float, viscosity: float
) -> dict[str, float]:
    """Window-free concentration barrier for growth of the critical determinant.

    ``A_spec`` is the squared norm of the actual Euler velocity projected onto
    the state's curl-cyclic subspace.  ``V_spec=Delta/A_spec`` therefore scales
    as a physical interaction volume.  Completing the square in the exact
    Lanczos determinant law gives a universal upper bound on ``(log Delta)'``.
    """

    A = float(spectral_velocity_norm_squared)
    nu = float(viscosity)
    if not math.isfinite(A) or A < 0.0 or not math.isfinite(nu) or nu <= 0.0:
        raise ValueError("finite nonnegative spectral action and positive viscosity required")
    if state.defect_energy <= 0.0 or state.enstrophy <= 0.0:
        return {
            "interaction_volume": math.inf,
            "spectral_reynolds": 0.0,
            "fixed_state_log_delta_upper": -math.inf,
            "relaxed_log_delta_upper": -math.inf,
            "growth_requires_reynolds_above_two": True,
        }
    B = state.defect_energy
    N2 = state.enstrophy / state.energy
    delta = state.critical_determinant
    V = math.inf if A == 0.0 else delta / A
    re2 = A / (nu * nu * N2 * B)
    fixed = (
        2.0 * state.beta2 * math.sqrt(A / B)
        - 2.0 * nu * (
            state.center**2
            + 2.0 * state.beta1**2
            + state.alpha1**2
            + state.beta2**2
        )
    )
    relaxed = (
        A / (2.0 * nu * B)
        - 2.0 * nu * (
            state.center**2 + 2.0 * state.beta1**2 + state.alpha1**2
        )
    )
    if fixed > relaxed + 8.0e-12 * max(1.0, abs(fixed), abs(relaxed)):
        raise AssertionError("completed-square determinant upper bound failed")
    # Since center^2+2 beta1^2+alpha1^2 >= N^2, relaxed positivity
    # forces A/(nu^2 N^2 B)>4.  This is a necessary, not sufficient, condition.
    return {
        "interaction_volume": V,
        "spectral_reynolds": math.sqrt(max(0.0, re2)),
        "spectral_reynolds_squared": max(0.0, re2),
        "fixed_state_log_delta_upper": fixed,
        "relaxed_log_delta_upper": relaxed,
        "growth_requires_reynolds_above_two": True,
        "intrinsic_critical_mass": A / (N2 * B),
        "critical_mass_threshold": 4.0 * nu * nu,
    }

def critical_beltrami_split(
    signed_frequencies: Sequence[float], modal_energies: Sequence[float]
) -> dict[str, float]:
    """Split Delta into radial-magnitude variance plus helicity coexistence.

    K=sum |a|e is the positive critical H^{1/2} stock in the helical spectral
    normalization.  The identity

        EZ-H^2 = (EZ-K^2) + (K^2-H^2)

    separates radial diversity from coexistence of the two curl signs without
    introducing a branch label.
    """

    a, e = _spectral_data(signed_frequencies, modal_energies)
    E = math.fsum(e)
    H = math.fsum(ai * ei for ai, ei in zip(a, e))
    Z = math.fsum(ai * ai * ei for ai, ei in zip(a, e))
    K = math.fsum(abs(ai) * ei for ai, ei in zip(a, e))
    Kp = math.fsum(abs(ai) * ei for ai, ei in zip(a, e) if ai > 0.0)
    Km = math.fsum(abs(ai) * ei for ai, ei in zip(a, e) if ai < 0.0)
    delta = pairwise_curl_shear_capacity(a, e)
    radial = math.fsum(
        e[i] * e[j] * (abs(a[i]) - abs(a[j])) ** 2
        for i, j in combinations(range(len(a)), 2)
    )
    helicity_mix = 4.0 * Kp * Km
    scale = max(1.0, delta, radial, helicity_mix)
    if abs(delta - radial - helicity_mix) > 8.0e-12 * scale:
        raise AssertionError("critical Beltrami positive split failed")
    # Moment forms are recorded values, not the numerically preferred proof.
    moment_scale = max(1.0, E * Z, H * H, K * K)
    if abs((E * Z - K * K) - radial) > 2.0e-10 * moment_scale:
        raise AssertionError("radial moment reconstruction failed")
    return {
        "critical_stock": K,
        "positive_helicity_stock": Kp,
        "negative_helicity_stock": Km,
        "critical_determinant": delta,
        "radial_variance_component": max(0.0, radial),
        "helicity_coexistence_component": max(0.0, helicity_mix),
    }


def critical_tangent_correlation_geometry(
    signed_frequencies: Sequence[float], modal_energies: Sequence[float]
) -> dict[str, float]:
    """Exact covariance geometry behind the critical ``phi(a)=|a|`` source."""

    a, e = _spectral_data(signed_frequencies, modal_energies)
    split = critical_beltrami_split(a, e)
    E = math.fsum(e)
    H = math.fsum(ai * ei for ai, ei in zip(a, e))
    K = split["critical_stock"]
    L = math.fsum(ai * abs(ai) * ei for ai, ei in zip(a, e))
    delta = split["critical_determinant"]
    radial = split["radial_variance_component"]
    gram = observable_tangent_gram(a, e, tuple(abs(ai) for ai in a))["gram_determinant"]
    cross = math.fsum(
        e[i] * e[j] * (a[i] - a[j]) * (abs(a[i]) - abs(a[j]))
        for i, j in combinations(range(len(a)), 2)
    )
    lhs = cross * cross + E * gram
    rhs = delta * radial
    scale = max(1.0, lhs, rhs)
    if abs(lhs - rhs) > 2.0e-10 * scale:
        raise AssertionError("critical tangent covariance determinant failed")
    rho = 0.0
    decorrelation = 0.0
    if delta > 0.0 and radial > 0.0:
        denom = delta * radial
        # The positive Gram quotient is the stable definition of 1-rho^2 near
        # almost-affine states.  Recover only the sign of rho from the moment
        # cross covariance; do not form 1-rho^2 by catastrophic subtraction.
        decorrelation = max(0.0, min(1.0, E * gram / denom))
        rho_mag = math.sqrt(max(0.0, 1.0 - decorrelation))
        rho = math.copysign(rho_mag, cross) if cross != 0.0 else 0.0
    return {
        **split,
        "signed_radius_cross_covariance": cross,
        "signed_radius_correlation": rho,
        "critical_tangent_gram": gram,
        "critical_tangent_residual_energy": 0.0 if delta <= 0.0 else gram / delta,
        "decorrelation_factor": decorrelation,
    }

def rms_curl_scale_log_rate(
    state: CurlKrylovState,
    defect_strain_inner: float,
    viscosity: float,
) -> dict[str, float]:
    """Exact logarithmic speed of N=sqrt(Z/E).

    The nonlinear contribution is ``theta*sigma_r`` with
    ``theta=B/Z``.  The viscous contribution is minus the variance of ``a^2``
    under the normalized physical energy law, divided by its mean ``Z/E``.
    Hence viscosity can never increase the RMS curl scale.
    """

    s = float(defect_strain_inner)
    nu = float(viscosity)
    if not math.isfinite(s) or not math.isfinite(nu) or nu < 0.0:
        raise ValueError("finite strain and nonnegative viscosity required")
    if state.enstrophy <= 0.0:
        return {
            "log_rms_scale_rate": 0.0,
            "nonlinear_log_rate": 0.0,
            "viscous_log_rate": 0.0,
            "defect_fraction": 0.0,
            "squared_frequency_variance": 0.0,
        }
    if state.defect_energy <= 0.0:
        sigma = 0.0
    else:
        sigma = s / state.defect_energy
    R = state.enstrophy / state.energy
    theta = state.defect_energy / state.enstrophy
    var_a2 = state.fourth_curl_moment / state.energy - R * R
    vscale = max(1.0, state.fourth_curl_moment / state.energy, R * R)
    if var_a2 < -5.0e-12 * vscale:
        raise AssertionError("squared curl-frequency variance lost nonnegativity")
    var_a2 = max(0.0, var_a2)
    nonlinear = theta * sigma
    viscous = -nu * var_a2 / R
    return {
        "log_rms_scale_rate": nonlinear + viscous,
        "nonlinear_log_rate": nonlinear,
        "viscous_log_rate": viscous,
        "defect_fraction": theta,
        "squared_frequency_variance": var_a2,
    }

def _three_point_equality_weights(a: tuple[float, float, float]) -> tuple[tuple[float, float, float], int]:
    order = sorted(range(3), key=lambda i: a[i])
    lo, mid, hi = order
    total_gap = 2.0 * (a[hi] - a[lo])
    if total_gap <= 0.0:
        raise ValueError("three distinct signed frequencies required")
    out = [0.0, 0.0, 0.0]
    # Each energy is proportional to the gap opposite that node.
    out[lo] = (a[hi] - a[mid]) / total_gap
    out[mid] = (a[hi] - a[lo]) / total_gap
    out[hi] = (a[mid] - a[lo]) / total_gap
    return (out[0], out[1], out[2]), mid


def three_point_current_law(
    signed_frequencies: Sequence[float],
    modal_energies: Sequence[float],
    *,
    waleffe_magnitude: float,
    phase_cosine_abs: float = 1.0,
) -> ThreePointCurrentLaw:
    """Sharp amplitude law for one three-node closed-helical-triad state.

    The theorem is algebraic once the physical Waleffe magnitude and phase are
    supplied by the closed triad.  It does not assert that arbitrary three curl
    nodes form a Fourier-closed triad.
    """

    a0, e0 = _spectral_data(signed_frequencies, modal_energies)
    if len(a0) != 3:
        raise ValueError("exactly three signed frequencies required")
    if any(x <= 0.0 for x in e0):
        raise ValueError("strictly positive three-point energies required")
    a = (a0[0], a0[1], a0[2])
    e = (e0[0], e0[1], e0[2])
    gaps = (abs(a[0] - a[1]), abs(a[0] - a[2]), abs(a[1] - a[2]))
    scale_a = max(1.0, *(abs(x) for x in a))
    if min(gaps) <= 1.0e-14 * scale_a:
        raise ValueError("pairwise distinct signed frequencies required")
    g = float(waleffe_magnitude)
    pc = abs(float(phase_cosine_abs))
    if not math.isfinite(g) or not math.isfinite(pc) or g < 0.0 or g > 0.5 + 1.0e-12:
        raise ValueError("physical Waleffe magnitude must lie in [0,1/2]")
    if pc > 1.0 + 1.0e-12:
        raise ValueError("absolute phase cosine must lie in [0,1]")
    pc = min(1.0, pc)

    st = curl_krylov_state(a, e)
    E = st.energy
    B = st.defect_energy
    p = tuple(x / E for x in e)
    gap_sum = math.fsum(gaps)
    variance = B / E
    chi = gap_sum * math.sqrt(math.prod(p) / variance)
    if chi > 1.0 + 8.0e-12:
        raise AssertionError("barycentric current efficiency exceeded one")
    chi = min(1.0, max(0.0, chi))

    equality, median = _three_point_equality_weights(a)
    vand = (a[0] - a[1]) * (a[0] - a[2]) * (a[1] - a[2])
    common_abs = 4.0 * g * pc * math.sqrt(math.prod(e))
    curvature = common_abs * abs(vand)
    curvature_from_D2 = 4.0 * g * pc * math.sqrt(st.third_hankel_determinant)
    cscale = max(1.0, curvature, curvature_from_D2)
    if abs(curvature - curvature_from_D2) > 6.0e-12 * cscale:
        raise AssertionError("closed-triad curvature current lost D2 factorization")

    # The singleton median signed-curl slot has work magnitude equal to the
    # total gross positive (and negative) same-triad energy current.
    spread = max(a) - min(a)
    gross = common_abs * spread
    gross_upper = 2.0 * g * pc * E * math.sqrt(B) * chi
    gscale = max(1.0, gross, gross_upper)
    if abs(gross - gross_upper) > 6.0e-12 * gscale:
        raise AssertionError("gross-current barycentric factorization failed")
    global_upper = E * math.sqrt(B)
    if gross > global_upper + 8.0e-12 * max(1.0, global_upper):
        raise AssertionError("single-triad gross work exceeded E*sqrt(B)")

    return ThreePointCurrentLaw(
        total_energy=E,
        defect_energy=B,
        critical_determinant=st.critical_determinant,
        third_hankel_determinant=st.third_hankel_determinant,
        barycentric_efficiency=chi,
        sharp_efficiency_bound=1.0,
        waleffe_magnitude=g,
        phase_cosine_abs=pc,
        curvature_current_magnitude=curvature,
        gross_energy_current_magnitude=gross,
        gross_current_upper=gross_upper,
        global_gross_current_upper=global_upper,
        equality_weights=equality,
        median_index=median,
    )


def symmetric_heterochiral_upward_efficiency(radius_ratio: float) -> float:
    """Sharp amplitude-reduced log-upward factor on D=S=r, r in (1/2,1).

    This is a one-dimensional physical geometry slice, not a claim that the
    global two-parent variational problem has already been proved symmetric.
    """

    r = float(radius_ratio)
    if not math.isfinite(r) or not (0.5 < r < 1.0):
        raise ValueError("symmetric strict-upward ratio r must lie in (1/2,1)")
    return (
        math.sqrt(4.0 * r * r - 1.0)
        / (2.0 * math.sqrt(2.0) * r * (1.0 + r))
        * math.log(1.0 / r)
    )


def symmetric_heterochiral_stationarity_residual(radius_ratio: float) -> float:
    """Stationary equation for the preceding exact one-dimensional factor."""

    r = float(radius_ratio)
    if not math.isfinite(r) or not (0.5 < r < 1.0):
        raise ValueError("ratio r must lie in (1/2,1)")
    return (
        math.log(1.0 / r) * (1.0 + 2.0 * r - 4.0 * r**3)
        - (1.0 + r) * (4.0 * r * r - 1.0)
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "primitive_operator": "C=curl is self-adjoint; J_u v=P(u x v) is skew-adjoint and J_u u=0; NS is u_t=J_u(C-lambda)u-nu C^2u",
        "critical_determinant": "Delta=EZ-H^2=E||r||^2=sum_{a<b}(a-b)^2 E_a E_b is the exact pairwise signed-curl shear capacity",
        "pair_interaction": "for curl eigenstates Cu=a u, Cv=b v, the symmetrized Euler interaction is ((b-a)/2) P(u x v); equal curl eigenvalues do not interact",
        "second_residual": "h=Cr-(B/E)u-(<r,Cr>/B)r is orthogonal to u,r and D2=Delta||h||^2",
        "tangent_gradient": "h is the tangent projection of C^2u (equivalently of (C-lambda)^2u) to the common E,H level set",
        "commutator": "B'_NL=<r,[C,J_u]r>=2<Cr,u x r>=2 int r.S r",
        "critical_log_law": "(log Delta)'=2 int r.S r/B-2nu(Z/E+||Cr||^2/B)",
        "three_point_curvature": "on one closed three-node helical triad, |Q_curv|=4|g||cos Phi| sqrt(D2) <= 2 sqrt(D2)",
        "three_point_gross": "same-triad gross energy current is <=2|g| E sqrt(B)<=E sqrt(B), sharply over amplitudes; equality places one half of energy on the median signed-curl node",
        "hankel_meaning": "Dn is the squared curl-Krylov volume and equals a Vandermonde-squared average over n+1 signed-curl samples",
        "universal_observable_law": "for every spectral phi, production is the pairing of one Euler spectral velocity with the tangent residual of phi(C)u; on a three-node triad |M_phi'|=4|g||cosPhi|sqrt(det Gram(u,Cu,phi(C)u))",
        "martingale_grammar": "positive three-node spread is the unique signed-curl martingale transport; defect/enstrophy production is its quadratic variation and |C| production is its Tanaka reading",
        "fisher_speed": "the full curl-spectral Euler speed is A_spec=(1/4) sum S_a^2/E_a; on one triad A_spec=4|g|^2 cos^2(Phi) Delta",
        "phase_circle": "isolated triad total complex-amplitude speed is 4|g|^2 Delta and phase rotates it exactly between radial energy motion and quadrature shape motion; Im(A0 A1 A2 conj(g)) is invariant",
        "critical_correlation": "the |C| tangent volume equals Delta*Delta_rad*(1-rho(a,|a|)^2)/E; critical growth needs radial diversity and sign-radius decorrelation",
        "interaction_volume": "V_spec=Delta/A_spec is a PDE-defined physical-volume scale; Delta growth requires Re_spec=sqrt(A_spec)/(nu N sqrt(B))>2",
        "global_extension_guard": "isolated triads are self-contained one-dimensional martingale/phase systems, but the full PDE can coherently reorient overlapping triads through vertical shape motion; sustained collapse of V_spec is the remaining constitutive problem",
        "symmetric_transport_note": "the amplitude-reduced heterochiral symmetric slice has stationary ratio near 0.5981296, but global symmetry of the two-parent optimizer is not claimed here",
        "case_taxonomy_used": False,
        "temporal_matching_used": False,
        "global_regularity_claimed": False,
    }
