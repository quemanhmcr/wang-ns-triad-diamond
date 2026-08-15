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
    curvature = side * donor * child
    if progress > curvature + 8.0e-13 * max(1.0, progress, curvature):
        raise AssertionError("UV log progress exceeded signed-curl quadratic variation")
    return {
        "donor_work_per_common_current": donor,
        "child_work_per_common_current": child,
        "side_work_per_common_current": side,
        "child_log_progress_per_common_current": progress,
        "progress_upper_from_side": upper,
        "signed_curl_curvature_per_common_current_at_unit_child_radius": curvature,
        "log_progress_upper_from_curvature": curvature,
        "high_child_retained_fraction": child / donor,
        "side_fraction": side / donor,
    }


def modal_euler_action_decomposition(
    signed_frequencies: Sequence[float],
    modal_amplitudes: Sequence[complex],
    modal_sources: Sequence[complex],
) -> dict[str, object]:
    """Exact Pythagorean split of one complex Euler velocity.

    For ``z=sqrt(e) exp(i theta)`` and ``f=zdot``, each occupied helical mode obeys
    ``|f|^2=edot^2/(4e)+e*thetadot^2``.  Grouping the radial part by equal signed
    curl eigenvalue gives a second orthogonal split into total curl-spectral
    energy motion and redistribution inside a degenerate curl eigenspace.  An
    exactly zero amplitude carrying nonzero source is recorded as birth action.
    """
    if len(signed_frequencies) != len(modal_amplitudes) or len(modal_amplitudes) != len(modal_sources):
        raise ValueError("matching signed-frequency/amplitude/source data required")
    if not signed_frequencies:
        raise ValueError("nonempty modal data required")
    a = tuple(float(x) for x in signed_frequencies)
    z = tuple(complex(x) for x in modal_amplitudes)
    f = tuple(complex(x) for x in modal_sources)
    if not all(math.isfinite(x) for x in a):
        raise ValueError("finite signed frequencies required")
    if not all(math.isfinite(w.real) and math.isfinite(w.imag) for w in z + f):
        raise ValueError("finite complex amplitudes/sources required")

    total = math.fsum(abs(fi) ** 2 for fi in f)
    radial = 0.0
    phase = 0.0
    birth = 0.0
    groups: dict[float, list[float]] = {}
    for ai, zi, fi in zip(a, z, f):
        ei = abs(zi) ** 2
        if ei == 0.0:
            birth += abs(fi) ** 2
            continue
        q = zi.conjugate() * fi
        radial += q.real * q.real / ei
        phase += q.imag * q.imag / ei
        row = groups.setdefault(ai, [0.0, 0.0])
        row[0] += ei
        row[1] += 2.0 * q.real

    spectral = math.fsum(S * S / (4.0 * E) for E, S in groups.values() if E > 0.0)
    within = radial - spectral
    scale = max(1.0, total, radial, phase, birth, spectral, abs(within))
    if within < -2.0e-11 * scale:
        raise AssertionError("within-curl-eigenspace radial action lost nonnegativity")
    within = max(0.0, within)
    represented = spectral + within + phase + birth
    if abs(total - represented) > 3.0e-11 * scale:
        raise AssertionError("modal Euler Pythagorean action decomposition failed")
    keys = tuple(sorted(groups))
    return {
        "total_euler_action": total,
        "curl_spectral_action": spectral,
        "within_eigenspace_radial_action": within,
        "phase_rotation_action": phase,
        "new_amplitude_birth_action": birth,
        "represented_total_action": represented,
        "curl_spectral_fraction": 0.0 if total == 0.0 else spectral / total,
        "vertical_reconfiguration_action": within + phase + birth,
        "grouped_signed_frequencies": keys,
        "grouped_energies": tuple(groups[x][0] for x in keys),
        "grouped_energy_rates": tuple(groups[x][1] for x in keys),
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



def spectral_fitness_replicator_law(
    signed_frequencies: Sequence[float],
    modal_energies: Sequence[float],
    euler_energy_rates: Sequence[float],
    viscosity: float = 0.0,
) -> dict[str, object]:
    """One-scalar fitness form of the signed-curl energy law.

    On every occupied curl level ``a``, put ``f(a)=S_a/(2E_a)`` where ``S_a``
    is the actual Euler energy rate.  The two Euler null laws become
    ``int f d rho=int a f d rho=0`` and the Fisher action is ``int f^2 d rho``.
    Adding viscosity gives the exact unnormalized score ``2(f-nu a^2)`` and,
    after normalizing by total energy, the replicator score
    ``2(f-nu(a^2-N^2))``.
    """
    a, e = _spectral_data(signed_frequencies, modal_energies)
    if len(euler_energy_rates) != len(a):
        raise ValueError("one Euler energy rate per signed curl node required")
    rates = tuple(float(x) for x in euler_energy_rates)
    if not all(math.isfinite(x) for x in rates):
        raise ValueError("finite Euler energy rates required")
    nu = float(viscosity)
    if not math.isfinite(nu) or nu < 0.0:
        raise ValueError("finite nonnegative viscosity required")
    action = spectral_source_action(a, e, rates)["spectral_velocity_norm_squared"]
    E = math.fsum(e)
    Z = math.fsum(ai * ai * ei for ai, ei in zip(a, e))
    N2 = Z / E
    fit: list[float] = []
    raw: list[float] = []
    normalized: list[float] = []
    for ai, ei, si in zip(a, e, rates):
        if ei <= 0.0:
            if abs(si) > 5.0e-12 * max(1.0, *(abs(x) for x in rates)):
                raise ValueError("zero-energy node cannot carry first-order Euler energy rate")
            fi = 0.0
        else:
            fi = si / (2.0 * ei)
        fit.append(fi)
        raw.append(2.0 * (fi - nu * ai * ai))
        normalized.append(2.0 * (fi - nu * (ai * ai - N2)))
    mean_f = math.fsum(fi * ei for fi, ei in zip(fit, e))
    mean_af = math.fsum(ai * fi * ei for ai, fi, ei in zip(a, fit, e))
    scale = max(1.0, math.sqrt(max(0.0, action * E)), abs(mean_f), abs(mean_af))
    if abs(mean_f) > 8.0e-11 * scale or abs(mean_af) > 8.0e-11 * scale:
        raise AssertionError("Euler fitness lost the affine null laws")
    action_from_fit = math.fsum(ei * fi * fi for ei, fi in zip(e, fit))
    if abs(action-action_from_fit) > 5.0e-11 * max(1.0, action, action_from_fit):
        raise AssertionError("Euler fitness L2 norm lost Fisher action")
    normalized_mass_residual = math.fsum(ei * gi for ei, gi in zip(e, normalized)) / E
    if abs(normalized_mass_residual) > 8.0e-11 * max(1.0, *(abs(x) for x in normalized)):
        raise AssertionError("normalized spectral replicator lost unit mass")
    return {
        "fitness": tuple(fit),
        "unnormalized_log_energy_rates": tuple(raw),
        "normalized_log_probability_rates": tuple(normalized),
        "spectral_action": action,
        "fitness_action": action_from_fit,
        "fitness_energy_mean": mean_f,
        "fitness_helicity_mean": mean_af,
        "rms_curl_scale_squared": N2,
        "normalized_mass_residual": normalized_mass_residual,
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
    theta = state.beta1**2 / N2
    skew = state.alpha1**2 / N2
    sharp_threshold = 2.0 * math.sqrt(1.0 + theta + skew)
    return {
        "interaction_volume": V,
        "spectral_reynolds": math.sqrt(max(0.0, re2)),
        "spectral_reynolds_squared": max(0.0, re2),
        "fixed_state_log_delta_upper": fixed,
        "relaxed_log_delta_upper": relaxed,
        "defect_fraction": theta,
        "defect_skew_fraction": skew,
        "state_sharp_reynolds_threshold": sharp_threshold,
        "growth_requires_reynolds_above_two": True,
        "intrinsic_critical_mass": A / (N2 * B),
        "critical_mass_threshold": 4.0 * nu * nu,
    }


def critical_determinant_live_balance(
    state: CurlKrylovState,
    defect_strain_inner: float,
    spectral_velocity_norm_squared: float,
    viscosity: float,
) -> dict[str, float]:
    """Exact Krylov-impedance form of the critical determinant balance.

    ``defect_strain_inner=<h,F_E>`` is the actual productive Euler pairing and
    ``A_spec`` is the squared norm of the whole spectral Euler velocity.  Writing
    its q2 coefficient as ``c2`` exposes the exact dimensionless law

        (log Delta)'/(2 nu N^2)
          = gamma2 Re_spec y - (1+theta+chi+y^2),

    where ``y=beta2/N``.  Thus too little next-Krylov opening gives no nonlinear
    leverage, while too much is quadratically visible to viscosity.
    """
    s = float(defect_strain_inner)
    A = float(spectral_velocity_norm_squared)
    nu = float(viscosity)
    if not all(math.isfinite(x) for x in (s, A, nu)) or A < 0.0 or nu <= 0.0:
        raise ValueError("finite strain/action and positive viscosity required")
    if state.energy <= 0.0 or state.enstrophy <= 0.0 or state.defect_energy <= 0.0:
        return {
            "normalized_log_delta_rate": -math.inf,
            "q2_alignment": 0.0,
            "spectral_reynolds": 0.0,
            "krylov_opening": 0.0,
            "krylov_impedance": math.inf,
            "effective_drive": 0.0,
        }
    B = state.defect_energy
    N2 = state.enstrophy / state.energy
    N = math.sqrt(N2)
    beta2 = state.beta2
    hnorm = math.sqrt(B) * beta2
    if beta2 <= 0.0:
        if abs(s) > 5.0e-11 * max(1.0, abs(s)):
            raise ValueError("nonzero productive pairing with zero q2 Krylov opening")
        c2 = 0.0
    else:
        c2 = s / hnorm
    rootA = math.sqrt(A)
    gamma2 = 0.0 if rootA == 0.0 else c2 / rootA
    if abs(gamma2) > 1.0 + 3.0e-9:
        raise ValueError("productive q2 coefficient exceeds total spectral action")
    gamma2 = max(-1.0, min(1.0, gamma2))
    re_spec = rootA / (nu * N * math.sqrt(B))
    theta = state.beta1**2 / N2
    chi = state.alpha1**2 / N2
    y = beta2 / N
    baseline = 1.0 + theta + chi
    normalized = gamma2 * re_spec * y - (baseline + y * y)
    direct = critical_determinant_log_rate(state, s, nu)["log_delta_rate"] / (2.0 * nu * N2)
    if abs(normalized-direct) > 3.0e-10 * max(1.0, abs(normalized), abs(direct)):
        raise AssertionError("Krylov impedance factorization failed")
    impedance = math.inf if y <= 0.0 else y + baseline / y
    return {
        "normalized_log_delta_rate": normalized,
        "q2_coefficient": c2,
        "q2_alignment": gamma2,
        "spectral_reynolds": re_spec,
        "defect_fraction": theta,
        "defect_skew_fraction": chi,
        "krylov_opening": y,
        "krylov_impedance": impedance,
        "minimum_krylov_impedance": 2.0 * math.sqrt(baseline),
        "effective_drive": gamma2 * re_spec,
        "spectral_reconfiguration_action": max(0.0, A - c2 * c2),
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



def curl_nijenhuis_torsion_eigenfactor(
    parent_a: float, parent_b: float, child_c: float
) -> float:
    """Eigenvalue factor of the curl Nijenhuis torsion on one bracket output.

    If ``Cu=a u`` and ``Cv=b v``, then

        T_C(u,v)=[Cu,Cv]-C([Cu,v]+[u,Cv])+C^2[u,v]
               =(C-a)(C-b)[u,v].

    Hence the component at a child curl eigenvalue ``c`` carries exactly the
    factor returned here.  This is zero precisely when the bracket output stays
    on one of the two parent curl levels.
    """
    a, b, c = map(float, (parent_a, parent_b, child_c))
    if not all(math.isfinite(x) for x in (a, b, c)):
        raise ValueError("finite curl eigenvalues required")
    return (c - a) * (c - b)


def nijenhuis_curvature_from_root_work(
    parent_a: float,
    parent_b: float,
    child_c: float,
    child_energy_work: float,
) -> float:
    """Recover the full three-node quadratic-curvature current from one root.

    For the repository cyclic convention ``T_c=(a-b)R`` and

        Q=sum lambda_i^2 T_i,

    exact interpolation gives

        Q=(c-a)(c-b) T_c.

    The multiplier is exactly the child component of the curl Nijenhuis torsion.
    """
    t = float(child_energy_work)
    if not math.isfinite(t):
        raise ValueError("finite physical child work required")
    return curl_nijenhuis_torsion_eigenfactor(parent_a, parent_b, child_c) * t


def critical_escape_balance(
    signed_frequencies: Sequence[float],
    modal_energies: Sequence[float],
    modal_energy_rates: Sequence[float],
    viscosity: float,
) -> dict[str, float]:
    """Exact one-current balance for the critical ``K=sum |a|e`` stock.

    ``modal_energy_rates`` is the actual Euler signed-curl energy source and must
    satisfy the two affine null laws.  The source is compared with the one
    Fisher--Rao spectral speed, the tangent geometry of ``|a|``, and the exact
    cubic viscous moment.  No shell, owner, packet clock or temporal matching is
    introduced.
    """
    a, e = _spectral_data(signed_frequencies, modal_energies)
    if len(modal_energy_rates) != len(a):
        raise ValueError("one modal-energy rate per signed curl node required")
    rates = tuple(float(x) for x in modal_energy_rates)
    if not all(math.isfinite(x) for x in rates):
        raise ValueError("finite modal-energy rates required")
    nu = float(viscosity)
    if not math.isfinite(nu) or nu <= 0.0:
        raise ValueError("positive finite viscosity required")

    action = spectral_source_action(a, e, rates)["spectral_velocity_norm_squared"]
    split = critical_beltrami_split(a, e)
    corr = critical_tangent_correlation_geometry(a, e)
    gram = observable_tangent_gram(a, e, tuple(abs(ai) for ai in a))

    E = math.fsum(e)
    H = math.fsum(ai * ei for ai, ei in zip(a, e))
    Z = math.fsum(ai * ai * ei for ai, ei in zip(a, e))
    K = split["critical_stock"]
    B = Z - H * H / E
    delta = split["critical_determinant"]
    radial = split["radial_variance_component"]
    N2 = Z / E
    N = math.sqrt(max(0.0, N2))
    M3 = math.fsum(abs(ai) ** 3 * ei for ai, ei in zip(a, e))
    K_nl = math.fsum(abs(ai) * si for ai, si in zip(a, rates))
    # On either exact intrinsic boundary, |a| is affine on the occupied law:
    # one helicity sign gives |a|=+/-a, while zero radial variance gives
    # |a|=constant.  Euler annihilates both affine observables exactly.  Use
    # the positive spectral geometry rather than dividing two roundoff-scale
    # numbers to manufacture a spurious alignment cosine.
    if split["helicity_coexistence_component"] == 0.0 or radial == 0.0:
        K_nl = 0.0
    K_total = K_nl - 2.0 * nu * M3

    if B <= 0.0 or Z <= 0.0 or delta <= 0.0 or N <= 0.0:
        return {
            "critical_stock": K,
            "nonlinear_critical_rate": K_nl,
            "viscous_critical_rate": -2.0 * nu * M3,
            "total_critical_rate": K_total,
            "spectral_action": action,
            "spectral_reynolds": 0.0,
            "critical_alignment": 0.0,
            "defect_fraction": 0.0,
            "radial_variance_fraction": 0.0,
            "radial_coefficient_of_variation": 0.0,
            "decorrelation_factor": 0.0,
            "cubic_viscous_factor": math.inf,
            "live_escape_number": 0.0,
            "normalized_critical_rate": -math.inf,
            "radial_impedance_threshold": math.inf,
        }

    tangent_norm2 = gram["gram_determinant"] / delta
    source_envelope = 2.0 * math.sqrt(max(0.0, tangent_norm2 * action))
    if source_envelope == 0.0:
        gamma = 0.0
        if abs(K_nl) > 5.0e-10 * max(1.0, abs(K_nl)):
            raise AssertionError("critical source survived zero tangent/source envelope")
    else:
        gamma = K_nl / source_envelope
        if abs(gamma) > 1.0 + 2.0e-9:
            raise AssertionError("critical source alignment left the Hilbert cosine range")
        gamma = max(-1.0, min(1.0, gamma))

    theta = B / Z
    eta = radial / (E * Z)
    eta = max(0.0, min(1.0, eta))
    decor = corr["decorrelation_factor"]
    re_spec = math.sqrt(action) / (nu * N * math.sqrt(B))
    mu3 = M3 / (E * N**3)
    # Z^2 <= K M3, equivalently mu3 >= 1/sqrt(1-eta).
    impedance_lower = math.inf if eta >= 1.0 else 1.0 / math.sqrt(max(1.0e-300, 1.0 - eta))
    if mu3 + 2.0e-10 * max(1.0, mu3, impedance_lower) < impedance_lower:
        raise AssertionError("radial log-convex viscous impedance failed")
    radial_cv = math.inf if eta >= 1.0 else math.sqrt(max(0.0, eta / max(1.0e-300, 1.0 - eta)))
    effective_drive = gamma * re_spec * math.sqrt(max(0.0, theta * decor))
    radial_threshold = math.inf if radial_cv == 0.0 or not math.isfinite(radial_cv) else radial_cv + 1.0 / radial_cv
    radial_optimized_log_upper = 0.25 * effective_drive * effective_drive - 1.0

    live = gamma * re_spec * math.sqrt(max(0.0, theta * eta * decor))
    normalized = K_total / (2.0 * nu * E * N**3)
    represented = live - mu3
    one_minus_eta = max(0.0, 1.0 - eta)
    log_critical_per_parabolic_clock = (
        -math.inf if K <= 0.0 or one_minus_eta <= 0.0
        else normalized / math.sqrt(one_minus_eta)
    )
    if abs(normalized - represented) > 3.0e-9 * max(1.0, abs(normalized), abs(represented)):
        raise AssertionError("critical live-number factorization failed")
    if math.isfinite(log_critical_per_parabolic_clock):
        if log_critical_per_parabolic_clock > radial_optimized_log_upper + 5.0e-10 * max(1.0, abs(log_critical_per_parabolic_clock), abs(radial_optimized_log_upper)):
            raise AssertionError("radial-optimized critical log upper bound failed")

    return {
        "critical_stock": K,
        "nonlinear_critical_rate": K_nl,
        "nonlinear_source_envelope": source_envelope,
        "viscous_critical_rate": -2.0 * nu * M3,
        "total_critical_rate": K_total,
        "spectral_action": action,
        "critical_growth_action": gamma * gamma * action,
        "spectral_orthogonal_reconfiguration_action": max(0.0, (1.0 - gamma * gamma) * action),
        "spectral_reynolds": re_spec,
        "critical_alignment": gamma,
        "defect_fraction": theta,
        "radial_variance_fraction": eta,
        "radial_coefficient_of_variation": radial_cv,
        "decorrelation_factor": decor,
        "cubic_viscous_factor": mu3,
        "cubic_viscous_lower_from_radial_variance": impedance_lower,
        "effective_phase_geometry_drive": effective_drive,
        "radial_impedance_threshold": radial_threshold,
        "radial_optimized_log_critical_upper": radial_optimized_log_upper,
        "energy_clock_mean_abs_curl_log_upper": 0.25 * effective_drive * effective_drive,
        "live_escape_number": live,
        "normalized_critical_rate": normalized,
        "normalized_log_critical_rate_per_parabolic_clock": log_critical_per_parabolic_clock,
        "critical_growth_requires_reynolds_above_two": True,
    }








def critical_hilbert_square_balance(
    curvature_height: float,
    cubic_viscous_moment: float,
    radial_companion_norm_squared: float,
    viscosity: float,
) -> dict[str, float]:
    """Scalar audit of the exact critical Hilbert square completion.

    The operator identity uses

        A=|C|^(1/2) omega,
        B=|C|^(-1/2) P(u x |C|u),
        <A,B>=-kappa(0).

    Hence ``||A||^2=M3``, ``||B||^2=Q`` and

        K' = 2 kappa - 2 nu M3
           = -2nu ||A+B/(2nu)||^2 + Q/(2nu).

    This helper records only the scalar Gram data and refuses inputs which
    violate the Hilbert Cauchy constraint ``kappa^2<=M3 Q``.
    """

    kap, M3, Q, nu = map(
        float, (curvature_height, cubic_viscous_moment, radial_companion_norm_squared, viscosity)
    )
    if not all(math.isfinite(x) for x in (kap, M3, Q, nu)):
        raise ValueError("finite critical square data required")
    if M3 < 0.0 or Q < 0.0 or nu <= 0.0:
        raise ValueError("nonnegative square norms and positive viscosity required")
    if kap * kap > M3 * Q + 3.0e-12 * max(1.0, kap * kap, M3 * Q):
        raise ValueError("critical square data violate Hilbert Cauchy")
    square = M3 - kap / nu + Q / (4.0 * nu * nu)
    if square < -4.0e-12 * max(1.0, M3, abs(kap / nu), Q / (nu * nu)):
        raise AssertionError("completed critical Hilbert square became negative")
    square = max(0.0, square)
    direct = 2.0 * kap - 2.0 * nu * M3
    represented = -2.0 * nu * square + Q / (2.0 * nu)
    if abs(direct - represented) > 5.0e-12 * max(1.0, abs(direct), abs(represented)):
        raise AssertionError("critical Hilbert square completion failed")
    necessary_ratio = 0.0 if M3 == 0.0 else math.sqrt(Q / M3) / nu
    if direct > 0.0 and necessary_ratio <= 1.0:
        raise AssertionError("positive critical growth violated companion-norm necessity")
    return {
        "critical_rate": direct,
        "completed_square_norm_squared": square,
        "represented_rate": represented,
        "radial_companion_upper": Q / (2.0 * nu),
        "companion_to_viscous_norm_ratio": necessary_ratio,
    }


def critical_boost_logistic_bound(
    energy: float,
    critical_stock: float,
    curvature_height: float,
    cubic_viscous_moment: float,
    viscosity: float,
) -> dict[str, float]:
    """Exact critical rate plus its Jensen/Krein logistic upper.

    ``M3>=K^3/E^2`` and ``kappa=<y,S_u y>`` give

        K' <= 2 K (kappa/K - nu (K/E)^2).

    Consequently ``K'>0`` forces the actual critical boost Rayleigh rate
    ``kappa/K`` above ``nu (K/E)^2``.  No operator-norm replacement is made.
    """

    E, K, kap, M3, nu = map(
        float, (energy, critical_stock, curvature_height, cubic_viscous_moment, viscosity)
    )
    if not all(math.isfinite(x) for x in (E, K, kap, M3, nu)):
        raise ValueError("finite critical boost data required")
    if E <= 0.0 or K <= 0.0 or M3 < 0.0 or nu <= 0.0:
        raise ValueError("positive E,K,nu and nonnegative M3 required")
    lower = K**3 / (E * E)
    if M3 + 4.0e-12 * max(1.0, M3, lower) < lower:
        raise ValueError("cubic moment violates Jensen lower bound")
    direct = 2.0 * kap - 2.0 * nu * M3
    boost_rate = kap / K
    heat_rate = nu * (K / E) ** 2
    upper = 2.0 * K * (boost_rate - heat_rate)
    if direct > upper + 5.0e-12 * max(1.0, abs(direct), abs(upper)):
        raise AssertionError("critical boost logistic upper failed")
    return {
        "critical_rate": direct,
        "cubic_moment_lower": lower,
        "boost_rayleigh_rate": boost_rate,
        "quadratic_heat_rate": heat_rate,
        "logistic_upper": upper,
        "positive_growth_requires_boost_above_heat": direct <= 0.0 or boost_rate > heat_rate,
    }


def fixed_curl_cocycle_rhs(
    curl_eigenvalues: Sequence[float],
    amplitudes: Sequence[float],
    structure_triples: Sequence[tuple[int, int, int, float]],
    *,
    viscosity: float = 0.0,
) -> dict[str, object]:
    """Fixed-basis Cartan-tensor form of a finite curl-spectral NS system.

    In one real orthonormal curl eigenbasis ``Ce_i=lambda_i e_i`` put
    ``f_ijk=Omega(e_i,e_j,e_k)`` for ``i<j<k``.  The Euler RHS from one triple is

        F_i += -(lambda_k-lambda_j) f_ijk z_j z_k,
        F_j += +(lambda_k-lambda_i) f_ijk z_i z_k,
        F_k += -(lambda_j-lambda_i) f_ijk z_i z_j.

    This is exactly the contraction of the fixed alternating Cartan tensor with
    the signed-curl gaps.  Viscosity adds ``-nu lambda_i^2 z_i``.  The helper
    audits the two Euler affine invariants; Jacobi is an additional constraint
    on physical structure triples and is not manufactured here.
    """

    lam = tuple(float(x) for x in curl_eigenvalues)
    z = tuple(float(x) for x in amplitudes)
    if len(lam) != len(z) or not lam:
        raise ValueError("matching nonempty curl eigenvalues/amplitudes required")
    if not all(math.isfinite(x) for x in lam + z):
        raise ValueError("finite curl eigenvalues/amplitudes required")
    nu = float(viscosity)
    if not math.isfinite(nu) or nu < 0.0:
        raise ValueError("finite nonnegative viscosity required")
    n = len(z)
    euler = [0.0] * n
    for row in structure_triples:
        if len(row) != 4:
            raise ValueError("structure triple must be (i,j,k,f_ijk)")
        i, j, k, raw = row
        if not isinstance(i, int) or not isinstance(j, int) or not isinstance(k, int):
            raise ValueError("structure indices must be integers")
        if not (0 <= i < j < k < n):
            raise ValueError("structure indices must satisfy 0<=i<j<k<n")
        f = float(raw)
        if not math.isfinite(f):
            raise ValueError("finite Cartan structure coefficient required")
        li, lj, lk = lam[i], lam[j], lam[k]
        zi, zj, zk = z[i], z[j], z[k]
        euler[i] -= (lk - lj) * f * zj * zk
        euler[j] += (lk - li) * f * zi * zk
        euler[k] -= (lj - li) * f * zi * zj
    energy_rate = 2.0 * math.fsum(zi * fi for zi, fi in zip(z, euler))
    helicity_rate = 2.0 * math.fsum(li * zi * fi for li, zi, fi in zip(lam, z, euler))
    scale = max(1.0, math.sqrt(math.fsum(x * x for x in euler)))
    if abs(energy_rate) > 4.0e-12 * scale:
        raise AssertionError("fixed Cartan Euler RHS lost energy conservation")
    if abs(helicity_rate) > 4.0e-12 * max(scale, *(abs(x) for x in lam)):
        raise AssertionError("fixed Cartan Euler RHS lost helicity conservation")
    viscous = tuple(-nu * li * li * zi for li, zi in zip(lam, z))
    full = tuple(fi + vi for fi, vi in zip(euler, viscous))
    phase_divergence = -nu * math.fsum(li * li for li in lam)
    return {
        "euler_rhs": tuple(euler),
        "viscous_rhs": viscous,
        "full_rhs": full,
        "euler_energy_rate": energy_rate,
        "euler_helicity_rate": helicity_rate,
        "phase_space_divergence": phase_divergence,
        "euler_phase_space_divergence": 0.0,
    }


def sharp_helicity_flip_boost_geometry(
    advecting_radius: float, input_radius: float, output_radius: float
) -> dict[str, float]:
    """Sharp geometry bound for one critical helicity-flip commutator matrix element.

    Let ``q+l-k=0`` be one Fourier triangle.  In the critical coordinate
    ``y=|C|^(1/2)u``, the symmetric/Krein-boost matrix element from an input
    helicity ``s`` to output helicity ``-s`` generated by a unit helical velocity
    mode at radius ``q`` has magnitude at most

        2 sqrt(k l) |g|,

    with the low-mode helicity chosen for the larger Waleffe coupling.  Exact
    triangle geometry gives

        2 sqrt(k l)|g| <= (3 sqrt(6)/16) q.

    Same-helicity input/output is absent from the symmetric critical generator.
    """

    q, l, k = map(float, (advecting_radius, input_radius, output_radius))
    if not all(math.isfinite(x) for x in (q, l, k)) or min(q, l, k) <= 0.0:
        raise ValueError("positive finite Fourier radii required")
    if not (abs(l - k) < q < l + k):
        raise ValueError("strict nondegenerate Fourier triangle required")
    # Sum/difference coordinates are the numerically native positive chart.
    # Direct Heron factors catastrophically subtract when q << k,l.
    m = (k + l) / q
    x = (k - l) / q
    if abs(x) >= 1.0:
        raise AssertionError("strict triangle left |difference ratio|<1")
    area = 0.25 * q * q * math.sqrt(max(0.0, (m * m - 1.0) * (1.0 - x * x)))
    boost = (
        q
        * (1.0 + abs(x))
        * math.sqrt(max(0.0, (m * m - 1.0) * (1.0 - x * x) / (m * m - x * x)))
        / (2.0 * math.sqrt(2.0))
    )
    sharp_constant = 3.0 * math.sqrt(6.0) / 16.0
    upper = sharp_constant * q
    scale = max(1.0e-300, boost, upper)
    if boost > upper + 2.0e-12 * scale:
        raise AssertionError("critical helicity-flip Galilean null bound failed")
    return {
        "triangle_area": area,
        "sum_ratio": m,
        "difference_ratio": x,
        "max_helicity_flip_boost": boost,
        "sharp_upper": upper,
        "sharp_constant": sharp_constant,
        "normalized_boost": boost / q,
    }


def radial_fitness_selection_balance(
    signed_frequencies: Sequence[float],
    modal_energies: Sequence[float],
    modal_energy_rates: Sequence[float],
    viscosity: float,
) -> dict[str, float]:
    """Exact normalized mean-|curl| balance as Euler selection versus heat selection.

    Put ``p_a=E_a/E``, ``r_a=|a|``, ``m=E_p r=K/E`` and let
    ``f_a=S_a/(2E_a)`` be the actual Euler fitness on occupied curl levels.
    Then on the energy-loss clock ``d tau_E=2 nu N^2 dt``,

        d log m / d tau_E = R*x - I,

    where ``x=std_p(r)/m``, ``R=Cov_p(r,f)/(nu N^2 std_p(r))`` and
    ``I=Cov_p(r,r^2)/(N^2 m)``.  The positive moment inequality
    ``(E r^2)^2 <= (E r)(E r^3)`` gives ``I>=x^2`` and hence

        d log m / d tau_E <= R*x-x^2 <= R^2/4.

    This is an exact state/current law; ``R`` is not a new owner or event score.
    """

    a, e = _spectral_data(signed_frequencies, modal_energies)
    if len(modal_energy_rates) != len(a):
        raise ValueError("one Euler energy rate per signed-curl node required")
    rates = tuple(float(x) for x in modal_energy_rates)
    if not all(math.isfinite(x) for x in rates):
        raise ValueError("finite Euler energy rates required")
    nu = float(viscosity)
    if not math.isfinite(nu) or nu <= 0.0:
        raise ValueError("positive finite viscosity required")
    E = math.fsum(e)
    p = tuple(x / E for x in e)
    r = tuple(abs(x) for x in a)
    m = math.fsum(pi * ri for pi, ri in zip(p, r))
    N2 = math.fsum(pi * ri * ri for pi, ri in zip(p, r))
    if m <= 0.0 or N2 <= 0.0:
        return {
            "mean_absolute_curl": 0.0,
            "rms_curl_squared": N2,
            "radial_variance": 0.0,
            "radial_coefficient_of_variation": 0.0,
            "productive_fitness_score": 0.0,
            "viscous_radial_selection": 0.0,
            "normalized_log_mean_curl_rate": -math.inf,
            "quadratic_upper": 0.0,
        }
    # Exact Euler affine-null laws.  They are audited by spectral_source_action.
    spectral_source_action(a, e, rates)
    f = tuple(0.0 if ei == 0.0 else si / (2.0 * ei) for ei, si in zip(e, rates))
    fmean = math.fsum(pi * fi for pi, fi in zip(p, f))
    # Occupied zero-energy nodes have zero p and do not affect any covariance.
    var = math.fsum(pi * (ri - m) ** 2 for pi, ri in zip(p, r))
    sigma = math.sqrt(max(0.0, var))
    cov_rf = math.fsum(pi * (ri - m) * (fi - fmean) for pi, ri, fi in zip(p, r, f))
    M3 = math.fsum(pi * ri**3 for pi, ri in zip(p, r))
    cov_r_r2 = M3 - m * N2
    x = sigma / m
    R = 0.0 if sigma == 0.0 else cov_rf / (nu * N2 * sigma)
    I = cov_r_r2 / (N2 * m)
    lower = x * x
    scale = max(1.0, abs(I), abs(lower))
    if I + 3.0e-12 * scale < lower:
        raise AssertionError("quadratic radial viscous selection lower bound failed")
    exact = R * x - I
    first_upper = R * x - x * x
    quadratic_upper = 0.25 * R * R
    if exact > first_upper + 3.0e-12 * max(1.0, abs(exact), abs(first_upper)):
        raise AssertionError("radial selection exact rate exceeded moment upper")
    if first_upper > quadratic_upper + 3.0e-12 * max(1.0, abs(first_upper), abs(quadratic_upper)):
        raise AssertionError("completed-square radial selection upper failed")
    # Direct NS reconstruction for m=K/E.
    K_nl = math.fsum(ri * si for ri, si in zip(r, rates))
    if all(ai >= 0.0 for ai in a) or all(ai <= 0.0 for ai in a) or sigma == 0.0:
        # |a| is affine on the occupied support; remove roundoff from the exact null.
        K_nl = 0.0
    K = E * m
    Z = E * N2
    radial_delta = E * E * var
    curvature_height = 0.5 * K_nl
    curvature_score = 0.0 if radial_delta == 0.0 else curvature_height / (nu * N2 * math.sqrt(radial_delta))
    if abs(curvature_score - R) > 8.0e-10 * max(1.0, abs(curvature_score), abs(R)):
        raise AssertionError("productive fitness score disagrees with normalized curvature height")
    fisher_variance = math.fsum(pi * (fi - fmean) ** 2 for pi, fi in zip(p, f))
    productive_action = 0.0 if var == 0.0 else cov_rf * cov_rf / var
    curvature_action = 0.0 if radial_delta == 0.0 else curvature_height * curvature_height / radial_delta
    if abs(productive_action - curvature_action) > 8.0e-10 * max(1.0, productive_action, curvature_action):
        raise AssertionError("productive Fisher action disagrees with curvature-action quotient")
    if productive_action > fisher_variance + 5.0e-12 * max(1.0, productive_action, fisher_variance):
        raise AssertionError("radial productive Fisher action exceeded total Euler fitness variance")
    productivity_fraction = 0.0 if fisher_variance == 0.0 else productive_action / fisher_variance
    physical_log_upper = productive_action / (2.0 * nu * N2)
    Kdot = K_nl - 2.0 * nu * E * M3
    Edot = -2.0 * nu * Z
    direct = (Kdot / K - Edot / E) / (2.0 * nu * N2)
    if abs(direct - exact) > 8.0e-10 * max(1.0, abs(direct), abs(exact)):
        raise AssertionError("mean absolute curl selection balance reconstruction failed")
    return {
        "mean_absolute_curl": m,
        "rms_curl_squared": N2,
        "radial_variance": var,
        "radial_coefficient_of_variation": x,
        "euler_fitness_mean": fmean,
        "euler_radial_fitness_covariance": cov_rf,
        "nonlinear_critical_rate": K_nl,
        "curvature_height": curvature_height,
        "radial_critical_determinant": radial_delta,
        "productive_fitness_score": R,
        "normalized_curvature_height_score": curvature_score,
        "total_euler_fitness_variance": fisher_variance,
        "productive_fisher_action": productive_action,
        "productive_fisher_fraction": productivity_fraction,
        "physical_log_mean_curl_upper": physical_log_upper,
        "viscous_radial_selection": I,
        "viscous_selection_lower": lower,
        "normalized_log_mean_curl_rate": exact,
        "moment_upper": first_upper,
        "quadratic_upper": quadratic_upper,
        "direct_normalized_log_mean_curl_rate": direct,
    }


def parabolic_energy_clock_from_endpoints(
    initial_energy: float, final_energy: float, viscosity: float
) -> float:
    """Exact ``int N(t)^2 dt`` implied by the Navier--Stokes energy law.

    With ``N^2=Z/E`` and ``E'=-2 nu Z``, every interval on which both endpoint
    energies are positive satisfies

        int N^2 dt = log(E_initial/E_final)/(2 nu).

    This helper records the endpoint value; it does not manufacture a clock from
    an analysis partition.
    """

    E0, E1, nu = map(float, (initial_energy, final_energy, viscosity))
    if not all(math.isfinite(x) for x in (E0, E1, nu)):
        raise ValueError("finite energy endpoints and viscosity required")
    if E0 <= 0.0 or E1 <= 0.0 or E1 > E0 or nu <= 0.0:
        raise ValueError("require 0<E_final<=E_initial and positive viscosity")
    return math.log(E0 / E1) / (2.0 * nu)


def canonical_spectral_triple_source(
    left_frequency: float,
    median_frequency: float,
    right_frequency: float,
    triple_current: float,
) -> dict[str, float]:
    """One canonical curl-spectral three-current before any Fourier-edge split.

    For ordered signed-curl nodes ``a<m<b`` and the alternating physical current
    ``tau``, the induced Euler energy source is

        tau * (b-m, a-b, m-a).

    It annihilates the affine observables ``1`` and ``lambda`` identically.  If
    ``tau>0`` the median is the donor and the same vector is the barycentric
    martingale spread with donor mass ``tau*(b-a)``.
    """

    a, m, b, tau = map(
        float, (left_frequency, median_frequency, right_frequency, triple_current)
    )
    if not all(math.isfinite(x) for x in (a, m, b, tau)):
        raise ValueError("finite ordered curl nodes and triple current required")
    if not (a < m < b):
        raise ValueError("strictly ordered signed-curl nodes a<m<b required")
    source = (tau * (b - m), tau * (a - b), tau * (m - a))
    energy_residual = math.fsum(source)
    helicity_residual = math.fsum(x * y for x, y in zip((a, m, b), source))
    scale = max(1.0, *(abs(x) for x in source))
    if abs(energy_residual) > 4.0e-14 * scale:
        raise AssertionError("canonical triple source lost energy conservation")
    if abs(helicity_residual) > 4.0e-14 * max(scale, abs(a), abs(m), abs(b)):
        raise AssertionError("canonical triple source lost helicity conservation")
    donor = tau * (b - a)
    return {
        "left_source": source[0],
        "median_source": source[1],
        "right_source": source[2],
        "energy_residual": energy_residual,
        "helicity_residual": helicity_residual,
        "median_donor_mass": donor,
        "left_recipient_fraction": (b - m) / (b - a),
        "right_recipient_fraction": (m - a) / (b - a),
    }


def curl_spectral_bundle_base_velocity(
    signed_frequencies: Sequence[float],
    hellinger_amplitudes: Sequence[float],
    triple_coherences: Sequence[tuple[int, int, int, float]],
    *,
    viscosity: float = 0.0,
) -> dict[str, object]:
    """Exact base equation of the curl-spectral energy bundle.

    Write ``u=sqrt(E) sum_a q_a n_a`` with one unit state direction ``n_a`` in
    every occupied curl eigenspace and ``sum q_a^2=1``.  The supplied alternating
    coefficient on ``i<j<k`` is

        chi_ijk = sqrt(E) int n_i . (n_j x n_k).

    Holding the fiber directions fixed at the instant, Euler gives

        qdot_i += chi_ijk (a_k-a_j) q_j q_k,
        qdot_j += chi_ijk (a_i-a_k) q_i q_k,
        qdot_k += chi_ijk (a_j-a_i) q_i q_j.

    Viscosity adds exactly ``-nu(a_i^2-N^2)q_i`` after total-energy
    normalization.  No shell, packet, owner, or temporal matching is present.
    """

    a = tuple(float(x) for x in signed_frequencies)
    q = tuple(float(x) for x in hellinger_amplitudes)
    if len(a) != len(q) or not a:
        raise ValueError("matching nonempty curl/Hellinger coordinates required")
    if not all(math.isfinite(x) for x in a + q) or any(x < 0.0 for x in q):
        raise ValueError("finite curl nodes and nonnegative Hellinger amplitudes required")
    qnorm2 = math.fsum(x * x for x in q)
    if abs(qnorm2 - 1.0) > 2.0e-11:
        raise ValueError("Hellinger amplitudes must have unit squared norm")
    nu = float(viscosity)
    if not math.isfinite(nu) or nu < 0.0:
        raise ValueError("finite nonnegative viscosity required")

    euler = [0.0 for _ in q]
    n = len(q)
    for row in triple_coherences:
        if len(row) != 4:
            raise ValueError("each triple coherence must be (i,j,k,chi)")
        i, j, k, raw = row
        if not isinstance(i, int) or not isinstance(j, int) or not isinstance(k, int):
            raise ValueError("triple indices must be integers")
        if not (0 <= i < j < k < n):
            raise ValueError("triple indices must satisfy 0<=i<j<k<n")
        chi = float(raw)
        if not math.isfinite(chi):
            raise ValueError("finite triple coherence required")
        ai, aj, ak = a[i], a[j], a[k]
        qi, qj, qk = q[i], q[j], q[k]
        euler[i] += chi * (ak - aj) * qj * qk
        euler[j] += chi * (ai - ak) * qi * qk
        euler[k] += chi * (aj - ai) * qi * qj

    euler_energy_tangent = math.fsum(qi * vi for qi, vi in zip(q, euler))
    euler_helicity_tangent = math.fsum(ai * qi * vi for ai, qi, vi in zip(a, q, euler))
    escale = max(1.0, math.sqrt(math.fsum(x * x for x in euler)))
    if abs(euler_energy_tangent) > 3.0e-12 * escale:
        raise AssertionError("bundle Euler velocity left the energy sphere")
    if abs(euler_helicity_tangent) > 3.0e-12 * max(escale, *(abs(x) for x in a)):
        raise AssertionError("bundle Euler velocity left the helicity level")

    N2 = math.fsum(ai * ai * qi * qi for ai, qi in zip(a, q))
    viscous = [-nu * (ai * ai - N2) * qi for ai, qi in zip(a, q)]
    full = [x + y for x, y in zip(euler, viscous)]
    sphere_residual = math.fsum(qi * vi for qi, vi in zip(q, full))
    if abs(sphere_residual) > 3.0e-12 * max(1.0, math.sqrt(math.fsum(x * x for x in full))):
        raise AssertionError("normalized NS bundle velocity left the Hellinger sphere")
    return {
        "rms_curl_squared": N2,
        "euler_velocity": tuple(euler),
        "viscous_velocity": tuple(viscous),
        "full_velocity": tuple(full),
        "euler_hilbert_speed_squared": math.fsum(x * x for x in euler),
        "viscous_hilbert_speed_squared": math.fsum(x * x for x in viscous),
        "full_hilbert_speed_squared": math.fsum(x * x for x in full),
        "euler_energy_tangent_residual": euler_energy_tangent,
        "euler_helicity_tangent_residual": euler_helicity_tangent,
        "normalized_mass_residual": sphere_residual,
    }


def self_return_operator_geometry(
    energy: float,
    rms_curl_scale: float,
    horizontal_normalized_lamb_norm_squared: float,
    viscosity: float,
) -> dict[str, float]:
    """Operator form of the PDE-defined interaction volume and Reynolds number.

    For ``q0=u/sqrt(E)``, ``q1=(C-lambda)u/sqrt(B)`` and the curl-cyclic
    subspace ``K_u=closure{p(C)u}``, put

        g = P_K P(q0 x q1) = -C^{-1} P_K [q0,q1].

    The input ``horizontal_normalized_lamb_norm_squared`` is ``||g||_2^2``.
    Then exactly ``V_spec^{-1}=||g||^2`` and

        Re_spec = sqrt(E)||g||/(nu N).
    """

    E = float(energy)
    N = float(rms_curl_scale)
    g2 = float(horizontal_normalized_lamb_norm_squared)
    nu = float(viscosity)
    if not all(math.isfinite(x) for x in (E, N, g2, nu)):
        raise ValueError("finite self-return geometry required")
    if E <= 0.0 or N <= 0.0 or g2 < 0.0 or nu <= 0.0:
        raise ValueError("positive energy/scale/viscosity and nonnegative self-return norm required")
    V = math.inf if g2 == 0.0 else 1.0 / g2
    re = math.sqrt(E * g2) / (nu * N)
    return {
        "interaction_volume": V,
        "inverse_interaction_volume": g2,
        "spectral_reynolds": re,
        "spectral_reynolds_squared": re * re,
        "intrinsic_critical_mass": E * g2 / (N * N),
    }

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
        "krylov_impedance": "exactly (log Delta)'/(2nuN^2)=gamma2 Re_spec y-(1+theta+chi+y^2); Delta growth requires gamma2 Re_spec>y+(1+theta+chi)/y>=2sqrt(1+theta+chi)",
        "three_point_curvature": "on one closed three-node helical triad, |Q_curv|=4|g||cos Phi| sqrt(D2) <= 2 sqrt(D2)",
        "three_point_gross": "same-triad gross energy current is <=2|g| E sqrt(B)<=E sqrt(B), sharply over amplitudes; equality places one half of energy on the median signed-curl node",
        "hankel_meaning": "Dn is the squared curl-Krylov volume and equals a Vandermonde-squared average over n+1 signed-curl samples",
        "universal_observable_law": "for every spectral phi, production is the pairing of one Euler spectral velocity with the tangent residual of phi(C)u; on a three-node triad |M_phi'|=4|g||cosPhi|sqrt(det Gram(u,Cu,phi(C)u))",
        "martingale_grammar": "positive three-node spread is the unique signed-curl martingale transport; defect/enstrophy production is its quadratic variation and |C| production is its Tanaka reading",
        "fisher_speed": "the full curl-spectral Euler speed is A_spec=(1/4) sum S_a^2/E_a; on one triad A_spec=4|g|^2 cos^2(Phi) Delta",
        "spectral_fitness": "on occupied curl levels f(a)=S_a/(2E_a), so rho_t=2(f-nu a^2)rho, int f d rho=int a f d rho=0 and A_spec=int f^2 d rho",
        "lamb_ontology": "L=u x curl u is the one nonlinear activity field: Helmholtz projection gives pressure reaction plus Euler velocity; the latter splits into curl-spectral and isospectral shape motion, while curl L=-[u,curl u] and curl Nijenhuis curvature gives vortex stretching",
        "phase_circle": "isolated triad total complex-amplitude speed is 4|g|^2 Delta and phase rotates it exactly between radial energy motion and quadrature shape motion; Im(A0 A1 A2 conj(g)) is invariant",
        "critical_correlation": "the |C| tangent volume equals Delta*Delta_rad*(1-rho(a,|a|)^2)/E; critical growth needs radial diversity and sign-radius decorrelation",
        "interaction_volume": "V_spec=Delta/A_spec is a PDE-defined physical-volume scale; Delta growth requires Re_spec>2 and, state-sharply, Re_spec>2sqrt(1+theta+chi)",
        "nijenhuis_curvature": "with the divergence-free Lie bracket, C J_u v=-[u,v]; curl Nijenhuis torsion is T_C, K_C(x,y)=T_C(Cx,y)-T_C(x,Cy), and Z'_NL=B'_NL=(1/3)<C^-1 u,K_C(u,u)>=2 int r.Sr",
        "critical_live_balance": "for K=int|a|d rho, K'/(2nu E N^3)=gamma Re_spec sqrt(theta eta d)-mu3, with mu3>=1/sqrt(1-eta); hence K'>0 requires gamma Re_spec sqrt(theta d)>x+1/x>=2 where x=sqrt(eta/(1-eta))",
        "uv_curvature": "for every strict heterochiral UV spread, M^2 times true high-child log progress is <= the same positive Nijenhuis/enstrophy quadratic-variation current; constant 1 is sharp only at the degenerate boundary",
        "lamb_participation": "V_spec>=Delta/||u x curl u||_2^2>=||u x curl u||_1^2/||u x curl u||_2^2, so collapse of the PDE interaction volume requires physical concentration of the one Lamb activity field",
        "spectral_monotonicity_guard": "fixed curl-energy magnitudes with a pi triad-phase reversal flip every non-affine observable current; no spectral-only monotone beyond affine energy/helicity can contain the missing mechanism",
        "spectral_three_current": "before any Fourier-triad split, tau_abc=-2 Omega(u_a,u_b,u_c) reconstructs the full signed-curl node source tau(c-b,a-c,b-a); every spectral observable is its determinant/divided-difference readout",
        "spectral_bundle": "u=sqrt(E) sum q_a n_a gives qdot_a=sum chi_abc(c-b)q_bq_c-nu(a^2-N^2)q_a; the fiber direction uses the orthogonal remainder of the same projected cross product and viscosity does not rotate a curl eigenspace direction",
        "self_return_operator": "for q0=u/sqrt(E), q1=(C-lambda)u/sqrt(B), ell=P(q0 x q1)=-C^-1[q0,q1] and K_u=closure{p(C)u}, V_spec^-1=||P_K ell||^2 and Re_spec=sqrt(E)||P_K ell||/(nu N)",
        "cartan_jacobi": "B(x,y)=<C^-1x,y> is ad-invariant, Omega(x,y,z)=B(x,[y,z]) is a closed 3-cocycle, and the reduced K_u Jacobiator is exactly supplied by the vertical closure defect Q[x,y]",
        "lax_material_heat": "Euler obeys d/dt ad_omega=[ad_omega,ad_u]; after material pullback the NS vorticity 2-form obeys pure heat partial_t Phi^*varpi=nu Delta_(Phi^*g) Phi^*varpi while partial_t(Phi^*g)=2 Phi^*S and det(Phi^*g)=1",
        "parabolic_energy_clock": "-E'/E=2nuN^2, so int N^2 dt=log(E_initial/E_final)/(2nu) on positive-energy intervals; critical escape therefore requires unbounded persistence of the dimensionless productive self-return factor, not an event count",
        "master_commutator": "for every spectral observable M_phi=<u,phi(C)u>, M_phi'_NL=2<u,[phi(C),J_u]Cu>; spectrally [phi(C),J_u] is the divided-difference functional calculus of the one primitive commutator [C,J_u]",
        "critical_krein_boost": "with y=|C|^(1/2)u and J=sgn(C), L_u=|C|^(1/2)J_u J|C|^(1/2) is J-skew; its skew J-even part only rotates sectors, its self-adjoint J-odd part S_u=(1/2)|C|^(1/2)[J_u,J]|C|^(1/2) alone gives kappa(0)=<y,S_uy>",
        "normalized_curvature_score": "R_*=kappa(0)/(nu N^2 sqrt(EZ-K^2))=Cov_p(|a|,f)/(nu N^2 Std_p(|a|)); d log(K/E)/d tau_E=R_*x-I_r <=R_*x-x^2<=R_*^2/4",
        "galilean_null": "one helicity-flip critical boost matrix element obeys 2sqrt(k l)|g(q,l,-k)| <= (3sqrt(6)/16)|q| sharply; uniform/near-uniform sweeping therefore loses critical boost linearly in its advecting wave number",
        "critical_hilbert_square": "with R_u=P(u x |C|u), A=|C|^(1/2)curl u and B=|C|^(-1/2)R_u, K'=-2nu||A+B/(2nu)||^2+||B||^2/(2nu); the companion norm threshold is necessary only, while its Hilbert angle remains decisive",
        "critical_boost_logistic": "kappa(0)=<y,S_uy> and M3>=K^3/E^2 give K'<=2K(kappa(0)/K-nu(K/E)^2); any instant of positive K growth therefore requires the actual Krein boost Rayleigh rate above the quadratic mean-curl heat rate",
        "fixed_cartan_ode": "in one fixed curl eigenbasis the full Euler RHS is zI'=-(1/2)sum_JK(lambdaK-lambdaJ)f_IJK zJ zK with one time-independent alternating Cartan tensor f=Omega; NS adds only -nu lambdaI^2 zI and finite Galerkin Euler is phase-space divergence free",
        "productive_fisher_action": "A_prod=kappa(0)^2/(EZ-K^2)=Cov_p(|a|,f)^2/Var_p(|a|)<=A_spec/E; d log(K/E)/dt<=A_prod/(2nuN^2), so critical mean-curl escape requires divergence of the single curvature/Fisher action integral",
        "persistence_monotonicity_guard": "neither higher Krylov volumes nor the normalized curvature score R_* are universal monotones; explicit spectral and Galerkin falsifications force the remaining theorem to control persistence/action rather than instantaneous sign",
        "formalism_guard": "the physical alternating Euler triple generator is exact, but a direct-sum counterexample rejects promoting it to an unrestricted Nambu-Poisson bracket; Cartan/Jacobi compatibility is used instead",
        "hminus1_cartan_gradient": "in the curl H^-1 metric the Euler Cartan motion is exactly orthogonal to grad(E/2)=C^2u; ||u_t||_-1^2=||F_E||_-1^2+nu^2 Z and E/2 descends only through heat",
        "master_sobolev_strain_transfer": "for G_s=|C|^s J_u C |C|^-s and A=P S=(1/2)[C,J_u], Sym G_s is the hyperbolic functional calculus sinh((s-1/2)D)/sinh(D/2) on helicity-even A and -cosh((s-1/2)D)/cosh(D/2) on helicity-odd A, D=ad_log|C|",
        "critical_sech_strain": "at the unique midpoint s=1/2, all helicity-even strain vanishes and Sym G_1/2=-sech((1/2)ad_log|C|) A_odd; equivalently it is the two-sided Poisson/Sylvester average of physical odd strain",
        "critical_loggap_locality": "the critical sech map is a contraction and a cross block separated by log gap L obeys the dimension-free bound min(1,csch(L/2)); deep signed-curl separation is therefore continuously suppressed without shell thresholds",
        "poisson_scale_law": "dmu(t)=2||exp(-t|C|)|C|u||^2 dt/K is a canonical probability scale law with mean E/(2K)=1/(2 mean|curl|); the critical boost Rayleigh rate is the mu-average of Poisson-smoothed odd-strain Rayleigh quotients",
        "poisson_heat_subordination": "exp(-t|C|)=t/(2sqrt(pi)) int s^(-3/2) exp(-t^2/(4s)) exp(-sC^2) ds, so the scale filter intrinsic to critical Euler growth is subordinate to the same quadratic heat generator used by viscosity",
        "radial_resolvent_square": "for m=K/E and s_m=(|C|-m)u, m'=2/E(<s_m,F_spec>-nu<(s_m),(|C|+m)s_m>); completing the exact (|C|+m) Hilbert square gives a sharper resolvent productive-action upper than the variance-only Fisher bound",
        "viscous_orientation_guard": "pure viscosity rescales each curl eigenspace amplitude but does not rotate its normalized direction, so normalized Cartan triple orientation has zero viscous derivative; no proof may rely on fictitious heat dephasing",
        "resolvent_gap_guard": "the weighted radial resolvent Cauchy fraction has no universal gap below one; closed three-node currents can approach saturation on extreme amplitude states",
        "critical_operator_heat": "because the midpoint strain operator Sigma_c(u) is linear and translation covariant, its full NS evolution is Sigma_c,t=Sigma_c(F_E)-nu sum_j[D_j,[D_j,Sigma_c]]; the operator heat semigroup is Gaussian averaging by translation conjugations and is norm-contractive",
        "helicity_polar_modulus": "C=J|C| is the polar decomposition of the conserved helicity form; K=<u,|C|u> is its canonical positive modulus, frozen Euler G_u=J_uC is C-skew, and the critical cocycle is J-pseudo-unitary with reciprocal singular-value pairs in finite Galerkin dimension",
        "material_strain_metric": "on divergence-free L2, Sym(ad_u)=-P S(u); physical strain is exactly the metric-deformation part of material Lie transport, and the master Sobolev law is this same deformation read in interpolating curl metrics",
        "uniform_sweeping_conjugation": "for constant U, Sigma_c(-(U.grad)u)=-i[U.D,Sigma_c(u)]; uniform sweeping only unitarily conjugates the dangerous midpoint operator and cannot regenerate its norm",
        "objective_strain_regeneration": "Sigma_c(F_E) is the fixed midpoint filter of -(u.grad S+S^2+Omega^2+Hess p); pressure/corotation reorient strain but satisfy int S:Hess p=0 and S:[S,Omega]=0, while global strain amplitude obeys the Betchov cubic law",
        "critical_enstrophy_sign_guard": "actual dealiased Galerkin states realize all four sign quadrants of (K'_NL,Z'_NL), so curvature height kappa(0) and curvature area/int enstrophy current cannot be substituted for one another",
        "critical_spectral_hminus_half_action": "with F_spec the actual curl-spectral Euler velocity, the K balance completes an exact |C|^-1/2 spectral Hilbert square; critical blow-up therefore requires divergence of the Galilean-clean critical probability action E_pi[(f/|a|)^2]",
        "master_sobolev_hilbert_square": "for every s, C^2u is grad(K_s/2) in the shifted H^(s-1) metric and K_s completes the exact square with B_s=|C|^(s-1)F_spec; the normalized action is always E_pi_s[(f/|a|)^2], so the whole Sobolev hierarchy uses one local escape currency",
        "isolated_triad_scale_free_action": "every physical closed triad obeys the sharp root law 8|g|^2(a_j-a_k)^2<=(r_i)(r_j+r_k), hence |||C|^-1/2 F_spec,tri||^2 <= (1/2)E_tri K_tri and (log K_tri) <= E_tri/(4nu) independently of absolute Fourier scale; only coherent cross-triad action can evade this isolated ceiling",
        "critical_productive_volume": "Vprod=E K M3/(2 kappa(0)^2) is a representation-free productive physical volume, Vnu=E K/(2 nu^2 M3) is the matching viscous volume, and exactly K'/(2nuM3)=sgn(kappa)sqrt(Vnu/Vprod)-1; positive critical growth is precisely the race Vprod<Vnu",
        "critical_volume_scaling": "under NS dilation all of V_action, Vprod and Vnu scale as lambda^-3 while sqrt(Vnu/Vprod)=|kappa|/(nu M3) is invariant; finite-time critical escape requires divergence of the scale-invariant action integral int E/Vprod dt",
        "critical_volume_heat_guard": "neither total nor productive critical action volume is a heat Lyapunov scalar; pure heat can increase or decrease normalized productivity through amplitude reweighting even though the midpoint operator carrier itself obeys contractive operator heat",
        "critical_reynolds_operator": "with z=|C|^(3/2)u and Rc=nu^-1 |C|^-1 Sigma_c |C|^-1, Rc is self-adjoint/helicity-odd and the exact critical balance is 2nu times the Rc-minus-identity Rayleigh form; heat is identity in this critical dissipation coordinate",
        "reynolds_helicity_neutrality": "every nonzero eigenvector of the self-adjoint helicity-odd critical Reynolds operator has zero dissipation-coordinate helicity; the actual Rayleigh quotient is bounded by ||Rc|| sqrt(1-(H3/M3)^2)",
        "reynolds_heat_dirichlet": "at a top eigenvector v, <v,Delta_op Rc v>=2 sum_j <D_jv,(lambda_max I-Rc)D_jv>>=0; operator heat erodes the upper Reynolds edge exactly by translation leakage out of its top eigenspace",
        "continuum_midpoint_hs_isometry": "on continuum R3 with the repository unitary Fourier convention, Qc=|C|^-1 Sigma_c |C|^-1 is exactly a 1/8 Hilbert-Schmidt isometric embedding of Hdot^1/2: ||Qc(u)||_HS^2=K/64; it intertwines C^2 with Delta_op, so M3/64 is the operator Dirichlet form and kappa(0)/64 is the HS pairing with Qc(F_E)",
        "reynolds_spectral_capacity": "on continuum R3, K=128 nu^2 sum_{lambda_j>0} lambda_j(Rc)^2; hence ||Rc||op<=sqrt(K)/(8sqrt(2)nu), the number of positive eigenvalues above one is <=K/(128nu^2), and K<128nu^2 is an invariant no-critical-growth radius on smooth intervals",
        "continuum_midpoint_lossless": "the continuum midpoint transform T:u->Qc(u) satisfies T* T=(1/64)|C| and is injective on mean-zero critical states, so u=64|C|^-1 T*Qc; physical heat and nonlinear source preserve its Hilbert-Schmidt image",
        "midpoint_spectral_polynomial_guard": "Qc is self-adjoint/helicity-odd, hence every defined odd spectral trace Tr Qc^(2m+1) vanishes; kappa(0) is generally nonzero, so critical regeneration cannot be reduced to eigenvalues or a cubic Riccati trace of Qc alone",
        "continuum_midpoint_sobolev_scale": "the continuum midpoint transform intertwines |C|^alpha with Delta_op^(alpha/2), giving ||u||_Hdot^s^2=64||Delta_op^(s/2-1/4)Qc(u)||_HS^2; s=1/2 is uniquely the plain HS energy level, with E and Z at symmetric operator powers -1/4 and +1/4",
        "primitive_vacuum_state_chain": "the scalar unit is a canonical basepoint: Q*1=alpha and (Q*)^2 1=nu beta; with H=QQ*+Q*Q, [H,Q*]1=nu Q beta=nu e, so projected NS is partial_t(Q*1)=-(1/nu)P[H,Q*]1",
        "primitive_metric_current_law": "for every positive curl functional metric G=f(C), one-half d<alpha,G alpha>/dt=-(1/nu)||G^(1/2)(Q*)^2 1||^2-(1/nu)<[Q*,G]Q*1,(Q*)^2 1>; energy is one-way because [Q*,I]=0, while every two-way Sobolev motion is exactly the current/metric commutator",
        "primitive_critical_two_way_channel": "with Lambda=|C|, A=[Lambda^-1,Q*] is the normalized critical two-way channel and [Q*,Lambda]=Lambda A Lambda; its anticommutator with the same current creator is {A,Q*}=nu B, B=[Lambda^-1,beta wedge]",
        "primitive_critical_channel_isometry": "on continuum R3 with unitary Fourier convention and the full graded exterior HS norm, ||A||^2=(2/pi^2)K and ||B||^2=(1/pi^2)M3; both come from the exact raw translation integral 4pi with exterior creation multiplicities 4 and 2",
        "primitive_critical_channel_derivative": "B=i sum_j (dx^j wedge)[D_j,A], hence sum_j||[D_j,A]||^2=2||B||^2=(2/nu^2)||{A,Q*}||^2; critical heat is exactly the squared current-incompatibility derivative of the two-way channel",
        "primitive_critical_channel_rigidity": "{A,Q*}=0 forces B=0 and therefore M3=0 by the exact continuum isometry; on mean-zero/decaying states this implies u=0.  This is zero-self-frustration rigidity, not a uniform instantaneous gap",
        "primitive_critical_channel_dynamics": "A is linear and translation covariant, so A_t=A(F_E)-nu Delta_op A and one-half d||A||^2/dt=(2/pi^2)kappa(0)-2nu||B||^2=(2/pi^2)kappa(0)-(2/nu)||{A,Q*}||^2",
        "primitive_critical_pair_area_loop": "with E=alpha wedge and R=Lambda^-1, V=E R E is forced by E^2=0 and {A,E}=0; on continuum R3 ||V||_HS,gr^2=pi^-2||u(x) cross u(y)||_pair^2 and kappa(0)=-pi^2 Re<B,V>, equivalently the nonlinear critical work is one orientation/current cubic loop rather than an eigenvalue-only trace",
        "primitive_critical_gauss_residual": "Gc=V+2nu B is exactly twice the Hom-connection current nabla A between D0=nu d and D1/2=nu d+(1/2)E; K'=pi^2(||V||^2-||Gc||^2)/(2nu), equivalently K'=(||u(x) cross u(y)||_pair^2-||u(x) cross u(y)-2nu(omega(x)-omega(y))||_pair^2)/(2nu)",
        "primitive_critical_gauss_bianchi": "the Hom curvature is nabla^2 A=(nu/2)(beta wedge)A, hence the exact residual deciding K' obeys nabla Gc=nu(beta wedge)A; its physical kernel is omega(x).(u(y)-u(x))/(2pi^2|x-y|^2)",
        "primitive_critical_gauss_null_geometry": "exact vanishing of (beta wedge)A means omega(x) is orthogonal to the affine velocity-value span for every x; span dimension 3 is irrotational, dimension 2 is embedded 2D, dimension 1 is a shear u=phi b with (u.grad)u=0, and dimension 0 is uniform flow",
        "primitive_critical_six_dimensional_gauss_floor": "top-degree Hom tests give the quadratic-form normal operator L6=-nu^2(Delta_x+Delta_y)+|u(x)|^2/4 and ||Gc||^2>=<nu(beta wedge)A,L6^-1 nu(beta wedge)A>; the six variables are the two points already forced by |D|, not extra physical dimensions",
        "primitive_critical_gauss_stability_guard": "the Gauss floor gives exact control of the source in its native L6^-1 metric and the exact zero source is lower-dimensional, but quantitative near-zero proximity to the 2D/shear/null manifold is a separate stability/compactness theorem and is not assumed",
        "primitive_critical_carre_du_champ": "Gamma_u,ij=u_i Lambda u_j+u_j Lambda u_i-Lambda(u_i u_j)=pi^-2 int delta u_i delta u_j/|x-y|^4 dy is positive semidefinite and int tr Gamma_u=2K; critical H^1/2 size is the trace mass of the intrinsic Lambda increment metric",
        "primitive_critical_vorticity_metric": "||F A||_HS^2=(4pi^2)^-1 int omega^T Gamma_u omega; ker Gamma_u is the orthogonal complement of the affine velocity-value span, so the Gauss source is vorticity read in the same positive metric whose trace is critical size",
        "primitive_critical_metric_heat_law": "with L=Lambda^2=-Delta, (partial_t+nu L)Gamma_u=Gamma(F_E,u)+Gamma(u,F_E)-2nu sum_j Gamma_{partial_j u}; the sink is Loewner-positive and integrated trace is exactly K'=2kappa(0)-2nu M3",
        "primitive_critical_rank_persistence": "Cauchy-Binet makes e2(Gamma_u) the squared two-increment area and det Gamma_u the squared three-increment volume; exact null is fixed 2D/shear/irrotational geometry, finite-energy R3 has only vacuum, and periodic 2D/shear null classes are NS-invariant rather than reset surfaces",
        "primitive_critical_stretching_bridge_guard": "omega(x).(u(x+h)-u(x))=(S omega)(x).h+O(|h|^2), but the static Hdot^-1/2 bridge is withdrawn: near-parallel two-mode carriers k_N=(N,0,0), ell_N=(N,1,0) produce scale-growing exact Fourier ratios while localized projected packets retain the low-beat mechanism; the remaining bridge must use endogenous Gamma_u heat/turning history",
        "primitive_actual_current_operator": "the endogenous state operator is Q=nu delta+i_u with Q*=nu d+u^flat wedge; the physical electromotive current is e=Q beta, the projected momentum law is u_t=-P(Q beta), and therefore Q_t=-i_{P(Q beta)}",
        "primitive_current_curvature_square": "exterior Leibniz gives {delta,i_u}=(beta wedge)^*, hence Q^2=nu(beta wedge)^* on the whole graded exterior algebra and Q^2 beta=nu|beta|^2; vorticity is the failure of nilpotence of the actual NS current operator",
        "primitive_current_nilpotent_chords": "Q(u)-Q(v)=i_{u-v}, every chord/tangent contraction squares to zero, and any two contractions anticommute; the non-nilpotence of Q comes only from interaction of fixed delta with the physical contraction, not from nonlinear state-space chords",
        "primitive_full_hodge_gradient": "for H=Q*Q+QQ*, H^(1)=(-nu^2 Delta+|u|^2)I+2nu S and star^-1 H^(2) star=(-nu^2 Delta+|u|^2)I-2nu S; together with Q^2 this recovers the symmetric and antisymmetric parts of grad u from one Q algebra",
        "primitive_intertwining_lamb": "for every closed eta=*b^flat, H^(1)(Q eta)-Q(H^(2)eta)=[Q*,Q^2]eta=nu^2 d(omega.b)+nu[b cross (u cross omega)]^flat; at b=omega the non-exact term locally reconstructs the Lamb field away from omega=0",
        "primitive_finite_current_chains": "Q* alpha=nu beta gives alpha -> nu beta -> nu e -> nu^2|beta|^2 ->0 and energy dissipation is adjointness of its first two arrows; on top degree dV -> *alpha -> nu omega^flat -> nu(u.omega) ->0 and Q^4=0 in three dimensions",
        "primitive_midpoint_status": "M=nu delta+(1/2)i_u is the canonical midpoint between pure Hodge heat and the actual operator Q, selected by the critical/Poynting reflection; it is a reading rather than the fundamental state operator and satisfies M(M beta)=(nu/2)|beta|^2",
        "primitive_turning_frontier": "the previous persistent critical near-kernel frontier is sharpened by Gc=2 nabla A, nabla Gc=nu(beta wedge)A and the positive Gamma_u; static Sobolev coercion is excluded by the near-parallel low-beat referee, so blocking infinite productive regeneration leaves an unproved quantitative stability/history bridge from persistent small native L6^-1 residual through the endogenous Gamma_u heat/turning law",
        "graded_current_strain_parent": "with q=nu delta+(1/2)i_u, the positive Hodge-Dirac square q^*q+qq^* is (-nu^2 Delta+|u|^2/4)I+nu S on one-forms and the same scalar part minus nu S on two-forms; physical strain, hence the critical sech/Reynolds operator, is the degree imbalance of one midpoint current complex",
        "critical_gram_self_frustration": "for G_u b=u cross b-2nu curl b, the polarized Poynting law gives 4nu(PS)_odd=[P(U_u^*U_u-G_u^*G_u)P]_odd; the dangerous critical block is exactly a helicity-odd Gram imbalance, but odd projection destroys simple positivity so no free negative-square contraction follows",
        "critical_curvature_floor": "the mixed Gauss law (div-u/(2nu).)G_u b=omega.b implies G_u^*G_u>=M_omega^*(-Delta+|u|^2/(4nu^2))^-1M_omega in quadratic-form order; a persistent expanding Reynolds direction must continually avoid this state-generated curvature floor at the intrinsic Poisson depths",
        "critical_projected_gauss_tax": "the endpoint-current form of the s=1/2 Hilbert square is K'=(||Lambda^-1/2 F_E||^2-||Lambda^-1/2 P G||^2)/(2nu), while the Euler energy null forces ||Lambda^-1/2 P G||^2>=4nu^2 Z^2/K; positive action and compulsory residual tax are measured in the same critical metric",
        "critical_two_null_gauss_tax": "helicity conservation adds the second normal direction Lambda^(1-s)Cu to the endpoint-current residual; the resulting two-by-two Gram projection gives the exact sharpened two-null lower recorded in material Section 47 and introduces no new stock or branch",
        "graded_current_persistence_guard": "the new current algebra does not prove large-data regularity: static Schrödinger coercivity enters the exact/pressure sector and helicity-odd projection loses Gram order; the missing theorem remains dynamical persistence of the Reynolds/Krein alignment against current turning and material curvature-memory cost",
        "global_extension_guard": "isolated triads are self-contained martingale/phase systems, but the full PDE can reorient a high-dimensional fitness through shape motion and spatial concentration; persistence of productive Lamb/fitness alignment while V_spec collapses is the remaining constitutive problem",
        "symmetric_transport_note": "the amplitude-reduced heterochiral symmetric slice has stationary ratio near 0.5981296, but global symmetry of the two-parent optimizer is not claimed here",
        "case_taxonomy_used": False,
        "temporal_matching_used": False,
        "global_regularity_claimed": False,
    }


def cartan_hminus1_gradient_split(
    signed_frequencies: Sequence[float],
    modal_amplitudes: Sequence[complex],
    euler_sources: Sequence[complex],
    viscosity: float,
) -> dict[str, float]:
    """Exact H^{-1}_curl orthogonality of Euler motion and energy-gradient heat.

    On a nonzero curl eigenmode ``C e_i=a_i e_i``, the metric is
    ``||v||_{-1,C}^2=sum |v_i|^2/a_i^2``.  The H^{-1}_C gradient of E/2 is
    ``C^2 u``.  Euler energy conservation is exactly the metric orthogonality
    between its source and this gradient, hence the full NS speed satisfies a
    Pythagorean identity before any estimate.
    """
    if not (
        len(signed_frequencies) == len(modal_amplitudes) == len(euler_sources)
        and signed_frequencies
    ):
        raise ValueError("matching nonempty curl/amplitude/source data required")
    a = tuple(float(x) for x in signed_frequencies)
    z = tuple(complex(x) for x in modal_amplitudes)
    f = tuple(complex(x) for x in euler_sources)
    nu = float(viscosity)
    if not all(math.isfinite(x) and x != 0.0 for x in a):
        raise ValueError("finite nonzero curl eigenvalues required")
    if not all(math.isfinite(w.real) and math.isfinite(w.imag) for w in z + f):
        raise ValueError("finite complex modal data required")
    if not math.isfinite(nu) or nu < 0.0:
        raise ValueError("finite nonnegative viscosity required")

    energy = math.fsum(abs(zi) ** 2 for zi in z)
    enstrophy = math.fsum(ai * ai * abs(zi) ** 2 for ai, zi in zip(a, z))
    euler_energy_half_rate = math.fsum((zi.conjugate() * fi).real for zi, fi in zip(z, f))
    escale = max(1.0, math.sqrt(max(0.0, energy * math.fsum(abs(fi) ** 2 for fi in f))))
    if abs(euler_energy_half_rate) > 2.0e-11 * escale:
        raise ValueError("Euler source must be energy tangent")

    euler_action = math.fsum(abs(fi) ** 2 / (ai * ai) for ai, fi in zip(a, f))
    heat_action = enstrophy
    metric_cross = euler_energy_half_rate
    ns = tuple(fi - nu * ai * ai * zi for ai, zi, fi in zip(a, z, f))
    ns_action = math.fsum(abs(gi) ** 2 / (ai * ai) for ai, gi in zip(a, ns))
    represented = euler_action + nu * nu * heat_action
    scale = max(1.0, ns_action, represented)
    if abs(ns_action - represented) > 4.0e-11 * scale:
        raise AssertionError("H^-1 Cartan/gradient Pythagorean identity failed")
    return {
        "energy": energy,
        "enstrophy": enstrophy,
        "euler_hminus1_action": euler_action,
        "energy_gradient_hminus1_action": heat_action,
        "cartan_gradient_metric_cross": metric_cross,
        "ns_hminus1_action": ns_action,
        "represented_ns_hminus1_action": represented,
        "viscous_hminus1_action": nu * nu * heat_action,
    }


def normalized_triple_orientation_heat_law(
    curl_radii: Sequence[float],
    modal_energies: Sequence[float],
    triple_current: float,
    viscosity: float,
) -> dict[str, float]:
    """Heat damps raw Cartan current but leaves normalized triple orientation fixed.

    For one curl-spectral triple, ``tau`` is trilinear in the three state
    components.  Pure viscosity therefore gives
    ``tau'=-nu(r0^2+r1^2+r2^2)tau``.  Since each ``sqrt(E_i)`` has exactly the
    same individual heat rate ``-nu r_i^2``, the normalized orientation
    ``chi=tau/sqrt(E0 E1 E2)`` has zero viscous derivative.  Only Euler fiber
    motion can reorient it.
    """
    if len(curl_radii) != 3 or len(modal_energies) != 3:
        raise ValueError("exactly three radii and energies required")
    r = tuple(float(x) for x in curl_radii)
    e = tuple(float(x) for x in modal_energies)
    tau = float(triple_current)
    nu = float(viscosity)
    if not all(math.isfinite(x) and x > 0.0 for x in r + e):
        raise ValueError("positive finite triple radii/energies required")
    if not math.isfinite(tau) or not math.isfinite(nu) or nu < 0.0:
        raise ValueError("finite current and nonnegative viscosity required")
    root = math.sqrt(math.prod(e))
    chi = tau / root
    sigma2 = math.fsum(x * x for x in r)
    tau_dot_heat = -nu * sigma2 * tau
    root_dot_heat = -nu * sigma2 * root
    chi_dot_heat = (tau_dot_heat * root - tau * root_dot_heat) / (root * root)
    scale = max(1.0, abs(chi), abs(tau_dot_heat / root))
    if abs(chi_dot_heat) > 5.0e-14 * scale:
        raise AssertionError("normalized Cartan orientation acquired a viscous rotation")
    return {
        "normalized_orientation": chi,
        "raw_current_heat_rate": tau_dot_heat,
        "amplitude_root_heat_rate": root_dot_heat,
        "normalized_orientation_heat_rate": chi_dot_heat,
        "triple_heat_rate": nu * sigma2,
    }


def radial_mean_resolvent_balance(
    signed_frequencies: Sequence[float],
    modal_energies: Sequence[float],
    euler_energy_rates: Sequence[float],
    viscosity: float,
) -> dict[str, float]:
    """Exact radial-mean heat metric and resolvent Fisher-action upper.

    Put ``m=E_p |a|=K/E`` and ``s=(|C|-m)u``.  The normalized heat drift has the
    exact positive factor

        <s,|C|^2 u> = <s,(|C|+m)s>.

    The curl-spectral Euler velocity at level ``a`` is ``f(a)u_a`` with
    ``f=S/(2E_a)``.  Therefore the one radial cross term admits the exact
    weighted square in the metric ``|C|+m`` and the sharp Cauchy upper

        m' <= ||(|C|+m)^(-1/2) F_spec||^2 / (2 nu E).

    No shell/window is inserted; the inverse is the resolvent of the physical
    absolute-curl operator at its current mean scale.
    """
    a, e = _spectral_data(signed_frequencies, modal_energies)
    if len(euler_energy_rates) != len(a):
        raise ValueError("one Euler rate per curl node required")
    rates = tuple(float(x) for x in euler_energy_rates)
    if not all(math.isfinite(x) for x in rates):
        raise ValueError("finite Euler rates required")
    nu = float(viscosity)
    if not math.isfinite(nu) or nu <= 0.0:
        raise ValueError("positive finite viscosity required")
    spectral_source_action(a, e, rates)

    E = math.fsum(e)
    r = tuple(abs(x) for x in a)
    K = math.fsum(ri * ei for ri, ei in zip(r, e))
    if K <= 0.0:
        raise ValueError("positive mean absolute curl required")
    m = K / E
    Z = math.fsum(ri * ri * ei for ri, ei in zip(r, e))
    M3 = math.fsum(ri**3 * ei for ri, ei in zip(r, e))
    kappa = 0.5 * math.fsum(ri * si for ri, si in zip(r, rates))
    radial_metric = math.fsum(ei * (ri - m) ** 2 * (ri + m) for ri, ei in zip(r, e))
    direct_metric = M3 - m * Z
    mscale = max(1.0, radial_metric, abs(direct_metric))
    if abs(radial_metric - direct_metric) > 8.0e-11 * mscale:
        raise AssertionError("radial heat metric factorization failed")

    weighted_action_terms = []
    for ri, ei, si in zip(r, e, rates):
        if ei <= 0.0:
            if abs(si) > 5.0e-12 * max(1.0, *(abs(x) for x in rates)):
                raise ValueError("zero-energy node cannot carry interior spectral fitness")
            continue
        fi = si / (2.0 * ei)
        weighted_action_terms.append(ei * fi * fi / (ri + m))
    weighted_action = math.fsum(weighted_action_terms)

    cauchy_fraction = 0.0
    if radial_metric > 0.0 and weighted_action > 0.0:
        cauchy_fraction = kappa * kappa / (radial_metric * weighted_action)
        if cauchy_fraction > 1.0 + 8.0e-10:
            raise AssertionError("radial resolvent weighted Cauchy law failed")
        cauchy_fraction = min(1.0, max(0.0, cauchy_fraction))
    elif abs(kappa) > 8.0e-11 * max(1.0, abs(kappa)):
        raise AssertionError("radial current survived zero heat/action metric")

    direct_rate = (2.0 / E) * (kappa - nu * radial_metric)
    square_norm = radial_metric - kappa / nu + weighted_action / (4.0 * nu * nu)
    if square_norm < -8.0e-11 * max(1.0, radial_metric, weighted_action / (nu * nu)):
        raise AssertionError("radial resolvent square lost nonnegativity")
    square_norm = max(0.0, square_norm)
    represented_rate = -2.0 * nu * square_norm / E + weighted_action / (2.0 * nu * E)
    if abs(direct_rate - represented_rate) > 8.0e-10 * max(1.0, abs(direct_rate), abs(represented_rate)):
        raise AssertionError("radial resolvent Hilbert square failed")

    weighted_upper = weighted_action / (2.0 * nu * E)
    optimal_upper = 0.0
    if radial_metric > 0.0:
        optimal_upper = kappa * kappa / (2.0 * nu * E * radial_metric)
    if direct_rate > optimal_upper + 8.0e-10 * max(1.0, abs(direct_rate), optimal_upper):
        raise AssertionError("optimal radial scalar square upper failed")
    if optimal_upper > weighted_upper + 8.0e-10 * max(1.0, optimal_upper, weighted_upper):
        raise AssertionError("optimal radial upper exceeded resolvent action upper")

    return {
        "energy": E,
        "critical_stock": K,
        "mean_absolute_curl": m,
        "enstrophy": Z,
        "third_absolute_curl_moment": M3,
        "curvature_height": kappa,
        "radial_heat_metric": radial_metric,
        "resolvent_weighted_spectral_action": weighted_action,
        "weighted_cauchy_fraction": cauchy_fraction,
        "mean_curl_rate": direct_rate,
        "represented_mean_curl_rate": represented_rate,
        "mean_curl_resolvent_upper": weighted_upper,
        "mean_curl_optimal_upper": optimal_upper,
        "log_mean_curl_resolvent_upper": weighted_upper / m,
        "log_mean_curl_optimal_upper": optimal_upper / m,
    }


def sobolev_strain_transfer_multiplier(
    sobolev_exponent: float,
    output_radius: float,
    input_radius: float,
    *,
    same_helicity: bool,
) -> dict[str, float]:
    """Exact scalar block multiplier in the master Sobolev strain-transfer law.

    Put ``Lambda=|C|``, ``J=sgn(C)``, ``X=J_u`` and
    ``A=(1/2)[C,X]=P S(u)`` on divergence-free fields.  For

        G_s = Lambda^s X C Lambda^{-s},

    the symmetric part of ``G_s`` is obtained from the J-even and J-odd blocks
    of ``A`` by the hyperbolic multipliers returned here.  If ``D`` denotes the
    log-curl commutator, the operator formula is

        Sym G_s = sinh((s-1/2)D)/sinh(D/2) A_even
                  - cosh((s-1/2)D)/cosh(D/2) A_odd.

    The removable same-radius limit on the even block is ``2s-1``.  At the
    critical midpoint ``s=1/2`` the even multiplier is exactly zero while the
    odd multiplier is ``-sech(log(r_out/r_in)/2)``.
    """

    s = float(sobolev_exponent)
    ro = float(output_radius)
    ri = float(input_radius)
    if not all(math.isfinite(x) for x in (s, ro, ri)) or min(ro, ri) <= 0.0:
        raise ValueError("finite exponent and positive finite curl radii required")
    x = math.log(ro / ri)
    if same_helicity:
        if abs(x) <= 1.0e-10:
            multiplier = 2.0 * s - 1.0
        else:
            multiplier = math.sinh((s - 0.5) * x) / math.sinh(0.5 * x)
    else:
        multiplier = -math.cosh((s - 0.5) * x) / math.cosh(0.5 * x)
    midpoint = 0.0 if same_helicity else -1.0 / math.cosh(0.5 * x)
    return {
        "sobolev_exponent": s,
        "log_radius_ratio": x,
        "same_helicity": bool(same_helicity),
        "strain_multiplier": multiplier,
        "critical_midpoint_multiplier": midpoint,
        "energy_endpoint_multiplier": -1.0,
        "enstrophy_endpoint_multiplier": 1.0 if same_helicity else -1.0,
    }


def critical_logscale_strain_kernel(output_radius: float, input_radius: float) -> dict[str, float]:
    """Critical sech/Sylvester kernel relating physical odd strain to H^{1/2} boost.

    For radii ``r,s>0`` the critical midpoint multiplier is

        sech( (log r-log s)/2 ) = 2 sqrt(rs)/(r+s).

    Relative to the primitive master commutator ``[C,J_u]=2 P S(u)`` the
    multiplier is half as large, ``sqrt(rs)/(r+s)``.  The kernel is one on the
    diagonal and decays exponentially in log-scale separation.
    """

    r = float(output_radius)
    s = float(input_radius)
    if not all(math.isfinite(x) and x > 0.0 for x in (r, s)):
        raise ValueError("positive finite curl radii required")
    x = math.log(r / s)
    sech = 1.0 / math.cosh(0.5 * x)
    sylvester = 2.0 * math.sqrt(r * s) / (r + s)
    scale = max(1.0, abs(sech), abs(sylvester))
    if abs(sech - sylvester) > 8.0e-13 * scale:
        raise AssertionError("critical sech/Sylvester kernel identity failed")
    return {
        "log_radius_gap": abs(x),
        "strain_to_critical_multiplier": sech,
        "master_commutator_to_critical_multiplier": 0.5 * sech,
        "sylvester_multiplier": sylvester,
        "critical_strain_sign": -1.0,
    }


def critical_loggap_collective_bound(log_scale_gap: float) -> dict[str, float]:
    """Dimension-free cross-block bound for the critical log-scale strain filter.

    If one radius block lies below ``R`` and the other above ``exp(L)R``, the
    critical sech Schur/Sylvester map obeys

        ||B_cross|| <= min(1,csch(L/2)) ||A_odd,cross||.

    The ``1`` is the global contraction coming from the positive-definite sech
    kernel.  The ``csch`` term follows from the Neumann expansion in the
    separated block, and gives exponential decay for large log gap.  This is an
    operator-block statement rather than a per-edge shell threshold.
    """

    L = float(log_scale_gap)
    if not math.isfinite(L) or L < 0.0:
        raise ValueError("finite nonnegative log-scale gap required")
    if L == 0.0:
        separated = math.inf
        bound = 1.0
    else:
        separated = 1.0 / math.sinh(0.5 * L)
        bound = min(1.0, separated)
    return {
        "log_scale_gap": L,
        "global_contraction_bound": 1.0,
        "separated_block_bound": separated,
        "collective_multiplier_bound": bound,
    }


def poisson_critical_scale_measure_moments(
    signed_frequencies: Sequence[float], modal_energies: Sequence[float]
) -> dict[str, float]:
    """Canonical Poisson-scale probability attached to the critical stock.

    With ``Lambda=|C|`` and ``K=<u,Lambda u>``, put

        w_t = exp(-t Lambda) Lambda u,
        dmu(t) = 2 ||w_t||_2^2 dt / K.

    The measure is a probability.  It is equivalently obtained by first
    selecting a curl radius ``r`` with critical weight ``r E_r/K`` and then an
    exponential Poisson scale of rate ``2r``.  Consequently

        E_mu t = E/(2K) = 1/(2m),  m=K/E,
        E_mu t^2 = (sum E_r/r)/(2K).

    Thus the state itself supplies a cutoff-free scale law whose mean is exactly
    half the inverse mean absolute curl.
    """

    a, e = _spectral_data(signed_frequencies, modal_energies)
    r = tuple(abs(x) for x in a)
    if any(x <= 0.0 for x in r):
        raise ValueError("nonzero curl radii required for Poisson-scale moments")
    E = math.fsum(e)
    K = math.fsum(ri * ei for ri, ei in zip(r, e))
    if K <= 0.0:
        raise ValueError("positive critical stock required")
    m = K / E
    mass = math.fsum((ri * ei / K) for ri, ei in zip(r, e))
    mean_t = math.fsum((ri * ei / K) * (1.0 / (2.0 * ri)) for ri, ei in zip(r, e))
    second_t = math.fsum((ri * ei / K) * (1.0 / (2.0 * ri * ri)) for ri, ei in zip(r, e))
    represented_mean = 1.0 / (2.0 * m)
    represented_second = math.fsum(ei / ri for ri, ei in zip(r, e)) / (2.0 * K)
    scale = max(1.0, abs(mean_t), abs(represented_mean), abs(second_t), abs(represented_second))
    if abs(mass - 1.0) > 5.0e-13:
        raise AssertionError("Poisson critical scale lost unit mass")
    if abs(mean_t - represented_mean) > 5.0e-13 * scale:
        raise AssertionError("Poisson critical mean scale identity failed")
    if abs(second_t - represented_second) > 5.0e-13 * scale:
        raise AssertionError("Poisson critical second moment identity failed")
    return {
        "energy": E,
        "critical_stock": K,
        "mean_absolute_curl": m,
        "probability_mass": mass,
        "mean_poisson_scale": mean_t,
        "mean_poisson_scale_from_mean_curl": represented_mean,
        "second_poisson_scale_moment": second_t,
        "poisson_scale_variance": max(0.0, second_t - mean_t * mean_t),
    }


def critical_spectral_hminus_half_square(
    signed_frequencies: Sequence[float],
    modal_energies: Sequence[float],
    modal_energy_rates: Sequence[float],
    viscosity: float,
) -> dict[str, float]:
    """Exact critical Hilbert square using the actual curl-spectral Euler velocity.

    At an occupied signed-curl level ``a`` with energy ``E_a`` and nonlinear
    energy rate ``S_a``, the radial/curl-spectral Euler velocity has amplitude
    ``S_a/(2 sqrt(E_a))``.  Put ``r=|a|``.  Then

        A = Lambda^(3/2) u,
        B = Lambda^(-1/2) F_spec,

    satisfy ``<A,B>=kappa(0)`` and ``||A||^2=M3``.  Hence

        K' = -2 nu ||A-B/(2nu)||^2 + ||B||^2/(2nu).

    This quotient has already removed isospectral phase/shape motion and is
    therefore sharper than a square built from the full Lamb companion field.
    """

    a, e = _spectral_data(signed_frequencies, modal_energies)
    if len(modal_energy_rates) != len(a):
        raise ValueError("one nonlinear energy rate per curl level required")
    rates = tuple(float(x) for x in modal_energy_rates)
    if not all(math.isfinite(x) for x in rates):
        raise ValueError("finite nonlinear energy rates required")
    nu = float(viscosity)
    if not math.isfinite(nu) or nu <= 0.0:
        raise ValueError("positive finite viscosity required")
    spectral_source_action(a, e, rates)
    r = tuple(abs(x) for x in a)
    if any(x <= 0.0 for x in r):
        raise ValueError("nonzero curl radii required")

    E = math.fsum(e)
    K = math.fsum(ri * ei for ri, ei in zip(r, e))
    M3 = math.fsum(ri**3 * ei for ri, ei in zip(r, e))
    kappa = 0.5 * math.fsum(ri * si for ri, si in zip(r, rates))
    b2_terms = []
    for ri, ei, si in zip(r, e, rates):
        if ei <= 0.0:
            if abs(si) > 5.0e-12 * max(1.0, *(abs(x) for x in rates)):
                raise ValueError("zero-energy node cannot carry interior spectral velocity")
            continue
        b2_terms.append(si * si / (4.0 * ri * ei))
    B2 = math.fsum(b2_terms)
    cauchy = 0.0
    if M3 > 0.0 and B2 > 0.0:
        cauchy = kappa * kappa / (M3 * B2)
        if cauchy > 1.0 + 8.0e-10:
            raise AssertionError("critical spectral H^-1/2 Cauchy law failed")
        cauchy = min(1.0, max(0.0, cauchy))
    elif abs(kappa) > 8.0e-11 * max(1.0, abs(kappa)):
        raise AssertionError("critical current survived zero critical action")

    # The represented square norm is ||A-B/(2nu)||^2.
    square = M3 - kappa / nu + B2 / (4.0 * nu * nu)
    if square < -8.0e-10 * max(1.0, M3, B2 / (nu * nu)):
        raise AssertionError("critical spectral H^-1/2 square lost nonnegativity")
    square = max(0.0, square)
    K_rate = 2.0 * kappa - 2.0 * nu * M3
    represented = -2.0 * nu * square + B2 / (2.0 * nu)
    if abs(K_rate - represented) > 8.0e-10 * max(1.0, abs(K_rate), abs(represented)):
        raise AssertionError("critical spectral H^-1/2 square failed")

    log_upper = math.inf if K <= 0.0 else B2 / (2.0 * nu * K)
    scalar_optimal = 0.0 if M3 <= 0.0 else kappa * kappa / (2.0 * nu * M3)
    if K_rate > scalar_optimal + 8.0e-10 * max(1.0, abs(K_rate), scalar_optimal):
        raise AssertionError("critical scalar-optimal square upper failed")
    if scalar_optimal > B2 / (2.0 * nu) + 8.0e-10 * max(1.0, scalar_optimal, B2 / (2.0 * nu)):
        raise AssertionError("critical scalar action exceeded spectral action upper")

    # Two scale-covariant physical volumes.  The productive volume uses only the
    # current component aligned with the actual critical heat direction; the
    # broader action volume uses the whole curl-spectral Euler velocity.
    productive_action = 0.0 if M3 <= 0.0 else kappa * kappa / M3
    action_volume = math.inf if B2 <= 0.0 else E * K / (2.0 * B2)
    productive_volume = math.inf if productive_action <= 0.0 else E * K / (2.0 * productive_action)
    viscous_volume = math.inf if M3 <= 0.0 else E * K / (2.0 * nu * nu * M3)
    productive_reynolds = 0.0 if M3 <= 0.0 else abs(kappa) / (nu * M3)
    if math.isfinite(productive_volume) and math.isfinite(viscous_volume):
        volume_ratio = viscous_volume / productive_volume
        if abs(volume_ratio - productive_reynolds * productive_reynolds) > 2.0e-9 * max(1.0, volume_ratio):
            raise AssertionError("productive/viscous volume Reynolds identity failed")
    else:
        volume_ratio = 0.0

    # Critical probability pi_a=r E_a/K: B2/K=E_pi[(f/r)^2].
    critical_fitness_action = 0.0 if K <= 0.0 else B2 / K
    return {
        "energy": E,
        "critical_stock": K,
        "third_absolute_curl_moment": M3,
        "curvature_height": kappa,
        "critical_rate": K_rate,
        "spectral_hminus_half_action": B2,
        "critical_cauchy_fraction": cauchy,
        "represented_critical_rate": represented,
        "square_norm": square,
        "critical_rate_upper": B2 / (2.0 * nu),
        "critical_scalar_optimal_upper": scalar_optimal,
        "log_critical_upper": log_upper,
        "critical_probability_fitness_action": critical_fitness_action,
        "productive_critical_action": productive_action,
        "critical_action_volume": action_volume,
        "productive_action_volume": productive_volume,
        "critical_viscous_volume": viscous_volume,
        "productive_reynolds": productive_reynolds,
        "viscous_to_productive_volume_ratio": volume_ratio,
        "productive_growth_sign": 0.0 if kappa == 0.0 else math.copysign(1.0, kappa),
    }


def sobolev_spectral_hilbert_square(
    signed_frequencies: Sequence[float],
    modal_energies: Sequence[float],
    modal_energy_rates: Sequence[float],
    sobolev_exponent: float,
    viscosity: float,
) -> dict[str, float]:
    """Master Sobolev Hilbert square for one curl-spectral Euler velocity.

    For ``K_s=sum |a|^(2s) E_a`` and the physical Euler spectral velocity
    ``F_spec``, put

        A_s=|C|^(s+1)u,  B_s=|C|^(s-1)F_spec.

    Then the full NS rate obeys exactly

        K_s'=-2nu||A_s-B_s/(2nu)||^2+||B_s||^2/(2nu).

    The same heat operator is the gradient of ``K_s/2`` in the shifted
    ``H^(s-1)`` metric for every s.  Moreover ``||B_s||^2/K_s`` is the
    ``|a|^(2s)E``-weighted expectation of the single local score ``(f/|a|)^2``.
    """

    a, e = _spectral_data(signed_frequencies, modal_energies)
    if len(modal_energy_rates) != len(a):
        raise ValueError("one nonlinear energy rate per curl level required")
    rates = tuple(float(x) for x in modal_energy_rates)
    s = float(sobolev_exponent)
    nu = float(viscosity)
    if not all(math.isfinite(x) for x in rates + (s, nu)) or nu <= 0.0:
        raise ValueError("finite rates/exponent and positive viscosity required")
    spectral_source_action(a, e, rates)
    r = tuple(abs(x) for x in a)
    if any(x <= 0.0 for x in r):
        raise ValueError("nonzero curl radii required")

    K_s = math.fsum((ri ** (2.0 * s)) * ei for ri, ei in zip(r, e))
    D_s = math.fsum((ri ** (2.0 * s + 2.0)) * ei for ri, ei in zip(r, e))
    cross = 0.5 * math.fsum((ri ** (2.0 * s)) * si for ri, si in zip(r, rates))
    B2_terms = []
    for ri, ei, si in zip(r, e, rates):
        if ei <= 0.0:
            if abs(si) > 5.0e-12 * max(1.0, *(abs(x) for x in rates)):
                raise ValueError("zero-energy node cannot carry interior spectral velocity")
            continue
        B2_terms.append(0.25 * si * si * (ri ** (2.0 * s - 2.0)) / ei)
    B2 = math.fsum(B2_terms)
    square = D_s - cross / nu + B2 / (4.0 * nu * nu)
    if square < -1.0e-9 * max(1.0, D_s, B2 / (nu * nu)):
        raise AssertionError("master Sobolev square lost nonnegativity")
    square = max(0.0, square)
    rate = 2.0 * cross - 2.0 * nu * D_s
    represented = -2.0 * nu * square + B2 / (2.0 * nu)
    if abs(rate - represented) > 1.0e-9 * max(1.0, abs(rate), abs(represented)):
        raise AssertionError("master Sobolev Hilbert square failed")

    local_score = 0.0 if K_s <= 0.0 else B2 / K_s
    expected = 0.0
    if K_s > 0.0:
        for ri, ei, si in zip(r, e, rates):
            if ei <= 0.0:
                continue
            f = si / (2.0 * ei)
            pi = (ri ** (2.0 * s)) * ei / K_s
            expected += pi * (f / ri) ** 2
    if abs(local_score - expected) > 1.0e-10 * max(1.0, abs(local_score), abs(expected)):
        raise AssertionError("Sobolev local fitness expectation identity failed")
    return {
        "sobolev_exponent": s,
        "sobolev_stock": K_s,
        "sobolev_heat_moment": D_s,
        "nonlinear_half_rate": cross,
        "full_rate": rate,
        "spectral_shifted_action": B2,
        "square_norm": square,
        "represented_rate": represented,
        "rate_upper": B2 / (2.0 * nu),
        "log_rate_upper": math.inf if K_s <= 0.0 else B2 / (2.0 * nu * K_s),
        "sobolev_probability_local_fitness_action": local_score,
    }


def closed_triad_critical_action_bound(
    signed_frequencies: Sequence[float],
    modal_energies: Sequence[float],
    *,
    phase_cosine_abs: float = 1.0,
) -> dict[str, float | tuple[float, float, float]]:
    """Sharp scale-free H^-1/2 spectral-action bound for one physical closed triad.

    For three signed curl eigenvalues ``a_i=s_i r_i`` whose positive radii form
    a strict Fourier triangle, use the exact Waleffe magnitude determined by the
    three radii/helicity signs.  If the common physical phase has cosine ``c``,
    the closed-triad nonlinear energy rates are

        T_i = 4 |g| c sqrt(e0 e1 e2) (a_j-a_k)

    up to one common orientation sign.  The critical shifted action satisfies

        ||Lambda^-1/2 F_spec,tri||^2 <= (1/2) E_tri K_tri,

    and the constant 1/2 is a sharp supremum at the degenerate low-high-high
    boundary.  The theorem is independent of the absolute Fourier scale.
    """

    a, e = _spectral_data(signed_frequencies, modal_energies)
    if len(a) != 3 or any(x <= 0.0 for x in e):
        raise ValueError("three distinct occupied triad nodes required")
    r = tuple(abs(x) for x in a)
    scale = max(r)
    if not (r[0] + r[1] > r[2] and r[1] + r[2] > r[0] and r[2] + r[0] > r[1]):
        raise ValueError("strict physical Fourier triangle required")
    signs = tuple(1 if x > 0.0 else -1 for x in a)
    if any(x == 0.0 for x in a):
        raise ValueError("nonzero signed curl nodes required")
    from .helical import coupling_magnitude_closed

    g = coupling_magnitude_closed(r[0], r[1], r[2], signs[0], signs[1], signs[2])
    pc = abs(float(phase_cosine_abs))
    if not math.isfinite(pc) or pc > 1.0 + 1.0e-12:
        raise ValueError("absolute phase cosine must lie in [0,1]")
    pc = min(1.0, pc)
    common = 4.0 * g * pc * math.sqrt(e[0] * e[1] * e[2])
    rates = (
        common * (a[1] - a[2]),
        common * (a[2] - a[0]),
        common * (a[0] - a[1]),
    )
    # Affine conservation is exact algebraically; the generic helper also checks it.
    square = critical_spectral_hminus_half_square(a, e, rates, viscosity=1.0)
    B2 = square["spectral_hminus_half_action"]
    E = math.fsum(e)
    K = math.fsum(ri * ei for ri, ei in zip(r, e))
    ratio = B2 / (E * K)

    root_ratios = []
    pairs = ((1, 2), (2, 0), (0, 1))
    for i, (j, k) in enumerate(pairs):
        d = a[j] - a[k]
        root = 8.0 * g * g * d * d / (r[i] * (r[j] + r[k]))
        # Extreme near-degenerate triangles can lose a few ulps in Heron; this
        # guard is deliberately looser than theorem tests at conditioned states.
        if root > 1.0 + 2.0e-8:
            raise AssertionError("physical root action coefficient exceeded sharp unit bound")
        root_ratios.append(min(1.0, max(0.0, root)))
    if ratio > 0.5 + 2.0e-8:
        raise AssertionError("closed-triad critical action exceeded sharp one-half bound")
    return {
        "total_energy": E,
        "critical_stock": K,
        "waleffe_magnitude": g,
        "phase_cosine_abs": pc,
        "critical_shifted_action": B2,
        "action_to_energy_critical_ratio": min(0.5, max(0.0, ratio)),
        "sharp_action_ratio_upper": 0.5,
        "root_geometric_ratios": (root_ratios[0], root_ratios[1], root_ratios[2]),
        "log_critical_rate_upper_at_viscosity_one": B2 / (2.0 * K),
        "scale_free_log_rate_upper_at_viscosity_one": E / 4.0,
    }


def continuum_critical_operator_isometry_constant() -> dict[str, float]:
    """Exact R^3 unitary-Fourier constant for the critical midpoint HS isometry.

    For the viscosity-free operator

        Q_c(u)=|C|^-1 Sigma_c(u) |C|^-1
              =(1/2)|C|^-1/2 [J_u,sgn C] |C|^-1/2,

    the continuum helical kernel calculation gives

        ||Q_c(u)||_{HS}^2 = (1/64) ||u||_{Hdot^(1/2)}^2.

    The raw fixed-unit-wavevector geometric integral is ``pi^3/8`` and the
    repository unitary Fourier product contributes ``C_F^2=(2pi)^-3``.  Their
    product is exactly ``1/64``.  This is a continuum R^3 identity; discrete
    torus/lattice truncations have boundary/lattice weights and must not use the
    constant as an exact finite-cutoff equality.
    """

    raw = math.pi**3 / 8.0
    cf2 = (2.0 * math.pi) ** -3
    coeff = raw * cf2
    if abs(coeff - 1.0 / 64.0) > 5.0e-15:
        raise AssertionError("continuum critical operator isometry constant lost 1/64")
    return {
        "raw_fixed_wavevector_integral": raw,
        "unitary_fourier_factor_squared": cf2,
        "hilbert_schmidt_norm_squared_coefficient": coeff,
        "critical_norm_to_hs_isometry_factor": 8.0,
        "reynolds_hs_coefficient": coeff,
    }


def continuum_primitive_critical_channel_constants() -> dict[str, float]:
    """Exact continuum constants for the primitive critical two-way channel.

    On the full graded exterior algebra over R^3 define

        A=[|C|^-1,Q*]=[|C|^-1,u^flat wedge],
        B=(1/nu)[|C|^-1,(Q*)^2]=[|C|^-1,beta wedge].

    With the repository unitary Fourier convention, both kernels use the raw
    translation integral

        int_R3 (|p+e|^-1-|p|^-1)^2 dp = 4*pi.

    Exterior creation by a one-form has full graded Hilbert--Schmidt
    multiplicity 4, while creation by a two-form has multiplicity 2.  Hence

        ||A||_HS,gr^2 = (2/pi^2) ||u||_Hdot^(1/2)^2,
        ||B||_HS,gr^2 = (1/pi^2) ||u||_Hdot^(3/2)^2.

    Since B=i sum_j (dx^j wedge)[D_j,A], the exact constants also give

        sum_j ||[D_j,A]||_HS,gr^2 = 2 ||B||_HS,gr^2.

    Algebraically {A,Q*}=nu B.  The function records normalization only; it
    does not discretize the continuum operator or assert a finite-cutoff
    equality.
    """

    raw = 4.0 * math.pi
    cf2 = (2.0 * math.pi) ** -3
    one_form_creation_multiplicity = 4.0
    two_form_creation_multiplicity = 2.0
    channel = one_form_creation_multiplicity * cf2 * raw
    curvature = two_form_creation_multiplicity * cf2 * raw
    if abs(channel - 2.0 / math.pi**2) > 5.0e-15:
        raise AssertionError("primitive critical-channel constant lost 2/pi^2")
    if abs(curvature - 1.0 / math.pi**2) > 5.0e-15:
        raise AssertionError("primitive curvature-channel constant lost 1/pi^2")
    if abs(channel - 2.0 * curvature) > 5.0e-15:
        raise AssertionError("primitive channel Dirichlet/curvature ratio lost factor two")
    return {
        "raw_translation_integral": raw,
        "unitary_fourier_factor_squared": cf2,
        "one_form_creation_graded_multiplicity": one_form_creation_multiplicity,
        "two_form_creation_graded_multiplicity": two_form_creation_multiplicity,
        "critical_channel_hs_coefficient": channel,
        "curvature_channel_hs_coefficient": curvature,
        "translation_dirichlet_hs_coefficient": channel,
        "dirichlet_to_curvature_hs_ratio": channel / curvature,
    }



def continuum_critical_carre_du_champ_constants() -> dict[str, float]:
    """Exact R^3 constants for the critical ``Lambda=|D|`` carre-du-champ.

    For scalar ``f,g`` use ``Gamma(f,g)=f Lambda g+g Lambda f-Lambda(fg)``.
    The fractional-Laplacian kernel is ``pi^-2 |x-y|^-4``.  Consequently the
    velocity matrix ``Gamma_u=(Gamma(u_i,u_j))`` is positive, its trace mass is
    ``2K``, and the primitive Gauss source satisfies

        ||(beta wedge) A||_HS^2
          = (1/(4*pi^2)) int omega^T Gamma_u omega.

    Cauchy--Binet gives the universal ``e2`` and determinant prefactors below.
    The heat product law has sink ``-2 nu sum_j Gamma_{partial_j u}``.
    """

    kernel = 1.0 / math.pi**2
    gauss_metric = 1.0 / (4.0 * math.pi**2)
    e2 = 1.0 / (2.0 * math.pi**4)
    det = 1.0 / (6.0 * math.pi**6)
    return {
        "fractional_laplacian_kernel_coefficient": kernel,
        "critical_trace_integral_factor": 2.0,
        "gauss_hs_to_vorticity_metric_factor": gauss_metric,
        "second_elementary_symmetric_prefactor": e2,
        "determinant_prefactor": det,
        "heat_carre_du_champ_sink_coefficient": 2.0,
    }


def continuum_critical_gauss_bianchi_constants() -> dict[str, float]:
    """Exact continuum kernel constants for the primitive critical Gauss--Bianchi law.

    On R^3, R=|D|^-1 has kernel ``(2*pi^2*|x-y|^2)^-1``.  For

        E=alpha wedge,  F=beta wedge,
        A=[R,E], B=[R,F], V=E R E, Gc=V+2 nu B,

    full graded Hilbert--Schmidt multiplication by a two-form has multiplicity two,
    while the top-form source ``F A`` has multiplicity one.  Therefore

        ||V||_HS,gr^2
          = (1/(2*pi^4)) int int |u(x) cross u(y)|^2/|x-y|^4,

        ||F A||_HS^2
          = (1/(4*pi^4)) int int
              |omega(x).(u(y)-u(x))|^2/|x-y|^4.

    With the canonical pair weight ``1/(2*pi^2)``, the operator/pair norm factor is
    exactly ``1/pi^2``.  The Hom connection between ``D0=nu d`` and
    ``D1/2=nu d+(1/2)E`` has top-channel normal potential coefficient ``1/4``.
    These are continuum constants only; no hard-cutoff exactness is claimed.
    """

    riesz = 1.0 / (2.0 * math.pi**2)
    pair_weight = 1.0 / (2.0 * math.pi**2)
    pair_area_hs_raw = 2.0 * riesz * riesz
    gauss_source_hs_raw = riesz * riesz
    operator_to_pair = pair_area_hs_raw / pair_weight
    if abs(pair_area_hs_raw - 1.0 / (2.0 * math.pi**4)) > 5.0e-16:
        raise AssertionError("critical pair-area HS kernel constant changed")
    if abs(gauss_source_hs_raw - 1.0 / (4.0 * math.pi**4)) > 5.0e-16:
        raise AssertionError("critical Gauss source HS kernel constant changed")
    if abs(operator_to_pair - 1.0 / math.pi**2) > 5.0e-16:
        raise AssertionError("critical operator/pair normalization lost pi^-2")
    return {
        "riesz_kernel_coefficient": riesz,
        "canonical_pair_weight": pair_weight,
        "pair_area_hs_raw_coefficient": pair_area_hs_raw,
        "gauss_source_hs_raw_coefficient": gauss_source_hs_raw,
        "operator_to_pair_norm_squared_factor": operator_to_pair,
        "critical_square_hs_prefactor_times_viscosity": 0.5 * math.pi**2,
        "hom_midpoint_connection_coefficient": 0.5,
        "top_channel_normal_potential_coefficient": 0.25,
    }

def continuum_midpoint_operator_sobolev_dictionary(sobolev_exponent: float) -> dict[str, float]:
    """Exponent dictionary for the continuum midpoint Hilbert-scale isometry.

    If ``T u=Q_c(u)`` is the continuum midpoint transform, then

        T(|C|^alpha u)=Delta_op^(alpha/2) T(u)

    and therefore

        ||u||_{Hdot^s}^2
        =64 ||Delta_op^(s/2-1/4) Q_c(u)||_{HS}^2.

    This helper records only the exact exponent/normalization.  It does not
    discretize the continuum operator or claim the 1/64 identity for a hard
    torus/Galerkin cutoff.
    """

    s = float(sobolev_exponent)
    if not math.isfinite(s):
        raise ValueError("finite Sobolev exponent required")
    return {
        "physical_sobolev_exponent": s,
        "operator_laplacian_power": 0.5 * s - 0.25,
        "operator_norm_squared_multiplier": 64.0,
        "critical_midpoint_exponent": 0.5,
    }
