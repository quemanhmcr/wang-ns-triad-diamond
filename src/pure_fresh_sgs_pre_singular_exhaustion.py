from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from src.descending_fresh_sgs_scale_epoch_telescope import FreshSGSScaleOwnerCertificate
from src.resolved_objective_strain_collision import l32_derivative_to_linf_constant


STATUS = (
    "DRAFT_PURE_FRESH_SGS_PRE_SINGULAR_EXHAUSTION__"
    "WEIGHTED_SCALE_MOMENT_FORBIDS_ARBITRARILY_HIGH_DOUBLING__"
    "FIXED_SMOOTH_FILTER_ENERGY_ENVELOPE_GIVES_POSITIVE_FIRST_STOP_TIME__"
    "NO_ENTROPY_COST_NO_FINITE_RESET_NO_GLOBAL_H1_ASSUMPTION"
)

FRESH_OWNER_FRACTION = 0.25
INCREMENT_TO_BAND_ENERGY = 4.0
FRESH_HARD_SHELL_DENOMINATOR = 24.0
SGS_STRESS_SUPPORT_RATIO = 0.5


def _positive_finite(*values: float) -> bool:
    return all(math.isfinite(float(x)) and float(x) > 0.0 for x in values)


def fresh_service_weighted_scale_moment_upper(
    parent_frequency: float,
    scaled_lifetime: float,
    global_energy_upper: float,
) -> float:
    """Upper bound sum_(j<=0) 2^(-j) F_j by the square LP frame.

    The certified fresh band law has

        F_j <= 4 int M_j ||u_j||_2^2 d tau,  M_j=2^j N.

    Multiplying by 2^(-j), summing, and using sum ||u_j||_2^2=||u||_2^2
    for the fixed square-normalized frame gives

        sum 2^(-j) F_j <= 4 c N E_global.

    This is a local capacity inequality for one already-existing positive service
    law.  It is not a causal reweighting and no normalized scale law is used as a
    probability of child energy.
    """
    N = float(parent_frequency)
    c = float(scaled_lifetime)
    E = float(global_energy_upper)
    if not _positive_finite(N, c, E):
        raise ValueError("positive finite parent frequency/lifetime/energy required")
    return INCREMENT_TO_BAND_ENERGY * c * N * E


def fresh_scale_exponential_moment_upper_from_owner_floor(
    parent_frequency: float,
    *,
    global_energy_upper: float,
    forced_square_service_floor: float,
    scaled_lifetime_upper: float,
) -> float:
    """Bound E_p[2^(-j)] for every fresh owner with F>=Y_*/4 and c<=c_*.

    If p_j=F_j/F, then

        sum 2^(-j) p_j
        <= 4 c N E/F
        <= 16 c_* N E/Y_*.
    """
    N = float(parent_frequency)
    E = float(global_energy_upper)
    Y = float(forced_square_service_floor)
    c = float(scaled_lifetime_upper)
    if not _positive_finite(N, E, Y, c):
        raise ValueError("positive finite scale/energy/service/lifetime data required")
    K = 16.0 * c * N * E / Y
    # Any probability on j<=0 has sum 2^(-j)p_j >= 1.  If the coarse upper is
    # below one, the stated owner data are mutually inconsistent.
    if K < 1.0 - 2.0e-13:
        raise ValueError("fresh owner/service floor is incompatible with the global energy envelope")
    return max(1.0, K)


def max_scale_atom_lower_from_exponential_moment(moment_upper: float) -> float:
    """A countable j<=0 law with E[2^(-j)]<=K has p_max>=1/(2L).

    Put n=-j>=0.  Choose L=ceil(log_2(2K)).  Markov gives

        P(n>=L) <= K/2^L <= 1/2.

    Hence at least half the mass lies on the L atoms n=0,...,L-1, so one atom
    carries at least 1/(2L).  This is a deterministic concentration consequence
    of the service moment, not an entropy payment.
    """
    K = float(moment_upper)
    if not math.isfinite(K) or K < 1.0:
        raise ValueError("finite exponential-moment upper K>=1 required")
    L = max(1, int(math.ceil(math.log2(2.0 * K) - 2.0e-14)))
    return 1.0 / (2.0 * L)


def fresh_scale_active_band_count_upper(moment_upper: float) -> int:
    K = float(moment_upper)
    if not math.isfinite(K) or K < 1.0:
        raise ValueError("finite exponential-moment upper K>=1 required")
    return max(1, int(math.ceil(math.log2(2.0 * K) - 2.0e-14)))


def doubling_fresh_hard_shell_mass_lower(
    parent_frequency: float,
    *,
    global_energy_upper: float,
    forced_square_service_floor: float,
    scaled_lifetime_upper: float,
) -> float:
    """Uniform lower for a genuine fresh child at frequency 2N.

    A child 2N can occur only from selected canonical band j=0.  Because the
    selected band is a max atom of the actual fresh scale law, its mass fraction
    is p_max.  The certified two-hard-shell route gives

        mu_(2N) >= p_max Y/(24 c).

    Replacing p_max by the moment lower, Y by Y_* and c by c_* is legitimate and
    introduces no scale entropy as a causal currency.
    """
    K = fresh_scale_exponential_moment_upper_from_owner_floor(
        parent_frequency,
        global_energy_upper=global_energy_upper,
        forced_square_service_floor=forced_square_service_floor,
        scaled_lifetime_upper=scaled_lifetime_upper,
    )
    p = max_scale_atom_lower_from_exponential_moment(K)
    return p * float(forced_square_service_floor) / (
        FRESH_HARD_SHELL_DENOMINATOR * float(scaled_lifetime_upper)
    )


def pre_singular_h1_doubling_shell_mass_upper(
    parent_frequency: float,
    pre_singular_h1_seminorm_sq_upper: float,
) -> float:
    """H^1 upper for the exact fresh doubling shell {N<|xi|<=2N}.

    On a fixed compact pre-singular interval a classical NS solution has finite
    H^1 norm.  For the hard shell read at child frequency 2N,

        mu_(2N)=2N ||P_{N<|xi|<=2N}u||_2^2
               <= 2 ||grad u||_2^2/N.

    No bound uniform up to a hypothetical singular time is assumed.
    """
    N = float(parent_frequency)
    H = float(pre_singular_h1_seminorm_sq_upper)
    if not _positive_finite(N, H):
        raise ValueError("positive finite parent frequency and pre-singular H1 bound required")
    return 2.0 * H / N


def fresh_doubling_is_excluded(
    parent_frequency: float,
    *,
    global_energy_upper: float,
    forced_square_service_floor: float,
    scaled_lifetime_upper: float,
    pre_singular_h1_seminorm_sq_upper: float,
) -> bool:
    lower = doubling_fresh_hard_shell_mass_lower(
        parent_frequency,
        global_energy_upper=global_energy_upper,
        forced_square_service_floor=forced_square_service_floor,
        scaled_lifetime_upper=scaled_lifetime_upper,
    )
    upper = pre_singular_h1_doubling_shell_mass_upper(
        parent_frequency,
        pre_singular_h1_seminorm_sq_upper,
    )
    return lower > upper + 8.0e-13 * max(1.0, lower, upper)


def fresh_doubling_parent_ceiling(
    reference_frequency: float,
    *,
    global_energy_upper: float,
    forced_square_service_floor: float,
    scaled_lifetime_upper: float,
    pre_singular_h1_seminorm_sq_upper: float,
) -> float:
    """Return a dyadic parent scale at/above which a fresh doubling is impossible.

    Along N_q=2^q N_0, the moment-band count increases by exactly one per dyadic
    step while the H1 upper decays like 2^(-q).  Hence lower/upper grows like
    2^q/(L_0+q) and the loop must terminate for finite positive inputs.
    """
    N = float(reference_frequency)
    if not _positive_finite(
        N,
        global_energy_upper,
        forced_square_service_floor,
        scaled_lifetime_upper,
        pre_singular_h1_seminorm_sq_upper,
    ):
        raise ValueError("positive finite doubling-ceiling data required")
    for _ in range(4096):
        if fresh_doubling_is_excluded(
            N,
            global_energy_upper=global_energy_upper,
            forced_square_service_floor=forced_square_service_floor,
            scaled_lifetime_upper=scaled_lifetime_upper,
            pre_singular_h1_seminorm_sq_upper=pre_singular_h1_seminorm_sq_upper,
        ):
            return N
        N *= 2.0
        if not math.isfinite(N):
            break
    raise AssertionError("finite-input dyadic doubling exclusion failed to appear")


@dataclass(frozen=True)
class FixedSmoothSGSFilterEnvelope:
    """Scale-independent analytic constants of the already chosen smooth cutoff.

    If S_N has convolution kernel K^S_N=N^3 K^S(Nx), set
    `filter_kernel_l1=||K^S||_1`.  Choose one fixed smooth reproducing multiplier
    equal to one on B_(1/2); if its unit-scale kernel is K^rep, set
    `stress_reproducing_kernel_l32=||K^rep||_(3/2)`.

    Smooth compact Fourier support makes both constants finite.  They are analysis
    constants of the fixed resolved representation, not physical reset currencies.
    """

    filter_kernel_l1: float
    stress_reproducing_kernel_l32: float
    lowpass_l2_contraction: bool = True
    stress_support_ratio: float = SGS_STRESS_SUPPORT_RATIO

    def __post_init__(self) -> None:
        if not _positive_finite(self.filter_kernel_l1, self.stress_reproducing_kernel_l32):
            raise ValueError("positive finite smooth-filter kernel norms required")
        if not self.lowpass_l2_contraction:
            raise ValueError("canonical strict low-pass must retain |S_N(xi)|<=1 L2 contraction")
        if abs(float(self.stress_support_ratio) - SGS_STRESS_SUPPORT_RATIO) > 1.0e-14:
            raise ValueError("resolved SGS stress theorem requires support B_(N/2)")


def resolved_sgs_stress_l1_energy_upper(
    global_energy_upper: float,
    filter_envelope: FixedSmoothSGSFilterEnvelope,
) -> float:
    """||R||_1 <= (||K^S||_1+1) E for R=S(u tensor u)-V tensor V."""
    E = float(global_energy_upper)
    if not _positive_finite(E):
        raise ValueError("positive finite energy upper required")
    if not isinstance(filter_envelope, FixedSmoothSGSFilterEnvelope):
        raise TypeError("typed fixed smooth SGS filter envelope required")
    return (filter_envelope.filter_kernel_l1 + 1.0) * E


def resolved_sgs_source_rate_energy_upper(
    parent_frequency: float,
    global_energy_upper: float,
    filter_envelope: FixedSmoothSGSFilterEnvelope,
) -> float:
    """Energy-only pointwise upper for rho_R=N^-4||sym grad div R||_inf.

    Since supp Rhat subset B_(N/2), reproduce R by a fixed smooth kernel:

        ||R||_(3/2) <= C_rep N ||R||_1.

    The certified order-two L^(3/2)->L^inf multiplier constant C_2 on B_(N/2)
    then yields

        rho_R <= C_2 C_rep (C_S+1) N E_global.
    """
    N = float(parent_frequency)
    E = float(global_energy_upper)
    if not _positive_finite(N, E):
        raise ValueError("positive finite parent frequency/energy required")
    if not isinstance(filter_envelope, FixedSmoothSGSFilterEnvelope):
        raise TypeError("typed fixed smooth SGS filter envelope required")
    C2 = l32_derivative_to_linf_constant(2, SGS_STRESS_SUPPORT_RATIO)
    return (
        C2
        * filter_envelope.stress_reproducing_kernel_l32
        * N
        * resolved_sgs_stress_l1_energy_upper(E, filter_envelope)
    )


def fresh_sgs_physical_slab_duration_lower(
    parent_frequency: float,
    *,
    sgs_source_weight_floor: float,
    global_energy_upper: float,
    filter_envelope: FixedSmoothSGSFilterEnvelope,
) -> float:
    """Positive physical-time price of one actual SGS-source first-stop slab.

    `Sigma_R=int rho_R d tau` and `d tau=N^2 dt`.  If Sigma_R>=sigma_*, then

        Delta t >= sigma_*/[C_rate N^3 E_global],

    where C_rate is the fixed-filter constant above.  This is derived from the
    local filtered NS source itself; no observer-chosen event gap is postulated.
    """
    N = float(parent_frequency)
    sigma = float(sgs_source_weight_floor)
    if not _positive_finite(N, sigma, global_energy_upper):
        raise ValueError("positive finite scale/source-floor/energy required")
    rho_upper = resolved_sgs_source_rate_energy_upper(
        N,
        global_energy_upper,
        filter_envelope,
    )
    return sigma / (N * N * rho_upper)


@dataclass(frozen=True)
class PureFreshSGSFirstStopStep:
    owner: FreshSGSScaleOwnerCertificate
    child_frequency: float
    child_critical_mass: float
    sgs_source_weight: float
    earlier_time: float
    later_time: float

    def __post_init__(self) -> None:
        if not isinstance(self.owner, FreshSGSScaleOwnerCertificate):
            raise TypeError("typed fresh SGS owner required")
        vals = (
            self.child_frequency,
            self.child_critical_mass,
            self.sgs_source_weight,
        )
        if not _positive_finite(*vals):
            raise ValueError("positive finite fresh first-stop data required")
        if not all(math.isfinite(float(t)) for t in (self.earlier_time, self.later_time)):
            raise ValueError("finite physical first-stop times required")
        if not self.earlier_time < self.later_time:
            raise ValueError("fresh recursive slab must have positive physical duration")
        a, b = self.owner.hard_shell_candidates
        tol = 8.0e-13 * max(1.0, self.child_frequency, a, b)
        if min(abs(self.child_frequency - a), abs(self.child_frequency - b)) > tol:
            raise ValueError("fresh child must be one of the certified physical hard-shell witnesses")
        if self.child_critical_mass + tol < self.owner.selected_hard_shell_mass_lower:
            raise ValueError("fresh child lost its certified hard-shell lower")
        if self.child_frequency > 2.0 * self.owner.parent_frequency + tol:
            raise ValueError("fresh scale route escaped the certified <=2N geometry")
        if self.child_frequency > self.owner.parent_frequency + tol:
            if abs(self.child_frequency - 2.0 * self.owner.parent_frequency) > tol:
                raise ValueError("fresh scale can increase only through the j=0 child 2N")
            if self.owner.selected_band_index != 0:
                raise ValueError("a fresh doubling child must come from selected canonical band j=0")

    @property
    def physical_duration(self) -> float:
        return self.later_time - self.earlier_time


@dataclass(frozen=True)
class PureFreshSGSPreSingularExhaustionCertificate:
    event_count: int
    root_frequency: float
    maximum_parent_frequency: float
    physical_span: float
    global_energy_upper: float
    pre_singular_h1_seminorm_sq_upper: float
    forced_square_service_floor: float
    sgs_source_weight_floor: float
    scaled_lifetime_upper: float
    doubling_parent_ceiling: float
    reachable_frequency_upper: float
    minimum_event_duration_lower: float
    maximum_event_count_from_physical_time: int
    scale_entropy_used_as_cost: bool = False
    global_h1_to_singular_time_assumed: bool = False
    additive_reset_used: bool = False
    observer_event_gap_assumed: bool = False

    def __post_init__(self) -> None:
        if self.event_count < 1:
            raise ValueError("nonempty pure fresh SGS epoch required")
        if self.maximum_event_count_from_physical_time < self.event_count:
            raise ValueError("physical-time/source-rate exhaustion does not cover observed epoch")
        if not _positive_finite(
            self.root_frequency,
            self.maximum_parent_frequency,
            self.physical_span,
            self.global_energy_upper,
            self.pre_singular_h1_seminorm_sq_upper,
            self.forced_square_service_floor,
            self.sgs_source_weight_floor,
            self.scaled_lifetime_upper,
            self.doubling_parent_ceiling,
            self.reachable_frequency_upper,
            self.minimum_event_duration_lower,
        ):
            raise ValueError("positive finite pure-fresh exhaustion certificate required")
        if (
            self.scale_entropy_used_as_cost
            or self.global_h1_to_singular_time_assumed
            or self.additive_reset_used
            or self.observer_event_gap_assumed
        ):
            raise ValueError("pure fresh SGS exhaustion used a forbidden shortcut")


def pure_fresh_sgs_pre_singular_exhaustion(
    steps: Sequence[PureFreshSGSFirstStopStep],
    *,
    global_energy_upper: float,
    pre_singular_h1_seminorm_sq_upper: float,
    forced_square_service_floor: float,
    sgs_source_weight_floor: float,
    scaled_lifetime_upper: float,
    filter_envelope: FixedSmoothSGSFilterEnvelope,
) -> PureFreshSGSPreSingularExhaustionCertificate:
    """Close an eventually-pure fresh-SGS first-stop word on one smooth interval.

    The proof has two local NS pieces.

    1. The weighted service-scale moment plus the pre-singular H1 bound forbids a
       j=0 fresh doubling above a finite scale.  Since fresh children are always
       <=2N, all scales on the word are therefore uniformly bounded.
    2. On that bounded scale range, the filtered-stress energy envelope gives a
       uniform positive *derived* duration for every SGS-source first-stop slab.
       Consecutive first-stop slabs are disjoint except for endpoints, so only
       finitely many fit in the finite physical interval.

    No finite-reset count, entropy payment, global H1 control at a hypothetical
    singular time, or observer-chosen minimum event gap is used.
    """
    rows = tuple(steps)
    if not rows:
        raise ValueError("nonempty pure fresh SGS epoch required")
    E = float(global_energy_upper)
    H = float(pre_singular_h1_seminorm_sq_upper)
    Yfloor = float(forced_square_service_floor)
    sigma_floor = float(sgs_source_weight_floor)
    cmax = float(scaled_lifetime_upper)
    if not _positive_finite(E, H, Yfloor, sigma_floor, cmax):
        raise ValueError("positive finite pure fresh SGS exhaustion inputs required")
    if not isinstance(filter_envelope, FixedSmoothSGSFilterEnvelope):
        raise TypeError("typed fixed smooth SGS filter envelope required")

    root = rows[0].owner.parent_frequency
    ceiling = fresh_doubling_parent_ceiling(
        root,
        global_energy_upper=E,
        forced_square_service_floor=Yfloor,
        scaled_lifetime_upper=cmax,
        pre_singular_h1_seminorm_sq_upper=H,
    )
    # If a doubling parent must lie strictly below `ceiling`, its child is below
    # 2*ceiling.  If the root already lies above the exclusion scale it cannot
    # increase at all.  This crude bound is intentionally representation-stable.
    reachable_upper = max(root, 2.0 * ceiling)

    max_parent = 0.0
    for i, row in enumerate(rows):
        owner = row.owner
        N = owner.parent_frequency
        max_parent = max(max_parent, N)
        tol = 8.0e-13 * max(1.0, N, reachable_upper)
        if N > reachable_upper + tol:
            raise ValueError("fresh recurrence exceeded the H1/service-imposed reachable scale ceiling")
        if owner.forced_square_service_threshold + 8.0e-13 * max(1.0, Yfloor) < Yfloor:
            raise ValueError("fresh owner fell below the uniform forced square-service floor")
        if owner.scaled_lifetime > cmax + 8.0e-13 * max(1.0, cmax, owner.scaled_lifetime):
            raise ValueError("fresh owner exceeded the registered scaled-lifetime upper")
        if row.sgs_source_weight + 8.0e-13 * max(1.0, sigma_floor) < sigma_floor:
            raise ValueError("fresh source owner fell below its uniform positive source-weight floor")
        if i > 0:
            prev = rows[i - 1]
            ftol = 8.0e-13 * max(1.0, prev.child_frequency, N)
            if abs(N - prev.child_frequency) > ftol:
                raise ValueError("pure fresh recurrence must recurse through the selected physical child scale")
            ttol = 8.0e-13 * max(1.0, abs(prev.earlier_time), abs(row.later_time))
            if abs(row.later_time - prev.earlier_time) > ttol:
                raise ValueError("consecutive fresh first-stop slabs must meet at their physical endpoint")
        if row.child_frequency > N + tol:
            # Any actually realized doubling must remain below the local H1 upper.
            lower = doubling_fresh_hard_shell_mass_lower(
                N,
                global_energy_upper=E,
                forced_square_service_floor=Yfloor,
                scaled_lifetime_upper=cmax,
            )
            upper = pre_singular_h1_doubling_shell_mass_upper(N, H)
            if lower > upper + 8.0e-13 * max(1.0, lower, upper):
                raise ValueError("fresh doubling violates the compact pre-singular H1 shell upper")

        dt_floor = fresh_sgs_physical_slab_duration_lower(
            N,
            sgs_source_weight_floor=sigma_floor,
            global_energy_upper=E,
            filter_envelope=filter_envelope,
        )
        if row.physical_duration + 8.0e-13 * max(1.0, row.physical_duration, dt_floor) < dt_floor:
            raise ValueError("fresh SGS source slab is shorter than its local filtered-NS source-rate lower")

    span = rows[0].later_time - rows[-1].earlier_time
    if span <= 0.0 or not math.isfinite(span):
        raise ValueError("positive finite physical span required")

    # Since every parent on the word is <=reachable_upper and rho_R upper grows
    # linearly in N, use the worst allowed scale to get one common duration floor.
    dt_min = fresh_sgs_physical_slab_duration_lower(
        reachable_upper,
        sgs_source_weight_floor=sigma_floor,
        global_energy_upper=E,
        filter_envelope=filter_envelope,
    )
    count_upper = int(math.floor(span / dt_min + 8.0e-12))
    if len(rows) > count_upper:
        raise ValueError("observed pure fresh SGS word exceeds the physical-time exhaustion bound")

    return PureFreshSGSPreSingularExhaustionCertificate(
        event_count=len(rows),
        root_frequency=root,
        maximum_parent_frequency=max_parent,
        physical_span=span,
        global_energy_upper=E,
        pre_singular_h1_seminorm_sq_upper=H,
        forced_square_service_floor=Yfloor,
        sgs_source_weight_floor=sigma_floor,
        scaled_lifetime_upper=cmax,
        doubling_parent_ceiling=ceiling,
        reachable_frequency_upper=reachable_upper,
        minimum_event_duration_lower=dt_min,
        maximum_event_count_from_physical_time=count_upper,
    )


def theorem_certificate() -> dict[str, object]:
    K = fresh_scale_exponential_moment_upper_from_owner_floor(
        1.0,
        global_energy_upper=1.0,
        forced_square_service_floor=1.0,
        scaled_lifetime_upper=1.0,
    )
    p = max_scale_atom_lower_from_exponential_moment(K)
    if K != 16.0 or p <= 0.0:
        raise AssertionError("fresh weighted-moment registration failed")
    C2 = l32_derivative_to_linf_constant(2, 0.5)
    return {
        "status": STATUS,
        "weighted_scale_moment": "sum_(j<=0) 2^(-j) F_j <= 4 c N E_global; with F>=Y/4 this gives E_p[2^(-j)]<=16 c N E_global/Y",
        "deterministic_concentration": "for K=E_p[2^(-j)], L=ceil(log2(2K)); Markov leaves >=1/2 mass on L top bands, hence p_max>=1/(2L); H_inf is not charged as a cost",
        "doubling_exclusion": "a child 2N must be the j=0 max-scale atom, so mu_(2N)>=p_max Y/(24c), whereas on any fixed compact pre-singular interval mu_(2N)<=2 sup||grad u||_2^2/N; N/log N beats 1 and high-scale fresh doubling is impossible",
        "source_rate": f"for fixed smooth strict filter, ||R||_1<=(C_S+1)E and ||R||_(3/2)<=C_rep N||R||_1; certified order-two constant C2={C2:.12g} then gives rho_R<=C2 C_rep(C_S+1) N E",
        "physical_time": "Sigma_R=int rho_R d tau, d tau=N^2dt; once the doubling exclusion bounds N, every fresh SGS first-stop has a uniform derived Delta t>=Sigma_*/[C2 C_rep(C_S+1) N_max^3 E]",
        "topology": "consecutive recursive first-stop slabs meet only at physical endpoints, so their derived positive durations add inside the finite pre-singular physical interval",
        "forbidden_shortcuts": "no scale entropy is causal, no fresh service is summed as a global reset, no generic shell-registration ratio is called progress, no global H1 bound at a hypothetical singular time is assumed, and no minimum event gap is postulated",
        "scope": "draft closure of an eventually-pure fresh-SGS source word on any fixed compact pre-singular smooth interval; cross-family alternation with strain/dissipation and actual HH/high-tail work remains separate",
    }
