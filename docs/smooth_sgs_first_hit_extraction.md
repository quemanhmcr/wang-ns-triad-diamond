# Smooth-SGS measurable first-hit extraction: the physical corridor

## Status

**EXACT_SMOOTH_SGS_MEASURABLE_FIRST_EXIT_AND_LOCAL_NO_HIT_EXHAUSTION__RECURSIVE_REENTRY_REMAINS**

This note closes a local continuum gap without replacing physics by a synthetic state machine.  Its scope is one already-selected efficient smooth-SGS block on an interval where the Navier--Stokes solution is smooth.  Universal recursive re-entry is deliberately left open.

## 1. The clocks were the wrong object

The master ledger contains quantities which become physically meaningful at different causal stages.  Putting all of them into one simultaneous vector of normalized "clocks" creates structure which the PDE does not possess.

The natural order is instead a filtration of information:

1. **event/support layer:** the physical transfer selection, cross-cell `Xi` excision, and any fixed transfer-loss certificate are already decided at the selected event;
2. **smooth-SGS evolution layer:** resolved strain, coherent deformation, objective source action, aspect/radius and gauge-invariant phase evolve on the physical slab;
3. **generation layer:** the child energy balance creates the actual positive HH work measure `dT_HH`; this is a measure on generated events, not another state clock;
4. **backward registration layer:** for each generated parent mark, the exact adjoint Kelvin equation is read backward to the common slice and may encounter residual/source, HH regeneration, material relink or the initial boundary;
5. **ancestry-information layer:** Shannon/Renyi reuse is formed only after the physical parent law exists.

This is causal ordering by information available in the dynamics, not priority by theorem name.

## 2. First hit is first exit from a physical corridor

Let `I=[a,b]`.  For each already-quotiented physical cause `r`, keep its observable `f_r:I->R` in its native units and its own certified threshold `theta_r`.  There are two elementary topologies.

For a certificate valid at equality,

`tau_r^cl = inf { t : f_r(t) >= theta_r }`.

For a state whose safe corridor includes the boundary,

`tau_r^op = inf { t : f_r(t) > theta_r }`.

No source action is divided by a strain action, no Moyal energy is compared numerically with phase, and no common-unit Radon--Nikodym weight is invented.

If `f_r` is continuous, `tau_r^cl` is Borel on path space `C(I)`: for rational `q`,

`{tau_r^cl <= q} = { sup_[a,q] f_r >= theta_r }`,

and the supremum functional is Lipschitz in the uniform norm.  For the strict debut,

`{tau_r^op < q} = union_{p in Q, p<q} { sup_[a,p] f_r > theta_r }`,

so it is Borel as well.  When a strict state starts safe and crosses, continuity gives `f_r(tau_r^op)=theta_r`; the cause is the open-superlevel crossing germ immediately to its right.  A `dt`-absolutely-continuous physical work law assigns the single boundary time zero mass.

For finitely many causes set

`tau = min_r tau_r`,

`J = { r : tau_r = tau }`.

Both are measurable; `J` is a finite measurable set-valued mark.  Exact simultaneous causes are therefore not a singularity of the construction.  They are simply several faces of the physical safe corridor touched at the same time, exactly matching the unsplit joint-stop theorem.

The construction is invariant under an independent strictly increasing continuous change of units `f_r -> Phi_r(f_r)`, `theta_r -> Phi_r(theta_r)`.  This invariance is the formal expression of a physical point: heterogeneous causes need no artificial scalarization.

## 3. Why the smooth-SGS observables have the needed regularity

Fix any pre-singular smooth interval.  The strict low-pass field `V=S_{N/4}u` and all spatial derivatives used by the existing theorems are smooth.  The coherent affine variables solve finite-dimensional ODEs of the form

`Xdot = Vbar`,

`Ldot = Abar L`.

Their coefficients are continuous, hence `X,L` are `C^1`, `L` remains invertible while the selected block exists, and singular values, aspect and physical radius are continuous functions of time.

The natural action variables are integrals of nonnegative continuous or locally integrable densities.  In particular,

`K(t)=int_a^t ||S_V(s)||_op ds`,

`I_K(t)=int_a^t K_coh(s) ds`,

and the objective source/action integrals are absolutely continuous.  The selected-role equation already proved in the repository gives continuous child coefficients and forcing on the smooth interval.  Therefore

`2 [ Re <c(t),F(t)> ]_+`

is continuous, and its physical positive work is absolutely continuous in time.  Backward adjoint integrals `I_R[s,t]` and `I_HH[s,t]` are absolutely continuous in the backward endpoint `s`.

At the high-strain boundary there is no topological hole: the exact analytic collision estimate is

`D_V >= 384 pi^2 K^2/c`.

Thus at `K=1/30` itself one already has the non-strict critical lower bound `D_V >= 32 pi^2/(75c)`; the old strict wording only made the branch disjoint from the low-strain display.

## 4. Material coherence without selector chattering

A best-cell or maximum-energy selector is not a physical object and can chatter when two nearly equal cells exchange order.  It is unnecessary.

At the physical event, anchor the material coherent cell in the canonical intrinsic coordinate

`zeta=(L^-1 X/2, L^T k)`.

Common affine/Kelvin transport leaves `zeta` exactly invariant.  Hence the material cell itself is fixed in intrinsic coordinates.  Let

`F(t)=V_{g(t)}u(t)`

be the coherent transform with the moving normalized Gaussian window and let `C` be the anchored measurable cell.  Define

`E_C(t)=||1_C F(t)||_2^2`.

On a smooth interval,

`Fdot = V_g u_t + V_{g_t}u`.

Therefore, in the phase-space Hilbert space,

`E_C'(t)=2 Re <1_C F,1_C Fdot>`

for a.e. `t`, and polarized Moyal gives

`|E_C'| <= 2 sqrt(E_C) ||Fdot||_2`

`          <= 2 sqrt(E_C) ( ||g||_2 ||u_t||_2 + ||g_t||_2 ||u||_2 ).`

So the actual Moyal content of the transported material cell is absolutely continuous.  There is no derivative of a discrete cell name.

Moreover the Moyal energy and polarized work are represented by phase-space densities.  Dyadic cell boundaries are Lebesgue-null and consequently carry zero Moyal energy/work.  Choosing half-open representatives only names null sets; it cannot create physical relink mass.  A genuinely new material address, cross-cell nonlinear work, or symmetric-difference service is therefore physical provenance, whereas movement of an arbitrary optimizer is not.

## 5. Branch-free helical phase

A principal `arg` has an artificial jump on the negative real axis.  The physical holonomy is a point `h` on `S^1`.  Use instead

`d_S1(1,h)=acos(Re h) in [0,pi]`.

This is globally continuous on the circle and represents the shortest gauge-invariant phase mismatch.  The existing `1/5` holonomy threshold can therefore be monitored without an angle branch cut.

## 6. Generated no-hit events are exactly registered survivors

For one generated parent coefficient the already-proved exact identity is

`z(t)=z(s)+I_HH[s,t]+I_R[s,t]`.

As `s` moves backward inside the common natural window, the two impulses are AC.  Material-cell provenance is measurable by the preceding section, and `t=0` is an absorbing boundary.  Hence residual/source, HH regeneration, material relink and boundary obstruction have measurable first debuts and an unsplit joint cause set if they tie.

If none occurs, at the common slice

`|I_R| < |z(t)|/4`,

`|I_HH| < |z(t)|/2`,

with no material relink.  The exact triangle inequality then forces

`|z(s)| >= |z(t)|/4`.

Thus the no-hit generated event is not a packet-persistence assumption.  It is exactly the previously certified **registered generated survivor**, and two continuing parents retain the existing `1/16` product factor.

## 7. Whole-block no-hit is the Kelvin-flat alternative

The event-level transfer deficit is checked before dynamic evolution.  On the retained low-transfer block, the remaining whole-eddy service observables are the existing physical ones: objective strain/source action, total strain, coherent deformation, aspect, affine radius and branch-free phase holonomy.  Their regularity above makes their exits measurable.  If none exits, the hypotheses of the existing coherent service-or-flat theorem remain in the safe corridor and its conclusion is exactly the certified coherent Kelvin-flat alternative.

Hence after support-level `Xi` excision the only locally uncharged objects are the two natural survivors already present in the physics:

`registered generated survivor`

or

`certified Kelvin-flat continuation`.

They need not be manufactured as states of a common clock, and the theorem does not assert that the two descriptions are mutually exclusive on every block.

## 8. Scope and downstream status

The theorem begins with **an already recursively selected efficient smooth-SGS block**.  It proves that once such a block is supplied, its local first-hit cause set is measurable and that local no-hit behavior is exhausted by the generated-survivor / Kelvin-flat alternatives.

The generated-survivor part of recursive re-entry is now supplied by the companion material-carrier/cutoff relay theorems: the common-slice coefficient stays in the same smooth carrier, hard roles are created only by actual HH work, and changing to the parent resolved cutoff creates no new source.

Downstream work has now supplied explicit reentry for several of the recursive owners that were left abstract here: generic critical shells, high-strain dissipation, fresh SGS scale service, and high-tail regeneration through Fourier locality and sliding natural time.  The current continuum problem is to compose those supplied routes with the remaining material/source owners into one measurable master exhaustion theorem.  That composition still cannot be replaced by a dataclass constructor or by choosing a new packet by fiat.

No global-regularity claim is made.
