# High-strain descending-epoch physical dissipation telescope

## Status

**Independently audited theorem.**  The author's implementation on exact SHA
`774c702a692e67f5ccdf3a7028c16e437a0c5cc1` passed its original dedicated and
integration gates.  The trailing audit then exposed a scale-covariance failure in
the executable guards, repaired it without changing the continuum telescope, and
added a direct Fourier--Galerkin Navier--Stokes falsification lane.  The exact
audited PDE head is `70bb2e4a9ec5b7d8826dc1016da5157cbe5fb1ac`.

Theorem status:

`EXACT_HIGH_STRAIN_DESCENDING_EPOCH_TELESCOPE__PHYSICAL_GLOBAL_GRADIENT_RESERVOIR__THREE_SIXTEENTHS_RENEWAL_SCALE_DESCENT__ARBITRARY_TIME_OVERLAP_WEIGHTED_BY_SCALE__NO_EVENT_COUNT_RESET`.

The theorem addresses a genuine recursive owner rather than quotienting it away.  High strain is physical and remains event-facing.  The question is whether the **same high-strain route can recur consecutively forever**.

The answer encoded here is no, for a reason already contained in the Navier--Stokes viscous structure.

## 1. Every high-strain event pays actual resolved dissipation

At a first high-strain contact for the strict transporter

\[
V_j=S_{N_j/4}u,
\]

the already certified collision gives

\[
D_j
:=N_j\int_{I_j}\|\nabla V_j\|_2^2dt
\ge
D_*:=\frac{32\pi^2}{75c}.
\]

`D_j` is the actual normalized resolved dissipation on that physical first-hit history.  It is not a theorem counter.

## 2. The same physical route forces scale descent

The dissipation-weighted ancestor theorem places at least half of the actual `D_j` law on critical resolved shell-time atoms with

\[
M_j\le \frac14N_j.
\]

The canonical critical-shell renewal registers the next smooth carrier at

\[
N_{j+1}=A_j=\frac34M_j.
\]

Therefore every **consecutive high-strain renewal** obeys

\[
\boxed{
\frac{N_{j+1}}{N_j}\le\frac{3}{16}.
}
\]

This is not an analysis scale choice.  `M_j` is an actual `D_V|_G` shell-time ancestor and `A_j=3M_j/4` is the already certified renewal geometry.

If another physical owner occurs instead, the consecutive high-strain epoch ends.  That owner may later initiate another epoch, but it cannot be silently absorbed into the same monotone descent.

## 3. Arbitrary overlap does not duplicate the viscous reservoir without bound

Let

\[
G_*:=\int_0^{t_*}\|\nabla u(t)\|_2^2dt
\]

be the actual global gradient spacetime reservoir on the root history.

No disjointness of the high-strain histories `I_j` is assumed.  Since each `I_j` lies inside the global interval and the strict low-pass multiplier is an `L^2` contraction,

\[
\int_{I_j}\|\nabla S_{N_j/4}u\|_2^2dt
\le G_*.
\]

Hence

\[
\boxed{D_j\le N_jG_*.}
\]

This is exactly where overlap/reuse is handled.  The same piece of viscous spacetime dissipation may be seen by several nested histories, but each reading carries its **native physical scale weight** `N_j`.

Because those weights descend geometrically,

\[
\sum_jN_j
\le
\frac{N_0}{1-3/16}
=
\frac{16}{13}N_0.
\]

Therefore, even under complete time overlap,

\[
\boxed{
\sum_jD_j
\le
G_*\sum_jN_j
\le
\frac{16}{13}N_0G_*.
}
\]

No time bin, event disjointification, packet selector, or overlap entropy is introduced.

## 4. The physical frequency floor

Combining the two inequalities for each genuine high-strain event,

\[
D_*\le D_j\le N_jG_*,
\]

gives

\[
\boxed{
N_j\ge N_{min}:=\frac{D_*}{G_*}.
}
\]

But consecutive high-strain renewal also gives

\[
N_j\le\left(\frac{3}{16}\right)^jN_0.
\]

Thus only finitely many consecutive high-strain events can fit between `N_0` and the physical dissipation floor `D_*/G_*`.

The whole-epoch capacity gives a second finite bound:

\[
\#\{\text{high-strain steps}\}\,D_*
\le
\sum_jD_j
\le
\frac{16}{13}N_0G_*.
\]

The theorem records the minimum of the scale-floor count and this weighted-capacity count.

## 5. Navier--Stokes kinetic energy may close the reservoir numerically

For unforced Navier--Stokes with `nu>0`, in the repository normalization

\[
\|u(t)\|_2^2
+2\nu\int_0^t\|\nabla u\|_2^2ds
\le
\|u_0\|_2^2.
\]

Hence

\[
G_*\le\frac{\|u_0\|_2^2}{2\nu},
\]

and therefore

\[
\boxed{
N_{min}
\ge
\frac{2\nu D_*}{\|u_0\|_2^2}.
}
\]

The theorem API also accepts the actual `G_*` directly, because a sharper physical dissipation bound should not be discarded in favor of a coarser energy estimate.

## 6. Why this is not the forbidden `D_V` reset

The earlier anti-theorem remains correct.  On an **increasing** artificial/geometric scale chain `N_j=N_0q^j`, a fixed normalized `D_V` unit has actual viscous cost `nu D_V/N_j`, and the sum can be finite.  Therefore `D_V` is not a scale-independent additive reset.

The present theorem uses the opposite, physically supplied geometry:

\[
N_{j+1}\le\frac3{16}N_j.
\]

The bound depends on the epoch root scale `N_0` and on the global physical gradient reservoir `G_*`.  It applies only while the actual recursion remains on the consecutive high-strain ancestor-renewal route.

Thus we do not say “every high-strain event costs one unit of a finite budget.”  We say:

> repeated readings of one physical viscous reservoir are summable **along this actual descending lineage because the PDE itself supplies geometrically decreasing scale weights**.

## 7. Master consequence

A maximal consecutive high-strain epoch is finite.

Therefore an infinite recursive event path cannot eventually consist only of

`high strain -> D_V|_G ancestor -> renewed carrier -> high strain -> ...`.

If an infinite event path contains infinitely many high-strain events, it must also contain infinitely many **non-high-strain epoch breakers**.

Those breakers may be actual HH generation, source/SGS, independent service, high-tail work, material/new-ancestry relink, strain/deformation, causal reuse, or another genuine physical owner.  Smooth conservative `K_phys` role relink is not in this list because it is already certified same-event donor provenance with zero recursive depth.

This theorem therefore shrinks the mixed-owner recurrence frontier but does not close it.

No global Navier--Stokes regularity claim is made.


## 8. Independent audit certification

The audit first froze the author's theorem and asked whether all comparisons
remain valid after changing the physical unit of frequency and dissipation.  On
exact anti-test SHA `30a42157fc126ccf1bd4755cb6825e54eb3efdc6`, GitHub Actions
run `31461890777` failed all `6/6` adversarial tests:

- a tiny native-scale ancestor with `M=2N` and renewal `A=3M/4` was accepted as
  satisfying `M<=N/4` and `A<=3N/16`;
- a dissipation equal to only `D_*/2` was promoted to high strain;
- an event with `D_j>N_jG_*` over-drew the physical global gradient reservoir;
- a foreign next child at twice the actual renewed carrier frequency crossed the
  consecutive-epoch binding;
- a typed epoch certificate accepted `NaN` as its observed scale ratio;
- `_frequency_floor_count_upper(1e300,1e-300)` formed the ratio first, underflowed
  to zero, and raised `math domain error` although the logarithmic count is finite.

All five fail-open cases came from tolerances of the form
`epsilon*max(1,physical data)`.  At tiny native scales the artificial unit `1`
was larger than the entire event and therefore erased its actual relative law.

The repair at exact source SHA
`422ab677e635159d82720a2af60f7900e7b3be9f` compares only dimensionless native
ratios such as `M/N`, `A/M`, `A/N`, `D/D_*`, and consecutive
`N_next/A_previous`.  The reservoir check is performed as a logarithmic product
inequality, and the frequency count subtracts logarithms before division.  There
is no dimensional `max(1,...)` tolerance floor.  Nonfinite certificate data fail
closed.

Exact audited PDE head:

`70bb2e4a9ec5b7d8826dc1016da5157cbe5fb1ac`.

Adversarial run `31462711590` was **success**:

- `6/6` native-scale anti-tests and `801/801` theorem tests;
- `200,000` descending epochs with arbitrary interval overlap;
- native child-frequency range `[3.170e-199,2.431e+202]`;
- native normalized-dissipation range `[4.331e-140,1.032e+141]`;
- minimum native-relative geometric-frequency, normalized-capacity, and
  frequency-floor margins `3.270e-3`, `2.012e-2`, and `2.453e-5`;
- maximum certified count sampled `6`;
- all `200,000/200,000` nonconsecutive/ascending restarts rejected;
- `100,000` collision, resolved-ancestor, critical-carrier, and continuum-master
  dependency states at each audit lane, with master high-strain failures `0`.

The same run evolved the actual unforced incompressible three-dimensional
Fourier--Galerkin Navier--Stokes system with Leray projection, viscosity, `2/3`
dealiasing and RK4.  On each evolved trajectory it integrated the actual
low-pass strain action `K_N`, normalized resolved dissipation `D_N`, global
gradient reservoir `G_*`, and the positive `D_V` density disintegrated over the
physical dyadic ancestor shells.  For `N=16`, `M=4`, `A=3`, `c=1`, `nu=0.05`,
amplitude `256`, and `128` steps:

| resolution | div/`||u_0||` | energy balance | `K_N` | `D_N` | collision margin | reservoir margin | retained critical fraction | half-law margin | descendant/root gradient |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 12 | `6.679e-17` | `4.148e-9` | `2.199` | `3.534e6` | `0.9948` | `0.3027` | `1.0` | `0.5` | `0.0` |
| 16 | `6.933e-17` | `5.367e-9` | `2.230` | `3.521e6` | `0.9946` | `0.3190` | `1.0` | `0.5` | `0.0` |
| 20 | `6.203e-17` | `6.305e-9` | `2.246` | `3.521e6` | `0.9946` | `0.3217` | `1.0` | `0.5` | `0.0` |

The root-dissipation resolution spread was `3.777e-3`.  The renewed cutoff
`A/4=0.75` is below the first nonzero Fourier mode of this mean-zero periodic
Galerkin system, so its measured descendant resolved gradient is exactly zero;
the actual fixture contains one high-strain event and then leaves the pure
high-strain route.  This particular spectral-gap termination is only a numerical
fixture.  The continuum proof remains the general inequalities
`D_*<=D_j<=N_jG_*` and `N_(j+1)<=3N_j/16`.

Stored independent audit artifact:

`recorded-results/31462711590/audit-high-strain-descending-epoch-results/`.

GitHub artifact digest:

`sha256:de5a450a55267a40f0b2b843f7f6d29a44d9ed4dc949352b7818a093b10c8972`.

The exact-head dedicated run `31462711516` was also **success** with all `801`
tests, `50,000` wide native-scale epochs, the three-resolution physical PDE lane,
and every direct dependency.  Its artifact digest is
`sha256:e663f217f76b36a9b648cce354ba3691813d3303a8a3f23c84bdfb121231d040`.

Full physical-energy causal integration run `31462711518` was **success** on the
same exact head.  It passed the `801`-test theorem suite, the high-strain epoch
and PDE lanes, every existing physical owner/service/reuse compiler, the
`20,000`-episode master stress with worst margin `0.0`, and all `58` main
workflow steps.  Its `100`-file artifact digest is
`sha256:49239ed32c57e1f7a54c8fd351ac463b889816bf3d61b22464b95244a0997977`.

The numerical lane is falsification evidence, not a continuum regularity proof.
The executable theorem also remains conditional on each supplied step being the
actual high-strain event and critical-shell renewal from one common PDE history;
finite scalar records cannot independently prove that continuum provenance.

### Author pre-audit certification

Exact certified implementation SHA:

`774c702a692e67f5ccdf3a7028c16e437a0c5cc1`.

Dedicated GitHub Actions run:

`31460525711` — **success**.

It passed:

- `725` theorem tests;
- `50,000` descending high-strain epochs with arbitrary interval overlap permitted by construction;
- minimum geometric-frequency capacity margin `0.010872081991569138`;
- minimum weighted normalized-dissipation capacity margin `0.08466728285959269`;
- minimum last-scale/frequency-floor margin `1.2235436031599045e-05`;
- maximum sampled certified epoch count `6`;
- arbitrary-overlap cases `50,000`;
- non-consecutive/ascending restart rejections `50,000`;
- companion continuum-master high-strain epoch telescope failures `0`;
- companion master physical-time telescope residual `0.0`;
- companion master log-scale telescope residual `8.881784197001252e-16`.

Stored artifact:

`recorded-results/31460525711/high-strain-descending-epoch-telescope-results/`.

GitHub artifact digest:

`sha256:6869f07693605cc5ff1576b74b26a9bca79ff5a7c2dd0057135f17b1b3d41b56`.

Full physical-energy causal integration:

`31460525687` — **success** on the same exact SHA.

It passed the same `725`-test suite and `58` successful job steps through the complete causal spine.  The final master episode stress checked `20,000` traces with worst margin `0.0`.

Integration artifact digest:

`sha256:9b16322218e71e464ce0b2f1c69090a260550a1e10550239562c7b63126abd1f`.

Within the author's original suite, no fixture or theorem-source correction was
required: the first implementation SHA passed both original gates.  The stronger
independent native-scale failures and repair are recorded above.

This theorem closes eventually-pure consecutive high-strain recurrence.  It does not close mixed-owner recurrence and makes no Navier--Stokes global-regularity claim.
