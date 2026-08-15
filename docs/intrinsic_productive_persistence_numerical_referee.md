# Numerical referee: intrinsic productive persistence experiments

This note records independent small-Galerkin numerical experiments motivated by
Sections 35--42 of `native_curl_krylov_current_law.md`.  The computations were
performed directly from the rotational Navier--Stokes grammar

\[
C=\operatorname{curl},\qquad J_uv=\mathbb P(u\times v),\qquad
u C^2=-\nu\Delta,
\]

\[
\boxed{u_t=J_uCu-\nu C^2u.}
\]

The purpose is falsification/refereeing, not evidence for global regularity.
All observations below are finite-dimensional numerical checks.  They should
be read as diagnostics for the proposed intrinsic persistence frontier, not as
theorems.

## 1. Primitive operator identities reproduce the observable grammar

For a random real divergence-free Fourier-Galerkin state with normalized
energy, the nonlinear observable derivative was computed in two independent
ways: directly from `F_E = J_u C u`, and from

\[
\dot M_\phi^{NL}=2\langle u,[\phi(C),J_u]Cu\rangle.
\]

Representative results were

| observable | direct | commutator | absolute mismatch |
|---|---:|---:|---:|
| `I` | `-2.17e-18` | `0` | `2.17e-18` |
| `C` | `+4.29e-17` | `-3.21e-17` | `7.50e-17` |
| `C^2` | `+1.801419113907e-1` | same | `2.22e-16` |
| `|C|` | `+4.782605535373e-2` | same | `1.46e-16` |

Thus the energy/helicity null laws and the nontrivial enstrophy/critical
currents were recovered from the same master commutator at machine precision.

## 2. Critical Krein split isolates the productive Euler direction

With

\[
y=\Lambda^{1/2}u,\qquad \Lambda=|C|,\qquad \mathsf J=\operatorname{sgn}(C),
\]

and

\[
L_u=\Lambda^{1/2}J_u\mathsf J\Lambda^{1/2}=A_u+S_u,
\]

the numerical decomposition satisfied

- `||L_u y-(A_u+S_u)y|| = 3.90e-16`,
- `<y,A_u y> = -5.31e-18`,
- `<y,S_u y> = +2.391302767686e-2`,
- `<Lambda u,F_E> = +2.391302767686e-2`,
- mismatch between the last two quantities: `1.39e-17`.

This directly confirms on the Galerkin referee that the skew part is critical-
norm invisible and the self-adjoint helicity-flip part carries the exact
instantaneous critical Euler current.

## 3. The Section 39 completed square is an identity, not a loose estimate

For

\[
A=\Lambda^{1/2}\omega,\qquad
B=\Lambda^{-1/2}\mathbb P(u\times\Lambda u),
\]

one sample gave

\[
\langle A,B\rangle=-\kappa(0)
\]

to machine precision, and

- direct PDE `Kdot = -1.844549104113`,
- completed-square `Kdot = -1.844549104113`,
- mismatch `8.88e-16`.

A separate 500-state test verified the exact normalized sign law

\[
\frac{\dot K}{2\nu\|A\|^2}=-1-r\cos\theta,
\qquad
r=\frac{\|B\|}{\nu\|A\|},
\]

with maximum residual `4.0e-15`.  The classifier `-r cos(theta) > 1` matched
the sign of `Kdot` in all 500 samples.  In contrast, the norm threshold
`r>1` alone produced many false positives.  This numerically supports the
Section 39 warning that activity size without Hilbert angle is not causal.

## 4. Productive Fisher action is the radial regression component

For the actual occupied signed-curl spectral fitness

\[
f_a=\frac{S_a}{2E_a},\qquad p_a=E_a/E,
\]

a representative state gave

- `E_p f = 3.47e-18`,
- `Cov_p(|a|,f) = kappa(0)/E` to `3.47e-18`,
- `A_prod = 1.244948279238e-3` by both the covariance and
  `kappa(0)^2/(EZ-K^2)` formulas,
- `A_spec/E = 2.033178938471e-2`,
- productive spectral fraction `Corr_p(|a|,f)^2 = 0.061232`.

This confirms the regression interpretation used in Section 42.

## 5. Exact orthogonal chain: full Euler -> spectral -> productive

The full Euler force was decomposed orthogonally into its curl-cyclic spectral
projection and the one-dimensional productive radial regression component.
One representative state gave

\[
\|F_E\|^2=3.91952512\times10^{-1},
\]
\[
\|F_{spec}\|^2=2.48120316\times10^{-2},
\]
\[
\|F_{prod}\|^2=3.02909144\times10^{-3}.
\]

Both Pythagorean residuals were at machine precision (`5.55e-17` and `0`),
and the productive component alone reproduced the full critical current.

In 300 random states,

- median `||F_prod||^2/||F_E||^2 = 0.0010`,
- 90th percentile `0.0062`,
- maximum `0.0144`,
- median `||F_prod||^2/||F_spec||^2 = 0.0173`.

Thus generic Euler activity in this small referee is mostly nonproductive for
instantaneous critical radial growth.  This is only an ensemble observation,
not a universal bound.

## 6. Full signed-curl spectrum does not determine the nonlinear current

The complete signed-curl energy spectrum was held fixed while Fourier phases
were randomized.  The maximum spectral drift was at most `2.8e-17`, yet

\[
\kappa(0)\in[-6.15\times10^{-2},\,+5.36\times10^{-2}]
\]

in one phase ensemble, with approximately balanced signs.  Productive
fractions varied by more than an order of magnitude at exactly the same
spectrum.

A one-mode phase dial produced a first-harmonic critical current

\[
\kappa(\theta)=a+b\cos\theta+c\sin\theta
\]

with relative fit residual `1.81e-15`.

This is a direct numerical falsification of any spectral-state-only closure
for instantaneous nonlinear orientation.

## 7. Near-perfect productive alignment can be engineered instantaneously

At one fixed exact spectrum, phase-only optimization produced

\[
\boxed{\operatorname{Corr}_p(|a|,f)^2=0.998592}
\]

with spectral drift `2.8e-17`.  Across four independent spectra, analogous
phase-only searches reached

`0.9985, 0.9978, 0.9994, 0.9986`.

Therefore a universal pointwise gap

\[
\operatorname{Corr}_p(|a|,f)^2\le 1-\varepsilon
\]

cannot be expected from instantaneous spectral geometry alone.

Likewise, optimizing the stronger ratio

\[
\frac{\|F_{prod}\|^2}{\|F_E\|^2}
\]

at a fixed spectrum reached `0.27565`, with

- `Corr^2 = 0.83979`,
- `||F_spec||^2/||F_E||^2 = 0.32824`.

So the numerics also falsify a naive universal claim that the productive
fraction of the full Euler force is pointwise tiny.

## 8. The engineered dangerous orientation loses persistence under NS flow

A phase-engineered state with initial

\[
\operatorname{Corr}^2=0.9986
\]

was released into the actual Galerkin Navier--Stokes ODE with viscosity chosen
so that the initial critical derivative was positive.  The measured alignment
was

| time | `Corr^2` |
|---:|---:|
| `0` | `0.9986` |
| `0.300` | `0.5363` |
| `0.600` | `0.2260` |
| `0.900` | `0.0786` |
| `1.200` | `0.0129` |

while `K/K0` rose only to about `1.0106` by the final sample.  On this run,

- integrated `int R_*^2 d tau_E = 0.88870`,
- `4 Delta log(K/E) = 0.11625`,
- inequality slack `0.77245`.

Across four independently engineered spectra, the median productive-alignment
retention after time `0.45` was about `50%`.

These experiments do **not** prove a persistence theorem, but they support the
choice of persistence/action rather than instantaneous admissibility as the
remaining frontier: dangerous states exist and can be engineered, while the
full dynamics strongly reorients them.

## 9. Sharp Galilean low-frequency suppression is visible directly

For the critical helicity-flip matrix element in the low-frequency catalyst
geometry, scaling `q -> eps q` gave approximately constant `boost/|q|`:

| `eps` | `boost/|q|` |
|---:|---:|
| `0.5` | `0.447918` |
| `0.25` | `0.456952` |
| `0.125` | `0.458760` |
| `0.0625` | `0.458712` |
| `0.03125` | `0.458412` |

A random geometry search reached

\[
0.45927923,
\]

while

\[
\frac{3\sqrt6}{16}=0.45927933,
\]

leaving a gap below `1e-7`.  Uniform sweeping itself was critical-work
invisible at machine precision in the earlier operator test.  This numerically
supports the Section 38 Galilean null law and its sharp constant.

## 10. What these referees appear to rule out

Within their finite-dimensional numerical scope, the experiments argue against
trying to close the final theorem by any of the following shortcuts:

1. a spectral-state-only law for the sign of critical current;
2. a universal pointwise gap away from perfect productive alignment;
3. a universal pointwise smallness of productive/full Euler activity;
4. a norm-only Section 39 activity threshold;
5. a scalar instantaneous monotonicity principle.

The data instead repeatedly point to the same dynamical question already
isolated in Sections 41--42:

\[
\boxed{\text{can productive Krein self-return accumulate infinite action on the finite energy-loss clock?}}
\]

A plausible mechanism suggested by the trajectories is that the large
nonproductive/vertical part of the same Euler force may continually reorient
the state and destroy an engineered productive alignment.  This is only a
numerical hint.  A useful next theorem would have to turn that reorientation
into a quantitative integrated coercive estimate, while respecting the fixed
Cartan/Jacobi constraints, Galilean commutator null structure, and viscous
material geometry.

## 11. Scope guard

Nothing in this note establishes global regularity, a uniform persistence
bound, or even a dimension-independent Galerkin estimate.  The contribution is
strictly a numerical referee of the intrinsic formulas and of several tempting
false shortcuts.  Its main positive conclusion is narrower:

> the latest intrinsic reduction passes independent machine-precision identity
> checks, and adversarial finite-dimensional experiments continue to isolate
> **persistence of productive action**, rather than instantaneous size, as the
> nontrivial unresolved mechanism.
