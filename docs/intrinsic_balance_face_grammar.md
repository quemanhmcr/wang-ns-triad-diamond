# Intrinsic balance-face grammar

Status: **RESEARCH THEOREM CANDIDATE — PR #7. Not certified.**

The current programme contains many named mechanisms. The candidate claim here is that most names are not primitive dynamics. They are different coordinate readings of a smaller structure already forced by Navier--Stokes:

> **No recursive physical owner exists unless it is a nonzero face of an exact PDE balance after observer freedom and internal conservative circulation have been quotiented.**

This is not a new formalism imposed on the PDE. It is an attempt to identify one algebra already present independently in the certified mode, carrier, relink, material, and checkpoint theorems.

## 1. Native control-volume identity

Let `v` denote physical state nodes for one exact balance law. Let `E_v(t)>=0` be native stock when the state variable has a between-time stock meaning, `F_t(v,w)>=0` an already-existing same-time conservative flow, `S_v` independently derived production/source, and `D_v>=0` native dissipation. The local balance form is

\[
\boxed{\dot E_v=\sum_wF_t(w,v)-\sum_wF_t(v,w)+S_v-D_v.}
\]

For any physical subset `A`,

\[
\boxed{\frac d{dt}E_A=F(A^c,A)-F(A,A^c)+S_A-D_A.}
\]

After integration,

\[
\boxed{E_A(t_1)+\int D_A+\int F(A,A^c)=E_A(t_0)+\int S_A+\int F(A^c,A).}
\]

Every `A--A` flow cancels identically. This is already exact for the helical-mode energy law and, with zero time-stock term at one event, for the antisymmetric `K_phys` smooth-role relink matrix.

The primitive types are therefore: **stock; boundary flux; production/deformation face; dissipation/physical boundary; observation/provenance.**

## 2. Internal-flow cancellation

The exact helical donor kernel gives

\[
P_A=I_A+\Phi_{in,A},\qquad N_A=I_A+\Phi_{out,A},
\]

hence

\[
P_A-N_A=\Phi_{in,A}-\Phi_{out,A}.
\]

The internal nonlinear traffic `I_A` is real NS dynamics but has zero divergence through the boundary of `A`.

The smooth `K_phys` theorem has the same algebra in a different state space:

\[
T_{ab}^{phys}=-T_{ba}^{phys},\qquad R_a=\sum_bT_{ab}^{phys},
\]

so for a role subset `C`,

\[
\sum_{a\in C}R_a=F(C^c,C)-F(C,C^c).
\]

Again internal cycles cancel. Thus conservative internal motion may rearrange state identity and provenance, but only boundary divergence can change control-volume stock. Internal circulation is never a new source.

## 3. Pushforward cannot create source

Let `pi:v->ell` be a deterministic/measurable state label map which does not change the NS field. Push forward stock and flow by summing over fibers. The same divergence law holds on the label nodes. Internal flow hidden inside one fiber cancels; coarsening may change visible cancellation, but it cannot create a new source term.

This contains several existing anti-theorems at once: deterministic hard cells inherit canonical `dW+`; OO/ON/NN are restrictions of an existing positive law; material rereading can change provenance while the physical law stays fixed; cutoff repartition, common transported partitions, checkpoints, and hard-shell rereadings cannot manufacture physical source or event depth.

A representation change may expose a boundary or attach provenance. It is not allowed to manufacture `S`.

## 4. State witnesses and flux restrictions commute

At a fixed event let the physical role state be `U_C=(P_a u,P_bu,P_cu)` and let `M_C=M(U_C)` be a theorem witness of that state, for example the complex Young/Christ Gaussian mark. Let the canonical positive work law be restricted by a physical fate set `G`, `mu_G=1_G mu+`.

The restriction acts on the interaction/work measure, not on `u(t)`. Thus `U_C` is unchanged. Attaching the state witness gives

`mu_tilde_G=(id,M_C)_# mu_G`,

and forgetting the mark gives identically

`pr_# mu_tilde_G = mu_G`.

The PR #7 full-signed Young state-mark theorem is therefore a corollary of state/flux separation, not an isolated mixed-fate trick.

## 5. Exact carrier-energy faces

For a smooth self-adjoint carrier role `A`, with `eta=A^2`, Navier--Stokes gives exactly

`d/dt ||Au||_2^2 + 2 nu ||grad Au||_2^2 = <u,dot eta u> - 2 Re <eta u,B(u,u)>`.

After the strict resolved/high split and low--low support exclusion,

`d/dt ||Au||_2^2 + 2 nu ||grad Au||_2^2 = W_HH + J_A`.

The common observer/Kelvin motion is quotiented first. The remaining interface term splits into residual skew `K_phys`, which is same-event conservative boundary flux, and the symmetric part, which is existing strain/deformation.

So the PDE itself supplies the carrier-energy grammar

`endpoint stock | actual HH work | conservative relink flux | strain/deformation | viscous dissipation`.

No material name, checkpoint, coefficient locator, Christ mark, hard-cell label, or theorem horizon appears in this exact energy identity.

## 6. Derived state variables obey the same rule

For the coherent averaged resolved gradient `Abar=<grad V>_gamma`, the exact NS equation is

`dot Abar = -Abar^2 - <a^2> - <Hess P> - <grad div R_SGS> + nu <Delta A> - <r.grad A>`.

The apparent new coherent Reynolds terms are controlled by the already-existing coherent deformation variance and route to `D_V`; pressure, SGS and viscosity are inherited physical source faces, not new currencies created by averaging.

Thus “source” is relative to the native state variable whose exact PDE balance is being advanced. It need not mean creation of total kinetic energy. The invariant rule is: **the owner must be an actual term of an exact PDE balance for the physical state being advanced.**

## 7. Audit of the current master vocabulary

After the certified quotients, the historical names collapse as follows.

| Manifestation | Intrinsic type | Recursive generation? |
|---|---|---|
| same-carrier inherited energy | endpoint stock continuation | no |
| hard/smooth conservative relink | conservative same-event flux | no by itself |
| canonical `dW+` HH generation | physical nonlinear boundary work into selected state | yes when the energy gate selects it |
| high strain / `D_V` | deformation/dissipation face | yes |
| native pressure/SGS/viscous source | term of exact resolved/coherent state equation, then certified state consequence | yes through its native owner |
| full-natural service witness | consequence on an already completed corridor | no second event |
| checkpoint / cover rereading | observation of the same continuing state | no |
| material membership / OO-ON-NN | provenance restriction of existing law | no |
| raw material-relink/new-ancestry name | locator with no native face | rejected / fail closed |
| Young/Christ Gaussian mark | state witness | no charge |
| coefficient impulse threshold | locator for physical-energy reentry | no owner until exact work/energy is read |
| transfer loss / reuse / sideband | terminal functional of an existing physical law | no recursive state creation |
| `t=0` | physical absorbing boundary | terminates |

The historical enum is therefore not the intrinsic grammar. After quotient, the recursive alphabet is much smaller.

## 8. No-owner-without-balance-face candidate theorem

A canonical recursive owner declaration is admissible only if all of the following hold:

1. identify the physical state variable whose change is being explained;
2. derive its evolution directly from Navier--Stokes or an exact projection/average of it;
3. quotient common observer/representation motion before ownership;
4. bind the proposed owner to a nonzero endpoint, boundary-flux, production/deformation, dissipation, or physical-boundary term of that exact balance;
5. allow downstream restrictions/pushforwards/marks to carry the face but never create another copy;
6. quotient conservative internal flow to zero new recursion depth.

A label which fails item 4 is not a weak owner. It is not yet an owner at all.

## 9. Fundamental nonlinear face: energy--helicity barycentric redistribution

The closed helical triad formula sharpens the abstract boundary-flux face. Put `lambda_i=s_i|k_i|`. The exact rooted energy-work vector is

`T_triangle = R_triangle (lambda x 1)`.

Hence its zeroth and first curl-spectral moments vanish simultaneously:

`sum T_i=0`, `sum lambda_i T_i=0`.

So the full-field nonlinear energy face is not an arbitrary conservative network flow. On every three-mode interaction it is a one-scalar barycentric redistribution of positive energy mass preserving both total energy and helicity. The same event is either a mean-preserving spread or the reverse contraction on the signed curl spectrum.

Several historical manifestations are therefore readings of one face rather than independent verbs:

- rooted HH child work is one coordinate of `T_triangle`;
- cyclic donor/recipient provenance is its positive transport reading;
- radial crossing is a control-volume boundary reading after applying `|lambda|=|k|`;
- helicity is its conserved first moment;
- enstrophy/vortex stretching is its quadratic convex moment;
- the signed-good side recipient is the compulsory barycentric companion;
- equiradial transfer is a degenerate boundary on which all `|lambda|^p` moments remain unchanged.

This is developed in `helical_energy_helicity_barycentric_rigidity.md`.

The consequence for owner ontology is important but limited: actual rooted work charges remain real physical energy gains/losses. What disappears is the idea that their various scale/helicity/strain manifestations are unrelated mechanisms.

## 10. Critical-state refinement: the native observable is curl-spectral curvature

The energy--helicity theorem supplies a preferred critical state. Put `K=||u||_(dot H^(1/2))^2`. On the signed curl coordinate `lambda=s|k|`, the nonlinear triad law preserves the affine observables `1` and `lambda`; every other spectral production is the same three-point second-divided-difference curvature functional.

For `phi(lambda)=|lambda|`, curvature is concentrated at the helicity-sign interface `lambda=0`. Hence homochiral triads contribute exactly zero to the critical nonlinear source, while heterochiral spread/contraction produces the signed Jensen/Tanaka defect. Equivalently the two positive helicity reservoirs `K_+` and `K_-` receive the same nonlinear source `G`.

The existing log-progress multiplier `J` is **not** this intrinsic critical source. On the canonical `J`-good sublaw, signed-good geometry does give the useful comparison `18/49 < dK_NL/(M dW_good) < 20/49`. But the converse is false: the sharp strict-UV geometry for critical production occurs near `D=0.4539303`, `S=0.8242110`, where `J/J_*` is only about `0.272`. Thus the current good/bad transfer-efficiency split is a proof observable for one log-progress architecture, not the fundamental partition of nonlinear critical danger.

This correction is structurally important. A physical work atom may be terminal for the existing multiplicative-transfer recursion while still making a large positive contribution to the critical state balance. Therefore intrinsic regularity must ultimately control the signed critical curvature source itself rather than declare `J`-bad work physically harmless. Canonical `dW+` remains the causal energy-transfer law; the curvature law is a distinct state-balance observable derived from the same signed physical work.

See `curl_spectral_curvature_balance.py` and `helical_energy_helicity_barycentric_rigidity.md`.

## Curl-spectral curvature refinement

The full-field nonlinear face can be compressed further.  On the signed curl spectrum `lambda=s|k|`, each closed helical triad is a one-scalar barycentric spread/contraction preserving the affine moments `1` and `lambda`.  Its signed energy source is the second `lambda` derivative of a triangular barycentric tent.  Summing triads defines an intrinsic signed curvature potential `kappa` with

`partial_t rho = partial_lambda^2 kappa - 2 nu lambda^2 rho`.

This single current has several exact readings: `2 kappa(0)` is nonlinear `Hdot^(1/2)` production and opposite-helicity pair creation; `2 int kappa` is nonlinear enstrophy/vortex-stretching production; `kappa'(-R)-kappa'(R)` is net radial exterior-energy crossing.  Thus HH transfer, critical pair creation, radial cascade and strain/enstrophy are not primitive owner types at the full-field energy level; they are coordinate readings of one signed state current.

For every critical-growing barycentric spread with total physical donor work `Q`, median-donor radius `N_d` and extreme-recipient frontier `N_c`, the deterministic Tanaka-scale cocycle gives

`eta_crit := P_crit/(N_c Q) <= 1-N_d/N_c <= log(N_c/N_d) = (1/2) log(T_d/T_c)`.

This uses neither the log-progress `J` classifier nor a synthetic cost.  The remaining primitive question is constitutive: the actual Waleffe quadratic law must control how fast `kappa` can transport the positive state relative to viscous killing.  A balance law alone is not a regularity proof.
