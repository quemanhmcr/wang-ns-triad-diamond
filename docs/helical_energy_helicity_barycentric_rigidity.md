# Helical energy--helicity barycentric rigidity

Status: **RESEARCH THEOREM CANDIDATE — PR #7. Derived from the already-certified cyclic Waleffe work formula; not yet certified as a separate theorem.**

The certified cyclic triad theorem already gives more than energy conservation. Once the helical curl eigenvalue is named explicitly, the entire three-mode nonlinear energy-transfer vector is one-dimensional.

For a closed physical helical triad `k0+k1+k2=0`, define

\[
\lambda_i:=s_i|k_i|.
\]

Because `i k x h_s(k)=s|k|h_s(k)`, `lambda_i` is the curl eigenvalue of that physical helical mode. If `E_i=|a_i|^2`, then the triad contributions to kinetic energy and helicity are respectively the zeroth and first moments

\[
E_\triangle=\sum_iE_i,
\qquad
H_\triangle=\sum_i\lambda_iE_i.
\]

## 1. One scalar drives all three rooted works

The certified cyclic Waleffe identity is

\[
\boxed{
(T_0,T_1,T_2)
=R_\triangle
(\lambda_1-\lambda_2,\lambda_2-\lambda_0,\lambda_0-\lambda_1).
}
\]

Equivalently, with `1=(1,1,1)`,

\[
\boxed{T_\triangle=R_\triangle(\lambda\times\mathbf1).}
\]

Therefore two exact nonlinear conservation laws hold on the same physical triad before Hahn:

\[
\boxed{\mathbf1\cdot T_\triangle=0},
\qquad
\boxed{\lambda\cdot T_\triangle=0}.
\]

The first is kinetic-energy conservation. The second is helicity conservation.

If the three `lambda_i` are not all equal, the intersection

\[
\{T:\mathbf1\cdot T=0,\ \lambda\cdot T=0\}
\]

is one-dimensional. Thus simultaneous energy and helicity conservation determine the direction of the physical transfer vector uniquely; the Waleffe phase/amplitude factor `R_triangle` supplies its one remaining signed scalar magnitude. If all three `lambda_i` are equal, the certified work formula gives `T_triangle=0`.

## 2. The median curl-eigenvalue mode is the unique singleton side

Assume first that the curl eigenvalues are distinct and relabel the physical modes by

\[
\lambda_-<\lambda_0<\lambda_+.
\]

Energy and helicity conservation alone give

\[
\boxed{
T_-=-\frac{\lambda_+-\lambda_0}{\lambda_+-\lambda_-}T_0,
\qquad
T_+=-\frac{\lambda_0-\lambda_-}{\lambda_+-\lambda_-}T_0.
}
\]

Hence the two extreme-`lambda` modes always have the same work sign and the median-`lambda` mode has the opposite sign. The median mode is therefore always the unique singleton side of the donor/recipient pattern:

- `T_0<0`: one median donor feeds two extreme recipients;
- `T_0>0`: two extreme donors feed one median recipient.

This sharpens the generic cyclic sign theorem. The singleton donor/recipient is not arbitrary; it is selected by the ordering of the physical curl eigenvalues `s|k|`.

If two `lambda` values coincide, the third mode has zero work and the other two exchange equal and opposite energy. This includes genuine equiradial/helicity-degenerate transfer. If all three coincide, all work vanishes.

## 3. Exact barycentric split

Suppose the median mode donates `q=-T_0>0`. Put

\[
\alpha=\frac{\lambda_+-\lambda_0}{\lambda_+-\lambda_-},
\qquad
1-\alpha=\frac{\lambda_0-\lambda_-}{\lambda_+-\lambda_-}.
\]

Then

\[
T_-=\alpha q,
\qquad
T_+=(1-\alpha)q,
\]

and exactly

\[
\boxed{\lambda_0=\alpha\lambda_-+(1-\alpha)\lambda_+.}
\]

Thus the nonlinear event replaces energy mass at the median curl eigenvalue by the unique two-point distribution on the extreme curl eigenvalues with the same barycenter. The reverse event is the exact contraction back to the barycenter.

## 4. Convex-order law: only spread or contraction

For any convex `phi:R->R`, define `M_phi=sum_i phi(lambda_i) E_i`. The exact barycentric weights give

`dot M_phi^NL = T_0 [phi(lambda_0)-alpha phi(lambda_-)-(1-alpha) phi(lambda_+)]`.

The bracket is nonpositive by convexity. Therefore

- `T_0<0` implies `dot M_phi^NL>=0` for every convex `phi`;
- `T_0>0` implies `dot M_phi^NL<=0` for every convex `phi`.

So a nondegenerate closed helical NS triad has only two intrinsic nonlinear verbs on the curl spectrum:

`mean-preserving spread | mean-preserving contraction`.

This is exact algebra of the three simultaneous physical works, not a statistical model.

## 5. Enstrophy is the quadratic convex moment

Take `phi(lambda)=lambda^2`. Since `lambda_i^2=|k_i|^2`,

`Z_triangle=sum_i lambda_i^2 E_i=sum_i |k_i|^2 E_i`

is the triad contribution to `||curl u||_2^2`. The exact nonlinear production is

`dot Z_triangle^NL = sum_i lambda_i^2 T_i = -T_0 (lambda_0-lambda_-)(lambda_+-lambda_0)`.

In original cyclic indexing,

`sum_i lambda_i^2 T_i = R_triangle (lambda_0-lambda_1)(lambda_0-lambda_2)(lambda_1-lambda_2)`.

Hence a median-donor spread has strictly positive enstrophy production, while a median-recipient contraction has strictly negative enstrophy production. Degenerate curl-eigenvalue gaps give zero quadratic-moment production even when real equiradial energy exchange remains possible.

For the full incompressible NS field, the vorticity equation gives

`d/dt ||omega||_2^2 = 2 int omega.S omega dx - 2 nu ||grad omega||_2^2`.

Thus the sum of the triadic quantities above is exactly the nonlinear vortex-stretching/enstrophy face. **HH spectral work and strain/enstrophy production are two moments of the same nonlinear redistribution, not independent full-NS mechanisms.**

## 6. Existing signed-good donor/side law is a corollary

On the certified positive signed-good forward core, opposite parent helicities are already forced. After simultaneous helicity reversal if needed, write the child curl eigenvalue as `1`, the child-helicity parent ratio as `D`, and the opposite-helicity side-parent ratio as `S`. Then

`lambda_-=-S < lambda_0=D < lambda_+=1`.

The median mode is therefore the unique energy donor in the forward spread. The barycentric formulas give immediately

`W_child+/W_donor- = (D+S)/(1+S)`,

`W_side+/W_donor- = (1-D)/(1+S)`,

and hence

`W_side+/W_child+ = (1-D)/(D+S)`.

The current `3/5<D,S<5/8` window then yields exactly the already-certified side/donor and child/donor ratio bounds. Thus the rigid side-recipient law is not an additional dynamical mechanism: once opposite-helicity signed-good geometry is known, it is forced by simultaneous energy and helicity conservation.

## 7. Equiradial transfer is the degenerate convex-order boundary

If distinct physical modes have equal `|k|` but different helicity signs, their `lambda=s|k|` values may be `(+r,+r,-r)` or its reversal. Two curl eigenvalues coincide. The barycentric theorem then forces one rooted work to vanish and the remaining two to exchange equal and opposite energy.

All moments `|lambda|^p=|k|^p`, `p>=1`, are unchanged because every mode has the same radius. Hence real nonlinear energy transfer can coexist with zero radial displacement and zero enstrophy change. This recovers the certified equiradial radial-crossing anti-theorem as the degenerate boundary of the same energy--helicity law.

## 8. Divided-difference hierarchy

For any scalar observable `phi`, the exact triad production is

`sum_i phi(lambda_i) T_i = R_triangle [phi(lambda_0)(lambda_1-lambda_2) + phi(lambda_1)(lambda_2-lambda_0) + phi(lambda_2)(lambda_0-lambda_1)]`.

For distinct ordered eigenvalues this is the Vandermonde factor times the second divided difference of `phi`, up to the fixed orientation sign. Consequently affine observables `phi=1` and `phi=lambda` vanish identically: these are exactly energy and helicity conservation. The first polynomial moment which survives is `phi=lambda^2`, giving vortex-stretching/enstrophy production.

Thus the nonlinear triad does not act independently on every spectral moment. All moment production is generated by one three-point curvature functional on the curl spectrum.

## 9. Energy-side intrinsic grammar

At the level of positive helical modal energy, full incompressible Navier--Stokes has an exceptionally small local grammar:

1. physical modes are curl eigenstates `lambda=s|k|`;
2. quadratic interactions occur only on closed Fourier triads;
3. on each triad the nonlinear energy vector is a one-scalar barycentric redistribution preserving energy and helicity;
4. the redistribution is either a mean-preserving spread or its reverse contraction in `lambda`;
5. viscosity removes modal energy at rate `2 nu |k|^2 E = 2 nu lambda^2 E`.

The Waleffe geometry and phase determine **which** closed triads act and the signed scalar rate `R_triangle`; they cannot violate the barycentric transfer direction imposed by the two invariants.

Pressure, material labels, checkpoints, Young marks, and analysis partitions do not appear as independent helical-energy generators. Filtered pressure/SGS/source terms may be genuine faces of derived-state equations, but any claim that they create new full-field energy ancestry must reenter this fundamental stock/redistribution/dissipation law.

## 10. Strict UV-frontier creation is necessarily a spread

Consider one positive rooted child work in a nondegenerate closed triad. If

`|k_child| > max(|k_parent1|,|k_parent2|)`,

then `|lambda_child|` is larger than both parent absolute curl eigenvalues. Therefore `lambda_child` is an extreme of the ordered triple, never the median. Since the two extreme modes have the same work sign, positive child work forces the spread orientation: the median mode is the unique donor and the other extreme mode is a simultaneous positive side recipient.

Consequently every strict new radial-frontier child satisfies

`strict UV frontier creation -> mean-preserving spread -> positive triad enstrophy/vortex-stretching production`.

Conversely, if a positive recipient is the median curl eigenvalue, the event is a contraction. Because a number lying between two real numbers has absolute value no larger than the larger endpoint absolute value,

`|lambda_recipient| <= max(|lambda_donor1|,|lambda_donor2|)`.

Thus a contraction may move some donor energy upward across a chosen radial boundary, but at least one same-event donor already lives at radius at least as large as the recipient. It cannot be the first creation of a new radial frontier. Its high-frequency ancestry is already present in the same closed triad.

This distinction is stronger than raw low-to-high crossing. A crossing atom can occur inside a contraction, but then it is partly a redistribution of already-high stock. Genuine frontier creation is spread-only.

## 11. Critical `H^{1/2}` growth has a compulsory nonforward side flux

Define the positive absolute-helicity / critical Sobolev moment

`K = sum_(k,s) |k| E_(k,s) = sum |lambda| E`.

On one triad this is the convex moment `phi(lambda)=|lambda|`. If all three helicities have the same sign, `|lambda|` is affine on the whole triad interval, so the nonlinear contribution to `K` vanishes exactly. Homochiral transfer can redistribute energy across radius but cannot change this critical moment.

Now consider a heterochiral spread. The interval `lambda_-<lambda_+` crosses zero. Let `o` be the recipient whose helicity sign is opposite to the median donor. Helicity conservation and the barycentric split give exactly

`dot K_triangle^NL = 2 |lambda_o| W_o^+ > 0`.

Let `side` be the smaller-radial of the two positive recipients. The rooted edge of this side recipient has the larger-radial recipient among its two interaction parents. Hence

`J_side=0`:

the side is actual positive nonforward work and lies on the already-existing canonical bad/transfer-loss route.

If the opposite-helicity recipient is itself the side, then

`dot K_triangle^NL = 2 |lambda_side| W_side^+`.

If the opposite-helicity recipient is the larger-radial recipient, the exact barycentric weights give

`|lambda_o| W_o^+ < |lambda_side| W_side^+`.

Therefore every heterochiral critical-norm-growing spread obeys the intrinsic side-tax law

`0 < dot K_triangle^NL <= 2 |lambda_side| W_side^+`, with `J_side=0`.

Thus nonlinear growth of the scale-critical `H^{1/2}` energy can never occur as an isolated forward child gain. Energy plus helicity conservation forces a simultaneous positive nonforward side flux on the same closed physical triad. This statement does not declare that side work is dissipated or globally finite; it only identifies the compulsory same-event physical companion of critical growth.

## 12. Full-field critical dynamics has one nonlinear source

Let `u=u^++u^-` be the exact orthogonal helical decomposition and define
`K_+=sum_k |k||u^+(k)|^2`, `K_-=sum_k |k||u^-(k)|^2`.
Then `H=K_+-K_-` is signed helicity and `K=K_++K_-=||u||_(dot H^(1/2))^2` is the positive critical Sobolev stock.

Put `D_+^(3)=sum_k |k|^3|u^+(k)|^2` and `D_-^(3)=sum_k |k|^3|u^-(k)|^2`. The exact triad law gives one scalar nonlinear source `G(t)` in both helicity sectors:

`dot K_+ + 2 nu D_+^(3) = G`,

`dot K_- + 2 nu D_-^(3) = G`.

Thus `dot H + 2 nu (D_+^(3)-D_-^(3))=0` and `dot K + 2 nu (D_+^(3)+D_-^(3))=2G`. Homochiral triads contribute exactly zero to `G`; every nonzero contribution is heterochiral.

With the velocity-form convention `partial_t u=P(u cross omega)+nu Delta u`, the common source is

`G = -2 Re int (Lambda u^+) . [u cross (Lambda u^-)] dx`

`  = 2 Re int u . [(Lambda u^+) cross (Lambda u^-)] dx`.

The repeated-vector scalar-triple-product terms vanish pointwise. Hence the scale-critical nonlinear source is purely cross-helicity coupling. If either helical population is absent, the nonlinear `H^(1/2)` source vanishes exactly. This two-reservoir system is not a closure theorem by itself; it identifies the one native critical nonlinear face that a closure theorem must control.

## 13. Canonical good work is uniformly equivalent to native critical production

A positive contraction recipient is median in `lambda`, so one extreme interaction parent has radius at least as large and its rooted edge has `J=0`. In a spread, the smaller-radial positive recipient has the larger-radial recipient as an interaction parent and also has `J=0`. Therefore a canonical geometry-good positive edge can only be the larger-radial forward recipient of a spread.

The existing single-edge stability theorem then forces opposite parent helicities and the signed-good corridor. Write the child radius as `M`, the child-helicity donor ratio as `D`, and the opposite-helicity side ratio as `S`, with `3/5<D,S<5/8`. The energy--helicity law gives

`dot K_triangle^NL = 2 S M W_side^+`,

while barycentric conservation gives

`W_side^+/W_child^+ = (1-D)/(D+S)`.

Hence

`dot K_triangle^NL /(M W_child^+) = 2S(1-D)/(D+S)`.

This expression increases with `S` and decreases with `D`, so the exact endpoint bounds are

`18/49 < dot K_triangle^NL /(M W_child^+) < 20/49`.

Thus, on the entire canonical good law, child-scale-weighted physical work is uniformly equivalent to the native production law of the critical state `K=||u||_(dot H^(1/2))^2`. This does **not** make `M dW` a causal probability. Canonical causality remains `dW+`; the weighted law is a different, state-balance observable supplied by the PDE itself.
