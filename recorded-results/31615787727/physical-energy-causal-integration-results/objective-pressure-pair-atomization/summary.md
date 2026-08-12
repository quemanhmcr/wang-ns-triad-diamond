# Objective pressure Hessian: direct hard-pair source atomization

Status: **EXACT_OBJECTIVE_PRESSURE_HESSIAN_DUAL_PAIR_ATOMIZATION__PAIR_OWNER_ALWAYS_TO_CRITICAL_SHELL_WITH_ENTROPY_TRADEOFF__QUARTER_SPLIT_ONLY_A_COROLLARY__AGGREGATE_MU_V_NOT_CANONICAL**.

The canonical pressure route no longer needs to coarse-grain the resolved Hessian source into aggregate `mu_V`.  Let

`H=-N^-4 <Hess P>_gamma`

be the actual averaged pressure matrix in the coherent corotational strain equation and choose the measurable Frobenius dual `Z=H/||H||_F` (`Z=0` at `H=0`).  With the exact filtered split `H=H_V+H_R`,

`rho_P=||H||_F = Z:H_V + Z:H_R`.

Decompose only the resolved transporter `V=S_(N/4)u` into hard orthogonal dyadic shells at this **physical pressure event**.  The bilinear resolved pressure tensor then expands exactly into unordered atoms `{a,b}`: diagonal once, off diagonal as both orientations together.  Thus

`rho_P <= [r_R]_+ + sum_(a<=b)[p_ab]_+`.

These are positive source/service atoms, not child-energy causal probabilities.

The already-certified order-2 pressure and shell Bernstein constants give the exact rational ordered-pair coefficient

`256/1425 < 1/5`.

Therefore every unordered hard pair obeys

`|p_ab| <= (kappa_ab/5)(Mmax/N)^4 sqrt(mu_a mu_b)`,

with `kappa=1` diagonal and `2` off diagonal.  Gaussian averaging is a probability contraction and costs nothing further.

For integrated pressure source weight `Sigma_P`, exact positivity gives the joint half split

`SGS-positive source >=Sigma_P/2`  OR  `resolved positive pair source >=Sigma_P/2`.

The SGS branch still yields `int||R||_(3/2)>=190 Sigma_P` and enters the existing coherent-service compiler.  On the resolved branch normalize the actual unordered pair source law.  Its native statement is the threshold-free tradeoff

`mu_child exp(H2_pair) >= 320 Sigma_P/c`.

Indeed a countable positive pair law has an attained maximal atom, `q_max>=sum q^2=exp(-H2_pair)`, and that actual pair exposes the stated hard `u`-shell lower after the strict-lowpass contraction.  Therefore **every resolved pressure-pair owner already enters the generic critical-shell theorem**.  The familiar `theta=1/4` split is only a diagnostic corollary:

- if a pair is theta-dominant, its integrated capacity forces at some time an actual hard child shell with

  `mu_child >= [5 Sigma_P/(16c)](N/Mmax)^4 >= 80 Sigma_P/c`,

  because every resolved pair has `Mmax<=N/4`; this is a genuine input to the generic critical-shell theorem;
- if no pair exceeds one quarter, the actual source law has `H2_pair>=log 4`.  This quantifies why the unconditional shell seed is weaker, but it does not create another master fate and is not a causal HH probability.

At exact quarter mass both corollaries hold.  They are not competing routes: the physical resolved-pair owner has already been sent once to the same critical-shell recursion.

Material/coherent labels are deliberately absent from the scale proof.  They may be attached after the hard event as sidecars; on a supplied signed-good low-strain material lineage the previously certified fixed objective-Hessian pair contraction `<1/5` remains an optional reuse refinement.

The old coarse inequality `rho_P<=mu_V/5700+||R||_(3/2)/380` remains true as a diagnostic, but aggregate `mu_V` is no longer the canonical pressure renewal state.

Stress: `50000` pressure tensor/pair/source states
- worst dual scalar residual: `1.232e-14`
- minimum positive source-cover margin: `-1.066e-14`
- worst unordered reconstruction residual: `1.540e-14`
- minimum sampled pair-capacity margin: `8.889e-19`
- minimum primary owner half-split margin: `3.091e-04`
- minimum dominant-shell margin over clean lower: `3.476e-02`
- minimum diffuse-entropy margin: `9.525e-02`
- minimum entropy-shell tradeoff margin: `8.585e-03`
- maximum joint primary owner count: `2`
- maximum simultaneous quarter-cut corollaries: `1`

No packet synchronization, no coherent-frequency support fiction, no aggregate pressure-mass reset, and no Navier--Stokes global-regularity conclusion are asserted.
