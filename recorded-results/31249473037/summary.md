# Adjoint Kelvin--Duhamel causal gate

Status: **CERTIFIED_ADJOINT_DUHAMEL_AND_CAUSAL_TIME_GEOMETRY**.

Put the common low--high affine/Kelvin transport **and bulk viscosity** into the linear objective generator `G(t)`.  Along a Kelvin characteristic let the selected transverse coefficient obey

`c_dot = G c + F_HH + R_class`,

where `F_HH` is actual high--high generation and `R_class` contains only already classified cross-cell / moving-projector / H1-H3 / window / profile residuals.  Pressure is absent because the role equation has already been Leray projected.

For the backward adjoint dual

`psi_dot = -G^* psi`,

one has exactly

`d <psi,c>/dt = <psi,F_HH+R_class>`.

Thus on a child slab

`z_1 = z_0 + I_HH + I_R`.

If the terminal dual is chosen along the terminal coefficient, `|z_1|=||c(t_1)||=:A`, and the exact triangle inequality gives the clean causal alternative

- inherited material coefficient `|z_0|>=A/4`; or
- classified residual `|I_R|>=A/4`; or
- genuine high--high generation `|I_HH|>=A/2`.

No common affine strain, pressure or bulk viscosity is counted again as generation.

For the generated branch, decompose the high--high Duhamel impulse into quadratic parent-pair atoms `z_alpha`.  A **single** phase aligned with the total impulse gives

`sum [Re(conj(phase) z_alpha)]_+ >= |I_HH|`.

After normalization this is a positive causal-generation law on same-time quadratic parent-pair events.  No pointwise persistence is required.  This is an amplitude-generation law; identification with the physical positive child-energy transfer law is made only on the already-certified near-extremal phase-locked core.

Signed-good scale geometry also gives

`64/25 < T_parent/T_child < 25/9`

for natural parabolic lifetimes.  A half-child-slab carrying at least half of the positive generation mass therefore has parent **natural** backward windows with common overlap longer than `103/50 T_child`.  This is geometric synchronization only; actual packet persistence on that whole common window is not claimed.

Stress: `50000` instantaneous/phase checks plus `min(samples,5000)` exact block-exponential Duhamel histories
- worst instantaneous adjoint-pairing residual: `3.332e-14`
- worst exact piecewise Duhamel residual: `1.461e-15`
- minimum high-high generation margin: `4.055e-01`
- minimum positive aligned-mass margin: `-3.553e-15`
- minimum half-slab pigeonhole margin: `3.518e-06`
- branches: `{'material_inheritance': 5001, 'classified_residual': 1, 'high_high_generation': 1}`
