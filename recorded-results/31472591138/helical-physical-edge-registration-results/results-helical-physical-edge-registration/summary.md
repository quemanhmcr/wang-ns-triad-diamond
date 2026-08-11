# Helical physical edge registration: native NS capacity before Young

Status: **EXACT_HELICAL_PHYSICAL_EDGE_REGISTRATION__DIRECT_LERAY_CURL_EQUALS_WALEFFE__UNORDERED_PARENT_PAIR_FACTOR_FOUR__NATIVE_MODAL_CAPACITY__SIGNED_UPPER_PROGRESS_EQUALS_A_J_C**.

For one actual unordered helical parent pair `x,y` feeding `z=x+y`, the nonlinear Navier--Stokes source is read directly as

`F_z=P_z(u_x x omega_y + u_y x omega_x)`.

The repository Waleffe convention gives exactly

`<h_z,F_z>=2(s_x|x|-s_y|y|)conj(g_e)a_xa_y`.

Because physical child energy differentiates as `2 Re(conj(a_z)<h_z,F_z>)`, the signed work is

`T_e=4(s_x|x|-s_y|y|) Re[conj(a_z)conj(g_e)a_xa_y]`.

Define the native modal interaction capacity

`A_e=4|z||a_xa_ya_z|`,

the geometric upper-progress multiplier

`J_e=log_+(|z|/p_top)|s_x|x|-s_y|y|| |g_e|/|z|`,

and the signed phase/orientation alignment `c_e in [-1,1]`.  Then on the **same physical Fourier event**,

`T_e log_+(|z|/p_top)=A_e J_e c_e`.

`A_e` is not a Young norm product.  It is the available modal interaction amplitude of this triad, including the exact factor two from the unordered parent orbit and factor two from physical energy differentiation.  Young/Christ enters only downstream when one asks whether a block of these real edges is near saturation.

Stress: `50000` random physical helical triads/amplitudes
- worst direct-Leray / Waleffe coefficient relative residual: `3.047e-14`
- worst Leray-free child pairing residual: `1.036e-14`
- worst `T log = A J c` relative residual: `1.265e-13`
- worst unordered-parent swap residual: `8.394e-12`
- worst helical phase-gauge alignment residual: `1.277e-15`
- worst uniform-wavevector scale invariance residual: `3.978e-12`
- maximum sampled `J/J_*`: `0.999042728071`
- minimum phase alignment on sampled positive forward-work edges: `2.135e-06`
- positive forward-work samples: `16543`
- nonforward samples (zero upper-progress multiplier, retained rather than mislabeled signed-good): `16950`

This theorem registers one physical edge only.  It does **not** yet construct the continuum edge measure or claim that generic/nonforward HH is signed-good.  No Navier--Stokes global-regularity conclusion is asserted.
