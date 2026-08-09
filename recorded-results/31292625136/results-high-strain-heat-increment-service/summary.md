# High strain becomes heat-increment coherent service

Status: **EXACT_HIGH_STRAIN_TO_HEAT_INCREMENT_COHERENT_SERVICE__CRITICAL_RESOLVED_ANCESTOR_FRACTION_RETAINED__OLD_POOL_ROUTING_REMAINS**.

Use no arbitrary spatial packet.  Let `H_N` be the three-dimensional Navier--Stokes heat kernel at the block's intrinsic parabolic time `theta_N=1/(2N^2)`.  Equivalently `r` is centered Gaussian with covariance `N^-2 I`, so

`E exp(-i xi.r)=exp(-|xi|^2/(2N^2))`.

Plancherel gives the exact heat-defect identity

`E_H ||delta_r V||_2^2 = int 2(1-exp(-|xi|^2/(2N^2))) |Vhat(xi)|^2 dxi`.

Because the strict transporter has `supp Vhat subset B_(N/4)`, put `x=|xi|^2/(2N^2)<=1/32`.  The elementary integral identity `1-e^(-x)=int_0^x e^(-s)ds` gives

`e^(-1/32) x <= 1-e^(-x) <= x`.

Hence pointwise on the full resolved support

`e^(-1/32)||grad V||_2^2 <= N^2 E_H||delta_r V||_2^2 <= ||grad V||_2^2`.

After integrating a child lifetime, the positive heat-increment service

`S_heat=N^3 int dt E_H||delta_r V||_2^2`

satisfies

`e^(-1/32) D_V <= S_heat <= D_V`.

Thus the existing high-strain gate `D_V>=32 pi^2/(75c)` immediately forces a fixed positive **physical-space increment service**.  The Gaussian is not a replacement causal law: it is the heat semigroup already intrinsic to the viscous PDE, used only to expose where resolved gradient activity lives.

The preceding resolved-ancestor theorem says at least half the actual `D_V` law lies on shell-time atoms with critical mass `mu_j>=32pi^2/(75c^2)`.  Since the heat multiplier is pointwise between `e^(-1/32)` and `1` times the gradient multiplier, the same good shell-time set carries at least

`(1/2)e^(-1/32) = 0.484616617238`

of the **entire heat-increment service law**.  The frequency ancestor and spatial service are therefore simultaneous physical marks, not two unrelated pigeonholes.

For each good shell, each heat displacement `r`, and any normalized affine coherent window, Moyal gives an exact positive phase-space disintegration of `||delta_r P_jV||_2^2`.  Translation covariance

`V_g(delta_r f)(X,k)=e^(-ik.r)V_gf(X-r,k)-V_gf(X,k)`

shows that every cell atom is a real coherent edge between neighborhoods separated by the actual Brownian/heat displacement.  No coherent cell is selected by argmax and no global shell mass is divided among a guessed packet count.

Stress: `50000` spectral heat-defect/ancestor states and representative coherent identities
- minimum spectral lower margin: `2.167e-08`
- minimum spectral upper margin: `2.977e-12`
- minimum retained-fraction margin above `e^(-1/32)/2`: `6.353e-06`
- worst relative discrete Moyal residual: `6.278e-16`
- worst relative coherent increment-covariance residual: `8.016e-15`

This supplies the missing **spatial/coherent entrance** for the high-strain resolved-ancestor law.  What remains is material routing of these positive coherent edges: old--old service must enter the existing reservoir half-life, old--new edges must be physical relink/interface, and genuinely new--new edges must create coherent ancestry/service.  That routing must be proved for this dissipation-seeded measure rather than assumed from the SGS source theorem.  `D_V` remains an `O(1/N)` physical cost, not an additive reset.  No global-regularity claim is made.
