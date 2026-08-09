# Coherent objective-source owner compiler

Status: **EXACT_COHERENT_OBJECTIVE_SOURCE_OWNER_COMPILER__LOCAL_DV_AND_VISCOSITY_TO_CRITICAL_SHELL__SGS_TO_COHERENT_SERVICE__PRESSURE_TO_SGS_OR_RESERVOIR__NO_PACKET_SYNCHRONIZATION**.

The source branch of the coherent service-or-flat gate does not need a new packet object.  Group the exact coherent averaged strain source by its physical owner:

`local_DV + pressure + SGS + viscosity`.

If `A_obj` is the objective variation action on a scaled lifetime `c`, their positive scaled weights satisfy `sum Sigma_r >= A_obj/c`, so at least one owner has `Sigma_r>=A_obj/(4c)`.  Exact ties are retained jointly; there is no theorem-name priority.

The local coherent quadratic/Reynolds owner is already bounded by `C_local D_V`, and the viscous owner obeys `rho_nu<=nu sqrt(d_V)/1500`; both therefore feed the generic critical-shell-to-own-scale-service theorem as **recursive** resolved dissipation, never as an additive reset.

For the objective SGS owner, the clean collision `||R||_(3/2)>=380 rho_R`, Germano `3/2` power and coherent-service `2/3` power cancel exactly:

`Y_R >= C_Y rho_R`,

`C_Y = 380/[g1(1+g1)(C_LP C_B)^2]`.

Thus integrated SGS source weight produces integrated coherent square service with no persistence hypothesis and no affine-radius packet.  The positive service law routes exactly to high-frequency dissipation, old-pool capacity, selected-interface `Xi`, a dominant fresh critical shell, ancestry entropy, or same-ancestry cycle.  High-frequency dissipation is **not** renamed resolved `D_V`.

For pressure,

`rho_P <= mu_V/5700 + ||R||_(3/2)/380`

integrates to the honest alternative

`int mu_V >= 2850 Sigma_P`  OR  `int ||R||_(3/2) >= 190 Sigma_P`.

The stress alternative is exactly an effective SGS source weight `Sigma_P/2` and enters the same coherent-service route.  The resolved low-pass mass alternative is **not** promoted to a generic critical shell.  It remains objective pressure-Hessian reservoir occupation.  On a supplied signed-good low-strain lineage, one fixed materially reused low-low pair has the derivative-correct coefficient ratio `4084101/20971520<1/5`, hence total future fixed-pair capacity `<5/4` of generation zero.  Persistent pressure-Hessian service must relink pairs, fragment into component entropy/cycle, leave low strain, or use the SGS branch.  The separate H1 pressure-third source retains its own `<1/3` theorem.

Stress: `50000` source-owner states
- worst SGS closed-form relative residual: `1.207e-15`
- minimum owner-pigeonhole margin: `2.066e-03`
- minimum pressure split identity margin: `-1.110e-16`
- minimum local-DV identity margin: `-1.110e-16`
- minimum viscous-Cauchy identity margin: `0.000e+00`
- maximum sampled joint owner count: `4`

The resulting architecture is source-native: `local/viscous -> D_V -> critical shell`, `SGS -> coherent service`, `pressure -> SGS service or low-frequency reservoir`.  No packet synchronization theorem and no uniform finite resource are inserted.  Final continuum master assembly and supplier-specific signed-good scale geometry remain separate.  No Navier--Stokes global-regularity conclusion is asserted.
