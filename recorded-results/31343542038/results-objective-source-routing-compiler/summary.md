# Coherent objective-source owner compiler

Status: **EXACT_COHERENT_OBJECTIVE_SOURCE_OWNER_COMPILER__LOCAL_DV_AND_VISCOSITY_TO_CRITICAL_SHELL__SGS_TO_COHERENT_SERVICE__PRESSURE_TO_SGS_OR_ENTROPY_WEIGHTED_CRITICAL_SHELL__AGGREGATE_MUV_DIAGNOSTIC_ONLY__NO_PACKET_SYNCHRONIZATION**.

The source branch of the coherent service-or-flat gate does not need a new packet object.  Group the exact coherent averaged strain source by its physical owner:

`local_DV + pressure + SGS + viscosity`.

If `A_obj` is the objective variation action on a scaled lifetime `c`, their positive scaled weights satisfy `sum Sigma_r >= A_obj/c`, so at least one owner has `Sigma_r>=A_obj/(4c)`.  Exact ties are retained jointly; there is no theorem-name priority.

The local coherent quadratic/Reynolds owner is already bounded by `C_local D_V`, and the viscous owner obeys `rho_nu<=nu sqrt(d_V)/1500`; both therefore feed the generic critical-shell-to-own-scale-service theorem as **recursive** resolved dissipation, never as an additive reset.

For the objective SGS owner, the clean collision `||R||_(3/2)>=380 rho_R`, Germano `3/2` power and coherent-service `2/3` power cancel exactly:

`Y_R >= C_Y rho_R`,

`C_Y = 380/[g1(1+g1)(C_LP C_B)^2]`.

Thus integrated SGS source weight produces integrated coherent square service with no persistence hypothesis and no affine-radius packet.  The positive service law routes exactly to high-frequency dissipation, old-pool capacity, selected-interface `Xi`, a dominant fresh critical shell, ancestry entropy, or same-ancestry cycle.  High-frequency dissipation is **not** renamed resolved `D_V`.

For pressure, the coarse estimate `rho_P<=mu_V/5700+||R||_(3/2)/380` is retained only as a diagnostic.  The canonical compiler now uses the actual Frobenius-dual positive source law

`rho_P <= [r_SGS]_+ + sum_(a<=b)[p_ab]_+`.

Thus positive SGS pressure source or the resolved unordered pair law carries at least `Sigma_P/2`, with exact ties joint.  The SGS owner uses its **actual** positive source weight `r` and gives `int||R||_(3/2)>=380r`, hence direct coherent service.  Every resolved pair owner satisfies

`mu_child exp(H2_pair) >= 320 Sigma_P/c`

and therefore enters the generic critical-shell first-stop theorem.  On a full no-hit natural survivor the compiler records the corresponding own-scale service, but it does not turn that conditional service into an unconditional event.  `H2_pair` is only the logarithmic weakening of the shell seed; it is neither a causal child-energy probability nor a separate stop.

The old fixed-material-pair ratio `4084101/20971520<1/5` remains a valid optional reuse refinement after material sidecars are attached; it is no longer the pressure renewal entrance.  The separate H1 pressure-third source retains its own `<1/3` theorem.

Stress: `50000` source-owner states
- worst SGS closed-form relative residual: `1.215e-15`
- minimum owner-pigeonhole margin: `2.748e-03`
- minimum pressure diagnostic split margin: `-1.110e-16`
- minimum pressure entropy-shell margin: `-4.547e-13`
- minimum pressure full-survivor service registration margin: `-4.441e-16`
- minimum local-DV identity margin: `-1.110e-16`
- minimum viscous-Cauchy identity margin: `0.000e+00`
- maximum sampled joint owner count: `4`
- maximum sampled joint pressure owner count: `1`

The resulting architecture is source-native: `local/viscous -> D_V -> critical shell`, `SGS -> coherent service`, `pressure -> actual SGS service OR entropy-weighted critical shell`.  No packet synchronization theorem and no uniform finite resource are inserted.  Final continuum master assembly and supplier-specific signed-good scale geometry remain separate.  No Navier--Stokes global-regularity conclusion is asserted.
