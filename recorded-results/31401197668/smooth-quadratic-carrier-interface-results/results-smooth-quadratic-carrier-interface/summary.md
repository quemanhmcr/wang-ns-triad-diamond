# Smooth quadratic-carrier interface and physical energy reentry

Status: **EXACT_SMOOTH_QUADRATIC_CARRIER_INTERFACE__Q2_ENERGY_LAW__SQUARE_PARTITION_CONSERVATIVE_RELINK__SYMMETRIC_WORK_EXISTING_STRAIN__COEFFICIENT_OBSTRUCTION_ENERGY_REENTRY**.

For a smooth self-adjoint PDE carrier `w=Q u`, the native physical energy is

`E_Q=||Q u||_2^2=<u,Q^2u>`.

Thus smooth carrier roles are completed by a quadratic analysis partition `sum_a A_a^2=I`, not by treating `Q` and `I-Q` as hard projectors.  The direct Navier--Stokes energy identity is

`d||Q u||_2^2/dt + 2 nu ||grad Q u||_2^2 = <u,dot(Q^2)u> - 2 Re<Q^2u,B(u,u)>`.

After the exact resolved repartition with `V` and `h=u-V`, a selected outer role has low--low work zero, physical `q^2`-weighted HH work, and native interface work

`J_Q=<u,dot(Q^2)u>-2 Re<Q^2u,L_Vu>`.

This is also exactly the work of `(dot Q+[L_V,Q])u` in the outer-role equation **minus** the diagonal `L_V` work of `Q u`.  A smooth commutator is therefore never interpreted in isolation.

For `L_V=K+S`, `K*=-K`, `S*=S`, the moving plus skew part of `J_Q` is conservative relinking.  Across the complete quadratic partition its total is zero, and its synthesis pair flux on `A_a^2u` is antisymmetric.  The symmetric rows are the same physical resolved strain/deformation work and reconstruct `-2 Re<u,S u>` exactly.  Neither part is a new source currency or representation `Xi`.

The linear complement `I-Q` is forbidden for non-idempotent `Q`: its total skew-interface defect is `4 Re<Q(I-Q)u,K u>`.  Choose a smooth angle partition `Q=cos(theta)`, `R=sin(theta)` so `Q^2+R^2=I`.  Hard event registration remains exact on the plateau `QP=P`; no hard boundary is propagated.

Finally, a Duhamel/interface coefficient obstruction is only an interval locator.  At its first hit, actual terminal carrier energy and positive native interface work reenter the existing physical-energy gate.  The gate returns inherited energy, high strain, actual HH generation, or physical interface work.  Only the last branch is Hahn-routed to conservative relink or existing strain.  The coefficient impulse magnitude is never used as a causal weight.

Stress: `50000` smooth-partition/interface/PDE/reentry states
- worst quadratic-partition residual: `5.788e-15`
- worst differentiated-partition residual: `3.474e-15`
- worst native/outer recombination residual: `8.934e-16`
- worst conservative-relink residual: `8.741e-16`
- worst strain reconstruction residual: `2.923e-15`
- worst pair antisymmetry/symmetry residual: `4.258e-16`
- worst synthesis-pair row residual: `2.923e-15`
- worst direct carrier-energy residual: `7.854e-14`
- worst resolved-repartition residual: `1.492e-13`
- worst hard-event plateau registration residual: `9.515e-16`
- linear-complement counterexample defect: `1.000000`
- worst quadratic-complement skew residual: `8.006e-16`
- minimum clean HH-generation margin: `7.828e-03`
- energy-reentry branches: `{'material_energy_inheritance': 12730, 'smooth_interface_physical_work': 12378, 'physical_high_high_transfer_generation': 12469, 'high_strain_critical_dissipation': 12423}`

This theorem closes the algebraic smooth-envelope/projector mismatch and the local physical-energy reentry of coefficient obstructions.  It is complementary to hard event-role donor/circulation quotienting.  It does not prove global owner termination, UV closure, or Navier--Stokes regularity.
