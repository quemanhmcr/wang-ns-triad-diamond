# Smooth quadratic-carrier energy with observer-gauge quotient

Status: **EXACT_SMOOTH_QUADRATIC_CARRIER_INTERFACE__Q2_ENERGY_LAW__COMMON_GAUGE_QUOTIENT_BEFORE_PHYSICAL_RELINK__SYMMETRIC_WORK_EXISTING_STRAIN__COEFFICIENT_OBSTRUCTION_ENERGY_REENTRY**.

For a smooth self-adjoint PDE carrier `w=Q u`, the native carrier energy is

`E_Q=||Q u||_2^2=<u,Q^2u>`.

Smooth roles are therefore completed by `sum_a A_a^2=I`; `Q` and `I-Q` are not treated as hard projectors.  Direct Navier--Stokes differentiation gives

`d||Q u||_2^2/dt + 2 nu ||grad Q u||_2^2 = <u,dot(Q^2)u> - 2 Re<Q^2u,B(u,u)>`.

After resolved repartition the native interface is

`J_Q=<u,dot(Q^2)u>-2 Re<Q^2u,L_Vu>`,

and exactly equals outer commutator work minus the diagonal `L_V` work of `Q u`.  The commutator is never interpreted alone.

The crucial quotient is prior to ownership.  A moving square partition is accepted as transported PDE gauge only when one common skew generator `G` satisfies

`dot A_a + [G,A_a] = 0`

for every role.  Writing the actual skew resolved operator as `K=G+K_phys`, the `dot A` work cancels the common-`G` work exactly.  Only `K_phys` remains as antisymmetric conservative physical relink; `S` remains the existing symmetric strain/deformation work.  Arbitrary observer motion of the analysis windows fails this transport identity and cannot be Hahn-routed to a causal owner.

A coefficient obstruction remains only an interval locator.  Actual `Q^2` carrier energy and a **gauge-quotiented** interface-work certificate must pass through the physical-energy gate before inheritance, high strain, HH generation, physical relink provenance or strain ownership is named.  The coefficient impulse magnitude never enters a work threshold.

Stress: `50000` transported square-partition/interface/PDE/reentry states
- worst quadratic-partition residual: `4.584e-15`
- worst differentiated-partition residual: `1.670e-14`
- worst common-gauge transport residual: `1.594e-16`
- worst gauge-work cancellation residual: `1.716e-16`
- worst native/outer recombination residual: `6.011e-16`
- worst physical-relink total residual: `8.889e-16`
- worst strain reconstruction residual: `1.894e-15`
- worst pair antisymmetry/symmetry residual: `3.394e-16`
- worst pair row-sum residual: `2.274e-15`
- worst direct carrier-energy residual: `7.105e-14`
- worst resolved-repartition residual: `2.451e-13`
- worst hard-event plateau registration residual: `9.366e-16`
- linear-complement counterexample defect: `1.000000`
- worst quadratic-complement skew residual: `6.895e-16`
- arbitrary observer-motion rejections: `1`
- minimum clean HH-generation margin: `9.244e-03`
- energy-reentry branches: `{'material_energy_inheritance': 12460, 'physical_high_high_transfer_generation': 16359, 'smooth_interface_physical_work': 8632, 'high_strain_critical_dissipation': 12549}`

This theorem closes the local Q-vs-Q^2 mismatch without charging observer-selected window motion as physics.  It remains complementary to the hard event-role donor/circulation quotient.  It does not prove global owner termination, UV closure, or Navier--Stokes regularity.
