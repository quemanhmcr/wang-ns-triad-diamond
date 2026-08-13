# Actual Fourier--Galerkin Navier--Stokes referee: resolved-contact native binding

Status: **ORTHOGONAL_FOURIER_GALERKIN_NS_RESOLVED_CONTACT_NATIVE_BINDING__SAME_PHYSICAL_DWPLUS_UNDER_TWO_ADMISSIBLE_CUTOFF_REPARTITIONS__ACTUAL_U_EQUALS_V_PLUS_H_AND_SAME_EDGE_VH_HH__MIXED_KS_ADJOINT_REFEREE**.

The referee evolves the same real divergence-free six-mode Navier--Stokes data independently of the analysis cutoff.  At every snapshot it reconstructs canonical cyclic `dW+`, then reads two different admissible smooth `u=V+h` decompositions of that same physical interaction.  Both decompositions must reconstruct the identical full NS nonlinear source.  On interior contact one profile leaves genuine mixed+HH transition while the other resolves the same low parent completely; canonical cause and donor provenance do not change.

- maximum representation spread: `9.843e-16`
- worst full PDE bilinear repartition residual: `3.359e-16`
- worst cutoff-repartition gauge residual: `3.231e-16`
- worst high-shell low-low leakage: `0.000e+00`
- worst same-edge signed V/h repartition residual: `3.863e-17`
- worst K/S signed identity residual: `1.288e-17`
- worst K skew-pair residual: `1.610e-18`
- worst S symmetric-pair residual: `6.439e-18`
- worst canonical K/S positive-cover defect: `1.288e-17`
- worst Galerkin energy-balance residual: `8.560e-16`

This is a referee for the physical identities and type barriers, not a substitute for the continuum proof and not a Navier--Stokes regularity claim.
