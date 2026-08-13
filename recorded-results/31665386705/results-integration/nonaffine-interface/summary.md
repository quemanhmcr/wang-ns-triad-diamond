# Nonaffine hard-projector role-interface work

Status: **EXACT_NONAFFINE_ROLE_INTERFACE_SPLIT__SKEW_TRANSFER_CONSERVATIVE__SYMMETRIC_WORK_IS_STRAIN_PROVENANCE**.

For the nonaffine resolved linearized operator write `L=K+S`, `K*=-K`, `S*=S`.  For any complete orthogonal **event-role** partition `u=sum_a w_a`, the Heisenberg interface work splits exactly.

The skew/advection+rotation piece has pair flux

`T_ab=-2 Re <w_a,K w_b>`,

with `T_ab=-T_ba`.  Its row sums are the individual role-interface works, so the total skew interface work is **exactly zero**.  It is conservative physical role redistribution/relinking, never energy generation.

The symmetric piece has

`D_ab=-2 Re <w_a,S w_b>=D_ba`.

Its row sums are the symmetric interface works and

`sum_a I_a^S = 2 sum_a Re<w_a,S w_a> - 2 Re<u,S u>`.

Thus the symmetric hard-role interface is precisely the off-diagonal part of the **same physical strain/deformation work** already present in the resolved transporter.  It delegates to coherent deformation / strain / critical `D_V`; it is not a new currency and not representation `Xi`.

Stress: `50000` random complex Hilbert-space partitions/operators
- worst binary identity residual: `4.801e-16`
- worst skew antisymmetry residual: `1.409e-16`
- worst symmetric symmetry residual: `1.684e-16`
- worst skew row-sum residual: `6.609e-16`
- worst symmetric row-sum residual: `4.230e-16`
- worst total skew-interface residual: `8.172e-16`
- worst symmetric global strain-balance residual: `1.174e-15`

This theorem is deliberately a hard-projector lemma.  It must not be applied directly to the non-idempotent smooth PDE envelope.  The companion smooth quadratic-carrier theorem reads the propagated energy at `Q^2`, recombines the outer commutator with diagonal resolved-role work, then quotients the certified common transport gauge, and only then assigns residual skew work to conservative relink while symmetric work remains existing strain.  Quantitative first stopping remains governed by actual energy work, never raw coefficient impulse.  No global-regularity claim is made.
