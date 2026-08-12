# Actual Galerkin NS mixed-fate signed-cell audit

Status: **EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_MIXED_FATE_HARD_CELL_AUDIT__ACTUAL_GOOD_BAD_NEGATIVE_EDGE_WORK__FULL_SIGNED_CELL_RESERVATION_IDENTITIES**.

A real divergence-free Fourier polynomial is evolved by the repository's 2/3-dealiased Galerkin Navier--Stokes RK4 path.  At the selected physical child `z=(7, 6, 5)`, the initial state contains a near-extremal signed-good parent pair plus positive nonforward and negative-work parent pairs.  They are deliberately compressed into the same deterministic hard product cell only to stress the already-fixed edge-space Hahn provenance.

- resolution/cutoff: `24` / `7`
- steps/snapshots: `24` / `4`
- mixed-fate snapshots: `4`
- snapshots with good/bad/negative work: `4` / `4` / `4`
- initial maximum good signed efficiency: `0.999996646308`
- initial unresolved mixed good dW+ mass: `68.0888945495`
- worst actual-NS signed reconstruction relative residual: `6.439e-16`
- worst `T=g+b-n` hard-cell relative residual: `8.032e-17`
- worst `T-b=g-n` reservation-certificate relative residual: `1.004e-16`
- stage-zero first-time failures: `0`
- geometry-good marking promotions: `0`

The total-variation upper reported by the probe is only a sign/provenance diagnostic and is **not** substituted for the sharp continuous Young bound.  The analytic handoff theorem remains the statement about the actual full signed hard-cell trilinear form and a separately certified Young norm upper.  No global-regularity claim is made.
