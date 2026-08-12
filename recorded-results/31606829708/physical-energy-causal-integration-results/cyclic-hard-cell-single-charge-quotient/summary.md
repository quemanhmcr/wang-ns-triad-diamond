# Cyclic hard-cell single-charge quotient

Status: **EXACT_CYCLIC_HARD_CELL_SINGLE_CHARGE_QUOTIENT__RESTRICTED_DW_MINUS_ROWS__CANONICAL_DW_PLUS_GOOD_BAD_COLUMNS__OVERLAPPING_RECIPIENT_SUPPORT_NO_DOUBLE_CHARGE__COARSE_SELF_LOOPS_ZERO_RECURSION_DEPTH**.

The certified cyclic donor kernel is pushed through a deterministic hard Fourier/helicity map without changing cause.  For donor cell `C`, recipient cell `D`, and inherited recipient fate `F`, the positive table `K(C,D,F)` has row marginal `pi_#dW-` and column/fate marginal the already-canonical `pi_#dW+`.  Incoming atoms from several donor cells therefore disintegrate one recipient charge; they are not several causal charges.

Stress: `50000` physical closed triads
- resolved / numerically unresolved: `50000` / `0`
- one-donor / two-donor resolved cases: `25007` / `24993`
- overlapping-recipient-charge cases: `16112`
- coarse-self-loop cases: `22120`
- worst balance native residual: `1.586e-16`
- worst donor marginal native residual: `9.085e-17`
- worst recipient marginal native residual: `9.129e-17`
- worst recipient fate-partition native residual: `0.000e+00`
- worst restricted-donor pushforward native residual: `9.522e-17`
- signed-good good/bad recipient masses: `0.540401435212` / `0.17237839362`
- maximal-coarsening self-loop fraction: `1`
- generic two-donor -> one recipient single-charge anti-theorem: `True`

A geometry-bad recipient remains on the existing `TRANSFER_WORK_LOSS` stage-zero recursion route, but its modal energy remains in Navier--Stokes and may participate later.  A geometry-good recipient remains only Young-eligible.  A coarse self-loop is same-time physical redistribution with zero recursion depth and no scale progress.  The theorem does not say negative work pays failed good work and does not introduce any between-time deposit/withdrawal matching.  No global-regularity claim is made.
