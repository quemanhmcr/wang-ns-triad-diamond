# Experimental and certification log

This file is the home for **execution evidence**, not theorem architecture. Keep GitHub Actions run IDs, stress counts, actual-NS probes, representation residuals, artifact digests, and failure/correction lineage here. `RESEARCH_LEDGER.md` contains only the compact physical/theorem state.

Historical detail through the radial theorem is preserved verbatim in `docs/history/RESEARCH_LEDGER_full_before_compaction_2026-08-13.md`; earlier development is in `docs/history/RESEARCH_LEDGER_history_through_2026-08-10.md`. Raw certified outputs live under `recorded-results/<run-id>/`.

## Current exact-SHA frontier index

- helical physical edge registration: `6d2b721be579e299c42f4ad370448e5336f3f43c`
- continuum helical edge-measure registration: `6bc9099190048a8796b54fb2f4782314b699bd4b`
- canonical positive edge-work routing: `d78a3f4223c16f541548e7628cd73f70e6fdee6c`
- mixed-fate reserved Young/Christ handoff: `2600bcbe7db56f54af58d1e7685d1246ad3d7bee`
- cyclic helical-triad donor kernel: `647065ab1562e15da84c97d9be5265c2e43dc84d`
- cyclic hard-cell single-charge quotient: `964b72b0a614bb44c88a42d41160f11bd550d97a`
- helical mode-set energy continuity: `a39d502d9312ac3bd6613a780d60b22e88790863`
- radial spectral crossing layer cake: `667e687cb740a77df944753c575f581abad14199`
- hard-tail true-upward supply: `d064bc4d780d9c90d36e64c7b84c3b771b74c896`

For complete historical run IDs and metrics of the earlier entries, use the archived full ledger above.

## Hard-tail true-upward supply — exact certification

Exact theorem SHA: `d064bc4d780d9c90d36e64c7b84c3b771b74c896`.

Serious gates on that exact SHA:

- dedicated `31615787742` — **success**; full suite `892` tests; `75,000` true-upward support states plus radial/high-tail dependencies and evolved Navier–Stokes probes;
- independent audit `31615787779` — **success**; `100,000` true-upward support states, `100,000` radial dependency states, and longer actual-NS audit;
- physical-energy causal integration `31615787727` — **success**; full `892`-test suite and `50,000` true-upward states inserted after radial crossing in the complete causal spine.

Stored artifact trees and deterministic digests:

- `recorded-results/31615787742/` — `sha256:95147ccdcd456a1b8362c3c62846b88100574a9710806b7f475f974f01099f59`
- `recorded-results/31615787779/` — `sha256:72dc78c0c7bae44fb3eae5bd5a2dd8875da100dba1b9cb9fb6e93b80f0871d6c`
- `recorded-results/31615787727/` — `sha256:1c38edfed2b05843bb0eff21629390c923dc45b6ca575bf9aa912cc11033ae3e`

Independent `100,000` support audit: `99,999` resolved triads and one fail-closed unresolved triad; `100,032` upward atoms; `67,003` pure-UV and `33,029` resolved-contact atoms; `26,414` deep atoms; zero pure-UV nonfirst/deep violations; maximum sampled recipient shell index `6`; maximum pure-UV parent/shell ratio `1.49628423384 < 3/2`; all upward/radial-binding/common-unit partition residuals were exactly `0` in the reported stress.

The evolved deep-contact fixture retained both physical phenomena simultaneously: initial upward work `0.251755958403`, of which deep-contact work was `0.164487241289` and first-shell pure-UV work `0.0872687171139`; all `33` audited output snapshots retained positive deep upward work. This is why the theorem is atomwise rather than a whole-triad binary fate.

The load-bearing six-mode orthogonal Fourier–Galerkin NS audit reconstructed the **full** cyclic boundary law of the finite PDE. On FFT representations `20,24`, worst `Phi_up-Phi_down` versus direct signed tail nonlinear-work residual was below `5.4e-17`, tail continuity residual below `1.8e-15`, while internal high-tail circulation remained separately nonzero. The equiradial amplitude-scaling anti-test changed internal circulation `0.0777635049765 -> 77.7635049765` under amplitude `×10`, while true upward supply stayed exactly `0 -> 0`.

Failure lineage: initial candidate `dd3e7c730e9973eecd0d9a7168db10094abd0f0d` had two focused fixture errors before theorem stress: a test passed the `(triad, certificate)` tuple rather than the triad object, and the engineered deep-PDE fixture incorrectly demanded an exclusive whole-triad contact fate. Actual PDE bookkeeping showed that deep-contact and first-shell pure submeasures can coexist. Repair `d064bc4...` corrected only those fixtures/ontology expectations; no theorem identity, support threshold, causal unit, sample count, PDE equation, or numerical tolerance was weakened.

## Logging rule going forward

When a theorem is certified, record here: exact theorem SHA; dedicated/audit/integration runs; artifact digests; only the few metrics that could falsify the theorem; and failure lineage that changed understanding or certification. Keep theorem statements, owner semantics, and the current frontier in `RESEARCH_LEDGER.md`. Do not duplicate large CI transcripts into the root ledger.
