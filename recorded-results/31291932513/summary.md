# Resolved-cutoff repartition relay 31291932513

- branch: `research/smooth-material-carrier-relay`
- SHA: `4fcdc8c6c3cbb5492ad6acce193dd8754f62c912`
- conclusion: `success`
- full pytest: `492 passed`
- cutoff/scale stress: `50,000`
- worst old/new cutoff identity residual: `1.099e-14`
- worst residual against full `-Q B(u,u)`: `8.563e-15`
- minimum renewed low-low support gap: `3.242566e-01`
- minimum sampled lifetime-window margin: `5.371e-06`
- companion smooth material-carrier relay and outer moving-role stresses passed on the same SHA.

Run `31291900197` on precursor SHA `def38675...` was a certificate-fixture failure only: all `492` tests passed, and execution stopped before stress because Arb object equality was incorrectly used to check the exact rational lifetime endpoints. Commit `4fcdc8c...` changed only that check to exact `Fraction` arithmetic; no theorem formula, constant, or stress inequality changed.

The theorem proves that changing the resolved cutoff at parent-scale renewal is an exact repartition of the same Navier--Stokes quadratic interaction and introduces no cutoff-switch source/currency. Universal renewal/exhaustion for all recursive routes remains open; no global-regularity claim is made.
