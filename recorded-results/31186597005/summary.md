# Smooth-SGS two-level affine packet equation

Status: **CERTIFIED** for the good-core transport separation.

- signed-good parent/child ratio is strictly above `3/5`;
- choose a smooth common transport low-pass supported in `|xi|<=N/4`;
- its low-low product is supported in `|xi|<=N/2`, so it cannot feed any selected role directly;
- the selected role equation therefore splits exactly into low-high transport and high-high transfer;
- a Fourier packet multiplier transported by the affine dual flow has zero Heisenberg residual
  `partial_t m + [(Ax).grad,m(D)] = 0`;
- microscopic full-velocity role forcing and macroscopic resolved-SGS boundary transport are recorded as disjoint ledgers, preventing `RU`/pressure double counting.

Stress checks: `50000`
- worst affine Heisenberg residual: `0.000e+00`
- worst transported-symbol finite-difference residual: `5.718e-08`
- minimum sampled low-low/role support gap: `0.101207`
- source taxonomies disjoint: `True`
