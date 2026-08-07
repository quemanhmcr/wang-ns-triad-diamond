# Scale holonomy certificate

This note records an analytic obstruction suggested by the numerical experiment.
It is a statement about the reduced near-extremal triad model, not yet about the
full Navier–Stokes PDE.

Let `ell_v = log |k_v|`. A single near-optimal forward triad is expected to obey

1. its two parent scales are nearly equal;
2. its child lies approximately `gamma_*` log-frequency units above the parents,
   where `gamma_* = log R_*` and `R_* ≈ 1/0.610904` in the symmetric extremal family.

For the reuse motif

- `a+b -> m`,
- `m+c -> d`,
- `b+c -> n`,

define residuals

```
r1 = ell_a - ell_b
r2 = ell_m - (ell_a+ell_b)/2 - gamma_*
r3 = ell_m - ell_c
r4 = ell_b - ell_c
```

They satisfy the exact algebraic identity

```
r2 - r3 + (1/2) r1 + r4 = -gamma_*.
```

Hence

```
gamma_* <= |r2| + |r3| + (1/2)|r1| + |r4|.
```

At least one residual is therefore at least `gamma_*/3.5`. With
`gamma_* ≈ log(1/0.610904) ≈ 0.493`, this lower bound is about `0.141`.

That stability step is now complete in the finite-dimensional edge model.  The
Arb certificate proves, in the log-scale Hodge coordinates,

\[
\operatorname{Def}_e\ge\frac1{50}|u|+v^2
\ge\frac12(r_p^2+r_q^2)
\]

inside the near-extremal rectangle and `Def_e>=1/100` outside it.  Therefore
any Hodge threshold `E_H>=h_H` produces the certified block cost `h_H/2`.
The remaining task is no longer compactness of the single-edge multiplier but
the PDE-to-transfer-weight bridge that supplies the edge measure and its
summable errors.
