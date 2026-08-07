# Base-spinor cross forcing is the existing trilinear cross-error

For normalized actual parents `f,g` and Gaussian roles `F,G,H` with distance at most `eps`, the degree-zero/base-child source error obeys

`|T(f,g,H)-T(F,G,H)| <= (2 eps+3 eps^2+eps^3) ||T||`.

The child representation mismatch obeys

`|T(f,g,h)-T(f,g,H)| <= eps ||T||`.

Their sum is exactly the one-shot replacement polynomial
`3 eps+3 eps^2+eps^3` already used in the profile ledger.
At `eps=1%` this splits as `0.020301 + 0.010000 = 0.030301`.

Thus the **work-level degree-zero forcing** produced by parent remainders/cross components is already an omitted trilinear cross interaction and belongs to the existing `eta_j` / `Xi` ledger.  This does not claim an `L^2` bound for the entire nonlinear residual; orthogonal Hermite sidebands are handled separately.

Stress checks: `50000` rank-one norm-one complex trilinear forms
- worst base-source/bound ratio: `0.277831854`
- worst child-mismatch/bound ratio: `0.402957848`
- worst replacement-split identity residual: `5.551e-17`
