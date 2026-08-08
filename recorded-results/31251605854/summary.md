# Divergence-free coherent Parseval frame

Status: **EXACT_DIVERGENCE_FREE_COHERENT_PARSEVAL_FRAME**.

Let `P` be the Leray projector and `g_z` any normalized affine coherent Gaussian family. For a divergence-free velocity `Pu=u`, self-adjointness gives exactly

`<u,P(g_z e_a)>=<u,g_z e_a>`.

Componentwise Moyal therefore survives projection:

`sum_(a=1)^3 int |<u,P(g_z e_a)>|^2 dmu(z)=||u||_2^2`.

More generally polarized Moyal gives the exact work pairing on the divergence-free subspace. Each probe `P(g_z e_a)` is divergence-free, hence

`<P(g_z e_a),grad p>=0`

with no compact spatial cutoff and no pressure boundary term. The Leray-projected probe need not remain exactly Gaussian; on the signed-good narrow Fourier cell its multiplier distortion is precisely the already-certified smooth-symbol/freezing representation error.

Thus the canonical coherent ancestry/master analysis may use a **global divergence-free coherent frame** rather than a compact moving spatial window. Compact windows remain useful for optional local/CKN diagnostics, but their moving-boundary commutator and localized pressure work need not enter the canonical master `Xi` ledger.

Stress: `20000` random orthogonal-projector/Parseval checks
- worst projected coefficient residual: `1.208e-14`
- worst projected Parseval residual: `3.201e-15`
- worst gradient-pressure pairing residual: `6.804e-15`
- worst projector idempotence residual: `2.756e-15`
