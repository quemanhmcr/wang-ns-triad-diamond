# Odd-Hermite triad selection and daughter-transfer cost

Status: **EXACT_GIVEN_DEGREE3_GAUSSIAN_HYPERCONTRACTIVITY**.

For a centered resonant affine Gaussian triad, simultaneous inversion of all centered frequency deviations leaves the Gaussian trilinear weight invariant.  Therefore a polynomial/Hermite perturbation has the exact selection rule

`T(P_n1 G1,P_n2 G2,P_n3 G3)=0` whenever `n1+n2+n3` is odd.

In particular a single H1 or H3 daughter cannot feed the base Gaussian triad at first order.  Any sideband rescue contains at least two odd sidebands and is bounded by

`A3 (rho1 rho2 b3 + rho1 rho3 b2 + rho2 rho3 b1 + rho1 rho2 rho3)`.

There is also a quantitative one-role loss.  Let `R=P G` be an odd degree<=3 sideband and let `sigma^2=E_mu |P|^2` in the critical `|G|^(3/2)` Gaussian measure.  Gaussian hypercontractivity gives `E|P|^4<=729 sigma^4`.  If `sigma<=1/80`, parity plus uniform convexity yields

`||G+R||_(3/2)^(3/2) >= ||G||_(3/2)^(3/2) (1+3 sigma^2/16)`

and since the single-sideband numerator is exactly zero,

`Def_transfer >= sigma^2/16`.

Thus a coherent odd daughter has only two ways to remain efficient: recruit a second odd sideband (a genuine daughter/cross interaction component) or pay a quadratic transfer deficit.

Stress: `50000`
- worst H1 single-sideband parity residual: `1.019e-15`
- worst H3 single-sideband parity residual: `4.527e-14`
- minimum clean-deficit margin: `-3.088e-17`
- minimum quadratic-rescue margin: `9.088e-06`
