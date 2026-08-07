# Shell concentration versus affine Gaussian aspect

Status: **CERTIFIED**.

- exact uncertainty matrix: `Sigma_x Gamma_(3/2)=I/3`
- one-percent profile distance and the certified log shell imply every physical
  standard axis is `> (2/3) N^-1`
- on the radius-two covariance ellipsoid of the Gaussian profile, Hausdorff--Young
  plus local Holder gives actual physical mass
  `N integral_E |check f|^2 >= (3/10) N(det Sigma_x)^(1/6)`
- if `A=N l_max`, then
  `N integral_E |check f|^2 > (1/5) A^(1/3)`
- affine Young symmetry remains exact: elongation alone is not a transfer cost
- random covariance checks: `50000`
- worst uncertainty residual: `7.541e-07`
- local ellipsoid coefficient (numerical): `0.311109826`
- tested free affine Young aspect: `1.0e+08`

The last bullet is a required correction to a naive replication argument.  A
common affine Gaussian can be arbitrarily anisotropic while remaining an exact
Young extremizer.  Therefore aspect ratio alone must not be inserted as a
Bellman deficit.  Dynamics must instead see the anisotropy through the
ellipsoid-localized mass ledger and the grain-normalized curvature tensor.
