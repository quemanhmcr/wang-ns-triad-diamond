# Good-core helical frame Lipschitz theorem

Status: **CERTIFIED**.

- certified good-core parent angle: `1/4 < cos(theta) < 2/5`
- hence `sin(theta) > 9/10`
- exact normal derivative bound: `||dn|| <= (10/9)(||da||+||db||)`
- parent triad-normal helical frame: `||dh|| <= (5/2)(||da||+||db||)`
- child frame: the same `5/2` constant on the sum of the three carrier-direction rates
- random derivative checks: `50000`
- worst normal-bound ratio: `0.962540964`
- worst helical-bound ratio: `0.427795976`
- worst finite-difference derivative residual: `1.757e-07`
- minimum parent sine seen: `0.930294406`

The Chern obstruction is global, but the signed-good extremal core stays uniformly
away from collinearity.  Therefore the triad-normal gauge has no local chart
singularity on a good packet block: helical-frame variation is linearly
subordinate to carrier-direction variation with a scale-free constant.
