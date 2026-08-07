# Sharp-cutoff log-scale flux bridge

For one ordered conservative triad `k<=p<=q`,

`Integral Pi_K dK/K = -dE_k log(p/k) + dE_q log(q/p)`.

At equal parent scales `k=p=r*`, the lower segment vanishes exactly and the
upper segment is the logarithmic progress factor used by the single-edge
functional.

- r*: `0.610904101586766`
- gamma*: `0.492815285342135`
- symmetric lower segment: `-0.000e+00` times the common triad factor
- symmetric upper segment: `0.602125758280` times the common triad factor
- deterministic local grid points: `10000`
- worst adverse lower/upper ratio on grid: `0.054942561`
- minimum full/upper retention on grid: `0.945057439`

The 90% retention theorem is certified separately by Arb in the single-edge
certificate. The grid values above are regression evidence only.
