# Resistance-to-Bellman Poisson stopping

Random exact checks: `3000`; violations: `0`.
Minimum edge-entropy margin: `7.653e-04`.
Minimum exhaustive-witness margin: `2.825e-04`.

## Median constants

F0: `0.316060279414`
Bellman entropy floor h0: `0.172011060757`
simultaneous-cut multiplier K0: `10.655813654955`

## Toy trees

| type | median R | F | rho | cut bound | component-Q bound | edge H | edge-H lower bound |
|:--|--:|--:|--:|--:|--:|--:|--:|
| star | 34 | 0.594937 | 0.500000 | 1.858513 | 0.702532 | 2.833213 | 0.635989 |
| path | 85 | 0.592199 | 0.200000 | 0.861842 | 0.703901 | 2.833213 | 0.208545 |
| bottleneck | 128.04 | 0.642868 | 0.125273 | 0.497276 | 0.678566 | 2.773212 | 0.124955 |

The Poisson collision identity, simultaneous stopping bound, and edge-conductance collision bound are exact finite-dimensional statements. Random/exhaustive checks only validate implementation.
