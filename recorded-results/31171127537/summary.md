# Helical phase holonomy / diamond obstruction

Status: **CERTIFIED** for the clean finite-packet phase branch.

- exact diamond identity: modal phases cancel from the signed four-edge residual sum
- sharp phase cost: `sum(1-cos delta_i) >= 4(1-cos(|H|/4))`
- clean branch: `|H|>=1/5`, each edge capacity weight `>= beta`, each multiplier `>=1-1e-4`
  implies total polarization deficit `>= beta/250`
- certified coefficient before multiplying by beta: `[0.00499845852429300024714190369837822860206254401 +/- 5.81e-48]`
- random nondegenerate diamonds: `50000`
- worst triad-normal/global coupling reconstruction residual: `1.821e-12`
- worst modal-phase cancellation residual: `3.553e-15`
- worst rigid-rotation holonomy residual: `1.930e-12`
- minimum numerical sharp-cost margin: `4.077e-03`

The observable is relative incidence holonomy, not the Berry phase of one mode or
one triad.  Rigid rotation is exactly free.  A nonzero diamond holonomy is a
phase-lock obstruction and therefore feeds the existing positive polarization
deficit.
