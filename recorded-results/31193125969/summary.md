# H3 sideband local no-escape theorem

Status: **EXACT_GIVEN_COHERENCE_AND_ODD_SIDEBAND_THEOREMS**.

Let `I3=int ||Sym B_tilde|| dt` in the affine-curvature interaction frame and `J3=int ||Sym S_tilde|| dt`, with `T` the packet lifetime.

The coherence theorem gives

- `J3 >= I3/T`, or
- first-Duhamel H3 daughter `delta1^2 >= 3 I3^2/32`.

On the coherent branch, either nonlinear sideband feedback has size at least `delta1/2`, or the surviving daughter has `delta^2 >= 3 I3^2/128`.  Let `sigma` be its second moment in the critical `|G|^(3/2)` Gaussian measure.  The variance change gives `sigma>=delta`.

If `sigma>=1/80`, this is already a definite daughter-capacity event.  If `sigma<1/80`, odd-Hermite convexity gives single-role transfer loss at least `sigma^2/16`.  Splitting that loss against possible pair-sideband rescue yields the clean alternative

`net transfer deficit >= 3 I3^2/4096`

or

`pair-sideband rescue >= 3 I3^2/4096`.

Thus H3 curvature has five and only five exits in this model: acceleration-Hessian dephasing source, nonlinear sideband feedback, a definite large daughter, a quadratic transfer deficit, or a quadratic pair-sideband interaction.

Stress: `50000`
- branch counts: `{'dephasing_source': 9898, 'nonlinear_sideband_feedback': 10072, 'large_daughter_capacity': 9997, 'pair_sideband_rescue': 10003, 'transfer_deficit': 10030}`
- minimum first-impulse margin: `2.436e-08`
- minimum transfer-cost margin: `3.321e-11`
- minimum pair-cost margin: `4.293e-11`
