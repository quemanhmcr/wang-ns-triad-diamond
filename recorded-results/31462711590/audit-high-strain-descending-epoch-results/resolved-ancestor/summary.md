# High strain is a dissipation-weighted resolved-ancestor event

Status: **EXACT_HIGH_STRAIN_DISSIPATION_WEIGHTED_RESOLVED_ANCESTOR__HALF_DV_ON_CRITICAL_LOW_SHELLS__MATERIAL_RENEWAL_REMAINS**.

The normalized dissipation `D_V=N int ||grad V||_2^2 dt` should not be treated as an abstract reset currency.  The strict transporter itself reveals where that dissipation lives.  Decompose its Fourier ball into deterministic dyadic annuli

`A_j={M_j/2<|xi|<=M_j},  M_j=(N/4)2^(-j)`.

They satisfy `sum_j M_j=N/2`.  Put `E_j(t)=||P_j u(t)||_2^2` and the actual critical shell mass `mu_j(t)=M_j E_j(t)`.  Since the standard low-pass multiplier is an L2 contraction and `|xi|<=M_j` on `A_j`,

`N ||grad P_j V||_2^2 <= N M_j mu_j(t)`.

On a natural lifetime `T=cN^-2`, the part of the **actual resolved dissipation measure** lying where `mu_j(t)<mu_*` is therefore at most

`N T mu_* sum_j M_j = c mu_*/2`.

The high-strain collision gives

`K>=1/30 => D_V>=D_* = 32 pi^2/(75 c)`.

Choose `mu_*=D_*/c=32 pi^2/(75 c^2)`.  Then the low-mass part is at most `D_*/2`, so for every `D_V>=D_*`,

`D_V({mu_j>=mu_*}) >= D_V-D_*/2 >= D_V/2`.

Thus at least **half of the actual normalized dissipation law** already carries a simultaneous low-frequency ancestor with fixed critical mass.  No shell is selected by argmax, and there is no logarithmic loss from the infinitely many low shells.  The law is the physical density `N|xi|^2|Vhat|^2 dt dxi` disintegrated by dyadic frequency.

Every such ancestor has `M_j<=N/4`; its parabolic natural lifetime is therefore at least `16` child lifetimes.  This is a much stronger scale separation than the signed-good HH parent ratio, but it is a different physical object: a resolved reservoir/ancestor state, not yet a transfer-generated hard parent.

Stress: `100000` random shell/time dissipation laws
- minimum low-mass dissipation upper margin: `2.509e-01`
- minimum half-law margin: `2.598e-01`
- minimum sampled clean critical-mass threshold: `7.713e-02`
- minimum ancestor/child lifetime ratio: `16.000`
- maximum retained dissipation fraction: `1.000000`

This does **not** promote `D_V` to a globally finite reset: its physical viscous cost is still `nu D_V/N` and remains summable on geometric high-frequency chains.  What changes is its recursive meaning.  A high-strain stop is now accompanied, on at least half of its own physical dissipation law, by a critical resolved-shell ancestor at a genuinely lower scale.  The remaining bridge is to attach these dissipation-seeded shell ancestors to the existing material/coherent reservoir/reuse or renewed-slab machinery without inventing a packet selector.  No global-regularity claim is made.
