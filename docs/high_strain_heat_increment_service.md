# High strain becomes coherent heat-increment service

## Status

**EXACT_HIGH_STRAIN_TO_HEAT_INCREMENT_COHERENT_SERVICE__CRITICAL_RESOLVED_ANCESTOR_FRACTION_RETAINED__OLD_POOL_ROUTING_REMAINS**

The resolved-shell ancestor theorem identifies the **frequency** content of a high-strain event but deliberately does not manufacture a spatial packet.  The missing spatial structure is already present in the gradient itself.  The physically natural way to expose it is through the heat semigroup at the block's own parabolic scale.

## 1. The displacement law is the NS heat kernel

Let

`theta_N=1/(2N^2)`

and let `H_N(r)` be the three-dimensional heat kernel at time `theta_N`:

`H_N(r)=(N/sqrt(2pi))^3 exp(-N^2|r|^2/2)`.

It is a probability density and

`int H_N(r)e^(-i xi.r)dr=exp(-|xi|^2/(2N^2))`.

Thus averaging over `r` is not the introduction of a competing stochastic model.  It is simply the positive physical-space representation of the heat-semigroup defect at the same parabolic scale `N^-2` already present in the Navier--Stokes equation.

For `delta_rV(x)=V(x-r)-V(x)`, Plancherel gives exactly

`int H_N(r)||delta_rV||_2^2dr`

`=int 2(1-exp(-|xi|^2/(2N^2)))|Vhat(xi)|^2dxi`.

Equivalently this is

`2 <V,(I-exp(Delta/(2N^2)))V>`.

## 2. On the strict transporter it is the gradient, up to three percent

The transporter is strict low-pass:

`supp Vhat subset B_(N/4)`.

Put `x=|xi|^2/(2N^2)`.  Then `0<=x<=1/32`.  Since

`1-e^(-x)=int_0^x e^(-s)ds`,

one has

`x e^(-1/32)<=1-e^(-x)<=x`.

Therefore

`e^(-1/32)|xi|^2`

`<=2N^2(1-e^(-|xi|^2/(2N^2)))`

`<=|xi|^2`.

After integration,

`e^(-1/32)||grad V||_2^2`

`<=N^2 int H_N(r)||delta_rV||_2^2dr`

`<=||grad V||_2^2`.

So the heat-increment observable differs from the resolved gradient energy by less than the fixed support factor `e^(-1/32)`, uniformly in scale and aspect.

## 3. High strain forces a positive spatial service law

Define

`S_heat=N^3 int_I dt int H_N(r)||delta_rV(t)||_2^2dr`.

Since

`D_V=N int_I ||grad V||_2^2dt`,

we obtain

`e^(-1/32)D_V<=S_heat<=D_V`.

The high-strain collision gives

`D_V>=D_*=32pi^2/(75c)`

on a lifetime `T=cN^-2`.  Hence high strain forces

`S_heat>=e^(-1/32)D_*`.

This is a fixed, positive, scale-invariant **increment service**, not an additive global resource.  Its physical viscous cost remains tied to the same `D_V/N` scaling.

## 4. Frequency ancestor and spatial service occur on the same law

Let `G` be the shell-time set from `high_strain_resolved_ancestor` on which

`mu_j(t)=M_j||P_ju(t)||_2^2>=mu_*=32pi^2/(75c^2)`.

That theorem gives

`D_V(G)>=D_V/2`.

The heat/gradient comparison is pointwise in frequency, not merely after total integration.  Hence

`S_heat(G)>=e^(-1/32)D_V(G)`

while

`S_heat(total)<=D_V(total)`.

Consequently

`S_heat(G)/S_heat(total)>=e^(-1/32)/2>0.48`.

Thus almost one half of the whole heat-increment service law simultaneously carries a critical lower-frequency ancestor.  There is no independent frequency pigeonhole followed by an unrelated spatial pigeonhole.

## 5. Moyal turns this positive law into coherent spatial edges exactly

Because `delta_r` is a Fourier multiplier, it preserves every resolved dyadic shell.  For each shell `j`, heat displacement `r`, time `t`, and normalized affine coherent window `g_j`, define the positive cell atom

`s_(j,C)(t,r)=N^3 H_N(r) int_C |V_(g_j)(delta_r P_jV)(X,k)|^2 dmu`.

Moyal gives exactly

`sum_C s_(j,C)(t,r)=N^3H_N(r)||delta_rP_jV||_2^2`.

No frame loss and no packet coefficient occur.

Translation covariance is also exact:

`V_g(delta_rf)(X,k)=e^(-ik.r)V_gf(X-r,k)-V_gf(X,k)`.

Each positive atom is therefore a physical edge between two coherent neighborhoods separated by the actual Brownian/heat displacement `r`.  This is the spatial/material geometry that a bare global shell mass did not contain.

## 6. What remains

The measure is now both frequency-ancestral and spatially coherent, but its **material ownership** still has to be routed.  The next physical dichotomy should use the existing old-reservoir pool, not invent a new selector:

- old--old heat-increment edges must be bounded by the material reservoir service capacity;
- old--new edges are genuine material/interface relink;
- new--new edges create fresh coherent ancestry or collision entropy;
- if a critical ancestor cannot be represented by the transported old pool, that failure itself must be a material relink/source event.

This is close in spirit to `coherent_increment_service`, but it must be proved for the dissipation-seeded **heat increment measure** rather than borrowed from the SGS cubic source theorem.

`D_V` remains a scale-critical `O(1/N)` physical cost and is not promoted to a finite reset.  No global-regularity claim is made.
