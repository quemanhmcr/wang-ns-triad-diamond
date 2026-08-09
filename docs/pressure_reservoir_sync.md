# Pressure reservoir synchronization: derivative order sets the pair lifetime

Two pressure source objects appear in the continuum architecture and must not be identified merely because both come from pressure.

## 1. Resolved objective pressure Hessian

The objective-strain source uses

\[
\rho_P=N^{-4}\|\nabla^2P\|_\infty.
\]

For a low--low pair `V_a,V_b`, order-two differentiation and the band-limited `L^{3/2}->L^\infty` Bernstein step give, up to one fixed pair-independent Riesz/Bernstein constant,

\[
\rho_{P,ab}
\lesssim
\left(\frac{M_{max}}N\right)^4
\sqrt{\mu_a\mu_b}
=
\frac{M_{max}^4\sqrt{M_aM_b}}{N^4}\sqrt{E_aE_b},
\]

where `M_max=max(M_a,M_b)` and `mu_i=M_iE_i`.

On a supplied signed-good low-strain lineage each material reservoir frequency grows by less than `21/20`, while the block scale grows by more than `8/5`.  Therefore the coefficient **per unit `sqrt(E_aE_b)`** contracts by

\[
\boxed{
\left(\frac{21}{20}\right)^5
\left(\frac58\right)^4
=
\frac{4084101}{20971520}
<\frac15.
}
\]

Even if both reservoirs are adversarially allowed the entire global energy cap at every future service generation, one fixed objective-Hessian pair has

\[
\boxed{
\sum_{q\ge0}\mathcal P^{(2)}_{ab,q}
<\frac54\mathcal P^{(2),cap}_{ab,0}.
}
\]

The fifth power of the reservoir-growth factor is physical: four powers come from the Hessian plus `L^{3/2}->L^\infty` band scaling, and one comes from `sqrt(M_aM_b)` when critical mass is converted to physical energy.

## 2. H1 pressure-third source

The older H1 covariant source contains a different pressure object.  Its `L^{3/2}` pressure-third coefficient is

\[
\mathcal P^{(3)}_{ab}(N)
\sim
\frac{M_{max}^3\sqrt{M_aM_b}}{N^3}\sqrt{E_aE_b}.
\]

For the same material lineage this contracts by

\[
\boxed{
\left(\frac{21}{20}\right)^4
\left(\frac58\right)^3
=
\frac{194481}{655360}<\frac13,
}
\]

and its total future fixed-pair capacity remains `<3/2` of generation zero.

## 3. Master meaning

The mechanism is the same but the lifetimes are not: **derivative order is provenance**.

A persistent objective pressure-Hessian branch must therefore relink to new low-frequency material pairs, fragment its capacity across many pairs and pay pair/component entropy or cycle structure, leave the supplied low-strain signed-good lineage, or use the SGS-stress alternative already routed to coherent increments.  The H1 pressure-third branch obeys the analogous statement with its own one-third lifetime.

No low-pass mass is promoted to a generic critical shell merely because it participates in pressure.  The pair erosion theorem prices reuse; it does not invent a scale-independent reset.
