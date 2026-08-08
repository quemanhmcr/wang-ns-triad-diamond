# Flat companion gate: eliminating the equal-marginal assumption

The spherical master episode is often written with equal parent marginals.  A physical transfer block need not satisfy that symmetry exactly.  On a sufficiently flat block the asymmetry itself produces a clean entropy-or-companion alternative.

Let `b_1,b_2` be the transfer-weighted barycenters of the two parent direction marginals and `b_c` the actual child barycenter.  On the signed-good core, the physical triad calculation from `physical_flat_episode` gives, **without using `b_1=b_2`**,

\[
\boxed{
\left|b_c-\frac{b_1+b_2}{2c_*}\right|
\le e_g,
\qquad
e_g\le2\sqrt{E_H^{phys}}+\frac12E_H^{phys}.
}
\]

Take the service-or-flat value

\[
\tau=1/100,
\qquad
\sqrt{E_H^{phys}}\le\tau/3,
\]

and set

\[
\beta=99/100.
\]

## Low barycenter: fixed entropy

If either parent marginal has barycenter norm at most `beta`, the atomic barycenter--collision inequality gives

\[
H_2\ge\log\frac{2}{1+\beta}
=\log\frac{200}{199}.
\]

Arb certifies

\[
\boxed{H_2>1/200.}
\]

## Both barycenters concentrated: separated companion cores

Assume instead

\[
|b_1|,|b_2|>99/100.
\]

Let `phi` be the angle between their directions.  Since `|b_c|<=1`, the flat geometric error gives

\[
|b_1+b_2|\le2c_*(1+e_g).
\]

But

\[
|b_1+b_2|^2
=(|b_1|-|b_2|)^2+4|b_1||b_2|\cos^2(\phi/2)
\ge4\beta^2\cos^2(\phi/2).
\]

Therefore

\[
\cos(\phi/2)
\le\frac{c_*(1+e_g)}{\beta}.
\]

At `tau=1/100`, Arb gives

\[
\boxed{\phi>1\ \text{radian}.}
\]

For a probability measure on unit directions with barycenter `b`,

\[
\mathbb E|X-\widehat b|^2=2(1-|b|).
\]

Thus `|b|>0.99` implies this second moment is below `1/50`.  Markov at chord radius `3/10` gives

\[
\boxed{\mu\{|X-\widehat b|<3/10\}\ge7/9.}
\]

Because the two barycenter directions are more than one radian apart, these two `3/10`-chord caps have a clean chord gap

\[
\boxed{>1/3.}
\]

Hence a physical `1%` Kelvin-flat block has the unconditional finite-atomic alternative

\[
\boxed{
H_2\ge1/200
\quad\lor\quad
\text{two separated parent cores, each carrying at least }7/9.
}
\]

Relative to a distinguished old lineage, the second concentrated core is a trackable companion.  It must be classified as fresh or reused by the existing nested-grain/ancestry machinery; equal parent marginals are no longer required as an input to reach that gate.
