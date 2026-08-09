# First-hit heat reservoir erosion: high strain touches the boundary of a low-strain history

## Status

**EXACT_FIRST_HIT_HEAT_RESERVOIR_EROSION__CLEAN_OO_RATIO_441_640__FINITE_NON_OO_FORCE_ON_SUPPLIED_SIGNED_GOOD_EPOCH**

The remaining old--old (`OO`) high-strain question appears paradoxical only if a first-hit event is described after the fact as a whole "high-strain slab".  The physical object is instead the first boundary contact of a continuous strain action.

## 1. A first high-strain contact has a low-strain transport history

For one selected smooth-SGS block define

\[
K(t)=\int_{t_0}^t\|S_V(s)\|_{op}\,ds.
\]

If `tau` is the first contact with the high-strain face, then

\[
K(\tau)=1/30,
\qquad
K(t)\le1/30\quad(t_0\le t\le\tau).
\]

Kelvin transport therefore gives, including the terminal boundary time,

\[
\frac{M(\tau)}{M(t_0)}\le e^{1/30}<\frac{21}{20}.
\]

The high-strain collision already holds at equality.  Thus the event does not invalidate the low-strain kinematics that led to it.

## 2. Heat service has the coefficient `M^2 E/N`

Let `f` be supported in the deterministic Fourier band `|xi|<=M`.  For the Navier--Stokes heat probe at time `1/(2N^2)`,

\[
2\left(1-e^{-|\xi|^2/(2N^2)}\right)
\le\frac{|\xi|^2}{N^2}
\le\frac{M^2}{N^2}.
\]

Hence

\[
\mathbb E_{H_N}\|\delta_r f\|_2^2
\le\frac{M^2}{N^2}\|f\|_2^2.
\]

The canonical natural lifetime is `T(N)=cN^-2` with the same fixed scaled coefficient `c` used by the signed-good parabolic geometry.  Therefore

\[
\boxed{
N^3\int_I\mathbb E_{H_N}\|\delta_r f\|_2^2dt
\le c\frac{M^2}{N}\sup_I\|f(t)\|_2^2.
}
\]

This differs from the old SGS increment capacity `M^3E/N^2`.  Heat service has its own dimension and therefore its own erosion ratio.

## 3. The estimate is shellwise, not an old-field decomposition

The high-strain heat theorem first disintegrates the resolved field into deterministic dyadic Fourier shells, then applies Moyal to each shell to obtain a positive edge law.  Material ownership is imposed only afterwards.

Thus, within every shell,

\[
S_{OO,j}\le S_{heat,j}
\]

simply because OO is a positive submeasure of that already-existing shell law.  Summing the orthogonal shell capacities gives a whole-old-pool upper bound.  No object `u_old` is constructed, no identity `u=u_old+u_new` is used, and no quadratic cross term appears.

As in the established reservoir-pool theorem, the upper bound may adversarially allow the old pool the entire frame-energy budget at every service time,

\[
\sum_a E_{a,q}\le P E_{global}.
\]

This is deliberately pessimistic and requires no packet-mass floor.

## 4. The physical contraction is stronger than the clean rational envelope

On a supplied signed-good material epoch, the block frequency advances by

\[
N_{q+1}/N_q>8/5.
\]

A reused material frequency crossing a first-hit block grows by at most `e^(1/30)`.  Therefore its heat-capacity coefficient contracts by

\[
\boxed{
\rho_{phys}
\le\frac58 e^{1/15}
<\frac{441}{640}
<\frac{7}{10}.
}
\]

The number `441/640` is not asserted to be the exact physical ratio.  It is the clean rational envelope obtained from `e^(1/30)<21/20`:

\[
\frac{(21/20)^2}{8/5}=\frac{441}{640}.
\]

Consequently, if all old material frequencies initially obey `M_{a,0}<=alpha N_0`, then

\[
C_{OO}(q)
\le
c\alpha^2N_0PE_{global}
\left(\frac{441}{640}\right)^q.
\]

The total future OO heat capacity is finite:

\[
\boxed{
\sum_{q\ge0}C_{OO}(q)
\le\frac{640}{199}C_{OO}(0)
<\frac{13}{4}C_{OO}(0).
}
\]

## 5. Sufficient material age forces interface-or-new heat service

Every first high-strain contact has the existing scale-independent normalized heat lower

\[
S_{heat}\ge S_*(c)>0.
\]

Fix `0<f<1` and let `q_*` be the first material age satisfying

\[
C_{OO}(q_*)\le(1-f)S_*.
\]

At every high-strain contact in the same supplied material epoch with age `q>=q_*`, exact positive ownership gives

\[
\boxed{
S_{ON}+S_{NN}
=S_{heat}-S_{OO}
\ge fS_*.
}
\]

For `f=1/2`, at least half of the heat service is forced out of OO.  This is a finite material-capacity statement, not a count of `D_V` events and not a new reset currency.

## 6. What this actually closes

This theorem closes the **quantitative repeated-OO capacity seam** once the PDE has already supplied a signed-good material recursive epoch with the canonical natural lifetime `T(N)=cN^-2`.

It does not manufacture that epoch.  Universal slab renewal remains open.  Nor does it identify the forced remainder with a final master destination:

- `ON` is a same-time material-interface edge mark, but still needs attachment to an actual temporal interface/relink work law;
- `NN` is genuinely outside the old material pool, but still needs a renewal/ancestry theorem that does not postulate a packet mass floor.

No global-regularity claim is made.
