# Physical pair-weighted productivity: KL replaces a false Duhamel pair identification

The amplitude--entropy theorem needs a lower bound for the geometric mean of the
two parent amplitudes **under the physical causal law**.  It is unnecessary, and
in general false, to obtain this by first choosing a parent pair with an adjoint
Duhamel weight and then identifying that pair law with physical child-energy
work.

There is a direct physical proof.

## 1. Work density already contains the parent product

On a retained hard parent-pair cell `e`, let

\[
a_c(t)=\|\widehat w_c(t)\|_{3/2},
\qquad
a_{1,e}(t),a_{2,e}(t)
\]

be the critical Fourier amplitudes.  The actual positive child-energy work has
density

\[
d\mathcal T_e=r_e(t)dt.
\]

For the symmetrized high--high source, sharp Young and `|k|<=R_kN` give

\[
\boxed{
r_e(t)
\le C_YN\,a_c(t)a_{1,e}(t)a_{2,e}(t),
}
\]

where one may take

\[
C_Y=4R_kA_3,
\qquad
A_3=(\sqrt3/2)^3.
\]

This inequality uses **physical work itself**.  No adjoint probability appears.

## 2. Normalize physical time and pair cells

Let the child slab have

\[
T=cN^{-2}
\]

and suppose the retained frequency/helicity partition has `M` hard parent-pair
cells.  Put

\[
W=\sum_e\int r_e(t)dt
\]

and normalize

\[
dw=\frac{d\mathcal T}{W}.
\]

Use as reference probability measure

\[
d\mu_0=\frac{dt}{T}\frac1M.
\]

If `p=dw/dmu_0`, then

\[
r_e(t)=\frac{W}{TM}p(t,e).
\]

Therefore the pointwise Young bound gives

\[
\log(a_1a_2)
\ge
\log\frac{W}{TMC_YN}
+\log p
-\log a_c.
\]

Average with the **actual physical law** `w`.  Relative entropy is nonnegative,

\[
D(w\|\mu_0)=\int\log p\,dw\ge0,
\]

so

\[
\boxed{
\mathbb E_w\log(a_1a_2)
\ge
\log\frac{W}{TMC_YN}
-
\mathbb E_w\log a_c.
}
\]

Concentration of the physical work in time or among pair cells only makes the KL
term positive.  There is no discretization penalty beyond the explicit finite
pair-cell count `M`.

## 3. Energy controls the child amplitude from above

On the generated branch of the exact physical-energy gate,

\[
K\le1/30,
\qquad
E_0\le E_1/5,
\qquad
W_R^+\le E_1/5,
\qquad
W\ge8E_1/15.
\]

At every intermediate time, the same energy Gronwall inequality gives

\[
E(t)
\le e^{2K}(E_0+W+W_R^+).
\]

Since

\[
E_0+W_R^+\le\frac25E_1\le\frac34W,
\]

one obtains

\[
\boxed{
\sup_tE(t)
\le C_EW,
\qquad
C_E=\frac74e^{1/15}<2.
}
\]

The selected child role stays in a fixed Fourier ball `|xi|<=R N`.  Holder on
that finite Fourier volume gives

\[
\boxed{
a_c(t)
\le C_\Omega\sqrt{NE(t)},
\qquad
C_\Omega=|B_R|^{1/6}.
}
\]

Consequently

\[
\mathbb E_w\log a_c
\le
\log(C_\Omega\sqrt{NC_EW}).
\]

Insert this into the KL inequality and use `T=cN^-2`:

\[
\mathbb E_w\log(a_1a_2)
\ge
\frac12\log(NW)
-\log(cMC_YC_\Omega\sqrt{C_E}).
\]

Finally `W>=8E_1/15`, while every L2-normalized terminal analysis coefficient
satisfies

\[
\alpha_c\le\sqrt{NE_1}.
\]

Hence

\[
\boxed{
\mathbb E_w\log(a_1a_2)
\ge
\log\alpha_c+\log\Lambda_{role,M},
}
\]

with

\[
\Lambda_{role,M}
=
\frac{\sqrt{8/15}}
{cMC_YC_\Omega\sqrt{C_E}}.
\]

All powers of `N` have disappeared.

## 4. Register the physical parents, still with physical weights

For every retained signed-good parent role, complex Young and the dual-Gaussian
theorem give at its actual event time

\[
\alpha_{p,event}\ge\sqrt{\eta_{dual}}\,a_p.
\]

If no earlier causal stop occurs, common-slice registration gives the conservative
factor `1/4`.  Therefore

\[
\boxed{
\mathbb E_w
\log(\alpha_{p_1}\alpha_{p_2})
\ge
\log\alpha_c+\log\Lambda_M,
}
\]

where

\[
\boxed{
\Lambda_M
=
\frac{\eta_{dual}}{16}
\frac{\sqrt{8/15}}
{cMC_YC_\Omega\sqrt{C_E}}.
}
\]

The expectation is still taken with respect to **actual positive child-energy
work**.  Duhamel has not selected, reweighted, or changed a parent pair.

## 5. Conditioning on the first-stop survivor set

Let `C` be the set of physical pair events which survive all earlier causal stops and the complex-Young/common-slice registration gates, and write

\[
q={\mathcal T(C)\over W}.
\]

Apply the same KL proof to the normalized restricted law `dT|C/(qW)`.  The child shell-energy upper still uses the full generated work `W`, while the retained density numerator is `qW`.  Therefore the only change is

\[
\boxed{\Lambda_{C,M}=q\,\Lambda_M.}
\]

No new entropy inequality is needed.  In particular, if `q>=1/2`, survival conditioning loses at most `log 2` in the logarithmic productivity offset.  If `q<1/2`, more than half the physical work has already hit an earlier named causal stop and must be routed there rather than treated as a free continuing packet.

## 6. Fine symbol cells do not damage the asymptotic reuse slope

The smooth-symbol theorem may refine pair cells with

\[
h_j=(j+3)^{-2}.
\]

The compact pair-frequency domain is six dimensional, so a product grid has at
most

\[
M_j\le M_0(j+3)^{12}
\]

cells after absorbing the finite helical choices into `M_0`.

The binary amplitude recursion with variable productivity is

\[
\ell_j
\ge
\frac12\log\Lambda_j+rac12\ell_{j+1}.
\]

Thus a depth-`j` cell factor is weighted only by `2^(-j-1)`.  Since

\[
j+3\le3\,2^j,
\]

\[
\boxed{
\sum_{j\ge0}2^{-j-1}\log M_j
\le
\log M_0+12\log6<\infty.
}
\]

More generally, every polynomial cell refinement contributes only a finite
amplitude offset.  It cannot alter the linear `log(48/25)` Shannon/Renyi reuse
slope.

This geometric discount is the key reason one does not need a uniform
pointwise parent-product constant on every arbitrarily fine frequency cell.

## 7. What this removes

The previous architecture used Duhamel to prove a product bound and physical work
to define causal weights.  That left a hidden question: is the Duhamel-important
parent pair also important under `dT`?

The present theorem removes the question entirely.  Parent productivity is
proved **after averaging under `dT` itself**.  Duhamel can remain an exact causal
support/adjoint identity, but no `dGamma -> dT` pair change of measure is needed.

The remaining continuum assembly is local rather than measure-theoretic:

1. apply the sharp-Young work-density inequality on the retained signed-good hard
   event cells supplied by the event-role theorem;
2. when a cell fails complex Young, phase/polarization, or common-slice
   registration, stop at the already existing transfer/backscatter/relink/source
   cause;
3. otherwise use the physical weighted logarithmic productivity above in the
   amplitude--entropy telescope.

No Navier--Stokes global-regularity conclusion is made here.
