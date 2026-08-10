# High-tail HH ultraviolet locality is a dissipation/mass tradeoff

## 1. Read locality before coherent refinement

After the common-unit high-tail theorem, suppose HH is a physical regeneration owner.  The relevant scale law is not the positive Hahn atom mass.  It is the actual hard-output-shell work

\[
H_j=N\int [r_{HH,j}(t)]_+\,dt,
\qquad
H=\sum_{j\ge1}H_j\ge \frac{\nu D_{tail}}2.
\]

Normalize

\[
p_j=H_j/H,
\qquad
H_\infty^{out}=-\log p_{max}.
\]

This law exists before any coherent-cell refinement and is therefore the correct object on which to ask whether HH transfer is local in Fourier scale.

Choose the maximal shell `M=2^jN` and let

\[
\mu_{peak}=\max_t M\|P_Mu(t)\|_2^2.
\]

On the smooth observed slab this is an actual hard-shell PDE observable: `t -> P_Mu(t)` is `L^2`-continuous, so the maximum on the compact slab is attained.  No coherent cell maximum is used.

## 2. Exact high-pass Fourier geometry

For `supp fhat subset {|xi|>K}`, Hölder gives

\[
\boxed{
\|\widehat f\|_{3/2}
\le
\left(\frac{4\pi}{3}\right)^{1/6}
K^{-1/2}\|\nabla f\|_2.
}
\]

Indeed

\[
\int_{|\xi|>K}|\xi|^{-6}d\xi=\frac{4\pi}{3K^3}.
\]

Now fix a locality radius `R>1`.  If a parent producing a child `|xi|<=M` lies above `RM`, the triangle relation `xi=eta+zeta` forces the other parent above `(R-1)M`.

The strict resolved transporter `S_(M/4)u` is supported inside `B_(M/4)`.  Therefore on the UV parent frequencies used here the unresolved field `h=u-S_(M/4)u` equals `u` **exactly**; no cutoff-contraction loss and no new `Xi` term appears.

## 3. The constants collapse to `3 sqrt(pi)`

The physical pair-work bound at child derivative scale `M` uses

\[
C_Y=4A_3=\frac{3\sqrt3}{2}.
\]

The child hard shell and each high-pass parent contribute

\[
C_{hp}=\left(\frac{4\pi}{3}\right)^{1/6}.
\]

Thus

\[
\boxed{
C_Y C_{hp}^3
=\frac{3\sqrt3}{2}\sqrt{\frac{4\pi}{3}}
=3\sqrt\pi.
}
\]

Let `W_UV` denote the absolute work of the selected child shell coming from triads with at least one parent above `RM`.  If the lower parent cutoff `(R-1)M` lies above the parent block scale `N`, then its gradient energy is contained in `D_tail`.  Hence

\[
\boxed{
N W_{UV}^{abs}
\le
\frac{3\sqrt\pi}{\sqrt{R(R-1)}}
\sqrt{\mu_{peak}}\,D_{tail}.
}
\]

No parent-pair count appears.

## 4. The theorem is continuous, not a case split

Write the selected signed HH source as

\[
r_{HH}=r_{comp,R}+r_{UV,R},
\]

where the comparable part has both parent frequencies at most `RM`.  Since

\[
[r_{HH}]_+\le [r_{comp,R}]_+ + |r_{UV,R}|,
\]

its positive comparable work obeys

\[
W_{comp,R}^+
\ge
H_* - W_{UV}^{abs}.
\]

But `H_*=p_*H>=p_* nu D_tail/2`.  Dividing by `p_*D_tail` gives the native theorem

\[
\boxed{
\frac{W_{comp,R}^+e^{H_\infty^{out}}}{D_{tail}}
+
\frac{3\sqrt\pi}{\sqrt{R(R-1)}}
\sqrt{\mu_{peak}}e^{H_\infty^{out}}
\ge\frac\nu2.
}
\]

This is a one-parameter physical locality-radius relation.  `R` is not a new stop class: larger `R` asks for a weaker notion of locality and correspondingly reduces the UV dissipation coefficient.

## 5. Dyadic `R=2` is only a readable corollary

For one dyadic octave,

\[
C_2=3\sqrt{\pi/2}.
\]

At least one nonnegative term in the continuous relation is at least `nu/4`.  Therefore

\[
\boxed{
\mu_{peak}e^{2H_\infty^{out}}
\ge\frac{\nu^2}{72\pi}
}
\]

or

\[
\boxed{
W_{comp,2}^+e^{H_\infty^{out}}
\ge\frac{\nu D_{tail}}4.
}
\]

Equality may register both.  This balanced split is compiler convenience, not the theorem's ontology.

The first branch is already an actual hard critical-shell event and enters the generic shell first-stop/service theorem.

The second branch has actual positive HH work whose parent frequencies are at most `2M`.  **Only after this Fourier locality statement has been proved** do we apply exact coherent localization to that restricted source.  Its positive binary Hahn mass dominates the comparable aggregate positive work exactly as in the upstream binary-work theorem.

## 6. Scope and downstream status

This theorem itself does not claim:

- temporal concentration on one `M`-natural interval;
- signed-good forward scale progress;
- Young near-extremality;
- the generated-energy productivity gate;
- a finite reset count.

The temporal seam identified when this note was written is now closed downstream by `high_tail_natural_window_reentry.md`.  That theorem uses the same positive comparable-work measure and a **sliding** window of exact length `cM^-2`; it proves the scale--time concentration relation and produces an actual critical hard-shell event without packet persistence or fixed time bins.  The present locality theorem remains the Fourier-scale input to that result.

No 3D Navier--Stokes global-regularity conclusion is asserted.
