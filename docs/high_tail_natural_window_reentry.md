# High-tail comparable HH work: sliding natural-time concentration


> **Current placement.** The pure-UV branch of `hard_tail_true_upward_supply.md` arrives already at the first dyadic shell `M=2N` with comparable interaction parents.  It is therefore geometrically ready for this natural-window mechanism once the exact donor-restricted canonical positive submeasure is bound into the window law.  That binding must preserve the common `N dW` causal unit; the resolved-contact branch remains separate.

## 1. Time must be read from the same physical work measure

The ultraviolet theorem already selects one hard output shell

\[
M=2^jN,\qquad j\ge1,
\]

and, on its comparable-parent owner, supplies actual positive comparable HH work.  The remaining question is temporal: does that work produce a hard-shell event on the shell's own natural time?

No packet persistence is introduced.  Let `mu_comp` be the positive comparable-work measure on the parent source block.  The selected shell's natural duration is

\[
T_M=cM^{-2}.
\]

Define the **sliding** concentration

\[
\boxed{
p_t=
\frac{\sup_s\mu_{comp}([s,s+T_M])}
{\mu_{comp}(I)}.
}
\]

There is no fixed time partition and no distinguished bin origin.  For smooth Navier--Stokes, the comparable-work density is continuous, so the sliding integral is continuous in `s` and attains its maximum on the compact admissible start interval.
The concentration datum is inseparable from this physical window length: a value measured on any other duration is not admissible for the selected shell.

The logarithmic coordinate

\[
H_\infty^{time}=-\log p_t
\]

is only a way to record concentration of this same physical measure.

## 2. Time-origin, time-unit, and representation gauges

The sliding fraction is intrinsic:

- translating all times by a common constant changes neither total nor maximal-window work;
- under `t'=lambda t`, the density transforms as `rho'=rho/lambda` and the window as `T'=lambda T`, leaving the measure and `p_t` unchanged;
- subdividing a representation of the same positive measure leaves `p_t` unchanged.

Thus the theorem cannot be altered by choosing a different clock origin, unit, or bookkeeping mesh.

## 3. Natural-window work capacity

The canonical strict resolved cutoff currently certifies only

\[
|S|\le1.
\]

Therefore the unresolved field obeys only the universal estimate

\[
\|h\|_2=\|(I-S)u\|_2\le2\|u\|_2.
\]

We do **not** assume the sharper nonnegative-cutoff contraction before it has its own certificate.

Suppose the localized Fourier theorem has restricted both HH parents to frequencies at most `RM`.  If

\[
\mu_{win}=\max_{t\in I_*}M\|P_Mu(t)\|_2^2,
\qquad
E_{global}=\sup_t\|u(t)\|_2^2,
\]

then the child `L^(3/2)` amplitude is bounded by the hard-shell Fourier-volume constant times `sqrt(mu_win)`, while the two comparable parent amplitudes are bounded by the radius-`RM` Fourier-volume constant and `2sqrt(E_global)` each.  With the certified physical work Young coefficient,

\[
r_{comp}(t)
\le
12\sqrt\pi\,R\,M^2E_{global}\sqrt{\mu_{win}}.
\]

Integrating over the exact natural duration `cM^-2` cancels `M^2`:

\[
W_{win}^{phys}
\le
12c\sqrt\pi\,R\,E_{global}\sqrt{\mu_{win}}.
\]

The high-tail causal law is stored in the common unit `N dW`, hence

\[
\boxed{
N W_{win}
\le
12c\sqrt\pi\,R\,N E_{global}\sqrt{\mu_{win}}.
}
\]

The denominator `N E_global` is the scale-critical global energy at the parent block scale.

## 4. Scale concentration and time concentration compose without a packet

Let

\[
p_s=e^{-H_\infty^{out}}
\]

be the selected hard-output-shell fraction from the Fourier locality theorem.  On the clean comparable owner,

\[
\frac{W_{comp}}{p_s}
\ge
\frac{\nu D_{tail}}4.
\]

By the sliding definition,

\[
W_{win}=p_tW_{comp}.
\]

Combining with the natural-window capacity gives

\[
\boxed{
\frac{\sqrt{\mu_{win}}}{p_sp_t}
\ge
\frac{\nu D_{tail}}
{48c\sqrt\pi\,R\,N E_{global}}.
}
\]

Equivalently,

\[
\boxed{
\sqrt{\mu_{win}}
\exp(H_\infty^{out}+H_\infty^{time})
\ge
\frac{\nu D_{tail}}
{48c\sqrt\pi\,R\,\mathcal E_N},
\qquad
\mathcal E_N=NE_{global}.
}
\]

For the dyadic locality corollary `R=2`, the denominator is `96c sqrt(pi) mathcal_E_N`.

This is not an entropy stop.  It is a strength--concentration relation of one positive physical work law projected first to scale and then to sliding natural time.

## 5. The resulting event has genuine hard-tail progress

The peak `mu_win` is an actual hard-shell event at `M=2^jN`.  Therefore it enters the existing generic critical-shell first-stop theorem without a coherent-cell selector or material label.

Moreover,

\[
\boxed{M/N=2^j\ge2}
\]

and

\[
\boxed{T_M/T_N=4^{-j}\le1/4}.
\]

This is genuine forward scale progress and natural-time shortening supplied by hard-tail support itself.  It does not use the signed-good Young ratio of the near-extremal generated lineage.

The generic shell theorem then gives the already-certified alternatives: named strain/interface/HH stop, `t=0`, or own-scale service on a full no-hit natural corridor.  The own-scale service lower is conditional exactly as before.

## 6. Scope

This theorem uses no:

- packet persistence;
- fixed time bins;
- time-origin convention;
- nonnegative cutoff assumption beyond the certified `|S|<=1`;
- signed-good parent ratio;
- Young near-extremality;
- generated-energy productivity gate;
- additive finite reset.

It closes temporal localization of the comparable high-tail HH work into an actual hard-shell recursive event.  No 3D Navier--Stokes global-regularity conclusion is asserted.
