# Asynchronous Duhamel slabs synchronize by parabolic causality

The one-step adjoint Kelvin--Duhamel gate gives an exact causal source law but its parent-pair atoms occur at different physical times. A full causal ancestry therefore needs a common time layer. The key point is that **packet persistence is stronger than necessary**. It is enough to synchronize the support of the positive Duhamel generation measure.

## 1. Natural parabolic geometry

On the signed-good core,

\[
\frac35<\frac{N_p}{N_c}<\frac58.
\]

For the natural lifetime \(T(N)=cN^{-2}\),

\[
\boxed{
\frac{64}{25}<\frac{T_p}{T_c}<\frac{25}{9}.
}
\]

The generated branch of the adjoint theorem carries a positive aligned parent-pair measure. Divide the terminal child slab into two halves and retain the heavier half. At least half of the positive generation mass remains. Its parent events occupy a physical interval of width at most \(T_c/2\). Since every such parent has lifetime \(>64T_c/25\), the normalized parent-event span obeys

\[
\boxed{
\alpha_1\le \frac12\frac{25}{64}=\frac{25}{128}.
}
\]

Only this first factor \(1/2\) is lost. Later layers use their complete generated Duhamel measure.

The geometry in this paragraph is in fact **measure agnostic**.  The heavy-half step uses only positivity and event support.  After the physical-energy causal bridge, the preferred master-facing choice is to place the actual positive high--high child-work measure

\[
 d\mathcal T_{HH}=2[\Re\langle c,F_{HH,\alpha}\rangle]_+dt
\]

on these same quadratic event times and choose its heavier half.  The bound `alpha_1<=25/128`, the recurrence below and the common-slice geometry are unchanged.  Raw adjoint `dGamma` remains an amplitude/support diagnostic and need not be identified with `dT_HH`.

## 2. A common reference slice without persistence

Suppose one generated layer has event support

\[
H_j=[a_j,b_j],
\qquad
\alpha_j=\frac{b_j-a_j}{T_j^{\min}}.
\]

Choose

\[
\boxed{
s_j=a_j-\frac25T_j^{\min}.
}
\]

If \(\alpha_j\le3/8\), then for every event \(t\in H_j\),

\[
t-s_j
\le \left(\frac38+\frac25\right)T_j^{\min}
=\frac{31}{40}T_j^{\min}
<T_j^{\min}.
\]

Hence the exact adjoint Duhamel gate may be applied from \(s_j\) to **every node in the layer**, regardless of whether a Gaussian packet would have remained a good static approximation over the whole interval.

Each node now has the already-certified alternative:

\[
\text{material inheritance}
\quad\lor\quad
\text{classified residual/interface}
\quad\lor\quad
\text{high--high generation}.
\]

The first two branches stop. Only the generated branch is expanded to the previous causal layer. Those new parent-pair atoms lie inside the common interval \([s_j,b_j]\), whose width is at most

\[
(\alpha_j+2/5)T_j^{\min}.
\]

Their own natural lifetimes are at least \(64/25\) times longer, so

\[
\boxed{
\alpha_{j+1}
\le
\frac{25}{64}\left(\alpha_j+\frac25\right).
}
\]

At the proposed cone boundary,

\[
\frac{25}{64}\left(\frac38+\frac25\right)
=\frac{155}{512}
<\frac38.
\]

The loose cone `alpha<=3/8` is invariant.  But after the first half-slab one gets a sharper statement.  Since\n\[\n\frac{25}{128}<\frac{10}{39}\n\]\nand `10/39` is the affine-recurrence fixed point, induction gives\n\[\n\boxed{\alpha_j\le\frac{10}{39}\qquad(j\ge1).}\n\]\nHence the common-slice margin is actually\n\[\n\boxed{\n1-\frac25-\frac{10}{39}=\frac{67}{195}.\n}\n\]

This is a **parabolic synchronization cone** for the causal Duhamel measure. It replaces the unproved requirement that every selected parent packet persist coherently through a common natural window.

## 3. The initial surface is reached after finite causal depth — corrected asynchronous geometry

A subtle point matters here.  The next generated support is contained in the
previous common Duhamel interval,

\[
H_{j+1}\subset[s_j,b_j],
\]

but its left endpoint need not equal \(s_j\).  Thus one must compare successive
**reference slices**, not simply add \(2T_j/5\) per generation.

Write

\[
s_j=a_j-\frac25T_j,
\qquad
b_j=a_j+\alpha_jT_j.
\]

The worst case for backward progress is \(a_{j+1}=b_j\).  Since
\(T_{j+1}\ge(64/25)T_j\) and \(\alpha_j\le10/39\),

\[
\begin{aligned}
s_j-s_{j+1}
&\ge
\frac25\frac{64}{25}T_j
-\left(\frac25+\frac{10}{39}\right)T_j\\
&=\boxed{\frac{1792}{4875}T_j}.
\end{aligned}
\]

Therefore after \(L\) interior synchronization steps,

\[
\begin{aligned}
\Delta s_L
&\ge
\frac{1792}{4875}T_0
\sum_{j=0}^{L-1}\left(\frac{64}{25}\right)^j\\
&=\boxed{
\frac{1792}{7605}T_0
\left[
\left(\frac{64}{25}\right)^L-1
\right].
}
\end{aligned}
\]

This corrects the earlier, overly strong coefficient \(10/39\) in the cumulative
backward-displacement formula.  The synchronization recurrence and common-slice
construction are unchanged.  Crucially, the corrected coefficient is still
strictly positive, so an ancestry starting at finite physical time cannot have
infinite interior causal depth.  It either stops earlier through
inheritance/residual/reuse, or its next common slice crosses \(t=0\).

At that point the adjoint identity is simply integrated on the truncated interval
\([0,t]\).  There is no missing negative-time packet and no interior fresh-grain
charge. The node is an **initial-boundary root**.

For an initial dyadic band \(M\),

\[
\|P_Mu_0\|_2^2
\le M^{-2m}\|u_0\|_{\dot H^m}^2.
\]

If every selected initial root has critical mass \(ME_a\ge\eta\), then

\[
\boxed{
\#\mathcal R_M(0)
\le
\eta^{-1}M^{1-2m}
\|u_0\|_{\dot H^m}^2.
}
\]

This is a boundary-band estimate, not by itself a scale-independent global reset
count.

## 4. Material labels at the common slice

A causal parent label is transported backwards with the same physical Kelvin/affine map used by its adjoint propagator. Common affine motion preserves the intrinsic coherent coordinate

\[
\zeta=(L^{-1}X/2,L^Tk),
\]

so it creates no interface cost.

At the common physical slice \(s_j\), all parent labels can be compared in one coherent ancestry graph. Two existing interface theorems give the only representation charges:

1. changing the selected material coherent set costs its symmetric-difference Moyal energy \(E(S_+\triangle S_-)\);
2. changing a Gaussian covariance representative by log-distance \(d_{\log}\) costs at most
   \[
   \frac{d_{\log}}{\sqrt2}\|u\|_2^2.
   \]

Thus one ordered layer boundary carries the single registration charge

\[
\boxed{
\Xi_j^{\rm reg}
=
E(S_+\triangle S_-)
+rac{d_{\log,j}}{\sqrt2}E_j.
}
\]

It is entered **once** at that boundary. Exact common transport is free; a non-small covariance change is not hidden here and remains in the previously named strain/source/fresh-relink alternatives.

## 5. Scope

This theorem closes the **time-slab synchronization geometry** of the adjoint causal construction and supplies a canonical initial-boundary rule. It does not by itself prove the final global telescoping theorem. The remaining continuum work is now narrower:

- identify the generated parent-pair weights after each adjoint stopping step with the same physical transfer law used by the Shannon/Renyi causal modules, up to the already-certified near-extremal change of measure;
- prove the sum of the ordered registration charges above, together with the already existing spatial/frequency/profile interfaces, fits inside one global summable \(\Xi\) ledger;
- telescope flat erosion, coherent stopping, causal reuse and source/dissipation costs to the master scale path.

No persistence hypothesis, raw packet count, or negative-time extension is used.
