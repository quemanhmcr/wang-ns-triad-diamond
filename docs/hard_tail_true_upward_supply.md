# Hard-tail true upward supply: radial boundary feeding before UV continuation

Status: **candidate theorem block; certification requires dedicated, independent adversarial/actual-NS, and full causal-integration GitHub Actions gates.**

Candidate status string:

`EXACT_HARD_TAIL_TRUE_UPWARD_SUPPLY__RADIAL_BOUNDARY_STOCK_FLOW_VISCOSITY__CANONICAL_LOW_TO_HIGH_SUPPLY_IN_COMMON_N_DW_UNIT__PURE_UV_HH_ONLY_FIRST_DYADIC_SHELL__DEEP_UPWARD_HAS_RESOLVED_SCALE_PARENT_CONTACT__NO_INTERNAL_CIRCULATION_REENTRY_OR_INTERFACE_OVERCLAIM`

This theorem is downstream of the certified radial spectral-crossing law and upstream of any new recurrence inference.  Its purpose is to replace the old phrase “positive nonlinear tail regeneration” by the more precise physical observable that actually feeds tail stock: **low-to-high boundary crossing**.

## 1. Start from the exact radial Navier--Stokes law

For

\[
H_N=\{(k,s):|k|>N\},
\qquad
E_{>N}(t)=\sum_{(k,s)\in H_N}E_{k,s}(t),
\]

the certified cyclic donor flow gives

\[
\Phi_\uparrow(N)
=\mathcal M\{|k_d|\le N<|k_r|\},
\qquad
\Phi_\downarrow(N)
=\mathcal M\{|k_r|\le N<|k_d|\}.
\]

The crucial point is that high→high cyclic traffic is **not** in either boundary flow.  It is real nonlinear redistribution inside the tail and cancels from tail divergence.

With the existing normalized hard-tail dissipation

\[
D_{tail}=N\int_I\|\nabla P_{>N}u\|_2^2dt,
\]

the exact radial continuity identity, multiplied by the parent-tail scale `N`, is

\[
\boxed{
N E_{>N}(t_1)
+2\nu D_{tail}
+N\int_I\Phi_\downarrow(N,t)dt
=
N E_{>N}(t_0)
+N\int_I\Phi_\uparrow(N,t)dt.
}
\]

Every term has the same physical energy×frequency unit.  No shell-dependent causal reweighting appears.

Because final stock and downward flux are nonnegative,

\[
N E_{>N}(t_0)+N\Phi_\uparrow\ge2\nu D_{tail}.
\]

Therefore, with exact ties retained,

\[
\boxed{
N E_{>N}(t_0)\ge\nu D_{tail}
\quad\text{or}\quad
N\Phi_\uparrow\ge\nu D_{tail}.
}
\]

The first is inherited physical tail stock.  The second is **true upward nonlinear supply**.  This owner cover does not use gross positive tail work.

## 2. Why this sharpens the old high-tail bridge

The older high-frequency bridge used

\[
N E_{>N}(t_0)+N W_{>N}^{+}\ge2\nu D_{tail}.
\]

That inequality remains valid, but `W_>^+` includes positive work on modes already inside the tail.  A closed high-frequency triad can carry arbitrarily large high→high traffic while crossing no radial boundary at all.

The radial theorem identifies the exact correction:

\[
W_{>N}^{+}=I_{high}+\Phi_\uparrow,
\qquad
W_{>N}^{-}=I_{high}+\Phi_\downarrow.
\]

The same internal `I_high` appears on both sides and is not tail supply.  Thus the new owner law does **not** discard physical circulation; it simply records its actual meaning and refuses to charge it as UV feeding.

## 3. Disintegrate `Phi_up`, not a new Hahn law

Take one already-certified cyclic donor atom

\[
(k_d,s_d)\longrightarrow(k_r,s_r),
\qquad |k_d|\le N<|k_r|,
\]

with physical mass inherited from the canonical `dW^- -> dW^+` donor kernel.

Assign the recipient to the unique boundary-anchored dyadic hard shell

\[
\frac{M}{2}<|k_r|\le M,
\qquad M=2^jN,\quad j\ge1.
\]

This is a positive restriction/pushforward of the same `Phi_up` measure.  The causal work unit remains

\[
\boxed{N\,dW.}
\]

The theorem forbids replacing it by `M dW`.  Recipient shell scale is geometry, not a new probability unit.

## 4. The energy donor is one of the recipient's interaction parents

For a closed helical triad, root `r` is the recipient child edge and root `d` is an energy donor.  Since the roots are distinct, the closed donor mode is one of the two interaction parents of the recipient edge.  The donor **child** mode is its reality partner, so radii are unchanged.

Thus the physical radial donor radius is also one actual interaction-parent radius of the recipient edge.

This fact is what makes the support theorem rigid.

## 5. Exact support split

At recipient shell `M`, define two positive restrictions of the same true-upward law.

### 5.1 Pure-UV HH-by-support

Both recipient interaction parents satisfy

\[
|k_{p_1}|>M/4,
\qquad
|k_{p_2}|>M/4.
\]

On this region the existing strict transporter `S_{M/4}` vanishes on both parent frequencies, so `h=u` there exactly.  This is the **pure UV HH support region**.

### 5.2 Resolved-scale parent contact

At least one interaction parent satisfies

\[
|k_p|\le M/4.
\]

This is only a Fourier support statement.  The current smooth cutoff theorem certifies `supp S_{M/4} subset B_{M/4}` and `|S|<=1`; it does **not** certify that `S=1` on every low mode.  Therefore this theorem deliberately does **not** rename resolved-scale parent contact as “interface work”.  Such a name would require a separate positive binding through the actual smooth `V/h` decomposition.

The two support restrictions partition `Phi_up` exactly.

## 6. Pure upward HH can only feed the first dyadic shell

Suppose a true-upward atom is pure-UV HH-by-support.  Its energy donor is an interaction parent, so

\[
|k_d|\le N
\quad\text{and}\quad
|k_d|>M/4.
\]

If `j>=2`, then `M>=4N`, hence `M/4>=N`, contradicting those two inequalities.

Therefore

\[
\boxed{M=2N.}
\]

This is not a norm estimate.  It is exact Fourier support geometry of the physical energy donor.

Moreover

\[
\boxed{M/4<|k_d|\le M/2.}
\]

Let `p` be the other interaction parent and `z` the recipient child.  Triad closure gives

\[
|p|\le |z|+|k_d|\le M+M/2=3M/2.
\]

Pure-UV support already gives `|p|>M/4`.  Thus **both** interaction parents obey

\[
\boxed{M/4<|k_{p_i}|\le3M/2.}
\]

So true-upward pure-UV HH supply is automatically comparable and lives only on the first shell `M=2N`.  No ultraviolet-locality norm estimate is needed to obtain this comparability.

## 7. Deep direct upward crossing has resolved-scale parent contact

If the recipient belongs to any deeper shell

\[
M\ge4N,
\]

then the energy donor itself satisfies

\[
|k_d|\le N\le M/4.
\]

Therefore

\[
\boxed{
\text{every direct upward atom into }M\ge4N
\text{ has resolved-scale parent contact.}
}
\]

Again, this does not yet say that its positive work is an existing smooth-interface owner.  It says only that deep direct boundary crossing cannot lie in the pure-UV HH support region.

This distinction is essential because the smooth cutoff may have a transition region below `M/4`.

## 8. Upward-owner support alternative

If the exact tail identity selects the true-upward owner,

\[
N\Phi_\uparrow\ge\nu D_{tail},
\]

and

\[
\Phi_\uparrow=\Phi_{pureUV}+\Phi_{contact},
\]

then one positive support restriction obeys

\[
\boxed{
N\Phi_{pureUV}\ge\frac{\nu D_{tail}}2
\quad\text{or}\quad
N\Phi_{contact}\ge\frac{\nu D_{tail}}2,
}
\]

with exact ties joint.

The first branch is first-shell and automatically comparable by §6.  The second branch is a genuine resolved-frequency-contact law whose correct PDE ownership remains a later theorem.

## 9. Anti-theorem: high-frequency circulation is not regeneration supply

Place one regular physical equiradial closed triad wholly above a lower radial boundary.  Its cyclic donor flow has

\[
I_{high}>0,
\qquad
\Phi_\uparrow=\Phi_\downarrow=0.
\]

Scaling actual amplitudes by `lambda` scales the trilinear internal work by `lambda^3` while both boundary flows stay zero.

Hence no amount of high→high positive work proves direct feeding of tail stock.  This is why `Phi_up`, rather than `W_>^+`, is the native input to the sharpened high-tail continuation.

## 10. What actual Navier--Stokes CI must test

The evolved PDE audit keeps distinct observables distinct.

1. On the existing cutoff-7 mixed-fate trajectory, it retains the full radial-tail stock/signed-work/viscosity reading on FFT grids `24,28`.  Separately, the selected actual cyclic triad at `N=8` is checked to be true upward **pure-UV first-shell** supply with automatic parent comparability.
2. A separate real divergence-free cutoff-2 Galerkin state embeds the closed physical triad
   \[
   (-2,-1,-1)+(1,0,0)+(1,1,1)=0.
   \]
   At `N=1`, the high recipient has radius `sqrt(6)` and lies in `M=4`, while the radius-one donor touches `M/4` exactly.  The state is evolved by the repository RK4 Navier--Stokes solver on two FFT representations, and the actual evolving triad is required to begin with positive deep upward work and zero pure-UV classification.
3. A second six-mode **orthogonal Fourier--Galerkin Navier--Stokes** audit keeps exactly one closed triad and its reality partner as the retained Galerkin subspace.  At every output time it enumerates both reality-partner closed triples and all eight helicity sectors, reconstructs the **full** cyclic donor boundary law of that finite PDE, and requires
   \[
   \Phi_\uparrow-\Phi_\downarrow
   =W_{>N}^{signed}
   \]
   against the direct Galerkin nonlinear term.  It then integrates the full `Phi_up/Phi_down`, stock and viscosity and applies the inherited-or-true-upward owner cover itself.  This is the load-bearing actual-PDE referee for the sharpened owner law; no selected triad stands in for the full finite system.
4. Representation comparisons are normalized only by native physical work/energy envelopes.  A vanishing sub-observable is never used as a denominator.

The selected triads are sub-observables of the evolved PDE.  They are not substituted for the full tail law.

## 11. Scope

This theorem would close the seam

`physical hard-tail dissipation -> inherited stock OR gross positive tail regeneration`

into the sharper physical statement

`physical hard-tail dissipation -> inherited stock OR true low-to-high boundary supply`.

On the upward branch it would additionally prove the exact support alternative

`pure-UV first-shell comparable HH OR resolved-scale parent contact`.

It would **not** yet:

- bind resolved-scale contact to the existing smooth `K/S` interface owner;
- identify the pure-UV submeasure with an aggregate hard-shell Hahn law formed later;
- apply Young/Christ to a donor-restricted submeasure without a proved compatible pushforward;
- convert radial crossing into an event count or minimum progress tax;
- make `D_tail`, `N Phi_up`, or shell mass an additive reset;
- prove mixed-owner recurrence termination or Navier--Stokes global regularity.

The next question after certification should be decided by which support branch survives actual work routing.  The pure-UV branch is already first-shell and comparable; the resolved-contact branch first needs a type-correct positive interface/repartition binding.  No shortcut should merge those two physical phenomena.
