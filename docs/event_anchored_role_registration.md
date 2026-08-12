# Event-anchored hard transfer roles inside a smooth moving PDE envelope

The outer-role equation and non-affine interface theorem leave one apparent
seam.  The transfer/Young ledger wants exact disjoint Fourier/helical cells,
while a differentiated moving PDE multiplier should be smooth.  Those two jobs
should not be forced onto the same operator.

The natural construction is **hard core, smooth envelope**.

## 1. Hard roles belong to the physical event

At a fixed physical interaction time partition the crossing frequency region
into deterministic disjoint Borel cells `C_a`.  On the divergence-free fiber one
may further use the exact pointwise helical projections `H_+(xi),H_-(xi)`.
Define

\[
\boxed{P_{a\sigma}=1_{C_a}(D)H_\sigma(D).}
\]

These are self-adjoint orthogonal projections and, on the covered high-frequency
subspace,

\[
\boxed{\sum_{a,\sigma}P_{a\sigma}=I.}
\]

Consequently the high field and its quadratic source disintegrate exactly:

\[
h=\sum_{a,\sigma}h_{a\sigma},
\qquad
\mathcal B(h,h)=
\sum_{a,b,\sigma,\tau}
\mathcal B(h_{a\sigma},h_{b\tau}).
\]

There is no wave-packet synthesis coefficient and no overlap constant.

A hard partition is also perfectly compatible with smooth-symbol freezing.  The
SGS/Leray/helical multiplier is smooth; its value may be frozen on a Borel
product cell just as on any other set of the same diameter.  The error is

\[
A_3L_Mh\prod_i\|f_i\|_{3/2},
\]

and comes solely from variation of the **physical multiplier**.  There is no
extra cutoff-smoothing error.

## 2. The smooth envelope belongs only to the PDE between slices

Do not differentiate the hard event boundary.  Given one selected hard frequency
cell `P`, choose a scalar smooth Fourier multiplier `Q` with

\[
\boxed{QP=P},
\]

so `Q=1` on the hard cell.  The hard signed-good role starts at radial frequency
at least `3N/5`.  Choose the smooth envelope with radial support beginning at

\[
\frac{11}{20}N.
\]

There is an anchor plateau of width `N/20`.  If the coherent affine strain action
is at most `1/30`, dual-affine transport contracts frequency radius by at most
`e^{-1/30}`.  Arb certifies

\[
\boxed{
\frac{11}{20}e^{-1/30}N>\frac12N.
}
\]

Since `V=S_(N/4)u` has `V tensor V` supported in `B_(N/2)`, the smooth envelope
still excludes low--low forcing throughout the whole low-strain slab.

Thus the smooth `Q` may be used in the exact moving-role PDE with ordinary
Egorov calculus, while the hard `P` remains only an event observable.

## 3. Exact coefficient registration: no contamination from the envelope

For any field `u` and terminal probe `phi`, self-adjointness and `QP=P` give

\[
\boxed{
\langle Pu,\phi\rangle
=
\langle Qu,P\phi\rangle.
}
\]

This is the key registration identity.  Christ/Young is applied to the actual
hard physical-transfer role `Pu`.  The resulting dual Gaussian probe is simply
projected by the same hard role when inserted into the smooth outer PDE.
No field component in the overlap part of the smooth envelope can cancel the
coefficient.

Moreover, in the frequency variable both operations are pointwise contractions:

\[
|1_C(\xi)H_\sigma(\xi)v|\le |v|.
\]

Hence for every `1<=p<=infinity`, in particular `p=3` and `p=2`,

\[
\boxed{\|P\phi\|_p\le\|\phi\|_p.}
\]

So hard registration cannot worsen either the `L^3` duality constant or the
`L^2` normalization used by the amplitude-productivity theorem.

## 4. Helicity is a terminal fiber mark, not a persistent projector

The outer PDE envelope needs only the scalar frequency carrier.  At the physical
event the hard role may be split helically to expose the scalar complex Young
interaction.  Its selected polarization enters the terminal dual vector.
Between the event and the common causal slice, the already exact adjoint Kelvin
fiber equation transports that vector.

Therefore one never assumes

\[
H_\sigma(t)u(t)
\]

persists as a fixed helical packet.  Rotation is gauge; symmetric strain remains
in the Kelvin generator and in the existing polarization/strain currencies.

## 5. Complete high--high source, designation only at work level

With `h=u-S_(N/4)u`, the outer PDE source should be the complete

\[
-Q\mathbb P\nabla\cdot(h\otimes h).
\]

Do not remove all but one parent pair at the forcing-norm level.  The hard roles
above decompose the **actual work** of this full source exactly.  Positive atoms
are then selected by the physical transfer law; bad/cancelling atoms are already
transfer/backscatter cost.

This removes the need for an unproved `L^2` bound on “all other high--high
interactions”.  A designated near-extremal parent pair is a work atom, not a
separate PDE forcing term.

## 6. Single-charge interpretation

The two operators have deliberately different jobs:

- `P`: physical event identity, exact transfer, Young role, dual coefficient,
  orthogonal energy accounting;
- `Q`: smooth PDE envelope only.

The equality `QP=P` registers them without a new measure.  Overlap of smooth
PDE envelopes is therefore **not** another physical transfer partition and must
not be entered into `Xi`.

The only `Xi` here remains the pre-existing variation of the physical
SGS/Leray/helical symbol on the hard cell.

## 7. Remaining scope

This closes the frequency-role / transfer-cell alignment seam at the exact event
and coefficient levels.  Combined with the exact outer moving-role equation, the
smooth-envelope continuation is read through the quadratic carrier energy
`<u,Q^2u>`.  Completing `Q` by a smooth square partition preserves `QP=P` on the
event plateau while making moving/skew work conservative and symmetric work the
same resolved strain provenance.  No idempotence of `Q` is asserted.

A large interface coefficient impulse is only a first-hit locator.  Actual
terminal carrier energy and native interface work must reenter the
physical-energy gate before inheritance, HH generation, strain or relink becomes
a causal owner.  Global assembly of the resulting physical owners remains open.
No global-regularity claim is made.

## Canonical positive edge-work inheritance

When a hard event role is downstream of the continuum helical edge law, the role map is a measurable label map `pi` on already-registered physical Fourier/helicity edges.  Its causal mass is therefore inherited from the fixed Hahn law,
\[
\boxed{\nu=\pi_\#(dW^+),}
\]
while the signed hard-cell trilinear work is
\[
T_C=(\pi_\#dW)(C).
\]
These are deliberately different objects.  In general `(pi_#dW)^+ <= pi_#(dW^+)`, with strict inequality when physical cancellation is hidden by the hard cell.  The hard role may expose the signed `T_C` to Young/Christ for saturation diagnostics, but it may not promote `(pi_#dW)^+` into a second causal law.

The certified positive-edge router further separates inherited good and bad positive work.  Geometry-bad positive work has already terminated through the same physical transfer-loss interface.  Geometry-good positive work is only Young-eligible, and the full signed `T_C` may bind it to Young only when the hard cell is fate-pure, i.e. carries zero inherited bad-positive mass.  A mixed-fate hard cell remains an unresolved handoff rather than being refined by observer choice.  Thus event geometry, causal inheritance, signed cancellation and Young marking remain distinct physical layers.
