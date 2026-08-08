# Bargmann root-cell registration: parent identity by actual coherent energy

The dual-Gaussian theorem produces an actual large coherent coefficient near the
complex Young parent mark.  A remaining worry was that the **positive nonlinear
work atom** need not be localized in the very same coherent cell.  That is the
wrong requirement.

A causal parent slot needs two logically different objects:

1. its positive weight comes from actual quadratic child-energy work;
2. its label identifies which physical parent reservoir is being reused.

The label therefore only needs to be an actual energetic material anchor for the
parent role.  Gaussian coherent analysis supplies such an anchor quantitatively.

## 1. Bargmann submean turns a coefficient into local Moyal energy

Use a normalized Gaussian coherent family whose intrinsic coordinate is

\[
\zeta\in\mathbb R^6\simeq\mathbb C^3
\]

and whose overlap has magnitude

\[
|\langle g_\zeta,g_\omega\rangle|
=e^{-|\zeta-\omega|^2/2}.
\]

After a Weyl translation, put the marked point at zero.  The coherent coefficient
has the Bargmann form

\[
c(z)=e^{-|z|^2/2}F(z),
\]

where `F` is entire holomorphic in `C^3` and the Moyal measure is normalized by
`pi^-3 dz`.

For a ball `B_R`, the holomorphic mean-value inequality gives

\[
\int_{B_R}|F(z)|^2\frac{dz}{\pi^3}
\ge
\frac{R^6}{3!}|F(0)|^2.
\]

Since `e^{-|z|^2}>=e^{-R^2}` in the ball,

\[
\boxed{
E(B_R)
\ge
e^{-R^2}\frac{R^6}{3!}|c(0)|^2.
}
\]

This is a local statement about the **actual positive Moyal energy measure**.

## 2. The optimal fixed ball

The function

\[
e^{-R^2}R^6/3!
\]

is maximized at

\[
R^2=3.
\]

Hence

\[
\boxed{
E(B_{\sqrt3})
\ge
\frac{9}{2e^3}|c(0)|^2.
}
\]

The coefficient cannot be an arbitrarily narrow analysis spike.

## 3. Pass to one canonical material cell

Use the fixed unit dyadic grid in the intrinsic six-dimensional coordinate.  A
ball of radius `sqrt3` projects to an interval of length `<4` in each real
coordinate, so it intersects at most

\[
5^6
\]

unit cells.

Therefore one actual canonical cell `C` satisfies

\[
\boxed{
E_C
\ge
\frac{9}{2e^3\,5^6}|c(0)|^2.
}
\]

With the dual-Gaussian root quantum

\[
N|c(0)|^2\ge\eta_{dual}>1/5,
\]

this gives a fixed positive critical cell mass

\[
\boxed{
NE_C
\ge
\eta_{cell}
:=\frac{9}{2e^3\,5^6}\eta_{dual}>0.
}
\]

The constant is intentionally not optimized.  Causal reuse needs positivity
independent of scale, not a large numerical value.

Choose the maximum-energy cell among the finitely many cells meeting the ball;
use a fixed lexicographic rule only for exact energy ties.  This produces a
deterministic **energy anchor** from the physical Moyal measure itself.

## 4. Push physical causal weights to energy anchors

Each generated high--high event already has two physical parent roles and a
positive child-work weight.  For each parent role, perform the complex Young /
dual-Gaussian analysis and attach its energy-anchor cell.

Push the positive parent-slot law through this map.

This does not change total physical transfer mass.  If two distinct parent slots
map to the same energy anchor, their weights add: that is precisely parent reuse
at the chosen physical resolution.  If they map to distinct anchors, those roots
carry distinct positive cell-energy quanta.

No assertion is made that the nonlinear work itself occurred inside the anchor
cell.  The weight and the identity label are different physical observables:

- work says **how much child energy this parent pair generated**;
- the anchor says **which coherent parent reservoir carries the role**.

This distinction removes an unnecessarily strong transfer-localization demand.

## 5. Global energy budget for distinct anchor cells

For one exact outer Fourier/helicity subrole and one covariance representative,
positive coherent cell energies sum exactly to the subrole's `L2` energy.

The selected frozen frequency subroles inside one scale bin form an exact
disjoint partition.  Logarithmic scale bins are colored with four colors so that
same-color outer shells have disjoint Fourier support.  Thus they are orthogonal.
For covariance representatives one pays only the finite compact-cover count.

Consequently all distinct material energy anchors obey

\[
\sum_{root}E_{C(root)}
\le
P_{cell}E_{global},
\]

with

\[
\boxed{
P_{cell}
=4\,(\#\text{covariance representatives})<\infty.
}
\]

There is no phase-space coloring or Gaussian synthesis/Riesz factor here: within
one coherent frame the Moyal partition already has `P=1`.

If each root has `N_aE_a>=eta_cell` and

\[
N_a\le N_{max},
\]

then

\[
\boxed{
 n_0
\le
\frac{P_{cell}E_{global}N_{max}}{\eta_{cell}}.
}
\]

This is exactly the form needed by the binary/Rényi root telescope.  The large
constant changes only the finite logarithmic depth offset.

## 6. Canonical-label interpretation

The parent label can now be the energy-anchor material cell itself.  It is not a
Gaussian synthesis index.

At a common causal slice:

- the outer role and complex Young mark are extracted from the actual field;
- the dual coefficient locates a nearby physical coherent energy ball;
- the maximum-energy canonical cell in that ball is selected;
- Duhamel/physical-work parent slots are pushed to this same label;
- Shannon/Rényi reuse operates on that pushforward;
- Hodge/holonomy may use the same material node with their already controlled
  frequency/helical representatives.

Common affine motion remains gauge because the intrinsic coherent coordinate is
materially invariant.  Small representative covariance/symbol changes retain
the existing summable `Xi` ledger.

## 7. Remaining scope

This theorem closes the **energy-label part** of transfer-cell registration.
What remains for the final continuum assembly is operational rather than a new
currency: prove that every recursive efficient physical block is decomposed into
the exact disjoint/orthogonal frozen outer roles assumed above, apply the complex
Young/dual/Bargmann anchor construction measurably at the common causal slice,
and use that anchor as the canonical material label throughout the compiler.

No global-regularity claim is made.

## Superseding amplitude-entropy use

The canonical cell quantum here is used in the homogeneous form

\[
N E_C\ge\beta\alpha^2,
\qquad \alpha=\sqrt N|\langle u,\phi\rangle|,
\]

rather than as an absolute root mass.  Combined with the transfer-weighted
log-amplitude recursion, this is exactly what is needed for the root
energy--entropy inequality.  The amplitude imbalance mentioned above is
therefore closed by `amplitude_entropy_causal_reuse.md`; no new service currency
or amplitude floor is introduced.
