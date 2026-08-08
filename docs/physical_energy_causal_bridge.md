# Physical-energy causal bridge: Duhamel gives support, energy work gives the weights

The single-charge compiler isolated one apparent remaining measure bridge:
identify the positive adjoint Duhamel generation law with the positive physical
child-energy transfer law.  That identification is **false in general**, even in
the flattest possible coefficient dynamics.  The correct architecture is simpler
and more physical: use Duhamel to expose the same-time quadratic parent-pair
event support, and use the exact child-energy balance to weight those events.

## 1. Exact countermodel to raw kernel equality

Take a scalar phase-locked child coefficient on `[0,1]` with

\[
G=0,\qquad R=0,\qquad c(0)=0,\qquad \dot c=1.
\]

The terminal adjoint is `psi=1`, so normalized amplitude generation is

\[
\boxed{d\Gamma=dt.}
\]

But `c(t)=t`, and the actual positive child-energy work is

\[
 d\mathcal T=2c\dot c\,dt=2t\,dt,
\]

which is already normalized because `int_0^1 2t dt=1`.  Hence

\[
\boxed{\Gamma([0,t])=t,\qquad \mathcal T([0,t])=t^2.}
\]

At half-time the masses are `1/2` and `1/4`.  Thus even perfect phase lock, zero
strain, zero viscosity and zero residual do **not** make the two normalized laws
equal.  The difference is physical: amplitude is linear in accumulated forcing,
whereas energy is quadratic in amplitude.

This disproves the previous proposed bridge
`K_Gamma = K_phys` from the current packet-ODE hypotheses.  It is not a
Navier--Stokes counterexample; it is a countermodel to an unnecessarily strong
intermediate theorem.

## 2. The correct adjoint object is an energy lift

For the vector coefficient equation let `psi` be any nonzero adjoint response,

\[
z=\langle\psi,c\rangle,\qquad n=\|\psi\|^2,
\qquad c_\parallel=\psi z/n,
\qquad c_\perp=c-c_\parallel.
\]

For every individual high--high forcing atom `F_alpha`, exact Hilbert geometry
gives

\[
\boxed{
2\Re\langle c,F_\alpha\rangle
=
\frac{2}{n}\Re\!\left(\overline z\,\langle\psi,F_\alpha\rangle\right)
+2\Re\langle c_\perp,F_\alpha\rangle.
}
\]

The first term is the **energy-lifted adjoint response**.  Raw `dGamma` keeps only
the linear response atom `Re<psi,F_alpha>` and therefore cannot itself be the
energy measure.  The second term is actual child-role cross work; if one chooses
the adjoint response direction as the designated child profile, this is precisely
a third-role profile/cross contribution, not a second source.

In the scalar flat model the lift reduces to the exact size bias

\[
\boxed{d\mathcal T=2U\,d\Gamma=d(U^2),}
\qquad U(t)=\Gamma([0,t]).
\]

So the mismatch is not noise.  It is the quadratic passage from amplitude to
energy.

## 3. The selected-child energy balance supplies the physical causal measure

On one exact selected role write

\[
\dot c=G(t)c+F_{HH}(t)+R_{class}(t),
\]

where the objective Kelvin generator has

\[
G=-S_\perp-\nu|k|^2I.
\]

Put `E=||c||^2` and

\[
K=\int\|S_\perp\|_{op}dt.
\]

Viscosity is nonpositive in the energy identity, hence

\[
\frac{dE}{dt}
\le
2\|S_\perp\|_{op}E
+2[\Re\langle c,F_{HH}\rangle]_+
+2[\Re\langle c,R_{class}\rangle]_+.
\]

Define the integrated positive works `W_HH^+` and `W_R^+`.  Gronwall gives

\[
\boxed{
E_1\le e^{2K}\big(E_0+W_{HH}^++W_R^+\big),
}
\]

and therefore

\[
\boxed{
W_{HH}^+\ge e^{-2K}E_1-E_0-W_R^+.
}
\]

This is a physical energy statement; no identification with an adjoint
probability law occurs.

## 4. Clean low-strain causal gate

The existing service-or-flat architecture already stops at

\[
K>1/30
\]

and routes that event to critical `D_V`.  On the complementary branch,

\[
e^{-2K}\ge1-2K\ge14/15.
\]

Hence there is an exact three-way energy gate:

1. if `E_0 >= E_1/5`, the child has definite **material inherited energy**;
2. if `W_R^+ >= E_1/5`, an already-classified residual/source/interface carries
   definite **physical work**;
3. otherwise
   \[
   \boxed{W_{HH}^+\ge\frac{8}{15}E_1.}
   \]

The third branch is genuine physical high--high generation in the energy sense.
It is stronger for the causal master than the linear amplitude statement because
it directly produces the conserved quantity whose weights are used by the
transfer/Hodge ledger.

The first branch is sticky critical ancestry, not a uniform energy reset.  The
second delegates to the causal root already represented by `R_class`.  Neither
is double charged.

## 5. Same causal event space, physical weights

Decompose the designated high--high source into its exact quadratic parent-pair
atoms,

\[
F_{HH}=\int F_{HH,\alpha}\,d\alpha.
\]

Each atom lives at one physical interaction time and contains the two parent
roles.  Define

\[
\boxed{
d\mathcal T_{HH}(t,\alpha)
=2\big[\Re\langle c(t),F_{HH,\alpha}(t)\rangle\big]_+\,dt.
}
\]

This is **actual positive child-energy work**.  On the signed-good forward core
it is the same physical child-transfer quantity already used by the smooth SGS
midgap, Hodge and flat-episode theorems.  Negative work is real backscatter and
remains in the transfer/cancellation ledger.

The causal support is nevertheless the same one exposed by Duhamel: the atom is
the same `F_HH,alpha`, with the same two parents at the same time.  What changes
is only the weight assigned to that event.

## 6. Asynchronous synchronization is measure agnostic

The parabolic synchronization theorem uses only:

- positivity of the event weights;
- the fact that all events lie in one child slab;
- the signed-good scale law `3/5 < N_parent/N_child < 5/8`.

Its first step is merely the fact that one temporal half carries at least half of
**any positive measure**.  Therefore choose the heavier half using
`dT_HH`, not `dGamma`.  Its parent-event span still satisfies

\[
\boxed{\alpha_1\le25/128.}
\]

The later recurrence

\[
\alpha_{j+1}\le\frac{25}{64}(\alpha_j+2/5)
\]

and the sharp cone `alpha<=10/39` depend only on support geometry, so they are
unchanged.

This is the key redesign:

\[
\boxed{
\text{Duhamel / quadratic PDE source: causal parent-pair support},
\qquad
\text{child energy balance: physical causal weights}.
}
\]

There is no reason to force the two measures to be identical.

## 7. Shannon/Renyi now receives the physical law by construction

Normalize `dT_HH` on the selected generated layer and call the resulting law
`w`.  Duplicate each child event into its two structural parent slots with
weight `w/2`, then push through the canonical material parent label

\[
\zeta=(L^{-1}X/2,L^Tk).
\]

This is exactly the input assumed by the existing transfer-weighted
Shannon/Renyi theorems.  No change-of-measure from `dGamma` is needed.  The
`1/2` two-role baseline stays free, while actual parent merging creates the
existing hidden-pair / entropy / cycle currencies.

Thus the previous Duhamel-to-physical **kernel equality** bridge disappears.  It
is replaced by a physical-energy causal gate plus a measure-agnostic support
synchronization theorem.

## 8. Compiler consequence

For a low-transfer-cost recursive node the preferred causal order is now:

\[
\text{high strain}
\to D_V,
\]

or on low strain,

\[
\boxed{
\text{material energy inheritance}
\ \lor\
\text{classified residual physical work}
\ \lor\
\text{positive physical HH transfer }\ge8E_1/15.
}
\]

Only the last branch expands into a weighted parent-pair layer, and its weights
are already `dT_HH`.  Duhamel remains useful as an exact interaction-picture
identity and amplitude diagnostic, but **master readiness must not require raw
`dGamma` to equal physical transfer**.

## 9. Exact binary coherent atomization

The companion theorem `recursive_coherent_witness_extraction.md` removes a further packetization step.  Once the outer selected parent/child roles exist, the positive coherent resolution `sum_C A_C=I` and bilinearity of the quadratic source give exact work atoms

\[
W_{CDE}=2\Re\langle A_Ew_3,\mathcal N(A_Cw_1,A_Dw_2)\rangle,
\]

with `sum W_CDE` equal to the actual selected high--high child work.  Their positive parts therefore give a binary parent-pair/child material law directly, while negative parts are physical backscatter.  This means Duhamel is not needed even to *invent* a packet graph: quadratic Navier--Stokes work itself already supplies the binary event support after coherent disintegration.

The coherent pieces are not treated as compact Fourier projections; the outer selected Fourier/helical role continues to carry the interscale geometry, with representative errors in the existing symbol/covariance `Xi`.

## 10. Scope

This theorem closes a structural measure mismatch at the exact selected
packet/coefficient level.  It does **not** prove that every efficient continuum
Navier--Stokes block has already been extracted with the required exact moving
coherent coefficient equation, selected `F_HH`, canonical material labels and
summable rejected work.  That recursive PDE witness extraction remains the main
continuum bridge.

No global-regularity claim is made.
