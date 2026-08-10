# Common-slice coefficient registration is a first-stop theorem, not packet persistence

The amplitude--entropy causal theorem needs one scale-independent statement:
a parent coefficient observed at its quadratic interaction time must either be
registered on the synchronized common slice with a fixed fraction, or the
attempt to register it must expose an earlier obstruction.  A coefficient
obstruction receives a physical owner only after energy/work reentry.

No new persistence hypothesis is needed.

## 1. The common slice lies inside every natural parent window

After the first heavy-half step, asynchronous synchronization gives

\[
\alpha_j\le\frac{10}{39}.
\]

With

\[
s_j=a_j-\frac25T_j^{\min},
\]

every parent event `t in [a_j,b_j]` satisfies

\[
\frac{t-s_j}{T_j^{\min}}
\le\frac25+\frac{10}{39}
=\frac{128}{195}<1.
\]

The clean remaining natural-window fraction is therefore

\[
\boxed{1-\frac{128}{195}=\frac{67}{195}.}
\]

So the exact adjoint gate is available from the common slice to every event in
the layer.

## 2. Apply the same adjoint gate to the parent mark

For one selected parent coefficient in the common Kelvin interaction picture,

\[
z(t)=z(s)+I_{HH}[s,t]+I_R[s,t].
\]

Set `A=|z(t)|`.  The existing exact triangle gate gives

\[
|z(s)|\ge A/4
\quad\lor\quad
|I_R|\ge A/4
\quad\lor\quad
|I_{HH}|\ge A/2.
\]

Interpret this as **registration by first obstruction**.

- If `R_class` reaches its threshold first, registration stops and the same
  smooth carrier reenters its actual `Q^2` energy/work gate.
- If high--high amplitude generation reaches its threshold first, registration
  stops at a coefficient obstruction.  It becomes a generated node only if the
  physical-energy gate selects actual positive HH work; raw Duhamel mass never
  supplies the causal weight.
- If a genuine material coherent-cell change occurs first, it is a physical
  relink/fresh event.
- If none of these happens, the only remaining branch is
  \[
  \boxed{|z(s)|\ge A/4.}
  \]

Thus "the packet persists" is never assumed, and a coefficient obstruction is
never silently promoted to a physical owner.

## 3. Product consequence

If both parents of one continuing generated child reach the common slice without
an earlier stop,

\[
\boxed{
\alpha_{p_1}(s)\alpha_{p_2}(s)
\ge\frac1{16}
\alpha_{p_1}(t)\alpha_{p_2}(t).
}
\]

This is the `1/16` factor used in
`docs/amplitude_entropy_causal_reuse.md`.

The factor is conservative.  Its role is only to be strictly positive and
scale independent; improving it changes a finite logarithmic offset, not the
causal reuse slope.

## 4. Representation changes are not coefficient-loss causes

Common affine/Kelvin transport preserves the intrinsic material coordinate
exactly.  Nested refinement has zero Moyal charge.  Small symbol/covariance
representative changes belong to the already summable representation `Xi`.

A **physical** material-cell switch is different: it is a relink/backflow/fresh
causal event and stops the registration interval.  It must not be hidden by
calling it packet decoherence.

## 5. Scope

This theorem is exact once the moving selected parent-role coefficient equation
has been constructed.  That equation and its native `Q^2` energy interface are
now supplied by the outer-role and smooth quadratic-carrier theorems.  The
remaining frontier is global assembly/termination of the physical owners reached
after typed energy reentry.  No global-regularity claim is made.

### Simultaneous obstructions are not ordered lexicographically

A registration interval can meet more than one obstruction at the same first
time.  The registration theorem therefore returns the complete triggered set
`{HH coefficient obstruction, residual coefficient obstruction, material
relink}` and does **not** choose a primary by theorem name.  Coefficient faces
must first reenter actual physical energy/work; exact physical ties are handled
only after that typed routing.

## Outer-role update

The moving selected role equation assumed in the original scope is now supplied by `outer_moving_role_extraction.md` together with exact hard-event/smooth-envelope registration.  The native smooth-carrier energy handoff is supplied by `smooth_quadratic_carrier_interface.md`: common-slice coefficient failure locates energy reentry and never defines work by amplitude.  On the continuing subset, `physical_pair_weighted_productivity.md` uses the registered `1/4` factors under the actual physical transfer law.
