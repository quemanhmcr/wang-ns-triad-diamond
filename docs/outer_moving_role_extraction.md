# Exact outer moving-role extraction: low--high transport is algebra, not a packet hypothesis

At the stage when this theorem was introduced, the continuum frontier required an exact moving Fourier/helical role equation for
every recursively selected smooth-SGS block.  The main point of this note is that
the difficult-looking low--high decomposition is already an exact algebraic
identity of Navier--Stokes.  No Gaussian persistence theorem is needed to obtain
it.

## 1. Start from the actual resolved transporter

Write Leray Navier--Stokes as

\[
\partial_tu+\mathcal B(u,u)=\nu\Delta u,
\qquad
\mathcal B(a,b)=\mathbb P\nabla\!\cdot(a\otimes b).
\]

At a selected child scale `N` use the strict resolved field

\[
V=S_{N/4}u,
\qquad h=u-V,
\]

and define the exact resolved linearized transport

\[
\boxed{
\mathcal L_Vf=\mathcal B(V,f)+\mathcal B(f,V).
}
\]

Let `Q(t,D)` be a real scalar Fourier role multiplier, commuting with Leray and
spatial derivatives, and put

\[
w=Qu.
\]

Direct differentiation gives

\[
\begin{aligned}
(\partial_t+\mathcal L_V-\nu\Delta)w
={}&Q\mathcal B(V,V)-Q\mathcal B(h,h)\\
&+(\partial_tQ+[\mathcal L_V,Q])u,
\end{aligned}
\]

because the elementary quadratic identity is

\[
-\mathcal B(u,u)+\mathcal L_Vu
=\mathcal B(V,V)-\mathcal B(h,h).
\]

This is exact.  Pressure has already disappeared through Leray, and scalar
Fourier `Q` commutes with `Delta`, so viscosity creates no moving-role interface
term.

## 2. Low--low forcing is excluded for the whole low-strain slab

At the anchor event the signed-good selected role satisfies

\[
|\xi|\ge\frac35N.
\]

Transport its scalar symbol by the dual affine flow of the coherent averaged
jet.  If

\[
K=\int\|\operatorname{sym}\bar A\|_{op}dt,
\]

then singular-value control gives

\[
|\xi(t)|\ge e^{-K}|\xi(t_*)|.
\]

On the low-strain branch `K<=1/30`, therefore

\[
\boxed{
|\xi(t)|\ge\frac35e^{-1/30}N>\frac12N.
}
\]

On the other hand

\[
\operatorname{supp}\widehat V\subset B_{N/4}
\quad\Longrightarrow\quad
\operatorname{supp}\widehat{\mathcal B(V,V)}\subset B_{N/2}.
\]

Hence

\[
\boxed{Q(t,D)\mathcal B(V,V)=0}
\]

throughout the entire low-strain role interval.  The strict inequality is
Arb-certified.  This closes a support issue which had previously only been stated
at the frozen selection slice.

The exact moving selected-role equation is therefore

\[
\boxed{
(\partial_t+\mathcal L_V-\nu\Delta)w
=-Q\mathcal B(h,h)+R_Q,
}
\]

with

\[
\boxed{R_Q=(\partial_tQ+[\mathcal L_V,Q])u.}
\]

## 3. There is only one nonlinear generation term

The first term

\[
-Q\mathbb P\nabla\!\cdot(h\otimes h)
\]

is the complete quadratic non-resolved source seen by the selected role.  It need
not be split into an arbitrary list of packet forcings.  At each physical event,
the existing frozen relative-frequency/helicity cells disintegrate its actual
child-energy work.  The coherent resolution then gives the exact binary atoms

\[
W_{CDE}=2\Re\langle A_Ew_c,
\mathcal N(A_Cw_1,A_Dw_2)\rangle.
\]

Positive atoms are the physical causal law and negative atoms are backscatter.
Low-cost atoms are the ones on which complex Young / dual Gaussian / Bargmann
marking is applied.  Thus the high--high term in the equation above is already
exactly the source required by the physical-energy causal theorem.

## 4. Common affine motion creates zero role-interface forcing

For the scalar symbol choose

\[
\partial_tm-(\bar A^T\xi)\cdot\nabla_\xi m=0.
\]

Then the common affine advection part of the Heisenberg operator vanishes
exactly.  Constant affine stretching also commutes with scalar `Q`; Leray commutes
with it as a Fourier multiplier on the selected divergence-free role.

The remaining scalar advection commutator has the exact kernel form

\[
\boxed{
(\partial_tQ+[V\cdot\nabla,Q])f(x)
=
\int K_t(y)
[V(x)-V(x-y)-\bar A y]\cdot\nabla f(x-y)\,dy.
}
\]

The stretching commutator is likewise

\[
\boxed{
[(\nabla V),Q]f(x)
=
\int K_t(y)
[\nabla V(x)-\nabla V(x-y)]f(x-y)\,dy.
}
\]

Both vanish for genuinely affine resolved flow.  Thus `R_Q` is not an
approximation artifact.  It is exactly **non-affine resolved low--high
role-interface work**.

This is the operator-level version of the coherent deformation observable

\[
\mathcal K_C^2
=
\operatorname{Var}_\gamma
(L^{-1}\nabla V(X+Lz)L).
\]

## 5. Event anchoring preserves the original physical normalization

The role is not chosen by evolving some synthetic packet from the remote past.
At every recursively selected transfer event:

1. the existing physical smooth-SGS selector chooses the frozen transfer cell;
2. `Q(t_*)` is anchored to that cell;
3. only then is its scalar support transported over the short causal interval;
4. at the next common causal slice the field may be re-extracted again.

Therefore the physical transfer measure at the anchor event is unchanged.  The
moving role is a device for writing the exact PDE between two causal slices, not
a new transfer measure and not a persistent packet identity.

## 6. What is now closed

The following part of the old frontier is now exact:

\[
\boxed{
\text{selected frozen physical role}
\longrightarrow
\text{moving outer role PDE}
\longrightarrow
\text{one HH source + one non-affine interface}.
}
\]

In particular:

- strict low--low exclusion persists throughout the low-strain slab;
- common affine carrier/cell motion is exact gauge;
- pressure is absent;
- viscosity creates no cell-interface source;
- the complete nonlinear generation is one quadratic `h-h` term;
- no frozen Gaussian persistence is required.

## 7. Historical interface seam (closed at the smooth energy level downstream)

At the time of this theorem one term remained to be routed at **work level**:

\[
R_Q=(\partial_tQ+[\mathcal L_V,Q])u.
\]

Its affine part is already zero.  Its remainder is physical non-affine
low--high transfer across the moving outer-role interface.  The required downstream theorem had to
show an exact single-charge alternative:

\[
\boxed{
\text{small coherent deformation contribution}
\quad\lor\quad
\text{critical }D_V
\quad\lor\quad
\text{physical role relink/transfer loss}.
}
\]

The point is important: the commutator may **not** be hidden as representation
`Xi`, and it may not become a new currency merely because it is written with
brackets.  At the same time, for a non-idempotent smooth `Q`, its pairing with
`Q u` is not by itself the native carrier-interface work.

This note does not yet prove that final work-routing lemma, and therefore it does
not prove Navier--Stokes global regularity.

## Smooth-energy interface update

`smooth_quadratic_carrier_interface.md` now performs the missing type-correct
handoff.  The carrier energy is `<u,Q^2u>`, and the work of `R_Q` is recombined
with the diagonal `L_V` work already present in this outer equation.  The result
is the native `Q^2`-weighted interface law.  On a complete smooth square
partition its moving/skew part is conservative relink and its symmetric part is
the same resolved strain/deformation work.

`nonaffine_role_interface_work.md` remains the exact orthogonal hard-event
specialization, while `event_anchored_role_registration.md` supplies `QP=P` on
the event plateau.  Coefficient obstruction is only an interval locator; actual
ownership is assigned after physical-energy reentry.  The remaining task is
global recursive owner assembly, not construction of `Q` or invention of an
interface currency.
