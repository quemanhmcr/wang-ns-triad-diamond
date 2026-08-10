# Amplitude--entropy causal reuse: the root budget should respect trilinear homogeneity

The old synchronized causal theorem assumed that every distinct root cell carried
one fixed critical mass `N E >= eta`.  That is too strong for a theorem derived
only from Young near-extremality.  The trilinear functional is separately
homogeneous in its three roles, so shape rigidity cannot impose an absolute
amplitude floor.

The correct invariant is multiplicative.

## 1. Critical selected coefficient

For a selected coherent/packet probe `phi` at scale `N`, define

\[
\alpha=\sqrt N\,|\langle u,\phi\rangle|.
\]

This is invariant under the isotropic Navier--Stokes scaling.

For the L3-normalized Gaussian dual `h` associated with an `L^(3/2)` Gaussian
profile of physical geometric radius `r_g`, exact Gaussian integration gives

\[
\frac{\|h\|_{3/2}}{\|h\|_2}
=\frac43\pi^{1/4}r_g^{-1/2}.
\]

After covariance quantization by log-radius `delta` and the shell uncertainty
`N r_g>=2/3`,

\[
N^{-1/2}\|\phi\|_{3/2}
\le
\frac43\pi^{1/4}
\left[\frac23e^{-\sqrt3\delta/6}\right]^{-1/2}.
\]

Thus the packet test has exactly the `sqrt(N)` critical norm required to cancel
one derivative and one parabolic lifetime.

## 2. Duhamel is an amplitude diagnostic, not the causal weight

For the exact symmetrized high--high source

\[
F_{HH}=-Q_c\mathbb P\nabla\cdot
(w_1\otimes w_2+w_2\otimes w_1),
\]

sharp Young and `|k|<=R_kN` give

\[
|\langle\phi,F_{HH}\rangle|
\le
2R_kN A_3\|\phi\|_{3/2}
\|\widehat w_1\|_{3/2}\|\widehat w_2\|_{3/2}.
\]

The coefficient-space backward adjoint has

\[
\|\psi(t)\|
\le
\exp(K+\nu cR_k^2)
\]

on `T=cN^-2`.  Hence all powers of `N` cancel.  On the existing Duhamel-generated
branch `|I_HH|>=|c_1|/2`,

\[
\boxed{
\|\widehat w_1\|_{3/2}\|\widehat w_2\|_{3/2}
\ge \Lambda_{L^{3/2}}\alpha_c,
}
\]

where `Lambda_(L3/2)>0` depends only on the fixed clean constants, `c`, and `nu`,
not on `N`.

This use of Duhamel is deliberately different from the rejected
`dGamma=dT` identification.  Duhamel certifies **amplitude productivity and
same-time support**.  The causal probability weights remain the actual positive
child-energy work `dT_HH`.

## 3. Register the parents

The dual-Gaussian theorem gives, for a parent role with critical Young amplitude
`a=||what||_(3/2)`,

\[
\alpha_{parent,event}^2\ge\eta_{dual}a^2.
\]

To move that event mark to the synchronized common slice use the existing adjoint
first-stop rule.  If no earlier causal stop occurs, the inherited coefficient
keeps at least the clean `1/4` fraction.  If it does not, the first failure is
already an HH-generation, classified-residual/source, or genuine material
relink stop and the branch terminates there.

Therefore every **continuing** node obeys

\[
\boxed{
\alpha_{p_1}\alpha_{p_2}\ge\Lambda\alpha_c
}
\]

for one fixed `Lambda>0`.

The present module proves everything after this registered productivity law and
records the explicit candidate constant.  Proving the common-slice registration
for the full PDE selector is still the continuum bridge.

## 4. Log amplitude is the natural binary observable

Give each child event its two structural parent-role slots with probability
`1/2`, exactly as in the existing Shannon/Renyi causal law.  For one event,

\[
\frac12(\log\alpha_{p_1}+\log\alpha_{p_2})
\ge
\frac12\log\Lambda+rac12\log\alpha_c.
\]

After transfer-weighted averaging,

\[
\ell_j:=\mathbb E_{w_j}\log\alpha
\ge
\frac12\log\Lambda+rac12\ell_{j+1}.
\]

Thus

\[
\boxed{
\ell_0
\ge
(1-2^{-L})\log\Lambda
+2^{-L}\log\alpha_L.
}
\]

This is why amplitude imbalance is not a new escape.  If one parent is tiny and
the other enormous, their **average logarithm** still sees their exact product.
No threshold separating "small" and "large" parents is needed.

If several slots merge to one material label, assign that label the maximum
coefficient among its incoming slots.  This only increases `ell_j`; merging can
never spoil the inequality.

## 5. Energy controls entropy, not the number of roots

For a coherent coefficient `alpha_r`, Bargmann submean on the radius-`sqrt(3)`
ball and the canonical Moyal partition give

\[
\boxed{N_rE_r\ge\beta\alpha_r^2}
\]

with one fixed `beta>0`.  Distinct root anchors have a common-slice energy
budget

\[
\sum_rE_r\le P E_{global}.
\]

Since signed-good synchronization gives

\[
N_r\le N_{base}(25/24)^L,
\]

\[
\sum_r\alpha_r^2
\le
\frac{P E_{global}N_{base}}{\beta}(25/24)^L.
\]

Now let `w_r` be the **physical transfer-weighted root law**.  The log-sum
inequality gives exactly

\[
\boxed{
H(w_0)+2\mathbb E_{w_0}\log\alpha_r
\le
\log\sum_r\alpha_r^2.
}
\]

Consequently

\[
\boxed{
H(w_0)
\le
\log\frac{P E_{global}N_{base}(25/24)^L}{\beta}
-2\ell_0.
}
\]

This replaces `H(root)<=log n_0` plus a uniform per-root mass floor.  It is
strictly more compatible with the homogeneity of the PDE.

## 6. Shannon and Renyi telescopes keep the same linear slope

For one terminal causal atom,

\[
\sum_j\mathcal R_j=L\log2-H(w_0).
\]

Combining the preceding bounds,

\[
\boxed{
\sum_j\mathcal R_j
\ge
L\log\frac{48}{25}
-\log\frac{P E_{global}N_{base}}{\beta}
+2(1-2^{-L})\log\Lambda
+2^{1-L}\log\alpha_L.
}
\]

The amplitude terms are finite-depth offsets.  The coefficient of `L` is still
`log(48/25)`.

For Renyi-2,

\[
H_2(w_0)\le H(w_0),
\]

so

\[
L\log2+\log Q_0
\ge L\log2-H(w_0).
\]

The **same lower bound** therefore applies to the exact Renyi action

\[
\sum_j\log(1+\theta_j).
\]

Once it exceeds `L log(4/3)`, some layer has `theta_j>1/3`, and the existing
transfer-weighted parent-slot pair / component-entropy / same-ancestry-cycle
routing applies without change.

## 7. Physical interpretation

The old uniform-root-mass hypothesis discretized a scale-invariant trilinear law
too early.  The replacement is closer to the PDE:

\[
\text{quadratic productivity}
\longrightarrow
\text{log-amplitude transport}
\longrightarrow
\text{energy--entropy inequality}
\longrightarrow
\text{reuse/collision}.
\]

A highly unbalanced interaction does not need to be called a new currency.  Its
small and large amplitudes are conjugate through the product, and logarithms are
the natural additive coordinate of that multiplicative physics.

At this stage the remaining continuum task was the registration step stated above: turn
the eventwise coefficient into the coefficient on the synchronized common slice,
or stop at the first physical cause preventing that registration.  That interface is now supplied by the common-slice, moving-role, event-role and first-hit theorems.  No global
regularity conclusion is asserted.

## 8. Companion registration theorem closes the selected-role time interface

`common_slice_coefficient_registration.md` now proves the registration premise
used above at the selected-role model level.  The asynchronous cone keeps every
parent event inside its natural adjoint interval, and the exact first-stop gate
returns either an earlier physical obstruction or at least `1/4` of the event
coefficient on the common slice.  Hence the `1/16` product factor in `Lambda` is
no longer a persistence hypothesis.

That historical outer-role bridge is now closed by `outer_moving_role_extraction.md` together with the nonaffine-interface and event-role companion theorems.  The amplitude--entropy telescope is therefore a supplied causal-reuse module; the remaining burden lies in final continuum assembly, not construction of a persistent selected packet.  No global-regularity conclusion is asserted.

## Physical-weighted productivity update

The pointwise Duhamel-product discussion above is retained as provenance but is stronger than the causal telescope needs.  The preferred input is now `physical_pair_weighted_productivity.md`, which proves directly under actual positive child-energy work

`E_dT log(alpha_p1 alpha_p2) >= E_dT log(alpha_child) + log Lambda_j`.

For `M_j` retained hard pair cells, `Lambda_j=Lambda_1/M_j`; polynomial symbol refinement contributes only the finite geometrically discounted offset `sum 2^(-j-1) log M_j`.  The exact depth formula is `ell_0 >= sum_j 2^(-j-1) log Lambda_j + 2^(-L) log alpha_L`, and the Shannon/Renyi lower therefore receives `sum_j 2^(-j) log Lambda_j`.  Thus no Duhamel-to-physical parent-pair identification is part of the preferred theorem.  Outer roles and event registration are supplied by their companion exact theorems.
