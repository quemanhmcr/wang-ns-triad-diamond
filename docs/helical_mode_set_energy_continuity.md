# Helical mode-set energy continuity

Status: **candidate theorem block; certification requires dedicated, independent adversarial/actual-NS, and full causal-integration GitHub Actions gates.**

The cyclic donor theorem and its hard-cell single-charge quotient settle a **same-time** question: canonical negative nonlinear work is physical donor work whose positive recipient marginal is the already-canonical `dW+` law.  They deliberately do not manufacture a between-time matching of earlier deposits to later withdrawals.

The next native object is therefore not a cancellation budget.  It is the ordinary conservation law of energy stored on the physical Fourier--helical modes themselves.

Hard product cells do not appear as persistent state variables in this theorem.  They are interaction labels.  The persistent stock is modal energy.

---

## 1. Push the same-time donor law to physical helical-mode nodes

Let

\[
m=(k,s),\qquad s\in\{-1,+1\}
\]

be a physical helical mode identity.  The certified cyclic donor measure already carries donor and recipient child-mode identities.  Push it only to those nodes and write the resulting positive flow as

\[
\mathcal F_t(dm,dn).
\]

Its first marginal is gross canonical Hahn-negative child work and its second marginal is gross canonical Hahn-positive child work.  No new Hahn split is taken after modal summation.

For a finite/measurable mode set `A`, decompose the same-time flow into

\[
I_A(t)=\mathcal F_t(A\times A),
\]

\[
\Phi_{in,A}(t)=\mathcal F_t(A^c\times A),
\qquad
\Phi_{out,A}(t)=\mathcal F_t(A\times A^c).
\]

The gross positive and negative nonlinear works on modes in `A` are therefore

\[
\boxed{P_A=I_A+\Phi_{in,A}},
\qquad
\boxed{N_A=I_A+\Phi_{out,A}}.
\]

Subtracting gives the exact graph-divergence law

\[
\boxed{P_A-N_A=\Phi_{in,A}-\Phi_{out,A}.}
\]

The internal flow cancels because it is energy redistribution between simultaneous physical modes in the same set.  It is not deleted from the PDE.  It is simply divergence-free with respect to that set boundary.

---

## 2. Native helical modal stock law

Write the divergence-free Fourier coefficient in the repository's deterministic helical basis,

\[
\widehat u(k,t)=\sum_{s=\pm1}a_{k,s}(t)h_s(k),
\]

and define

\[
E_{k,s}(t)=|a_{k,s}(t)|^2.
\]

For one helical mode, Navier--Stokes gives

\[
\frac{d}{dt}E_{k,s}
=
W_{k,s}^{+}-W_{k,s}^{-}
-2\nu |k|^2 E_{k,s}.
\]

Here `W+` and `W-` are the **gross canonical edge Hahn marginals for that child helical mode**.  They are not obtained by first summing all edge work at the mode and then Hahn-splitting that modal net work.  This distinction is essential: gross positive and negative transfers may coexist even when their signed difference is small.

Integrating from `t0` to `t1` gives

\[
\boxed{
E_{k,s}(t_1)
+
\int_{t_0}^{t_1}W_{k,s}^{-}\,dt
+
2\nu |k|^2\int_{t_0}^{t_1}E_{k,s}\,dt
=
E_{k,s}(t_0)
+
\int_{t_0}^{t_1}W_{k,s}^{+}\,dt.
}
\]

This is an exact stock/flow identity.  It says what entered, what left, what remains, and what viscosity removed.  It does not say which earlier incoming atom funded which later outgoing atom.

---

## 3. Mode-set continuity law

Sum the modal stock identity over `A` and use the same-time donor graph divergence:

\[
E_A(t)=\sum_{m\in A}E_m(t),
\]

\[
D_A[t_0,t_1]
=
2\nu\int_{t_0}^{t_1}
\sum_{(k,s)\in A}|k|^2E_{k,s}(t)\,dt.
\]

Then

\[
\boxed{
E_A(t_1)
+D_A[t_0,t_1]
+\int_{t_0}^{t_1}\Phi_{out,A}(t)\,dt
=
E_A(t_0)
+\int_{t_0}^{t_1}\Phi_{in,A}(t)\,dt.
}
\]

This is the native between-time continuity equation supplied by the PDE plus the already-certified same-time cyclic donor law.

It has the familiar physical interpretation of a control volume, but the “volume” is a set of Fourier--helical mode nodes:

- stock at the final time;
- viscous dissipation inside the set;
- nonlinear energy carried outward across the mode-set boundary;
- initial stock;
- nonlinear energy carried inward across the boundary.

No scalar Bellman currency has been introduced.

---

## 4. Full-set recovery

For the full retained Galerkin mode set `A_all`, every cyclic donor/recipient atom is internal:

\[
\Phi_{in,A_{all}}=\Phi_{out,A_{all}}=0.
\]

The mode-set identity therefore reduces to

\[
\boxed{
E(t_1)+2\nu\int_{t_0}^{t_1}\|\nabla u(t)\|_2^2\,dt=E(t_0),
}
\]

up to the repository's fixed Fourier normalization.  Thus the new continuity law is not a competing energy theorem.  It refines the ordinary Navier--Stokes energy balance by resolving nonlinear transport across a chosen physical mode boundary.

---

## 5. Closed-triad anti-theorem: gross transfer is not a finite inventory resource

Take `A` to contain all three helical mode nodes of one nonzero closed triad.  Its cyclic donor flow is entirely internal:

\[
I_A>0,
\qquad
\Phi_{in,A}=\Phi_{out,A}=0,
\qquad
P_A=N_A=I_A.
\]

Now scale the actual physical closed-triad amplitudes by `lambda>0`,

\[
a_i\mapsto \lambda a_i.
\]

The trilinear Navier--Stokes work scales cubically,

\[
I_A\mapsto \lambda^3 I_A,
\]

while the boundary flux remains exactly zero.  The base measure is unchanged.

This is a physical anti-theorem against a tempting shortcut:

> modal energy stock does not bound total gross nonlinear transfer variation.

Conservative nonlinear energy may circulate or be recharged many times.  The continuity law controls **net boundary transport plus stock plus viscosity**, not the total amount of internal donor/recipient traffic.

Therefore neither

\[
\int W^+
\]

nor

\[
\int W^-
\]

may be declared a finite reset budget merely from energy conservation.

---

## 6. No FIFO/LIFO temporal provenance

The interval identity implies the aggregate inequality

\[
\int W^-_{k,s}
\le
E_{k,s}(t_0)+\int W^+_{k,s},
\]

but this is only a stock accounting consequence.  It does not canonically pair individual positive deposits with later withdrawals.

No physical law in the theorem selects:

- first-in-first-out;
- last-in-first-out;
- oldest-energy-first;
- newest-energy-first;
- a proportional temporal matching.

All such rules would be synthetic temporal couplings unless separately derived from Navier--Stokes itself.

The theorem therefore records aggregate continuity only.

---

## 7. Hard cells are not wallets across time

A hard product cell labels one nonlinear interaction configuration.  It is not a persistent state node.  In particular, a coarse hard-cell self-loop from the preceding theorem is a statement about same-time representation of donor and recipient interaction roots.  It cannot be carried forward as a stored “cell energy” account.

Persistent stock in the present theorem lives on

\[
(k,s)\mapsto E_{k,s}(t).
\]

Any later material/carrier theorem that wants a different persistent stock must prove its own physical transport law.

---

## 8. Reality covariance remains unquotiented

The real field obeys the repository gauge

\[
h_s(-k)=\overline{h_s(k)}.
\]

The Fourier coefficients at `k` and `-k` are conjugately related, but the canonical edge law has always retained both reality partners rather than quotienting them away.  The present singleton modal identity follows the same convention: one node is one `HelicalModeIdentity(k,s)`.  Reality covariance may be checked, but no extra factor of two is inserted into the node law by hand.

---

## 9. Actual Navier--Stokes audit

The companion PDE probe evolves the repository's real `2/3`-dealiased incompressible Fourier--Galerkin Navier--Stokes system.

At **every RK4 output time**, not merely at a few selected snapshots, it:

1. reads the actual evolving Fourier coefficient of one child wavevector;
2. decomposes that coefficient and the actual Leray-projected nonlinear source into the deterministic helical basis;
3. reconstructs the full unordered parent edge ledger for that child from the actual PDE state;
4. restricts the already-canonical edge atoms to one child helicity;
5. forms gross `W+` and `W-` before any modal summation;
6. checks that `W+ - W-` reconstructs the independently projected actual helical modal work;
7. trapezoid-integrates the gross work, helical modal stock, and viscous term over the same RK4 output grid;
8. checks the interval stock/work/viscosity balance with a tolerance owned by the finite RK4/trapezoid audit, not by the exact theorem object;
9. compares the same cutoff-7 Galerkin system represented on FFT grids `24` and `28`.

The exact theorem itself has no numerical tolerance.  In particular, a helicity component whose initial energy is exactly or nearly zero can have a finite-step relative quadrature residual that is much larger than machine epsilon even while the same PDE law converges correctly.  The probe therefore records the residual and enforces its own finite-step threshold on the native energy-throughput scale.

Cross-FFT representation checks use that same physical mode energy-throughput as an immutable numerical envelope.  They do **not** divide by the particular energy/work term being compared, because an opposite-helicity stock or a Hahn marginal can vanish by genuine physical selection/cancellation.  The throughput scale is used only to measure representation error; it is not a causal weight, recurrence currency, or gross-transfer budget.

The audit also runs the sign-reversed divergence-free initial state and the opposite child helicity.  These are real Navier--Stokes states/modes, not synthetic balance fixtures.

---

## 10. Scope after certification

The certified theorem closes the **between-time native energy-continuity registration** for cyclic donor provenance:

\[
\text{same-time donor flow on helical modes}
+
\text{modal NS stock/viscosity}
\Longrightarrow
\text{exact mode-set boundary continuity}.
\]

It would **not**:

- prove a finite bound on gross nonlinear transfer;
- give a canonical temporal matching of deposits and withdrawals;
- make hard interaction cells persistent energy stores;
- declare internal nonlinear flow dissipative;
- manufacture scale progress from an arbitrary mode set;
- solve degenerate Young/Christ margin;
- terminate generic mixed-owner recurrence;
- prove 3D Navier--Stokes global regularity.

The next question after this theorem should be asked only after seeing what additional rigidity real choices of the mode boundary supply.  In particular, a radial high-frequency set may turn the abstract boundary flux into genuine spectral crossing, but that scale statement must be proved from its own Fourier geometry rather than assumed from the continuity equation.


---

## 11. Certification evidence

Exact certified theorem/provenance SHA:

`a39d502d9312ac3bd6613a780d60b22e88790863`.

Dedicated GitHub Actions run `31598313644` completed successfully.  The full theorem suite passed `870` tests.  The analytic lane ran `75,000` physical closed triads and checked `150,000` physical mode sets, with worst boundary-divergence native residual `2.138e-16`.  The largest sampled internal/boundary-flow ratio with nonzero boundary flow was `38904.2004084`.  On the physical amplitude-scaling anti-theorem the closed signed-good triad's internal flow increased from `0.712779828831` to `712.779828831` under a tenfold amplitude scaling, exactly the expected cubic factor, while inward and outward boundary flux of the full three-mode set remained zero.

The same run evolved actual `2/3`-dealiased Fourier--Galerkin Navier--Stokes at common cutoff `7` on FFT grids `24` and `28`.  For the selected `+` helicity, the `64`-step interval continuity residual was `1.563e-11`; under global phase/sign reversal it was `2.960e-11`.  The deliberately difficult opposite-helicity mode began at energy about `1e-32`, was nonlinearly born to final energy `1.80191453799e-05`, received integrated gross positive/negative work `1.81542566415e-05` / `9.53908587922e-08`, and had viscous dissipation `3.97413682063e-08`.  Its interval continuity native residual was `1.155e-06`, while the worst 24/28 representation native residual was `4.890e-14`.  Amplitude adversaries `0.5` and `2.0` also passed without changing the physical law.

Independent audit `31598313628` also completed successfully.  Its algebra lane ran `100,000` physical closed triads and `200,000` mode sets, with worst boundary-divergence native residual `2.141e-16`; the companion cyclic donor stress also passed `100,000` states.  Its longer actual-NS lane used amplitude `1.3`, `96` RK4 steps, duration `0.0012`, common cutoff `7`, and FFT representations `24,28`.  The selected/sign-reversed interval residuals were `4.051e-11` / `6.485e-11`.  The longer opposite-helicity branch began at `6.91e-32` / `2.27e-32` on the two representations and ended at `7.38660111133e-05`; integrated positive work was `7.44536161723e-05`, negative work `3.92012387067e-07`, viscosity `1.9564433437e-07`, interval residual `6.939e-07`, and worst cross-FFT native residual `6.444e-14`.

Full physical-energy causal integration `31598313576` completed successfully on the same exact SHA, passed the same `870`-test suite, and ran `50,000` mode-continuity stress states immediately after cyclic donor + hard-cell single-charge provenance and before complex Young.  Its worst boundary-divergence native residual was `2.103e-16`; the entire downstream witness, role, donor-quotient, causal-reuse, branch-compiler and master spine remained green.

Stored artifact trees:

- `recorded-results/31598313644/`, deterministic tree digest `sha256:9b4bff2e5ce26ac1afb4afb7ebe8f25f4bd4f5c9a754b52e1d51812550d3cfae`;
- `recorded-results/31598313628/`, deterministic tree digest `sha256:7f34ef66faef1e4d132ec59622c10ab04c1eb8d5086318645a28828f83b22f74`;
- `recorded-results/31598313576/`, deterministic tree digest `sha256:6eb24c013d52c54ef63e4196e5e8df565756859bcae2271f4b014e796c587da7`.

Failure lineage is part of the certificate.  Initial candidate `8ec1efe924c1a7cabf2499ff232e09f060ced2d8` encoded a universal `5e-8` numerical residual gate inside the exact interval object; the opposite-helicity finite-step audit correctly rejected that category error.  Candidate `e11b3464d437b42e0474da1e016e67a584bcfe92` removed the universal theorem-object tolerance, but its cross-FFT audit still divided each representation gap by the particular realized observable.  Because an opposite-helicity stock/Hahn marginal may vanish by genuine selection or phase cancellation, the 24/28 audit again failed.  Final `a39d502...` measures cross-representation gaps on the already-native physical mode energy-throughput envelope.  That envelope is recorded only as a numerical error scale; it never enters the donor kernel, causal probability, recurrence currency or gross-transfer law.  No PDE case, sample count, sign pattern, or physical threshold was removed.

The theorem closes native **mode-set continuity**, not recurrence.  Its next physically meaningful specialization is a radial Fourier mode boundary, where inward/outward mode-set flow becomes genuine upward/downward spectral crossing.  Any scale direction, layer-cake identity, or recurrence consequence must be derived from that radial geometry itself; it is not supplied by an arbitrary mode-set continuity equation.  No Navier--Stokes global-regularity claim is made.
