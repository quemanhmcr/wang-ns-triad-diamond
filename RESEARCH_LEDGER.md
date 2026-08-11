# Research Ledger — current canonical architecture for the 3D Navier–Stokes rigidity programme
> **Purpose.** This file is the compact current-state map of the programme.
> It is organized by the physical structure of the Navier–Stokes dynamics, not by discovery date.
> The full chronological research history, including superseded routes, failed fixtures, and detailed CI provenance, is frozen in [`docs/history/RESEARCH_LEDGER_history_through_2026-08-10.md`](docs/history/RESEARCH_LEDGER_history_through_2026-08-10.md).
> The present ledger is intended to let another mathematical physicist understand the proof architecture quickly and identify exactly what is proved, what is conditional, and what remains open.

---

## 0. Executive state
We study the 3D incompressible Navier–Stokes system
\[ \partial_tu+\mathbb P\nabla\cdot(u\otimes u)=\nu\Delta u, \qquad \nabla\cdot u=0, \]
through scale-critical physical observables, exact Fourier/coherent decompositions, material transport, and causal positive-work laws.
The current architecture has reached the following point.
1. The physical child-energy cause is actual positive nonlinear work, not a normalized Duhamel amplitude.
2. Smooth moving roles and event-anchored hard roles have an exact PDE interface.
3. Objective source terms have been routed to native physical owners: resolved dissipation, coherent SGS service, hard pressure pairs, or viscosity.
4. Generic critical hard-shell mass has a certified own-scale first-stop/service reentry theorem.
5. High strain is converted to critical lower-frequency resolved ancestors and then to own-scale service.
6. Fresh SGS service is quotiented by coherent-cell refinement and converted to an actual hard-shell seed through a canonical LP scale law.
7. High-frequency coherent service is converted to physical hard-tail energy and then to inherited critical shells or actual regeneration work.
8. High-tail regeneration now has a complete physical continuation chain:
   - common-unit causal HH/interface ownership;
   - resolved interface quotient to same-event conservative donor provenance or existing strain/deformation;
   - exact binary HH work atomization;
   - ultraviolet locality paid by physical tail dissipation;
   - sliding natural-time concentration;
   - actual critical hard-shell reentry with forward scale progress.
9. Material labels are a sidecar quotient: pure label/family changes do not create a second carrier impulse.
10. Exact ties are retained as joint physical cause sets; no lexicographic theorem priority is canonical.
11. No critical energy or normalized dissipation event is treated as a scale-independent additive reset.
12. The continuum master now separates physical time from event topology: relay depth is quotiented, bounded-scale no-hit checkpoint continuation hits `t=0` by real corridor time, and full-natural horizon endpoints do not become event vertices merely because one natural window ended. Global termination remains open.
13. The resolved cross/interface branch is no longer an independent recursive generator: for the actual resolved low--high operator `L_V=K+S`, skew `K` is finite same-event conservative donor flux while symmetric `S` is existing strain/deformation ownership. Pure interface circulation is quotiented before recursion.
14. The propagated smooth-carrier interface is now type-correct at energy and ownership level: `Q^2` is the carrier energy effect; common affine/Kelvin role motion is quotiented by `dot A+[G,A]=0`; only residual skew `K_phys` is physical relink; raw HH/interface coefficient hits remain first-stop locators until actual energy/work reentry.
15. A completed full-natural critical-shell corridor carries its own-scale service as a same-interval physical witness; OO/ON/NN rereading adds zero depth, and the endpoint carrier has a same-time two-shell cover at ratios `3/4` or `3/2`, with the actual shell masses deciding the unique/joint witness set.
16. A complete no-hit natural horizon is now a certified analysis checkpoint, not a recursive physical event. The corridor time is real, but checkpoint cover geometry supplies neither a causal charge nor high-tail/directional scale provenance. Infinite genuine-event recurrence and event-free UV checkpoint continuation are therefore distinct remaining problems.
There is **no claim here of a proof of global regularity for 3D Navier–Stokes**.

---

## 1. Physical normalization and canonical scale variables
Fix a parent block frequency `N>0`.
Its natural parabolic time is
\[ T_N=cN^{-2}, \]
where `c>0` is a dimensionless lifetime coefficient.
For a hard shell
\[ A_M=\{M/2<|\xi|\le M\}, \]
write
\[ E_M(t)=\|P_Mu(t)\|_2^2, \qquad \mu_M(t)=M E_M(t). \]
`mu_M` is the basic scale-critical hard-shell mass.
For a block-scale resolved transporter
\[ V=S_{N/4}u, \qquad h=u-V, \]
use normalized resolved dissipation
\[ D_V=N\int_I\|\nabla V\|_2^2dt. \]
For a physical hard tail
\[ w=P_{>N}u, \]
use
\[ D_{tail}=N\int_I\|\nabla w\|_2^2dt. \]
The conserved/global physical kinetic-energy scale is
\[ E_{global}=\sup_t\|u(t)\|_2^2. \]
The scale-critical version at block scale `N` is
\[ \mathcal E_N=N E_{global}. \]
This quantity may appear as a denominator in scale-critical capacity inequalities.
It is **not** a finite reset budget.
A critical event at frequency `M` costs physical energy on the order `mu/M`, so geometrically increasing frequencies can have summable physical-energy cost.
Likewise normalized dissipations such as `D_V` have physical viscous cost `nu D_V/N` and are not scale-independent reset currencies.
The repository therefore distinguishes carefully between:
- scale-critical event size;
- physical energy/dissipation cost;
- causal probability law;
- deterministic concentration coordinate;
- representation/interface error;
- genuinely globally bounded scale-independent resources.
These objects are never interchangeable merely because they have comparable numerical size.

---

## 2. The governing methodological principle
The architecture is built around one rule:
> **Quotient observer freedom before charging physics.**
A mathematical representation may introduce cells, frames, cutoffs, covariance charts, packet labels, shell orientations, or time coordinates.
Such choices are not physical owners unless changing them changes an actual Navier–Stokes observable.
The current canonical quotients are:
### 2.1 Common affine/Kelvin motion is free
Common affine motion is a gauge.
For affine flow `Ldot=A L` and Kelvin covector `kdot=-A^T k`, the intrinsic material coordinate
\[ \zeta=(L^{-1}X/2,L^Tk) \]
is invariant.
No charge is assigned merely because the observer follows the common affine/Kelvin frame.
### 2.2 Resolved cutoff change is a repartition gauge
For any resolved field `V`, `h=u-V`, and smooth scalar role `Q`, the cutoff-dependent nonlinear expression
\[ G_V=-L_V(Qu)+Q B(V,V)-Q B(h,h)+(L_VQ-QL_V)u \]
satisfies exactly
\[ G_V=-Q B(u,u). \]
Changing from `S_{N/4}u` to a renewed cutoff `S_{M/4}u` does not create a new source.
It repartitions the same nonlinearity between transport, HH source, and interface work.
### 2.3 Hard interaction roles are event-anchored; hard shell state may also be reread at checkpoints
Hard Fourier/helical projectors that resolve actual nonlinear interaction roles are used at physical events and are not propagated as packets. A hard shell energy observable may also be reread at a no-hit analysis checkpoint; that state reading does not make the checkpoint an event. Smooth scalar PDE envelopes are the objects propagated through the intervening corridors.
If `P` is the event hard projector and `Q` is smooth with `QP=P`, then
\[ \langle Pu,\phi\rangle=\langle Qu,P\phi\rangle. \]
Hard projection does not need to be differentiated in time.
The propagated smooth carrier is not a hard projector. Its physical energy is
\[ \boxed{\|Qu\|_2^2=\langle u,Q^2u\rangle.} \]
A complete smooth analysis may use a square partition
\[ \sum_a A_a^2=I. \]
That identity is analysis bookkeeping, not by itself a physical transfer law.
Before any moving-window exchange is charged, the common affine/Kelvin transport generator `G` must be identified and the roles must satisfy
\[ \boxed{\dot A_a+[G,A_a]=0\quad\text{for every }a.} \]
This common motion is observer gauge and is quotiented. Only residual skew work after that quotient may be called physical relink.
### 2.4 Material labels are sidecars
Pure old/new label changes do not create a second carrier impulse when the smooth carrier `Q` and dual `psi` are unchanged.
The coefficient identity remains
\[ z(t)=z(s)+I_{HH}+I_{interface}. \]
OO/ON/NN is reread from actual positive service whenever needed.
A selected-family change is controlled by its symmetric-difference physical energy.
### 2.5 Exact ties are joint causes
If multiple physical causes hit at the same first time, the preferred master stores the unsplit joint cause set.
There is no theorem-name ordering.
Fine Radon–Nikodym splitting is optional bookkeeping only.
### 2.6 `t=0` is absorbing
Backward recursion stops exactly at the initial surface.
The initial boundary is not a fresh interior event.
### 2.7 Only true globally bounded resources may reset additively
A normalized critical event such as `D_V~1` or `NE~1` is not enough.
Its physical cost can decay like `1/N` along geometric high-frequency chains.
An additive reset is permitted only when the resource is genuinely scale-independent and globally bounded in the physical problem.
### 2.8 Natural horizons are checkpoints unless physics fires
A theorem may inspect one native parabolic interval without creating an event at its endpoint.
For a generic shell at incoming frequency `M`, the actual no-hit corridor runs at `A=3M/4` for time `cA^-2`. If no physical first stop occurs and `t=0` is not reached, that elapsed time is real Navier--Stokes evolution but the earlier endpoint is only an analysis checkpoint.
Hard-shell rereading there is state observation. The actual pair `(mu_A,mu_2A)` decides which endpoint shell is realized; the observer may not choose the `3M/2` cover branch to manufacture UV progress.

---

## 3. Canonical observables and their meanings
The following objects play distinct roles.
### 3.1 Physical positive child-energy work
For a carrier/coefficient `c` and physical source `F`, the positive work element is
\[ dT=2[\operatorname{Re}\langle c,F\rangle]_+dt. \]
This is the canonical causal weight when the event represents actual child-energy generation.
### 3.2 Signed nonlinear work
Before taking positive parts, keep the signed work.
Backscatter is physically real.
Positive/negative Hahn splitting is applied only after exact signed reconstruction.
### 3.3 Duhamel amplitudes
Duhamel identities provide exact support and adjoint transport information.
They are **not** causal probability laws by themselves.
The old identification `dGamma=dT` is false.
### 3.4 Critical shell mass
\[ \mu_M=M\|P_Mu\|_2^2. \]
This is a recursive shell-state strength.
It is not a globally finite event count.
### 3.5 Coherent increment service
Square increment/heat service is a positive physical measure useful for SGS routing, high-strain routing, and material ownership.
Its positive Moyal disintegration may be read by shell, material endpoint, or LP band.
### 3.6 Deterministic concentration coordinates
Examples:
\[ H_\infty^{scale}=-\log p_{max}, \qquad H_2=-\log\sum_jp_j^2. \]
These are logarithmic coordinates of a specified positive physical measure.
They are not automatically causal child-energy entropy.
The measure must always be named.
### 3.7 Causal Shannon/Rényi information
Shannon/Rényi reuse applies only to a genuine causal law built from actual physical child-work.
Do not substitute a source/service concentration law unless a theorem explicitly identifies the two.
### 3.8 `Xi`
`Xi` records controlled representation/interface losses that are not native physical generation.
Common affine motion and exact cutoff repartition do not belong to `Xi`.
Genuine covariance mismatch, symbol freezing, selected-family switch, or nontrivial role interface may contribute only through certified interfaces.
### 3.9 Physical owner sets
A master event may simultaneously own several native causes.
The code should return the realized set, not the first true branch in source order.

---

## 4. Physical energy causal bridge
The first decisive correction to the old ancestry picture is that normalized Duhamel mass is not physical child-energy work.
Let the carried energy satisfy an energy inequality of the form
\[ E_1\le e^{2K}(E_0+W_{HH}^++W_R^+). \]
On the low-strain branch
\[ K\le1/30, \]
with
\[ E_0<E_1/5, \qquad W_R^+<E_1/5, \]
one obtains
\[ \boxed{W_{HH}^+\ge\frac{8}{15}E_1.} \]
This theorem has one precise purpose:
- it compares actual positive HH work with terminal child energy.
It is **not** required merely to define a causal HH work law.
This distinction matters later in high-tail regeneration.
A positive HH law exists before the `8/15` generated-energy gate.
The gate is needed only when the subsequent theorem requires productivity relative to child energy.

---

## 5. Exact coherent binary work atomization
Let `A_C` be an exact coherent analysis/synthesis partition resolving identity.
Bilinearity gives
\[ \mathcal N(w_1,w_2)=\sum_{C,D}\mathcal N(A_Cw_1,A_Dw_2). \]
Resolving the child as well gives signed atoms
\[ W_{CDE} =2\operatorname{Re}\langle A_Ew_3, \mathcal N(A_Cw_1,A_Dw_2)\rangle. \]
Exactly,
\[ \sum_{C,D,E}W_{CDE}=W_{HH}. \]
Define
\[ P=\sum[W_{CDE}]_+, \qquad N=\sum[-W_{CDE}]_+. \]
Then
\[ P-N=W_{HH}, \qquad P\ge[W_{HH}]_+. \]
Thus every positive HH event admits an exact binary physical parent-work law.
This atomization does not require:
- Duhamel pair weights;
- packet persistence;
- Young near-extremality;
- a child-energy gate;
- a coherent cell mass floor.
The coherent cells answer **which parent/child roles participate in work**.
A separate energetic material anchor may later answer **which material reservoir carries a parent role**.
Do not identify those two questions.

---

## 6. Smooth moving carriers and the exact outer-role PDE
Set
\[ V=S_{N/4}u, \qquad h=u-V, \]
and
\[ L_Vf=\mathbb P\nabla\cdot(V\otimes f+f\otimes V). \]
For a time-dependent scalar divergence-free Fourier role
\[ w=Q(t,D)u, \]
one has exactly
\[ (\partial_t+L_V-\nu\Delta)w =Q B(V,V)-Q B(h,h)+(\partial_tQ+[L_V,Q])u. \]
On the low-strain affine-dual-transport branch, the moving support of `Q` stays above `N/2`.
Since
\[ \operatorname{supp}\widehat{B(V,V)}\subset B_{N/2}, \]
one gets persistent low-low exclusion:
\[ Q B(V,V)=0. \]
Therefore
\[ \boxed{ (\partial_t+L_V-\nu\Delta)w =-Q\mathbb P\nabla\cdot(h\otimes h)+R_Q, } \]
with
\[ R_Q=(\partial_tQ+[L_V,Q])u. \]
The first term is the genuine HH quadratic source.
The second is an outer-equation forcing term; its work is not the commutator pairing by itself.
For the smooth carrier, with `eta=Q^2`, the native carrier-energy identity is
\[ \boxed{\frac d{dt}\|Qu\|_2^2+2\nu\|\nabla Qu\|_2^2
=\langle u,\partial_t(Q^2)u\rangle-2\operatorname{Re}\langle Q^2u,B(u,u)\rangle.} \]
The native resolved interface work is
\[ J_Q=\langle u,\partial_t(Q^2)u\rangle-2\operatorname{Re}\langle Q^2u,L_Vu\rangle, \]
and the outer-role equation reproduces it exactly only after diagonal resolved work is restored:
\[ \boxed{J_Q=2\operatorname{Re}\langle Qu,R_Q\rangle-2\operatorname{Re}\langle Qu,L_VQu\rangle.} \]
There is no separate pressure forcing in this Leray-projected role equation.
## 7. Exact Egorov structure and coherent deformation
For convolution kernel `K_N`, after subtracting the affine jet `A=\nabla V(X)`, the exact scalar commutator is
\[ (\dot Q+[V\cdot\nabla,Q])f(x) = \int K_N(y) [V(x)-V(x-y)-Ay]\cdot\nabla f(x-y)dy. \]
A center Hessian sampled only at `X` does **not** control all nonaffine effects.
There are smooth strict-lowpass divergence-free counterexamples with
\[ V(0)=0, \quad \nabla V(0)=0, \quad \nabla^2V(0)=0, \quad \nabla^3V(0)\ne0. \]
The correct affine-invariant observable is coherent deformation variance.
With Gaussian coordinate `z`, define
\[ F(z)=L^{-1}\nabla V(X+Lz)L, \qquad \bar A=\mathbb E_\gamma F, \]
and
\[ K_{coh}^2 = \mathbb E_\gamma\|F-\bar A\|^2. \]
Common affine flow gives `K_coh=0`.
Gaussian Poincaré controls it by averaged Hessian information across the coherent eddy.
The certified Hermite residual bounds include
\[ \mathbb E|R|^2\le K_{coh}^2/2, \qquad \mathbb E|z|^2|R|^2\le7K_{coh}^2. \]
The nonaffine interface therefore has physical coherent-deformation provenance.
It is not arbitrary packet error.

---

## 8. Nonaffine role-interface work
Write the resolved linear role operator as
\[ L_V=K+S, \qquad K^*=-K, \quad S^*=S. \]
For a complete orthogonal **hard event-role** partition, the `K` pair work is antisymmetric conservative same-event flux and the `S` pair work is symmetric strain/deformation. This is the hard donor/circulation theorem.
For the propagated non-idempotent smooth carrier, do not import hard-projector algebra. Complete the carrier by `sum_a A_a^2=I`, then first quotient the actual common role transport:
\[ \dot A_a+[G,A_a]=0. \]
After this has been verified, write
\[ K=G+K_{phys}. \]
The moving partition work cancels the common `G` work exactly. Only `K_phys` remains as antisymmetric conservative **physical** relink, while `S` remains the existing strain/deformation work.
Thus conservation of channel exchange is not enough to make an observer-selected window motion physical. A role change caused only by observer gauge or cutoff repartition is free; a genuine residual physical relink is not.
The hard event-role and smooth propagated measures share `K/S` provenance but are not identified as the same measure without an explicit physical-work pushforward.
## 9. Coherent averaged strain source and service-or-flat gate
The affine/Kelvin carrier uses a coherent averaged resolved gradient.
The exact moving-average source equation has the form
\[ \dot{\bar A} = -\bar A^2 -\langle a^2\rangle -\langle\nabla^2P\rangle -\langle\nabla\operatorname{div}R_{SGS}\rangle +\nu\langle\Delta A\rangle -\langle r\cdot\nabla A\rangle. \]
The final term is handled exactly by integration by parts.
The coherent nonaffine error obeys a resolved dissipation bound
\[ I_K^2\le0.275568824559\,cD_V. \]
For the canonical `q_max=4.71207563594`, the coherent service functional is
\[ F_{coh} = \sqrt{E_H^{phys}} +(dT)_{nonconf} +C(q_{max})I_K, \]
with
\[ C(q_{max})=5.65481629117. \]
At `tau=0.01`, sufficiently small `I_K` makes the nonaffine action less than one third of the service tolerance; otherwise the event routes to `D_V`.
The gate returns **all** realized physical roots.
It does not choose a lexicographic primary theorem.

---

## 10. Generic critical hard-shell reentry
This theorem is the central recursive landing zone for many supplier mechanisms.
Assume at an event time
\[ M\|P_Mu(t)\|_2^2\ge\mu_0>0. \]
Set
\[ A=\frac34M. \]
Choose a smooth scalar envelope `Q_A` equal to one on the whole hard shell.
With canonical dual
\[ \psi=P_Mu/\|P_Mu\|_2, \]
one has terminal coefficient mass
\[ A|z(t)|^2\ge\frac34\mu_0. \]
Run the backward natural first-stop corridor with only native monitors:
- renewed strain;
- nonaffine/interface coefficient obstruction;
- HH coefficient obstruction.
No material boundary monitor is needed in the generic theorem.
If a named stop occurs, route to its owner.
If `t=0` is reached, stop at the initial boundary.
If the full natural corridor survives, then every prefix satisfies
\[ A|z(s)|^2\ge\frac{3\mu_0}{64}. \]
Including registered affine/Kelvin/viscous dual growth `J`,
\[ A\|Q_Au(s)\|_2^2 \ge \frac{3\mu_0}{64J^2}. \]
A bounded heat-defect fraction then gives a uniform own-scale increment-service lower.
The theorem records service only on a full no-hit survivor.
It does not turn a shell seed into unconditional persistent service.
Material identity is deferred until actual renewed service is reread.

---

## 11. High-strain route: from resolved strain to critical service
For the strict transporter `V=S_{N/4}u`, Fourier Cauchy–Schwarz gives a scale-correct estimate of `||grad V||_infty`.
On `T=cN^-2`,
\[ K:=\int\|S_V\|_{op}dt \le \frac{\sqrt{cD_V}}{8\sqrt6\pi}. \]
Hence
\[ D_V\ge\frac{384\pi^2}{c}K^2. \]
At the low-strain boundary
\[ K=1/30, \]
one obtains
\[ \boxed{D_V\ge\frac{32\pi^2}{75c}.} \]
Now decompose `V` into deterministic dyadic resolved annuli
\[ M_j=(N/4)2^{-j}. \]
Let
\[ \mu_j=M_j\|P_ju\|_2^2. \]
Choose
\[ \mu_*=\frac{32\pi^2}{75c^2}. \]
The low-mass portion of actual normalized resolved dissipation is at most half the high-strain lower.
Therefore at least half of the actual `D_V` law lies on shell-time atoms satisfying
\[ \mu_j\ge\mu_*. \]
Every such ancestor has
\[ M_j\le N/4, \]
so its natural lifetime is at least `16` parent-block lifetimes.
This is a resolved critical ancestor law, not a selected transfer parent law.
The high-strain route can enter generic critical-shell service directly without an NN material entrance.

---

## 12. Heat increment service and material ownership
The resolved heat-increment observable uses a Gaussian heat kernel at the block scale.
For `V` supported below `N/4`,
\[ e^{-1/32}\|\nabla V\|_2^2 \le N^2\int H_N(r)\|\delta_rV\|_2^2dr \le \|\nabla V\|_2^2. \]
Thus normalized heat service satisfies
\[ e^{-1/32}D_V\le S_{heat}\le D_V. \]
The same critical shell-time set retains a fixed positive fraction of this service.
Material endpoint labels are assigned in intrinsic coordinates.
For a heat edge with endpoints
\[ \zeta_0=(L^{-1}X/2,L^Tk), \qquad \zeta_1=(L^{-1}(X-r)/2,L^Tk), \]
OO/ON/NN is an exact positive partition.
Common affine/Kelvin transport preserves the labels.
A membership change under continuous nonaffine evolution requires boundary contact.
The local heat-edge capacity satisfies
\[ |e^{-ik\cdot r}A_1-A_0|^2 \le2(|A_0|^2+|A_1|^2). \]
Material ownership is therefore read from positive physical service, not imposed as a packet decomposition of `u`.

---

## 13. Old-pool erosion and NN critical seeds
At a first high-strain contact, pre-hit strain remains below the boundary.
Kelvin frequency growth is therefore bounded.
On a supplied signed-good epoch, block scales advance by more than `8/5` while reused material frequencies grow only by the low-strain Kelvin factor.
For old incident heat service, the clean one-generation capacity contraction is
\[ \rho_{old}<\frac{441}{640}<0.7. \]
This bound applies to OO plus ON, not only OO.
Therefore mixed ON service cannot form a neutral infinite regime.
Eventually a fixed positive fraction of high-strain heat service is NN.
Combining this with the critical shell-time heat fraction by inclusion–exclusion gives a positive NN-critical heat sublaw.
Push this law to its deterministic hard shell rather than selecting a largest coherent cell.
If
\[ M\|P_Mu(t)\|_2^2\ge\frac{32\pi^2}{75c^2}, \]
renew at
\[ A=3M/4. \]
Then
\[ A\|P_Mu(t)\|_2^2 \ge \frac{8\pi^2}{25c^2}. \]
The canonical smooth whole-shell carrier therefore has a scale-independent critical coefficient.
The NN endpoint mark remains material provenance.
The whole `u` shell is **not** declared new material.

---

## 14. Fresh SGS service: quotient coherent-cell refinement first
The coherent SGS service law is positive.
Coherent-cell entropy depends on observer refinement and therefore is not canonical for renewal entrance.
The intrinsic fresh/old ownership indicator is material.
Push fresh positive service only to the fixed canonical LP band index.
Let `F_j` be the integrated fresh service in band `j` and
\[ F=\sum_jF_j. \]
Coherent-cell refinement subdivides the same positive measure and leaves every `F_j` unchanged.
Use a square-normalized dyadic annular LP frame
\[ \sum_j\phi_j^2=1, \qquad |\phi_j|\le1, \qquad \operatorname{supp}\phi_j\subset\{M_j/2<|\xi|<2M_j\}. \]
For one fresh band on scaled lifetime `c`,
\[ F_j\le4\int M_j\|u_j\|_2^2d\tau_N. \]
Hence at some time
\[ M_j\|u_j\|_2^2\ge\frac{F_j}{4c}. \]
The annular band touches at most two exact hard shells.
Therefore one actual hard shell satisfies
\[ \mu_{hard}\ge\frac{F_j}{6c}. \]
If the total fresh law obeys the certified lower
\[ F\ge Y/4, \]
and
\[ p_j=F_j/F, \qquad p_{max}=e^{-H_\infty^{scale}}, \]
then
\[ \boxed{ \mu_{hard}e^{H_\infty^{scale}} \ge \frac{Y}{24c}. } \]
The collision-entropy corollary follows because `p_max>=sum p_j^2`.
`H_inf^scale` is a service-scale concentration coordinate.
It is not a child-energy causal entropy.
Cell-dominance/ancestry entropy remains optional sideledger metadata only.

---

## 15. Objective source owner compiler
Large coherent objective variation is grouped by physical source owner:
\[ local\ D_V + pressure + SGS + viscosity. \]
If the scaled objective source action is `A_obj`, then
\[ \Sigma_{local}+\Sigma_P+\Sigma_R+\Sigma_\nu \ge A_{obj}/c. \]
At least one physical owner carries at least one quarter of that total.
Exact ties are joint.
### 15.1 Local resolved owner
The local coherent quadratic/Reynolds source satisfies
\[ \Sigma_{local}\le C_{local}D_V. \]
Thus it routes to resolved dissipation and then generic critical-shell reentry.
### 15.2 Viscous owner
The pointwise viscous source estimate has the form
\[ \rho_\nu\le\nu\sqrt{d_V}/1500. \]
Integrated source therefore forces a positive resolved dissipation lower and enters the generic shell theorem.
### 15.3 Direct SGS owner
The clean objective SGS collision gives
\[ \|R\|_{3/2}\ge380\rho_R. \]
Germano exponent `3/2` and coherent square-service exponent `2/3` cancel exactly.
Thus integrated SGS source yields coherent increment service **linearly** in source weight.
No packet persistence is needed.
The resulting service routes jointly to any realized native owners:
- high-frequency dissipation;
- old-pool capacity;
- selected-interface `Xi`;
- fresh-scale critical shell.
Source-code branch order has no theorem meaning.
### 15.4 Pressure owner
Pressure is handled by direct objective-Hessian pair atomization, not by aggregate low-pass mass.
The old `mu_V/5700` pressure reservoir inequality remains diagnostic only.

---

## 16. Objective pressure: exact dual split and unordered hard pairs
Let the actual coherent objective pressure Hessian matrix be `H`.
Use the measurable Frobenius dual
\[ Z=H/\|H\|_F \]
when `H!=0`, otherwise `Z=0`.
Then
\[ \rho_P=Z:H. \]
Split pressure exactly into resolved `V\otimes V` and SGS-pressure parts.
After positive parts,
\[ \rho_P \le [r_{SGS}]_+ + \sum_{a\le b}[p_{ab}]_+. \]
This is source positivity, not a child-energy causal law.
Decompose only the resolved transporter into hard orthogonal dyadic shells.
Use **unordered** physical pairs `{a,b}`:
- diagonal once;
- off-diagonal with both orientations combined.
This removes fake observer-orientation entropy.
The pressure stress-to-Hessian symbol has Frobenius operator norm at most `|xi|^2`.
No hidden Riesz `L^(3/2)` constant is needed.
The clean pair capacity is
\[ |p_{ab}| \le \frac{\kappa_{ab}}5 \left(\frac{M_{max}}N\right)^4 \sqrt{\mu_a\mu_b}, \]
with
\[ \kappa_{ab}=1 \]
on the diagonal and
\[ \kappa_{ab}=2 \]
off diagonal.
The positive pressure owner splits:
- SGS pressure source at least `Sigma_P/2`, or
- resolved positive pair source at least `Sigma_P/2`;
with ties joint.
For the pair law, its collision entropy `H_2^P` and shell mass obey the continuous relation
\[ \boxed{ \mu_{child}e^{H_2^P} \ge 320\frac{\Sigma_P}{c}. } \]
Thus every resolved pressure-pair owner already enters generic critical-shell reentry.
The familiar quarter split is only a corollary:
- `q_max>=1/4` gives `mu_child>=80 Sigma_P/c`;
- `q_max<=1/4` gives `H_2^P>=log 4`.
There is no separate diffuse-pressure fate.
The coarse aggregate `mu_V` pressure route is not canonical.

---

## 17. Pressure derivative orders must remain separate
Objective pressure Hessian and H1 pressure-third terms have different derivative order.
They must not share one reuse coefficient theorem.
For a fixed material objective-Hessian pair, the signed-good low-strain generation ratio is bounded by
\[ (21/20)^5(5/8)^4 = 4084101/20971520 <1/5. \]
Thus future fixed-pair Hessian capacity is geometrically summable with factor below `1/5`.
For pressure-third,
\[ (21/20)^4(5/8)^3 = 194481/655360 <1/3. \]
The derivative distinction is structural, not bookkeeping.
The fixed-material reuse theorem is optional after material sidecars are attached.
It is not the pressure renewal entrance.

---

## 18. High-frequency coherent service to physical hard-tail owners
The coherent `D_high` observable is a smooth LP quantity.
It is not identical to physical hard-tail dissipation.
Use a square-normalized Calderón analysis–synthesis frame
\[ u_j=\phi_j(D)u, \qquad \sum_j|\phi_j|^2=1. \]
With high-band lower support `aM_j` and square-Bessel upper `B`, Plancherel gives
\[ D_{tail}\ge c_{LP}D_{high}, \qquad c_{LP}=a^2/B. \]
Canonical values `a=1/2`, `B=1` give
\[ \boxed{D_{tail}\ge D_{high}/4.} \]
Do not relabel `D_high` as resolved `D_V`.
For
\[ w=P_{>N}u, \]
physical tail energy obeys
\[ N\|w(s)\|_2^2+NW_>^+ \ge 2\nu D_{tail}. \]
Thus one of two native owners carries at least `nu D_tail`:
1. inherited tail energy;
2. positive nonlinear regeneration.
### 18.1 Inherited tail energy
Hard-shell disintegration gives an actual shell `M_j>=2N` satisfying
\[ M_j\|P_ju(s)\|_2^2 \ge \nu D_{tail}. \]
This is an immediate generic critical-shell seed.
### 18.2 Regeneration
Orthogonal shell work satisfies
\[ \sum_jW_j^+\ge W_>^+. \]
At each shell choose `V=S_{M_j/4}u`.
Low-low output is below the child hard shell.
Therefore signed shell work is exactly HH plus resolved mixed/interface work.
Positive parts give
\[ W_{shell}^+ \le W_{HH}^++W_{interface}^+. \]
This regeneration owner is **not automatically** the generated-energy productivity branch.
Actual HH child-work and interface work continue through their own physical routes.

---

## 19. High-tail causal unit: `N dW`, not own-scale reweighting
For high-tail regeneration, cross-shell causal comparison must use a common physical work unit.
The canonical unit is
\[ N\,dW, \]
where `N` is the parent block scale whose tail energy identity generated the owner.
From the hard-tail energy branch,
\[ NW_>^+\ge\nu D_{tail}. \]
Hard-shell disintegration gives
\[ \sum_jNW_j^+ \ge NW_>^+. \]
Low-low exclusion yields
\[ H_N^++I_N^+ \ge \sum_jNW_j^+. \]
Hence actual positive HH or interface work carries at least
\[ \boxed{\nu D_{tail}/2.} \]
with exact ties joint.
If HH owns the event, exact coherent Hahn atomization gives a binary positive physical child-work law with total mass at least the same clean lower.
Do **not** normalize by `M_jdW_j` across shells.
That changes relative causal probabilities by the observer's scale-dependent unit.
Own-scale weights `M_jdW_j` may be read only after the common-unit causal law exists, as a strength diagnostic.
Causality and scale strength are different observables.

---

## 20. High-tail Fourier locality without a local-triad assumption
Push the actual positive HH common-work law to hard output shell level:
\[ H_j=N\int[r_{HH,j}(t)]_+dt. \]
Normalize
\[ p_j=H_j/H, \qquad H=\sum_jH_j. \]
Choose the maximal output shell `M=2^jN` and set
\[ H_\infty^{out}=-\log p_{max}. \]
This scale law is read **before coherent Hahn refinement**.
Atomic positive Hahn mass is not a substitute because cancellation can enlarge it and change `p_max`.
For any field Fourier-supported above `K`,
\[ \|\widehat f\|_{3/2} \le \left(\frac{4\pi}{3}\right)^{1/6} K^{-1/2}\|\nabla f\|_2. \]
If one parent feeding child shell `M` lies above `RM`, triad closure forces the other above `(R-1)M`.
On these UV parent frequencies, the strict transporter `S_{M/4}u` vanishes, so `h=u` exactly.
Combining both high-pass parent bounds, the child shell bound, and sharp Young collapses the constant to
\[ \boxed{3\sqrt\pi}. \]
If
\[ \mu_{peak}=\max_t M\|P_Mu(t)\|_2^2, \]
then
\[ NW_{UV}^{abs} \le \frac{3\sqrt\pi}{\sqrt{R(R-1)}} \sqrt{\mu_{peak}}D_{tail}. \]
The canonical theorem is the **continuous** locality-radius relation
\[ \boxed{ \frac{W_{comp,R}e^{H_\infty^{out}}}{D_{tail}} + \frac{3\sqrt\pi}{\sqrt{R(R-1)}} \sqrt{\mu_{peak}}e^{H_\infty^{out}} \ge \frac\nu2. } \]
`R` is a physical locality-radius parameter, not a stop class.
For `R=2`, a readable balanced corollary is
\[ \mu_{peak}e^{2H_\infty^{out}} \ge \frac{\nu^2}{72\pi} \]
or
\[ W_{comp,2}e^{H_\infty^{out}} \ge \frac{\nu D_{tail}}4. \]
The first owner is already a critical shell.
The second leaves actual positive comparable-parent HH work.
Only after this Fourier theorem is proved may the comparable source be coherently atomized.
Locality is a PDE triad property, not a coherent-cell property.

---

## 21. High-tail natural-time concentration and shell reentry
This closes the current high-tail temporal seam.
Suppose the certified Fourier locality route carries the clean comparable-HH owner on selected shell
\[ M=2^jN, \qquad j\ge1. \]
Then
\[ T_N=cN^{-2}, \qquad T_M=cM^{-2}=4^{-j}T_N. \]
Let `mu_comp` be the actual positive comparable-work measure on the parent block.
Define the sliding natural-window fraction
\[ \boxed{ p_t = \frac{\sup_s\mu_{comp}([s,s+T_M])} {\mu_{comp}(I)}. } \]
Set
\[ H_\infty^{time}=-\log p_t. \]
There is no fixed time grid.
For smooth Navier–Stokes the positive comparable-work density is continuous, so the sliding integral is continuous and its maximum is attained on the compact admissible start interval.
The fraction `p_t` is invariant under:
- translation of the time origin;
- simultaneous rescaling of time units and density;
- subdivision/refinement of a representation of the same positive measure.
The theorem API requires the concentration window to equal the selected shell's exact `T_M`.
A concentration number measured on another scale cannot be attached by hand.
The canonical strict cutoff currently certifies only
\[ |S|\le1, \]
hence
\[ \|h\|_2\le2\|u\|_2. \]
For comparable parents `<=RM`, sharp physical Young gives
\[ r_{comp} \le 12\sqrt\pi\,R\,M^2E_{global}\sqrt{\mu_{win}}, \]
where
\[ \mu_{win} = \max_{t\in I_*}M\|P_Mu(t)\|_2^2. \]
Integrating on `T_M=cM^-2` cancels `M^2`:
\[ W_{win}^{phys} \le 12c\sqrt\pi\,R\,E_{global}\sqrt{\mu_{win}}. \]
In the high-tail common causal unit,
\[ NW_{win} \le 12c\sqrt\pi\,R\,N E_{global}\sqrt{\mu_{win}}. \]
The locality owner gives
\[ W_{comp}/p_{scale} \ge \nu D_{tail}/4. \]
Since
\[ W_{win}=p_tW_{comp}, \]
one obtains
\[ \boxed{ \frac{\sqrt{\mu_{win}}}{p_{scale}p_t} \ge \frac{\nu D_{tail}} {48c\sqrt\pi\,R\,N E_{global}}. } \]
Equivalently,
\[ \boxed{ \sqrt{\mu_{win}} \exp(H_\infty^{out}+H_\infty^{time}) \ge \frac{\nu D_{tail}} {48c\sqrt\pi\,R\,\mathcal E_N}, \qquad \mathcal E_N=NE_{global}. } \]
For `R=2`, the denominator is
\[ 96c\sqrt\pi\,\mathcal E_N. \]
The peak `mu_win` is an actual hard-shell event and enters generic critical-shell first stopping.
Hard-tail support also supplies genuine scale geometry:
\[ M/N=2^j\ge2, \]
and
\[ T_M/T_N=4^{-j}\le1/4. \]
This progress does not use signed-good Young near-extremality.
Full-survivor own-scale service remains conditional on the generic shell corridor.

---

## 22. Physical pair productivity: what the energy gate is actually for
When a generated child-energy branch satisfies the physical energy gate, use actual positive child-work as the probability law.
For one physical pair-work density,
\[ r_e(t) \le C_Y N a_c a_1a_2. \]
KL positivity against normalized time × hard pair-cell reference gives
\[ \mathbb E_{dT}\log(a_1a_2) \ge \mathbb E_{dT}\log a_c + \log\Lambda_j. \]
This is the correct multiplicative productivity statement.
It does not identify Duhamel pair weights with physical work weights.
Amplitude imbalance between the two parents is allowed.
Only their transfer-weighted logarithmic product matters.
Common-slice registration and dual-Gaussian marking contribute finite multiplicative constants.
The energy gate is therefore a **productivity gate**, not a **causality gate**.
High-tail HH work can be causal before this gate.

---

## 23. Complex Young and dual-Gaussian parent marking
Near-extremal physical pair work can mark parent roles without asserting absolute parent mass.
For a frozen selected cell, if weighted transfer is large relative to the normalized parent shape and symbol error is small, complex Young gives a lower bound on the normalized parent transfer coefficient.
The Christ near-extremizer theorem is used at the modulus/shape level.
Restoring actual amplitude preserves trilinear homogeneity.
Therefore:
- normalized shape rigidity can force a scale-critical coefficient;
- it cannot force an absolute physical root mass independently of amplitude.
The dual-Gaussian mark gives a nearby energetic parent identity anchor.
Bargmann submean plus Moyal yields a physical cell energy proportional to the **square of the actual parent amplitude**.
The canonical normalized Bargmann cell constant is small but scale independent.
It is not an absolute root-energy quantum.

---

## 24. Common-slice registration and asynchronous support
Duhamel remains useful as an exact support/adjoint identity.
It provides common backward reference slices for generated roles.
On the signed-good generated lineage, parent/child scale ratios give compatible natural-time windows.
At a common slice the exact coefficient identity is
\[ z(t)=z(s)+I_{HH}+I_R. \]
A triangle split gives the alternatives
\[ |z(s)|\ge A/4, \]
or
\[ |I_R|\ge A/4, \]
or
\[ |I_{HH}|\ge A/2. \]
The first is retained carrier amplitude.
The second and third are coefficient **first-stop locators**, not physical work owners.
A large coefficient impulse alone is **not** automatically large physical work.
At a coefficient hit the same smooth carrier must reenter the actual `Q^2` energy law. For an interface hit, the common observer gauge must additionally be quotiented before gauge-quotiented native work is Hahn-routed to residual physical relink or existing strain. HH generation is named only when the physical-energy gate selects actual positive HH work.
Raw coefficient-obstruction labels are forbidden from entering the canonical master physical-owner state directly.
## 25. Recursive first stopping and joint master projection
The recursive physical witness constructor acts on retained **actual positive work** after the relevant physical gate.
`Xi` is excised once before causal routing.
At each physical event, the preferred master records the first actual cause set.
Coefficient threshold debuts are recorded one layer earlier as typed first-stop locators. They are not admitted to the canonical physical-owner set until actual `Q^2` energy/work reentry resolves them.
Terminal semantics are:
1. `t=0` → absorbing initial boundary;
2. fixed multiplicative transfer/reuse cost → `TRANSFER_COST`;
3. genuine scale-independent globally bounded resource → `ADDITIVE_RESET`;
4. source/SGS, critical dissipation, new ancestry/material relink, HH regeneration, shell reentry → `RECURSE_CRITICAL`.
Exact simultaneous causes remain one joint set.
No lexicographic priority is canonical.
No dummy RN tie weight may alter the coarse physical fate.
If a fraction `q` of generated physical work survives all earlier first stops and registration, then the physical productivity constant becomes
\[ \Lambda_{survivor}=q\Lambda_{full}. \]
If
\[ q\ge1/2, \]
the logarithmic conditioning penalty is at most `log 2`.
If
\[ q<1/2, \]
a majority of actual current work has already exited through a named cause or earlier regeneration.
There is no fourth free continuation branch.

---

## 26. Shannon/Rényi reuse: only under genuine causal weights
For generated child events, each event has two structural parent slots.
Under the actual positive child-work law, physical pair productivity gives a recursive lower bound on expected log parent amplitude.
Bargmann/Moyal gives an energy budget for distinct material anchors proportional to parent amplitude squared.
The log-sum inequality then bounds root entropy without assuming a uniform root mass floor.
The old hypothesis
\[ NE_{root}\ge\eta \]
for every root is unnecessary and generally incompatible with Young homogeneity.
The correct entropy bound uses
\[ H(w_{root})+2\mathbb E\log\alpha_r \le \log\sum_r\alpha_r^2. \]
The Shannon telescope retains a positive linear depth slope after finite logarithmic productivity losses.
Because
\[ H_2\le H_1, \]
Rényi reuse obtains the same structural lower.
A rich layer then routes to existing parent-pair, component-entropy, or same-ancestry-cycle structures.
This machinery applies to a genuine causal law.
Do not feed it fresh-SGS scale entropy, pressure-pair source entropy, or high-tail output-scale concentration merely because all are written with `H` symbols.

---

## 26.5 Continuum master event quotient: physical time, not theorem depth
The first final-master assembly theorem is now certified.
Its central correction is that theorem depth is not recursion depth.
Same-law owner relays preserve one unsplit physical measure, while certified witness relays may change observable and units without creating a second causal charge.
Thus source/SGS/pressure/service/shell compiler layers are quotiented before recursive depth is charged.
The canonical recursive state contains only:
- actual physical event time;
- a supplied physical shell/carrier frequency;
- the named physical owner measure;
- the unsplit joint physical cause set;
- optional material or diagnostic sidecars.
There is no common normalized master clock and no theorem-priority coordinate.
For any ordered backward physical event or checkpoint times,
\[ \boxed{\sum_j(t_j-t_{j+1})=t_0-t_L.} \]
A no-hit physical corridor at its own corridor frequency `M` consumes exactly
\[ cM^{-2} \]
of backward Navier--Stokes time unless that corridor reaches `t=0`, which is absorbing.
An interior no-hit horizon is an analysis checkpoint and contributes zero physical event vertices.
Hence if a no-hit checkpoint continuation has corridor frequencies bounded by `Mbar`, then after at most
\[ \boxed{\left\lceil t_*Mbar^2/c\right\rceil} \]
such physical corridors the initial boundary is forced.
The later checkpoint quotient sharpens the topology: these no-hit horizons add zero event vertices. Thus any infinite recursive **event** path avoiding `t=0` must contain infinitely many genuine physical owner events, while UV-unbounded no-hit checkpoint continuation is a separate event-free PDE seam rather than a second branch of the same event path.
The UV alternative is genuine, not a clock artifact, because for `M_j=M_0r^j`, `r>1`,
\[ \sum_{j\ge0}cM_j^{-2}=\frac{cM_0^{-2}}{1-r^{-2}}<\infty. \]
Scale progress remains supplier-specific:
- generated signed-good HH: `3/5<N_next/N<5/8`;
- resolved `D_V` ancestor: `N_next/N<=1/4`;
- resolved pressure pair: `N_next/N<=1/4`;
- high tail: `N_next/N>=2`;
- fresh SGS: only `N_next/N<=2`, with no directional progress theorem;
- generic shell, material/reuse, and unresolved HH regeneration: no synthetic scale progress.
The natural master Bellman object is therefore typed/direct-product rather than scalar.
Physical time and actual log shell scale telescope kinematically; actual multiplicative transfer cost, work-weighted causal reuse, `Xi`, and each genuinely globally bounded resource telescope only in their own native ledgers.
Own-scale service generated by a completed no-hit corridor is a same-interval witness and adds zero event depth. Independent source/service/reuse events remain physical owners. No service quantity is promoted to an additive finite resource.
Fresh/high-tail `H_inf/H2` remain diagnostic concentration coordinates and are forbidden from entering causal Shannon/Rényi action.
This theorem is an assembly/quotient theorem only; it does not prove global no-escape or Navier–Stokes regularity.

---

## 26.6 Resolved interface is flux/strain, not a third generator
The resolved cross/interface branch is now structurally quotiented before master recursion.
At a shell-time event, fix the actual resolved low--high linearized operator on the selected high field and take its adjoint split
\[
L_V=K+S,\qquad K^*=-K,\qquad S^*=S.
\]
For a complete orthogonal event-role partition `w_a=P_a h`, the signed resolved role work splits exactly as
\[
R_a=R_a^K+R_a^S.
\]
Hence atom by atom
\[
[R]_+\le [R_K]_+ + [R_S]_+.
\]
The two pieces have different but already native physical meanings.
The `S` part is the existing symmetric strain/deformation work and delegates once to the coherent-deformation, high-strain, objective-source, or critical-`D_V` owners.
The `K` part is conservative role redistribution.  Its pair work satisfies
\[
T_{ab}=-T_{ba},\qquad R_a^K=\sum_bT_{ab}.
\]
With `F_{b->a}=[T_{ab}]_+`, every role set `C` obeys the exact finite-role divergence identity
\[
\boxed{\sum_{a\in C}R_a^K=F_{C^c\to C}-F_{C\to C^c}.}
\]
Therefore any role with positive skew gain has a finite backward donor closure containing an actual negative-net donor.  After cycles are deleted, a donor path has at most `#roles-1` edges.
All such donor tracing occurs at the same physical event time.  Internal role circulation is real redistribution but is neither energy generation nor a new recursive PDE generation; it creates no second causal charge and no scale progress.
The full resolved mixed-work observable and the moving-projector commutator observable are companion readings of the same `K/S` operator provenance; they are not identified as the same measure.
On the clean high-tail interface owner, the unchanged common `N dW` law gives
\[
W_{interface}^+\ge\nu D_{tail}/2
\quad\Longrightarrow\quad
W_K^+\ge\nu D_{tail}/4
\quad\text{or}\quad
W_S^+\ge\nu D_{tail}/4,
\]
with exact ties joint.
No `M/N` causal reweighting, new entropy, synthetic clock, or interface reset is introduced.
This closes resolved interface as an independent recursive-generation loophole; it does not prove that the genuine donor/strain owner reached afterward globally terminates.

---

## 26.7 Smooth propagated interface: `Q^2` energy and observer-gauge quotient
The hard event-role donor theorem does not directly apply to the propagated smooth non-idempotent envelope.
For `w=Qu`, the physical carrier energy is
\[ \boxed{E_Q=\|Qu\|_2^2=\langle u,Q^2u\rangle.} \]
A complete smooth analysis uses `eta_a=A_a^2` with
\[ \sum_a\eta_a=I. \]
The native smooth-interface row is
\[ J_a=\langle u,\partial_t\eta_a u\rangle-2\operatorname{Re}\langle\eta_a u,L_Vu\rangle. \]
It is also exactly the outer Heisenberg work minus the diagonal resolved work, so the commutator is never interpreted alone.

The decisive quotient comes before ownership. Let `G` be the common skew generator of the affine/Kelvin transport actually used to propagate all smooth roles. Require
\[ \boxed{\dot A_a+[G,A_a]=0\quad\forall a.} \]
Then
\[ \langle u,\partial_t\eta_a u\rangle-2\operatorname{Re}\langle\eta_a u,Gu\rangle=0. \]
Writing the actual skew resolved operator as
\[ K=G+K_{phys}, \]
the common observer motion disappears and only
\[ -2\operatorname{Re}\langle\eta_a u,K_{phys}u\rangle \]
remains as smooth physical relink. Its synthesis-pair matrix is antisymmetric and has zero total work. The symmetric `S` rows reconstruct the existing resolved strain/deformation work.

A square partition whose motion is chosen independently by the observer may still exchange channel energy with zero total. That fact is **not** sufficient for physical ownership. Such motion is rejected unless it satisfies the certified common transport equation above.

The common-slice coefficient thresholds `|I_R|>=A/4` and `|I_HH|>=A/2` are therefore typed first-stop locators only. They cannot enter `RecursiveEventState` or `PhysicalOwnerBundle` as physical owners. The same carrier must first reenter actual `Q^2` energy/work causality. Interface work is owner-eligible only after the observer gauge has been quotiented; HH generation is owner-eligible only after actual positive HH work is selected.

This closes the smooth-envelope/projector mismatch, the arbitrary-moving-window loophole, and the coefficient-locator/master type hole. It creates no new interface currency and does not prove termination of the genuine physical owners reached afterward.

---

## 26.8 Full-natural service is a corridor witness, not a second event
A generic critical shell which survives its full backward natural interval has already traversed the physical corridor
\[
I=[t-cA^{-2},t],\qquad A=\frac34M.
\]
The own-scale bounded heat/increment service proved on that interval is a genuine positive Navier--Stokes observable, but it is supported on the **same completed corridor**. Reading that service does not create another event time, causal charge, or recursion edge.

The same is true of its material reading. Exact Moyal disintegration gives
\[
S_{OO}+S_{ON}+S_{NN}=S_{service}.
\]
OO/ON/NN are positive submeasures of the same service law. They may certify downstream physical states, but their mere rereading adds zero recursion depth.

The surviving smooth carrier is already present at the earlier corridor endpoint. Its transported support lies inside
\[
(A/2,2A).
\]
Splitting into the exact hard shells `(A/2,A]` and `(A,2A]`, with
\[
\mu_A=A\|P_Au\|_2^2,\qquad \mu_{2A}=2A\|P_{2A}u\|_2^2,
\]
gives
\[
A\|Q_Au\|_2^2\le \mu_A+\frac12\mu_{2A}
\le\frac32\max(\mu_A,\mu_{2A}).
\]
Hence
\[
\boxed{\max(\mu_A,\mu_{2A})\ge\frac23A\|Q_Au\|_2^2.}
\]
This is a same-endpoint hard-shell witness set. Relative to the incoming shell `M`, its frequencies are `3M/4` and `3M/2`. This is comparable-scale geometry, not monotone progress. Exact ties remain joint.

Therefore a chain whose apparent extra layers are only

`full natural corridor -> own-scale service -> Moyal/material rereading -> endpoint survivor`

contains one real physical corridor, but if no first stop fires its natural-horizon endpoint contributes zero physical event vertices. After this quotient, bounded-frequency checkpoint continuation still hits `t=0` by physical time. If it avoids `t=0`, the remaining no-hit possibility is an event-free UV continuation seam rather than a recursive-event chain.

This closes the service-theorem-depth and endpoint-service-attachment seams. The subsequent checkpoint theorem further removes the no-hit horizon itself from event depth. Independent source/service/reuse events remain physical, and UV-unbounded no-hit continuation remains open.

---

## 26.9 Full-natural horizons are analysis checkpoints, not events
For a generic critical shell at incoming frequency `M`, the actual no-hit theorem runs at the renewal scale
\[
A=\frac34M,\qquad T_A=cA^{-2}.
\]
A physical first stop before `T_A` remains event-facing after its required routing, and `t=0` remains absorbing. If neither occurs and the whole interval survives, the elapsed `T_A` is genuine Navier--Stokes time but the earlier endpoint was selected only because the theorem chose one natural horizon. It is therefore an **analysis checkpoint**, not a new `RecursiveEventState`.

The endpoint smooth carrier may be reread through the two exact hard shells at `A` and `2A`. Their ratios to the incoming shell are `3/4` and `3/2`, but this is checkpoint-cover provenance rather than high-tail provenance. Production re-registration accepts the actual endpoint critical masses `(mu_A,mu_2A)`, invokes the exact hard-shell realization internally, returns the unique physical maximum or the full exact tie, and never accepts an observer-chosen frequency branch.

Thus the upper `3M/2` witness is a genuine state observable when its mass wins, but its appearance does not by itself prove nonlinear UV generation, directional progress, or a high-tail owner. The independently certified high-tail route keeps its own physical `D_tail`/work provenance.

Checkpoint chains telescope their real physical corridor time exactly while adding zero event vertices and zero causal charges. A bounded-scale checkpoint continuation therefore reaches `t=0`. A UV-growing sequence of no-hit checkpoints may still have finite total physical duration; after quotienting horizon segmentation this is an **event-free PDE continuation seam**, not an infinite recursive event path. Closing that seam remains open.

---

## 27. Finite-dimensional rigidity modules: supporting geometry, not current PDE bottleneck
Earlier layers of the programme established a large library of exact or certified finite-dimensional geometric facts.
They remain important as local rigidity/registration tools but are no longer the main continuum seam.
The essential surviving modules include:
### 27.1 Helical triad geometry
The helical decomposition identifies signed interaction geometry, extremal polarization, and phase constraints.
Single-edge stability exposes a narrow signed-good region and log-scale progress window.
### 27.2 Smooth log-scale cocycle
Signed-good transfer produces controlled forward scale progress and bounded geometric defect.
This remains the natural scale geometry on near-extremal generated lineages.
It is distinct from hard-tail support progress `M/N>=2`.
### 27.3 Scale holonomy and Hodge structure
Reuse curvature, cycle defects, and Bellman/Hodge formulations quantify how repeated ancestry must either reuse structure or pay physical defects.
### 27.4 Flat-network erosion
Flat/cycle countermodels showed why purely combinatorial tree growth is insufficient.
Spherical erosion, component entropy, multicommodity Hodge routing, and resistance stopping remove several finite-dimensional escape modes.
### 27.5 Gaussian/affine grain modules
Affine Gaussian dynamics, critical-grain energy, material phase lock, spin transport, curvature, and sideband structures motivated the later coherent/Kelvin framework.
The current canonical continuum architecture uses only the pieces that survive exact gauge quotienting and physical ownership tests.
Historical derivations and superseded packet formulations are preserved in the history note.

---

## 28. Anti-theorems and forbidden inferences
The following statements are either false or non-canonical and must not silently re-enter the proof.
### 28.1 Raw Duhamel mass is not physical work
`dGamma=dT` is false.
Use Duhamel for support/adjoint identities, physical positive work for causality.
### 28.2 Large coefficient impulse is not automatically large physical work
A Duhamel or interface coefficient impulse requires an independent energy/work theorem before being promoted to physical generation.
### 28.3 Young shape does not create an absolute root mass floor
Trilinear homogeneity allows arbitrary amplitude scaling.
Any physical cell-energy quantum scales with the parent amplitude squared.
### 28.4 Small-parent / large-reservoir is not a new currency
Amplitude imbalance belongs inside multiplicative productivity or existing reservoir/service routing.
### 28.5 Center Hessian alone does not close nonaffine transport
Coherent deformation variance is the correct affine-invariant observable.
### 28.6 Higher Hermite towers are not the canonical closure mechanism
Use exact outer-role PDE plus native nonaffine/interface/source ownership.
Do not create an infinite observer-chosen Hermite hierarchy merely to describe the same physics.
### 28.7 Frozen packet/profile persistence is not assumed
The architecture uses event hard roles, smooth moving carriers, first stopping, and sliding physical measures.
### 28.8 Coherent-cell dominance is not canonical for fresh SGS renewal
Fresh renewal uses the refinement-invariant LP band pushforward.
Cell entropy/cycle is ancestry sideledger only.
### 28.9 Aggregate `mu_V` is not the canonical pressure renewal entrance
Pressure uses exact Frobenius dual source splitting and unordered hard pairs.
The coarse `mu_V/5700` estimate is diagnostic only.
### 28.10 `D_high` is not `D_V`
Smooth LP high-frequency service must first be compared with physical hard-tail dissipation.
### 28.11 Own-scale reweighting is not high-tail causal probability
Across shells, use common `N dW`.
`M_jdW_j` is a post-causal strength diagnostic.
### 28.12 Atomic Hahn mass is not the high-tail output-scale law
Locality is read from aggregate positive HH hard-shell work before coherent refinement.
### 28.13 A chosen time partition is not natural-time concentration
Use the sliding physical window `cM^-2`.
The result must be invariant under time-origin and unit gauges.
### 28.14 A concentration coordinate is not automatically causal entropy
Always name the underlying positive measure.
Pressure `H_2^P`, fresh `H_inf^scale`, high-tail `H_inf^out`, and `H_inf^time` are deterministic source/service/work concentration coordinates.
Only the work-derived ancestry law is causal Shannon/Rényi input.
### 28.15 Pure material relabeling is not a carrier impulse
If `Q` and `psi` are unchanged, label-sidecar changes do not enter the coefficient identity a second time.
### 28.16 Cutoff changes are not new sources
The exact cutoff repartition identity forbids a cutoff-switch currency.
### 28.17 Critical normalized energy/dissipation is not an additive reset
Physical cost decays with scale.
Only a globally bounded scale-independent resource may reset additively.
### 28.18 Exact ties must not be broken by theorem priority
The first actual physical causal root owns the event.
Joint ties remain joint.
### 28.19 Pressure Hessian and pressure-third are not one derivative theorem
Their scale erosion coefficients differ and must remain separate.
### 28.20 High-tail comparable work is not automatically productivity-good
The natural-window theorem creates an actual critical shell from positive comparable work.
It does not claim the HH work satisfies the generated-energy `8/15` gate.
### 28.21 Conditional full-survivor service is not unconditional service
A shell event enters a first-stop corridor.
Service is recorded only if the full no-hit natural corridor survives.
### 28.22 A critical shell is not automatically fresh material
Material provenance is attached only by actual physical service/ownership information.

---

### 28.23 Conservative analysis-window exchange is not automatically physical relink
A time-dependent square partition can move channel energies while their sum stays fixed even when the Navier–Stokes state is unchanged.
Before charging relink, prove `dot A_a+[G,A_a]=0` for the common physical transport and quotient that gauge. Only residual `K_phys` work is physical relink.
### 28.24 A coefficient first-stop locator is not a physical owner
The labels produced when `|I_R|` or `|I_HH|` crosses its registration face describe where clean coefficient continuation failed.
They are forbidden from the canonical physical-owner state until actual `Q^2` energy/work reentry resolves a physical owner.

---

### 28.25 A positive observable layer is not automatically a new physical event
A full-natural own-scale service law is physically real, but when it is proved on the corridor just traversed it is a witness of that interval, not a second recursive event. Moyal/material disintegration of the same service measure likewise adds no event time or causal charge. A downstream theorem must supply an actual new state/time if it creates a new recursion edge.

### 28.26 A theorem horizon and its cover branch are not physical dynamics
Completing one natural no-hit window does not make its endpoint a physical event. The physical corridor time is retained, but the horizon is an analysis checkpoint unless a first stop or `t=0` occurs. Likewise, exposing hard shells at `A` and `2A` does not allow the analyst to choose a desired scale branch. The actual endpoint masses decide the unique/joint witness set, and checkpoint-cover provenance is not high-tail provenance.

---

## 29. Structural constants worth remembering
These constants are not all fundamental; many are clean certified envelopes.
The important point is where each belongs.
### Low strain
\[ K_{low}=1/30. \]
### High-strain normalized dissipation lower
\[ D_*=\frac{32\pi^2}{75c}. \]
### High-strain critical ancestor mass
\[ \mu_*=\frac{32\pi^2}{75c^2}. \]
### Renewed high-strain shell mass
\[ \frac{8\pi^2}{25c^2}. \]
### Energy-generated HH work fraction
\[ 8/15. \]
### Common-slice retained coefficient factor
`1/4` on the clean survivor.
### Heat-shell lower comparison
\[ e^{-1/32}D_V\le S_{heat}\le D_V. \]
### Old-incident heat erosion
\[ 441/640<0.7. \]
### Objective SGS clean source coefficient
`380`.
### Objective pressure resolved pair coefficient
\[ \frac{\kappa_{ab}}5(M_{max}/N)^4\sqrt{\mu_a\mu_b}. \]
### Pressure pair entropy-shell conjugacy
\[ \mu_{child}e^{H_2^P}\ge320\Sigma_P/c. \]
### Fresh SGS scale-shell conjugacy
\[ \mu_{hard}e^{H_\infty^{scale}}\ge Y/(24c). \]
### Smooth LP high-tail comparison
\[ D_{tail}\ge D_{high}/4. \]
### High-tail common-unit HH/interface owner
\[ \ge\nu D_{tail}/2. \]
### High-tail resolved-interface K/S component
\[ \ge\nu D_{tail}/4. \]
on the clean interface owner, in the same common `N dW` unit.
### High-tail UV geometric constant
\[ 3\sqrt\pi. \]
### High-tail continuous locality theorem
\[ \frac{W_{comp,R}e^{H_\infty^{out}}}{D_{tail}} + \frac{3\sqrt\pi}{\sqrt{R(R-1)}} \sqrt{\mu_{peak}}e^{H_\infty^{out}} \ge\nu/2. \]
### Dyadic `R=2` locality corollary
\[ \mu_{peak}e^{2H_\infty^{out}}\ge\nu^2/(72\pi) \]
or
\[ W_{comp,2}e^{H_\infty^{out}}\ge\nu D_{tail}/4. \]
### High-tail natural-window capacity
\[ NW_{win} \le12c\sqrt\pi\,R\,NE_{global}\sqrt{\mu_{win}}. \]
### High-tail scale-time shell conjugacy
\[ \frac{\sqrt{\mu_{win}}}{p_{scale}p_t} \ge \frac{\nu D_{tail}} {48c\sqrt\pi\,R\,NE_{global}}. \]
### Hard-tail forward scale geometry
\[ M/N\ge2, \qquad T_M/T_N\le1/4. \]
These constants should not be promoted into stop taxonomies unless the theorem explicitly assigns that semantics.

---

## 30. Current canonical route map
A useful compressed view is:
### 30.1 Generated physical work route
`terminal child energy`
→ physical energy gate
→ actual positive HH child-work
→ exact coherent binary parent-work atoms
→ physical pair productivity when needed
→ common-slice first stopping
→ material energy anchors
→ Shannon/Rényi reuse or named physical stop.
### 30.2 Objective source route
`coherent objective variation`
→ local / pressure / SGS / viscosity owner set
→ local or viscosity: `D_V`
→ generic critical shell;
→ pressure: SGS or unordered hard pair
→ coherent service or entropy-weighted critical shell;
→ SGS: coherent increment service
→ high-tail / old capacity / Xi / fresh scale shell jointly.
### 30.3 High-strain route
`renewed high strain`
→ `D_V` lower
→ actual dissipation-weighted critical resolved ancestors
→ heat increment service
→ material ownership if needed
→ generic critical shell / own-scale service.
### 30.4 Fresh SGS route
`positive fresh coherent service`
→ quotient coherent-cell refinement
→ canonical LP band law
→ one of two hard shells
→ critical shell first stopping.
### 30.5 High-frequency route
`smooth LP D_high`
→ physical `D_tail`
→ inherited hard-tail shell
or
→ actual nonlinear regeneration
→ common-unit HH/interface owner.
If the interface owner wins, it is immediately quotiented by `L_V=K+S`: `K` is same-event conservative donor provenance and `S` is existing strain/deformation ownership. Only genuine owners reached after that quotient enter recursive depth.
### 30.6 High-tail HH continuation
`positive HH common-unit work`
→ actual positive HH output-shell law
→ UV dissipation/locality tradeoff
→ critical shell immediately
or
→ comparable HH work
→ sliding `M`-natural window
→ actual critical shell
→ generic first-stop/service reentry.
This route now supplies its own forward scale progress `M/N>=2`.
### 30.7 Material route
`actual positive heat/increment service`
→ intrinsic endpoint labels
→ OO/ON/NN positive ownership
→ old-incident geometric erosion
→ NN-critical seed when old capacity is exhausted
→ whole-shell carrier
→ first stopping.
Material labels do not alter the smooth carrier coefficient identity by themselves.

---

### 30.8 Smooth propagated carrier route
`hard physical event P`
→ exact registration into smooth `Q` with `QP=P`
→ carrier energy read as `<u,Q^2u>`
→ common transported gauge `dot A+[G,A]=0` quotiented
→ coefficient obstruction, if any, remains only a locator
→ actual physical-energy reentry
→ inheritance / high strain / actual HH work / gauge-quotiented interface work
→ interface work only: residual physical relink `K_phys` or existing strain `S`
→ canonical physical-owner state.
No arbitrary motion of the smooth analysis windows and no raw coefficient impulse becomes a recursive generation.

---

### 30.9 Full-natural service corridor route
`actual critical hard shell at (t,M)`
→ backward first-stop corridor at `A=3M/4`
→ named physical stop / absorbing `t=0` / completed full-natural corridor
→ on the full survivor: own-scale service is attached to that same interval
→ OO/ON/NN is a same-measure material witness partition
→ at the same endpoint: hard-shell witness set at `A` or `2A`, with critical mass at least `(2/3)A||Q_Au||_2^2`
→ only a downstream theorem supplying a genuinely new physical state/time creates another recursion edge.
The endpoint ratios `3/4` and `3/2` are comparable geometry only; no directional scale progress is inferred.

---

### 30.10 Full-natural checkpoint route
`actual critical hard shell (t,M)`
→ run the native corridor at `A=3M/4`
→ physical first stop, if one occurs, remains event-facing after routing
→ `t=0`, if reached, absorbs
→ otherwise complete real corridor time `cA^-2`
→ endpoint is `full_natural_analysis_checkpoint`, not `RecursiveEventState`
→ read actual endpoint masses `(mu_A,mu_2A)`
→ exact realization returns the physical lower-shell witness, upper-shell witness, or joint exact tie
→ re-register analysis with zero event depth, zero causal charge, and no directional/high-tail provenance
→ continue until a genuine physical stop/owner occurs or an event-free UV checkpoint continuation must be closed by a separate PDE theorem.
No observer-selected cover branch can manufacture an increasing recursive path.

---

## 31. What remains open
The programme is no longer missing a generic “packet persistence theorem”, a common clock, or a definition of recursive depth.
The first continuum master quotient has reduced the final seam to genuinely physical recurrence.
### 31.1 Exhaustive measurable owner assembly
Every certified supplier must now be wired into the quotient event state on smooth pre-singular intervals, with each transition proved to be either:
- a zero-charge relay/witness map;
- a named non-free physical owner event;
- a real no-hit physical corridor ending in a zero-event analysis checkpoint;
- a true terminal cost/resource event; or
- absorbing `t=0`.
The wiring must preserve actual owner sets, exact ties, conditional first-stop semantics, and no double counting.
### 31.2 Observer-clock seam is structurally closed
The master uses the actual Navier--Stokes time coordinate for both physical events and theorem checkpoints, but only genuine physical stops/owners enter the event state. Supplier natural times remain local corridor horizons and are never normalized into a synthetic common clock. The horizon itself is now quotiented from event topology; what remains is physical continuation, not a clock-definition problem.
### 31.3 Hard and smooth interface-owner seams are structurally closed
Resolved mixed/cross-interface work no longer supplies an independent recursive owner chain.
For a complete orthogonal hard event-role partition of the actual resolved operator `L_V=K+S`, the skew part is conservative same-event role flux and the symmetric part is existing strain/deformation work. The skew donor closure is finite, internal circulation cancels in the subset divergence identity, and donor traversal creates no recursion depth or scale progress.

For the propagated smooth envelope, the measure is different and remains separately typed. Its energy is read at `Q^2`; common transported role motion is certified by `dot A_a+[G,A_a]=0` and quotiented as observer gauge; only `K_phys=K-G` is physical skew relink. Arbitrary square-partition motion is rejected even if its channel exchanges sum to zero. The symmetric branch remains the existing strain/deformation work.

Raw HH/interface coefficient first stops are also excluded from physical ownership: they locate a failed continuation interval and must pass through actual `Q^2` energy/work reentry before the master receives inheritance, high strain, HH work, relink, or strain ownership.

What remains open is termination/telescoping of those genuine physical owners after these quotients, not an additional hard-interface, smooth-window, or coefficient-amplitude mechanism.
### 31.4 Full-natural service attachment is structurally closed
Service remains conditional: it is recorded only after a critical shell survives the complete no-hit natural corridor. But once that corridor is complete, the own-scale service is a positive law **on the corridor already traversed**, not a second event after it.
Exact Moyal OO/ON/NN reading of the same service law is a zero-depth witness disintegration. The surviving smooth carrier is already present at the corridor endpoint and has a same-time hard-shell witness set at `A` or `2A` with
\[
\max(\mu_A,\mu_{2A})\ge\frac23A\|Q_Au\|_2^2.
\]
Thus the master must not attach a second recursive service owner merely because a service theorem is invoked. The no-hit horizon endpoint is likewise only a checkpoint. A genuinely new event vertex requires an actual first stop/owner law, not merely a new state reading.
What remains open is not endpoint-service attachment, but the continuation/telescoping of genuine first-hit, work, source, reuse, or independent service events and the separate event-free UV checkpoint seam.
### 31.5 Supplier-specific scale geometry is now registered, not scalarized
Hard-tail gives forward ratio at least `2`.
Signed-good generated transfer gives `3/5<N_next/N<5/8`.
Resolved dissipation and pressure-pair shells give lower-frequency ratios at most `1/4`.
Fresh SGS gives only an upper ratio `<=2` and no directional progress.
Generic shell/material/reuse routes get no invented scale progress. The full-natural checkpoint cover at `3/4` or `3/2` is state geometry only; its actual endpoint masses choose the witness and it is not a supplier-progress theorem.
The remaining termination argument must use these branch facts exactly as supplied.
### 31.6 Global termination has split into two sharp physical frontiers
Pure theorem/representation depth has now been removed from both owner recurrence and no-hit horizon continuation. Hard-interface circulation, smooth observer motion, raw coefficient locators, same-corridor service layers, and natural-horizon checkpoints do not manufacture event vertices.

The remaining problems are therefore **not one artificial path dichotomy** but two physically distinct questions.

**A. Genuine event recurrence.** Any infinite recursive event path avoiding `t=0` must contain infinitely many genuine non-free physical owner events after all quotients. These include actual first-hit/work/source/reuse events and independent service events. Their recurrence must telescope only through the native typed laws they actually supply: physical work, causal reuse, independent source/service laws, `Xi`, and genuinely globally bounded resources where available. No critical `NE`, `D_V`, shell mass or scale-critical service may be promoted to a finite reset.

**B. Event-free UV checkpoint continuation.** A sequence of no-hit checkpoints can consume real corridor time while adding zero event vertices. If its corridor scales remain bounded, physical time forces `t=0`. If the checkpoint scales grow, the total corridor time can be finite. The two-shell checkpoint cover does not close this: even a realized upper `3M/2` shell has checkpoint-cover provenance, not physical high-tail provenance. A future PDE theorem must show whether such UV continuation necessarily forces actual tail dissipation/work, another physical first stop, or some other native obstruction.

The independently certified high-tail route remains available when its physical hypotheses are actually met; the checkpoint theorem deliberately does not infer those hypotheses from cover geometry. Closing A and B, then connecting the result to the initial-data and singular-time interfaces, is the present global frontier.
### 31.7 Initial data interface
Backward causal recursion reaching `t=0` is already absorbing.
For regular initial data, band-limited root counts/energies have scale-decaying bounds.
The eventual complete continuum theorem must state the exact initial-data hypothesis and connect it to the existing initial-boundary root estimates.
### 31.8 Singular-time conclusion
Even after global master termination is proved, one must state precisely what contradiction or a priori estimate is obtained as a hypothetical singular time is approached.
No such global-regularity conclusion is claimed in the current ledger.

---

## 32. The deepest structural simplifications discovered so far
Several major reductions repeatedly followed the same pattern.
### 32.1 Cause before representation
Physical work exists before coherent atomization.
The coherent frame refines cause; it does not create it.
### 32.2 Scale before cells
Fresh SGS renewal is controlled by canonical LP band pushforward, not coherent-cell dominance.
High-tail locality is controlled by hard output-shell work, not atomic Hahn mass.
### 32.3 Natural time before time bins
High-tail temporal concentration is a sliding measure on `cM^-2`, invariant under clock origin and units.
### 32.4 Pressure source before reservoir narrative
Direct Frobenius dual pair atomization exposes hard shell mass without aggregate low-pass reservoir synchronization.
### 32.5 Material identity as a sidecar
The same smooth carrier may be relabeled materially without inventing a second PDE impulse.
### 32.6 Cutoff as gauge
Changing the resolved transporter repartitions one exact Navier–Stokes interaction.
It does not create source.
### 32.7 Energy gate as productivity, not causality
Positive HH work already defines a physical cause law.
The child-energy gate is needed only when one wants productivity relative to child energy.
### 32.8 Nonlocality paid by dissipation
Instead of assuming local triads, UV high-high→low work is bounded by the physical high-tail dissipation that actually supports it.
### 32.9 Entropy as logarithmic coordinate of a named measure
Pressure, fresh service, high-tail scale, high-tail time, and causal ancestry all produce entropy-like quantities.
They are not one ontology.
They become meaningful only after the underlying physical measure is fixed.
### 32.10 Theorem depth is not recursion depth
A chain of certified source/pressure/service/shell consequence maps does not by itself create new physical generations.
Same-law relay duplication is quotiented, and different-unit witness relays create state but no second causal charge.
### 32.11 Physical time is universal; eventhood is not automatic
Events and checkpoints live in the same Navier--Stokes time coordinate, so physical corridor time telescopes without a synthetic clock. But a theorem-selected horizon does not become an event merely because time was spent reaching it. Bounded-scale checkpoint continuation reaches `t=0`; UV-growing checkpoint continuation can remain finite-time and must be treated as a separate PDE seam.
### 32.12 Circulation is not generation
The skew part of resolved low--high interface work is an antisymmetric finite-role flux. Internal role cycles cancel from every subset balance and remain at one physical event time. Quotient those cycles and donor relays before recursion; the symmetric remainder is the already existing strain/deformation work.
These simplifications are central to the programme's current direction.

---

### 32.13 Conservation of representation is not yet physics
A decomposition may conserve total energy while merely moving energy between observer-selected channels. Physical ownership begins only after the representation motion is tied to the actual PDE transport and that common gauge has been quotiented.
### 32.14 First stop is not necessarily first cause
A measurable coefficient threshold can correctly locate the earliest failure of clean continuation without itself carrying a physical work law. The master stores physical owners only after the required energy/work reentry.
### 32.15 Observable layer is not event depth
A physical observable may be discovered after traversing an interval without occurring after that interval in physical time. Full-natural own-scale service and its OO/ON/NN disintegration live on the corridor already counted. Preserve the measure; quotient only the duplicate theorem depth.
### 32.16 Horizon segmentation is not event topology
A natural-window endpoint is an analysis checkpoint unless a physical stop fires there. Preserve the real corridor time, but do not create a recursive vertex. At that checkpoint let the actual shell masses determine the hard-shell witness; never choose the upper cover branch to manufacture scale ascent.

## 33. Certified recent theorem blocks
### 33.1 High-tail sliding natural-window reentry
The latest completed high-tail bridge is the high-tail sliding natural-window reentry theorem.
Status:
`EXACT_HIGH_TAIL_COMPARABLE_HH_TO_SLIDING_NATURAL_WINDOW__TIME_ORIGIN_AND_UNIT_INVARIANT__SCALE_TIME_CONCENTRATION_TO_CRITICAL_SHELL__NO_PACKET_PERSISTENCE_OR_TIME_BINNING`.
Dedicated GitHub Actions run:
`31354438956`.
Exact theorem SHA:
`6c42d53f4f3903d141986d85a5b45954c27c18c8`.
Results:
- `641` tests passed;
- `50,000` sliding-measure / scale-time / shell-reentry states;
- worst time-origin invariance residual `7.861072903736499e-14`;
- worst time-unit invariance residual `1.3600232051658168e-14`;
- worst representation-refinement residual `1.5126788710517758e-15`;
- minimum scale-time tradeoff relative margin `0.942049661829323`;
- exact sharp boundary `M/N=2` attained in sampled geometry;
- exact sharp boundary `T_M/T_N=1/4` attained in sampled geometry;
- minimum conditional full-survivor service lower `9.309745536273158e-17`.
Stored artifact:
`recorded-results/31354438956/`.
Full causal integration:
`31354509984`.
It completed successfully on the same exact SHA through:
- source compiler;
- pressure pair atomization;
- fresh SGS reentry;
- high-frequency tail bridge;
- common-unit high-tail causality;
- ultraviolet locality;
- sliding natural-window reentry;
- event roles and material relay;
- physical productivity;
- recursive witness construction;
- Shannon reuse;
- Rényi reuse;
- physical branch compiler;
- master episode stress.
The initial candidate `ea9e89e...` had only two fixture-construction failures after `639` passing tests; no theorem stress ran there and no theorem equation changed in the correction.

### 33.2 Continuum master event quotient
Status:
`EXACT_CONTINUUM_MASTER_EVENT_QUOTIENT__ZERO_CHARGE_RELAYS_COLLAPSED__NATIVE_PHYSICAL_TIME_RECURSION__SUPPLIER_SPECIFIC_SCALE_PROGRESS__COMPACT_SCALE_FULL_SURVIVOR_NO_ESCAPE__NO_COMMON_CLOCK_OR_CAUSAL_REWEIGHTING`.
Dedicated GitHub Actions run:
`31369437763`.
Exact theorem SHA:
`37d43e189de6f8b0294cf7d18f2c672bcb419f87`.
Results:
- `652` tests passed;
- `50,000` quotient/path states;
- worst zero-charge owner-mass residual `0.0`;
- worst physical-time telescope residual `0.0`;
- worst log-scale telescope residual `8.881784197001252e-16`;
- bounded-scale boundary failures `0`;
- supplier-scale failures `0`;
- largest sampled relayed joint-owner set `8`.
Stored artifact:
`recorded-results/31369437763/`.
Full causal integration:
`31369437677`.
It completed successfully on the same exact SHA through all source, pressure, material, high-tail, role, joint-stop, causal-reuse, physical-branch, and master-episode stages.
The theorem does not assert global termination.  It certifies the quotient state and the compact-scale free-survivor no-escape dichotomy, reducing the remaining infinite-path problem to named non-free owner recurrence or UV-unbounded free survival.


### 33.3 Resolved interface donor/circulation quotient
Status:
`EXACT_RESOLVED_INTERFACE_DONOR_QUOTIENT__POSITIVE_INTERFACE_TO_CONSERVATIVE_SKEW_DONOR_OR_EXISTING_STRAIN__FINITE_SAME_EVENT_DONOR_EXHAUSTION__CIRCULATION_ZERO_RECURSION_DEPTH`.
Dedicated GitHub Actions run:
`31398210897`.
Exact theorem SHA:
`c51846914109abf9d881d0a4ef5545fa023677f3`.
Results:
- `660` tests passed;
- `50,000` split/flux/donor/high-tail states;
- worst signed `R=R_K+R_S` residual `8.881784197001252e-16`;
- minimum sampled positive-cover margin `-1.7763568394002505e-15` (floating roundoff around the exact analytic nonnegative cover);
- worst role-divergence residual `2.100962885439726e-15`;
- worst total skew-work residual `2.6968038962353156e-15`;
- worst donor-closure balance residual `2.7511442661221775e-15`;
- donor-existence failures `0`;
- largest sampled shortest donor path `4`;
- high-tail component failures `0`.
Stored artifact:
`recorded-results/31398210897/resolved-interface-donor-quotient-results/`.
Full causal integration:
`31398211279`.
It completed successfully on the same exact SHA through the full source, pressure, material, high-tail, resolved-interface, role, joint-stop, causal-reuse, physical-branch, and master-episode chain.
The initial implementation SHA `6ae71cba...` had one brittle certificate-wording assertion after `659` passing tests; no theorem identity or physical bound failed.  The assertion was relaxed to test semantics rather than exact prose, and the later notation-only `I_a^K -> R_a^K` cleanup changed no equation or routing.
The theorem closes resolved interface as an independent recursive-generation loophole.  It does not prove global termination of the donor/strain owners and makes no Navier--Stokes regularity claim.


---

### 33.4 Smooth quadratic-carrier observer-gauge quotient
Status:
`EXACT_SMOOTH_QUADRATIC_CARRIER_INTERFACE__Q2_ENERGY_LAW__COMMON_GAUGE_QUOTIENT_BEFORE_PHYSICAL_RELINK__SYMMETRIC_WORK_EXISTING_STRAIN__COEFFICIENT_OBSTRUCTION_ENERGY_REENTRY`.
Dedicated GitHub Actions run:
`31444417439`.
Exact implementation SHA:
`309d6ee83e39a96d8efd2a0ddcfcbf6839d9264b`.
Results:
- `676` tests passed;
- `50,000` transported square-partition/interface/PDE/reentry states;
- worst quadratic-partition residual `4.583749441582483e-15`;
- worst differentiated-partition residual `1.669620150224271e-14`;
- worst common-gauge transport residual `1.594436429147036e-16`;
- worst gauge-work cancellation residual `1.7157190746562474e-16`;
- worst native/outer recombination residual `6.010633403939076e-16`;
- worst physical-relink conservation residual `8.888567800392867e-16`;
- worst strain reconstruction residual `1.894411407492663e-15`;
- forbidden linear-complement counterexample defect `1.0`;
- arbitrary observer-motion rejection count `1`, with no admission;
- coefficient-obstruction master-barrier failures `0` in the companion `50,000`-state master stress.
Stored artifact:
`recorded-results/31444417439/smooth-quadratic-carrier-interface-results/`.
Full causal integration:
`31444417546`.
It completed successfully on the same exact SHA with the same `676`-test suite, all source/pressure/material/high-tail/master/role/reuse/compiler stages, and `20,000` master episode traces with worst margin `0.0`.
The theorem says neither that every genuine relink/strain/HH owner terminates nor that every UV-unbounded survivor is closed. It makes no Navier–Stokes regularity claim.


### 33.5 Full-natural service corridor quotient
Status:
`EXACT_FULL_NATURAL_SERVICE_CORRIDOR_QUOTIENT__OWN_SCALE_SERVICE_IS_SAME_INTERVAL_WITNESS_NOT_NEW_EVENT__MATERIAL_DISINTEGRATION_ZERO_RECURSION_DEPTH__ENDPOINT_SMOOTH_CARRIER_HAS_COMPARABLE_HARD_SHELL_WITNESS_SET`.
Dedicated GitHub Actions run:
`31448743219`.
Exact implementation SHA:
`e351d0d6bef5a6bd6275083e1d2e706acf717a18`.
Results:
- `684` tests passed;
- `50,000` corridor/service/material/endpoint-shell states;
- worst natural-time identity residual `4.440892098500626e-16`;
- worst OO/ON/NN same-measure partition residual `1.0658141036401503e-14`;
- minimum sampled two-hard-shell cover margin `2.840630292262874e-08`;
- exact joint hard-shell witness tie retained;
- companion continuum-master service-witness barrier failures `0`.
Stored artifact:
`recorded-results/31448743219/full-natural-service-corridor-quotient-results/`.
GitHub artifact digest:
`sha256:47a7296117a1dc966c47d325aecaab6754bb73bb9d9f8d59f5fc5b1dd49b0c5e`.
Full causal integration:
`31448763557`.
It completed successfully on the same exact SHA with the same `684`-test suite through source, pressure, material, critical-shell, service-corridor quotient, high-tail, hard/smooth interface, first-stop, causal-reuse, physical-branch, and master stages. The master checked `20,000` episode traces with worst margin `0.0`.
The initial SHA `4d19bf616c88e634698039c098e2c774be01a669` had only two brittle certificate-wording assertion failures after `682` passing tests; no theorem stress ran there and no theorem equation, physical bound, or routing changed in the correction.
The theorem removes service-theorem depth and closes endpoint-service attachment. It does not terminate genuine first-hit/work/source/reuse/independent-service recurrence or UV-unbounded full-survivor chains, and makes no Navier--Stokes regularity claim.


### 33.6 Full-natural horizon checkpoint quotient
Status:
`EXACT_FULL_NATURAL_HORIZON_CHECKPOINT_QUOTIENT__PHYSICAL_CORRIDOR_TIME_WITH_ZERO_EVENT_DEPTH__ENDPOINT_HARD_SHELL_COVER_IS_ANALYSIS_REREGISTRATION_NOT_SCALE_PROGRESS__UV_CHECKPOINT_CONTINUATION_SEPARATED_FROM_RECURSIVE_EVENT_PATH`.
Dedicated GitHub Actions run:
`31451492854`.
Exact implementation SHA:
`75ceff3481dccc41a9e915ce8c1400638e440820`.
Results:
- `696` tests passed;
- `50,000` checkpoint/corridor/cover states;
- worst physical-time telescope residual `0.0`;
- maximum sampled endpoint cover ratio `1.5000000000000002`;
- checkpoint-to-event failures `0`;
- cover-to-high-tail misclassification failures `0`;
- minimum sampled UV checkpoint time beyond the first corridor `8.874718028220728e-06`;
- companion `50,000`-state continuum-master checkpoint-barrier failures `0`;
- companion master physical-time telescope residual `0.0`;
- companion master log-scale telescope residual `8.881784197001252e-16`;
- the independent physical high-frequency dissipation dependency remained green.
Stored artifact:
`recorded-results/31451492854/full-natural-checkpoint-quotient-results/`.
GitHub artifact digest:
`sha256:212a85c8d90535c74dc4035b0e8372cfb4d51ca1b4282888a4852ed06d9e07cb`.
Full causal integration:
`31451492844`.
It completed successfully on the same exact SHA with the same `696`-test suite through the full physical-energy, source, material, generic-shell, service-corridor, checkpoint, high-tail, hard/smooth-interface, first-stop, causal-reuse, physical-branch and master chain. The final master checked `20,000` episode traces with worst margin `0.0`. Integration artifact digest:
`sha256:23cd2a39f9a19006008f4b29a99ef8dce2ee60b7f34952f4225e30b398027b89`.
The theorem removes natural-horizon segmentation and observer-selected cover ascent from recursive event topology. It does **not** prove that event-free UV checkpoint continuation forces high-tail work, and it makes no Navier--Stokes global-regularity claim.

---

## 34. Certification discipline
Serious theorem validation is performed only in GitHub Actions.
Local work is restricted to:
- reading;
- static inspection;
- editing;
- `git` operations;
- viewing/downloading CI results.
Do not run theorem stress, pytest, numerical experiments, or certificate execution locally.
For a new theorem:
1. branch from verified `origin/main`;
2. read the current ledger;
3. add source/test/doc/workflow wiring;
4. static-audit units, semantics, and anti-theorems locally;
5. push candidate;
6. run dedicated Actions;
7. classify failures as fixture/wiring versus theorem failure;
8. run full causal integration on the exact theorem SHA;
9. store the key artifact under `recorded-results/<runid>/`;
10. update history/current architecture only after Actions green;
11. run exact-SHA promotion validation;
12. fast-forward `main` only after the promotion sweep is green and `origin/main` is still an ancestor.
A passing numerical stress is not itself a proof.
The source module must encode the theorem algebra and guards; CI is regression/certificate evidence around that structure.

---

## 35. Recommended reading order for a new mathematical physicist
Start with this ledger.
Then read the following current-spine documents in order:
1. `docs/physical_energy_causal_bridge.md`
2. `docs/recursive_coherent_witness_extraction.md`
3. `docs/outer_moving_role_extraction.md`
4. `docs/smooth_quadratic_carrier_interface.md`
5. `docs/resolved_role_egorov.md`
6. `docs/nonaffine_role_interface_work.md`
7. `docs/resolved_interface_donor_quotient.md`
8. `docs/event_anchored_role_registration.md`
9. `docs/coherent_service_or_flat.md`
10. `docs/critical_shell_service_reentry.md`
11. `docs/full_natural_service_corridor_quotient.md`
12. `docs/full_natural_checkpoint_quotient.md`
13. `docs/high_strain_resolved_ancestor.md`
14. `docs/high_strain_heat_increment_service.md`
15. `docs/material_label_carrier_quotient.md`
16. `docs/objective_source_routing_compiler.md`
17. `docs/objective_pressure_pair_atomization.md`
18. `docs/fresh_service_scale_reentry.md`
19. `docs/high_frequency_dissipation_reentry.md`
20. `docs/high_tail_binary_work_reentry.md`
21. `docs/high_tail_ultraviolet_locality.md`
22. `docs/high_tail_natural_window_reentry.md`
23. `docs/physical_pair_weighted_productivity.md`
24. `docs/common_slice_coefficient_registration.md`
25. `docs/recursive_physical_witness_constructor.md`
26. `docs/joint_causal_stop_projection.md`
27. `docs/amplitude_entropy_causal_reuse.md`
28. `docs/weighted_causal_reuse.md`
29. `docs/renyi_causal_reuse.md`
30. `docs/physical_branch_compiler.md`
31. `docs/continuum_master_event_quotient.md`
32. `docs/master_no_escape.md`
For historical development, finite-dimensional precursor geometry, superseded packet formulations, and detailed CI chronology, read:
`docs/history/RESEARCH_LEDGER_history_through_2026-08-10.md`.

---

## 36. Compact one-page mental model
If only one picture is retained, use this one.
Navier–Stokes gives actual signed nonlinear work.
Positive child-energy work defines cause.
Coherent analysis refines that work into physical parent/child roles.
Smooth moving carriers transport roles between events, but their energy is `<u,Q^2u>`, not `<u,Qu>`.
Hard interaction projectors are event-anchored and are not differentiated through the slab; hard shell energy may also be reread at a zero-event checkpoint as a state observable.
Common affine/Kelvin role transport and cutoff repartition are gauges.
For a smooth square partition, require `dot A+[G,A]=0` and quotient that common observer motion before ownership.
Only residual skew `K_phys=K-G` is conservative physical relink; symmetric `S` is existing strain/deformation.
Conservative motion of observer-selected channels is not itself physics.
Hard resolved cross/interface work has its separate same-event donor quotient; the hard and smooth measures share provenance but are not identified.
Objective source variation is split by physical owner.
Resolved dissipation and pressure/fresh source laws produce actual critical hard shells.
Generic critical shells enter first-stop corridors at renewal scale `A=3M/4`. A full no-hit branch completes one real physical corridor which already carries its own-scale service as an attached witness; service/material rereading adds no event depth.
The natural-horizon endpoint is an analysis checkpoint, not an event. At that checkpoint the actual pair `(mu_A,mu_2A)` determines the unique lower/upper shell witness or exact joint tie. Ratios `3/4` and `3/2` are checkpoint-cover geometry, not directional or high-tail progress.
High strain produces critical resolved ancestors.
Fresh SGS service is quotiented to a scale law before shell extraction.
High-frequency service is converted to physical tail dissipation before any causal conclusion.
High-tail regeneration uses the common causal unit `N dW`.
Its HH part is already a binary physical work law.
Its UV nonlocality is paid by physical `D_tail`.
Its comparable remainder is localized in time by a sliding `M`-natural window.
That window produces an actual critical hard-shell event with `M/N>=2` and `T_M/T_N<=1/4`.
Material labels are attached as sidecars to actual service and do not duplicate carrier impulses.
Causal Shannon/Rényi reuse is reserved for actual positive child-work ancestry.
Source/service concentration entropies remain deterministic coordinates of their own measures.
All routed physical owners and simultaneous physical causes feed one joint recursive event master. Raw HH/interface coefficient threshold hits remain typed first-stop locators until actual `Q^2` energy/work reentry. Full-natural checkpoints and same-corridor service witnesses are explicitly excluded from that owner state.
Relay depth is quotiented before recursion depth is counted.
Actual physical time is the universal coordinate shared by events and checkpoints; each supplier keeps its own native corridor length, while eventhood is assigned only by physical stops/owners.
Bounded-scale no-hit checkpoint continuation hits `t=0`.
An infinite recursive **event** path avoiding `t=0` must therefore recur through infinitely many genuine named physical owners after all quotients. Separately, an event-free UV-growing checkpoint continuation can have finite physical duration and remains an open PDE seam.
The independently certified high-tail route applies only when its physical tail hypotheses are met; checkpoint cover ascent does not manufacture those hypotheses.
Only true globally bounded scale-independent resources may terminate additively.
`t=0` is absorbing.
The present frontier is twofold: telescope genuine event recurrence through native typed physical ledgers, and prove whether event-free UV checkpoint continuation necessarily enters actual tail work/dissipation or another native physical obstruction, without reintroducing observer-made structure.
