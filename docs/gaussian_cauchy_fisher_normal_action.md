# Gaussian–Cauchy–Fisher normal-action checkpoint

Status: **analytic research checkpoint, not a Navier–Stokes regularity proof**.

This note records a new physical/geometric route developed after the existing triad/diamond program. The guiding rule is to quotient the motions that the PDE itself treats as geometry — pressure constraint, translation, incompressible affine deformation, covariance transport, and Brownian ancestry — before estimating anything. Once this is done, the remaining nonlinear object is genuine non-affinity, and several apparently unrelated positive structures collapse onto the same residual.

## 1. Gaussian affine projection and exact scale sink

For `R ~ N(0,2a I)`, define

\[
U_a=e^{a\Delta}u,\qquad A_a=\nabla U_a,
\]

and the orthogonal affine residual

\[
N_a(R)=u(x+R)-U_a-A_aR,
\]

with

\[
\mathbb E N_a=0,\qquad \mathbb E[N_a\otimes R]=0.
\]

The subscale energy is

\[
e_a=\frac12\big((|u|^2)_a-|U_a|^2\big)
=a|A_a|^2+\frac12\mathbb E|N_a|^2.
\]

Using incompressibility,

\[
|A_a|^2=|S_a|^2+\frac12|\omega_a|^2,
\]

so the scale-critical density

\[
\Phi_a:=\frac{e_a}{a}
\]

has the exact positive decomposition

\[
\boxed{
\Phi_a
=|S_a|^2+\frac12|\omega_a|^2
+\frac{\mathbb E|N_a|^2}{2a}.
}
\]

The exact Gaussian scale law is

\[
\boxed{
(\partial_a-\Delta)\Phi_a
=-\frac{\mathbb E|N_a|^2}{2a^2}.
}
\]

Thus the true non-affinity component

\[
\mathfrak n_a:=\frac{\mathbb E|N_a|^2}{2a}
\]

is itself the positive scale sink after multiplication by `1/a`.

The same structure extends to full covariance `C`. For `R ~ N(0,C)`, with `U_C`, `A_C`, residual `N_C`, and residual covariance

\[
\mathcal R_C=\mathbb E[N_C\otimes N_C]\succeq0,
\]

one has

\[
(u\otimes u)_C-U_C\otimes U_C=A_CCA_C^T+\mathcal R_C.
\]

On every covariance ray `C=\lambda C_0`,

\[
\boxed{
\left(\partial_\lambda-\frac12 C_0:\nabla^2\right)
\left(\frac{R_{0,\lambda C_0}}{\lambda}\right)
=-\frac{\mathcal R_{\lambda C_0}}{\lambda^2}\preceq0.
}
\]

## 2. Vorticity area flux is already inside the affine Gaussian reservoir

For arbitrary covariance `C`, let `W_C=curl U_C` and `B=C^{1/2}`. Gaussian Stein algebra gives

\[
\boxed{
W_C^T\operatorname{cof}(C)W_C
=2\left|\operatorname{skew}(B A_C B)\right|_F^2.
}
\]

Hence the scale-critical vorticity flux through the Gaussian ancestry area is carried entirely by the affine Hermite mode. The true non-affine residual contributes no direct term at this rung.

Moreover,

\[
\boxed{
\operatorname{tr}(C\tau_C)
=\frac12W_C^T\operatorname{cof}(C)W_C
+\left|\operatorname{sym}(B A_C B)\right|_F^2
+\operatorname{tr}(C\mathcal R_C),
}
\]

where `\tau_C` is the full subscale covariance. Thus rotational area flux, affine strain, and true non-affinity are three orthogonal positive reservoirs of one Gaussian scale quantity.

A useful rigorous dichotomy follows. For any `0<\theta<1`, define

\[
\mathcal B_\theta
=\{(x,t,a):\mathfrak n_a\ge\theta\Phi_a\}.
\]

Then energy dissipation implies the spacetime–log-scale budget

\[
\boxed{
\int_{\mathcal B_\theta}\Phi_a\,dx\,dt\,\frac{da}{a}
\le \frac{K_0(0)}{\theta\nu}.
}
\]

On the complement the Gaussian velocity is quantitatively near affine:

\[
\boxed{
\frac{\mathbb E|N_a|^2}{\mathbb E|A_aR|^2}
<\frac{\theta}{1-\theta}.
}
\]

So persistent bad behavior across many logarithmic scales must either repeatedly pay the non-affinity sink or become asymptotically affine on the relevant Gaussian clouds.

## 3. Stochastic Cauchy geometry and 2-form area precision

For backward stochastic ancestry

\[
dY_s=-u(Y_s,T-s)\,ds+\sqrt{2\nu}\,dB_s,
\]

let `J_s=\nabla_yY_s`, `K_s=J_s^{-1}`. Incompressibility gives `det J_s=1`.

The Cauchy-transformed vorticity

\[
q_s=K_s\omega(Y_s,T-s)
\]

satisfies the exact martingale equation

\[
\boxed{
dq_s=\sqrt{2\nu}\,K_s\nabla\omega_s\,dB_s.
}
\]

Thus vortex stretching disappears from the drift.

Define accumulated material Malliavin covariance

\[
\mathcal G_s=2\nu\int_0^sK_rK_r^T\,dr,
\]

ordinary precision `P_s=\mathcal G_s^{-1}`, and the natural precision for a 2-form

\[
\boxed{
\Pi_s=\operatorname{cof}(P_s)=\frac{\mathcal G_s}{\det\mathcal G_s}.
}
\]

Then

\[
\boxed{\Pi_s'\preceq0}
\]

pathwise. More explicitly, with `B=\mathcal G^{-1/2}KK^T\mathcal G^{-1/2}=\sum b_\alpha b_\alpha^T`,

\[
-q^T\Pi' q
=\frac{2\nu}{\det\mathcal G}
\sum_\alpha|b_\alpha\times \mathcal G^{1/2}q|^2\ge0.
\]

This is the correct area-resolution Fisher sink for vorticity viewed as a transported 2-form.

## 4. Future-variance reservoir closes the vertical next-jet producer

Fix a terminal horizon `S` and define the conditional terminal second moment

\[
Z_s=\mathbb E[q_Sq_S^T\mid\mathcal F_s]
=q_sq_s^T+W_s,
\]

where

\[
W_s=\operatorname{Cov}(q_S\mid\mathcal F_s)\succeq0.
\]

The positive quadratic variation of `q_sq_s^T` is exactly the negative drift of `W_s`. Therefore the completed vertical functional

\[
\boxed{
\mathfrak V_s=\operatorname{tr}(\Pi_sZ_s)
}
\]

satisfies

\[
\boxed{
d\mathfrak V_s
+[-\operatorname{tr}(\Pi_s'Z_s)]\,ds=dM_s,
}
\]

so `\mathfrak V_s` is a supermartingale. The apparent `+\nabla\omega` producer is not a new source; it is internal transfer between coherent Cauchy vorticity and unresolved future variance.

In the affine sector, `W_s` is exactly the Gaussian/Hermite subscale covariance of pulled-back terminal vorticity. Hence the martingale future-variance ladder and the Gaussian Hermite ladder are the same reservoir viewed from opposite time directions.

The dual scale-critical area-flux tensor is

\[
A_s=\operatorname{cof}(\mathcal G_s),\qquad \Pi_s=A_s^{-1},
\]

and for every `Z\succeq0`,

\[
\boxed{
(\operatorname{tr}(A Z))(\operatorname{tr}(A^{-1}Z))
\ge(\operatorname{tr}Z)^2.
}
\]

The precision side is dissipative; the flux side is scale-critical and, by the Gaussian identities above, already belongs to the affine Gaussian scale reservoir.

## 5. Exact affine-quotiented ancestry normal form

Choose a self-consistent affine Gaussian reference centered at `X_s` with covariance `C_s`:

\[
\dot X_s=-U_{C_s}(X_s,T-s),
\]

\[
\dot C_s=-A_sC_s-C_sA_s^T+2\nu I,
\qquad A_s=\nabla U_{C_s}(X_s,T-s).
\]

Let `J_s` solve `\dot J_s=-A_sJ_s`, `K_s=J_s^{-1}`, and define the material coordinate

\[
Z_s=K_s(Y_s-X_s).
\]

With the Gaussian affine residual

\[
N_s(r)=u(X_s+r,T-s)-U_s-A_sr,
\]

and `\widehat N_s(z)=K_sN_s(J_sz)`, the true stochastic ancestry obeys the exact normal form

\[
\boxed{
dZ_s=-\widehat N_s(Z_s)\,ds+\sqrt{2\nu}\,K_s\,dB_s.
}
\]

The affine Gaussian reference obeys

\[
\boxed{
d\bar Z_s=\sqrt{2\nu}\,K_s\,dB_s.
}
\]

All affine stretching has disappeared. The only dynamical difference is genuine non-affinity `\widehat N`.

The reference covariance is

\[
G_s=2\nu\int_0^sK_rK_r^Tdr,
\]

and physical covariance is exactly

\[
\boxed{C_s=J_sG_sJ_s^T.}
\]

Thus Gaussian covariance and Malliavin covariance are the same object in physical and material frames.

## 6. Non-affinity is path-space information action

The Gaussian affine projection is the exact orthogonal best affine incompressible drift. For any affine incompressible `b+Br`,

\[
\boxed{
\mathbb E_C|u-b-Br|^2
=\mathbb E_C|N_C|^2
+|b-U_C|^2
+\operatorname{tr}[(B-A_C)C(B-A_C)^T].
}
\]

Let `\mathbb P` be the true path law of `Z` and `\mathbb Q` the Gaussian reference law of `\bar Z`. Girsanov gives

\[
\boxed{
D_{KL}(\mathbb P\|\mathbb Q)
=\frac12\mathbb E_{\mathbb P}\int
\widehat N^T\dot G^{-1}\widehat N\,ds
=\frac1{4\nu}\mathbb E_{\mathbb P}\int|N|^2ds.
}
\]

The reverse orientation is especially useful because it is purely Gaussian:

\[
\boxed{
D_{KL}(\mathbb Q\|\mathbb P)
=\frac1{4\nu}\int\mathbb E_{R\sim N(0,C_s)}|N_s(R)|^2ds.
}
\]

Thus the same non-affinity residual is simultaneously

- Gaussian affine projection error,
- positive Gaussian scale sink,
- and path-space relative-entropy rate away from affine-Gaussian dynamics.

If `p_s` and `\gamma_s=N(0,G_s)` are the true and reference material densities, with `r=\nabla\log(p/\gamma)`, then

\[
\boxed{
D_{KL}(\mathbb P\|\mathbb Q)
=D_{KL}(p_S\|\gamma_S)
+\frac12\int_0^S
\mathbb E_p\left|r+\dot G^{-1}\widehat N\right|_{\dot G}^2ds.
}
\]

Path non-affinity action therefore splits exactly into terminal non-Gaussianity plus an irreversible Fisher mismatch square.

## 7. Rigorous epsilon-stability without Gronwall in the affine strain

Couple the true and affine reference ancestries with the same Brownian motion. For `\delta_s=Y_s-\bar Y_s`,

\[
\boxed{
\frac d{ds}(K_s\delta_s)=-K_sN_s(Y_s-X_s).
}
\]

All affine stretching and all Brownian noise cancel exactly.

Using `G_S=2\nu\int_0^SK_sK_s^Tds`, the control-Gramian Cauchy–Schwarz inequality gives

\[
\boxed{
|K_S\delta_S|_{G_S^{-1}}^2
\le\frac1{2\nu}\int_0^S|N_s(Y_s-X_s)|^2ds.
}
\]

Hence

\[
\boxed{
\mathbb E_{\mathbb P}|K_S(Y_S-\bar Y_S)|_{G_S^{-1}}^2
\le2D_{KL}(\mathbb P\|\mathbb Q).
}
\]

Pinsker/data processing also gives terminal distributional stability from either path KL orientation. This is an affine-invariant perturbation theorem: the size of `A_s` never appears.

## 8. The next rung: deformation stability begins only at grad N

Position stability alone is insufficient for Cauchy vorticity because the true and affine Jacobians can differ. If `J_s^{NS}` is the true Jacobian and `J_s` the affine reference Jacobian, define

\[
H_s=K_sJ_s^{NS}.
\]

Then exactly

\[
\boxed{
\dot H_s
=-\mathcal E_sH_s,
\qquad
\mathcal E_s
=K_s(\nabla N_s)(Y_s-X_s)J_s.
}
\]

Equivalently, in material coordinates,

\[
\boxed{
\dot H_s=-(\nabla_z\widehat N_s)(Z_s)H_s.
}
\]

Again the affine gradient `A_s` cancels completely. Failure of affine Cauchy cancellation begins only at the **next rung**, `\nabla\widehat N`.

This same next rung already appears in the Gaussian/Hermite program:

- residual covariance production is generated by spatial variation of the affine deformation;
- Hermite feedback is nearest-neighbour;
- Gaussian Hodge spectral cost of a residual starting at chaos degree `2` controls its first derivative;
- second stochastic/Malliavin variation begins at `\nabla A=\nabla^2u`.

No new obstruction has been introduced.

## 9. Critical log-scale interpretation

In a Type-I near-affine regime with `|A|~1/s` and viscous scale `a~\nu s`, define the relative Gaussian non-affinity ratio

\[
\theta(s)=\frac{\mathbb E|N|^2}{\mathbb E|AR|^2}.
\]

Since `\mathbb E|AR|^2~\nu/s`, the path-information rate behaves like

\[
\frac1{4\nu}\mathbb E|N|^2
\sim\frac{\theta(s)}{s}.
\]

Thus

\[
\boxed{
\mathscr A_0\sim\int\theta(s)\,d\log(1/s).
}
\]

Uniformly small relative non-affinity over infinitely many dyadic scales is not enough; to remain genuinely affine-Gaussian in path information, the error must be **log-summable** toward a hypothetical singular endpoint.

This matches the exact spacetime–log-scale Gaussian sink budget. A persistent singular cascade is forced toward a stronger asymptotic-affinity condition, not merely a fixed small perturbation.

## 10. Current frontier

The nonlinear algebra has been reduced to a concentration/stability problem rather than a vortex-stretching estimate.

Define the rung-0 normal affine action along an ancestry path

\[
\boxed{
\mathfrak A_0
=\frac12\int\widehat N^T\dot G^{-1}\widehat N\,ds.
}
\]

Rung 0 is already controlled geometrically: small `\mathfrak A_0` gives path-law and trajectory closeness to the affine Gaussian reference without any Gronwall factor in the strain.

The decisive next object is the rung-1 deformation action

\[
\boxed{
\mathfrak A_1
\sim\int\|\nabla\widehat N\|_{\mathrm{natural\ covariance}}^2ds.
}
\]

The next target is an **epsilon-stability theorem for the relative deformation** `H` in the same covariance/Malliavin geometry, with error measured by the existing Gaussian/Hermite next-rung sink rather than by an external Sobolev norm.

The desired architecture is:

1. non-affine spacetime–log-scale intervals are paid for by the exact Gaussian scale sink;
2. on near-affine intervals, `\mathfrak A_0` controls ancestry/path departure;
3. `\mathfrak A_1` controls Jacobian/Cauchy departure;
4. the vertical Cauchy area-precision functional is a completed supermartingale after adding future variance;
5. therefore a hypothetical singular ancestry would have to concentrate the globally finite next-rung dissipation along one exceptional backward thread through infinitely many scales.

The unresolved problem is now sharply formulated:

\[
\boxed{
\text{Can a single stochastic ancestry repeatedly sample the exceptional part of a globally finite Gaussian/Hermite scale-dissipation measure strongly enough to make }\mathfrak A_1\text{ critical?}
}
\]

A successful answer should be formulated as a covariance-adapted Carleson/stopping-time or path-capacity statement, not as a crude pointwise Sobolev inequality. Until that bridge is proved, there is **no 3D Navier–Stokes global-regularity proof**.
