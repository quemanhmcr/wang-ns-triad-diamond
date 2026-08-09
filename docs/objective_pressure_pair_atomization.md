# Objective pressure Hessian: direct hard shell-pair atomization

## 1. Do not coarse-grain before reading the elliptic structure

The coherent corotational strain equation contains the actual averaged matrix source

\[
H_P(t)=-N^{-4}\langle \nabla^2P(t)\rangle_\gamma.
\]

The older bound

\[
\rho_P\le \mu_V/5700+\|R\|_{3/2}/380
\]

is valid, but the scalar `mu_V=N||V||_2^2` throws away the derivative suppression carried by the pressure Hessian.  It is therefore useful as a coarse diagnostic, not as the canonical renewal state.

The direct route keeps the tensor until the physical event.

## 2. Exact Frobenius dual: one scalar source law, no phase branch

Set

\[
Z(t)=\frac{H_P(t)}{\|H_P(t)\|_F}
\]

when `H_P!=0`, and `Z=0` at zero.  This is a measurable unit Frobenius dual.  The filtered pressure equation is linear in the resolved stress `V tensor V` and SGS stress `R`, hence

\[
H_P=H_V+H_R
\]

exactly and

\[
\rho_P:=\|H_P\|_F=Z:H_V+Z:H_R.
\]

Thus

\[
\boxed{
\rho_P\le [Z:H_R]_+ + [Z:H_V]_+.
}
\]

The positive parts here are **pressure source/service weights**.  They are not child-energy causal probabilities.

## 3. Hard only at the physical pressure event

Write the strict resolved transporter

\[
V=S_{N/4}u
\]

as an orthogonal hard dyadic shell decomposition at this event.  The nonzero shell frequencies satisfy

\[
M_a/2<|\xi|\le M_a,
\qquad M_a\le N/4.
\]

The zero mode, if present in a periodic representation, is pressure-null in the low--low cross term by incompressibility and does not create a pressure Hessian source.

Because the resolved pressure map is bilinear in `V tensor V`, its matrix source expands exactly into ordered shell pairs.  Orientation of a pair is not physical, so group them as **unordered atoms**:

\[
H_{aa}:=B(V_a,V_a),
\]

and, for `a<b`,

\[
H_{ab}:=B(V_a,V_b)+B(V_b,V_a).
\]

Then

\[
H_V=\sum_{a\le b}H_{ab}
\]

exactly.  With the same event dual,

\[
p_{ab}=Z:H_{ab},
\qquad
Z:H_V=\sum_{a\le b}p_{ab},
\]

so

\[
\boxed{
\rho_P\le [r_R]_+ + \sum_{a\le b}[p_{ab}]_+.
}
\]

No coherent localization is asked to preserve Fourier support.  Coherent/material labels may be attached only afterward as sidecars.  Before any entropy is read, identical unordered shell pairs are **quotiented and integrated once**.  A proof representation is not allowed to duplicate one physical pair into several records, because that would manufacture artificial collision entropy.  Likewise one hard shell label carries one physical frequency.

## 4. Pair capacity from already-certified constants

For one ordered pair with `M=max(M_a,M_b)`, the product is supported in the ball of radius at most `2M`.  The certified objective-source estimate says that at support radius `N/2`, the normalized order-two `L^(3/2)->L^infinity` constant is `<1/380`.  Exact homogeneity therefore gives

\[
N^{-4}\|D^2 B(V_a,V_b)\|_\infty
<\frac{256}{380}\left(\frac MN\right)^4
\|V_a\|_3\|V_b\|_3.
\]

The certified resolved-shell interpolation constant `<1/15` at support `(1/4)N'`, applied with `N'=4M_a`, gives

\[
\|V_a\|_3^2<\frac4{15}\mu_a,
\qquad
\mu_a=M_a\|V_a\|_2^2.
\]

Hence one ordered pair obeys

\[
|p_{a\to b}|
<\frac{256}{1425}
\left(\frac MN\right)^4
\sqrt{\mu_a\mu_b}
<\frac15
\left(\frac MN\right)^4
\sqrt{\mu_a\mu_b}.
\]

Grouping the two orientations yields the clean physical unordered bound

\[
\boxed{
|p_{ab}|
\le
\frac{\kappa_{ab}}5
\left(\frac{M_{\max}}N\right)^4
\sqrt{\mu_a\mu_b},
}
\]

where `kappa_aa=1` and `kappa_ab=2` for `a!=b`.

The coherent Gaussian average is a probability average, so it cannot increase the global `L^infinity` capacity and adds no constant.

The countable low-frequency shell tail causes no hidden divergence.  For the canonical hard dyadic family `M_j=(N/4)2^-j`, Cauchy gives

\[
\sum_j\sqrt{\mu_j}
\le
\left(\sum_jM_j\right)^{1/2}
\left(\sum_jE_j\right)^{1/2}
=\left(\frac N2\|V\|_2^2\right)^{1/2}.
\]

Using `M_max/N<=1/4` in the clean unordered capacity therefore yields

\[
\boxed{
\sum_{a\le b}\operatorname{cap}_{ab}
\le\frac{N\|V\|_2^2}{2560}.
}
\]

Thus the infinite pair expansion is absolutely source-summable.  Aggregate resolved energy appears here only as a convergence budget; it is still not the renewal state.

## 5. Positive pressure source has only two native owners

Let

\[
\Sigma_P=\int \rho_P\,d\tau,
\qquad d\tau=N^2dt,
\]

on the source episode.  Define actual positive integrated weights

\[
R_{SGS}=\int[r_R]_+d\tau,
\qquad
R_{ab}=\int[p_{ab}]_+d\tau.
\]

The exact dual law gives

\[
\Sigma_P\le R_{SGS}+\sum_{a\le b}R_{ab}.
\]

Therefore, with exact ties retained jointly,

\[
\boxed{
R_{SGS}\ge\Sigma_P/2
\quad\text{or}\quad
R_{pair}:=\sum_{a\le b}R_{ab}\ge\Sigma_P/2.
}
\]

The SGS owner satisfies

\[
\int\|R\|_{3/2}d\tau\ge190\Sigma_P
\]

and is exactly the effective objective-SGS source weight `Sigma_P/2` already routed to coherent square service.

## 6. Dominant pair gives a genuine critical shell

Normalize the positive pair law by

\[
q_{ab}=R_{ab}/R_{pair}.
\]

Take the physical dominance threshold `theta=1/4`.  If an unordered pair has `q_ab>=theta`, then on the pair owner branch

\[
R_{ab}\ge \frac{\theta}{2}\Sigma_P.
\]

The pair capacity implies

\[
\int\sqrt{\mu_a\mu_b}\,d\tau
\ge
\frac5{\kappa_{ab}}
\left(\frac N{M_{max}}\right)^4R_{ab}.
\]

If the scaled source horizon has length `c`, at some actual time

\[
\sqrt{\mu_a\mu_b}
\ge
\frac5{\kappa_{ab}c}
\left(\frac N{M_{max}}\right)^4R_{ab}.
\]

At least one of the two hard shells of the resolved transporter has critical mass at least this geometric mean.  This is not yet an identification with the full-velocity shell.  The canonical scalar low-pass satisfies `|S_(N/4)(xi)|<=1` and commutes with the hard shell projector, so exactly

\[
M\|P_MV\|_2^2\le M\|P_Mu\|_2^2.
\]

Therefore the same numerical lower is inherited by the corresponding hard shell of `u`, with no inverse low-pass estimate.  Uniformly over diagonal/off-diagonal pairs,

\[
\boxed{
\mu_{child}
\ge
\frac{5\theta}{4c}\Sigma_P
\left(\frac N{M_{max}}\right)^4.
}
\]

Since every resolved pair has `M_max<=N/4`, the canonical quarter threshold gives

\[
\boxed{
\mu_{child}\ge80\,\Sigma_P/c.
}
\]

This is an actual hard-shell event and enters the generic critical-shell first-stop/service theorem.  The child frequency is at most `N/4`, so its natural lifetime is at least `16` times the parent block lifetime scale.

## 7. Diffuse pair law is source entropy, not causal probability

If no unordered pair exceeds one quarter, then

\[
\sum_{a\le b}q_{ab}^2
\le
\max q_{ab}
\le\frac14.
\]

Therefore

\[
\boxed{H_2(q)\ge\log4.}
\]

This is the collision entropy of the **actual positive pressure source law**.  It is already a quantitative fragmentation certificate.  If ancestry labels are later attached, the exact atomic-to-component chain rule may further split it into component entropy or same-ancestry pair mass; however the weighted cycle/master conversion is a separate theorem interface.  Therefore this result does **not** by itself declare a terminal `TRANSFER_COST`, and it must not be relabeled a Shannon/Renyi child-energy causal law.

At exact `q_max=1/4`, both the dominant-pair shell route and the entropy route are valid and are kept jointly.

## 8. Material reuse is a sidecar refinement

The shell pair above is a hard physical event label, not a material ancestry label.  If a later argument attaches a material/coherent pair sidecar and supplies a signed-good low-strain lineage, the derivative-correct objective-Hessian reuse theorem gives the fixed-pair contraction

\[
\left(\frac{21}{20}\right)^5\left(\frac58\right)^4<\frac15,
\]

with total future fixed-pair capacity `<5/4` of generation zero.

The material-label quotient says a pure sidecar change need not destroy an unchanged smooth carrier.  None of this is needed to prove the dominant hard-shell event.

## 9. Architectural consequence

The canonical pressure route becomes

\[
\boxed{
\text{actual averaged pressure Hessian}
\to
\text{SGS positive source}
\quad\text{or}\quad
\text{resolved positive hard-pair law}
}
\]

and on the resolved law

\[
\boxed{
\text{dominant unordered pair}\to\text{critical shell},
\qquad
\text{diffuse pairs}\to H_2^{source}\ge\log4
\quad\text{(source fragmentation)}.
}
\]

Aggregate `mu_V` remains a correct coarse inequality but is no longer the physical state variable for pressure renewal.  The diffuse branch remains recursive until its independent ancestry/component/master conversion is supplied.  No packet synchronization, no coherent-frequency fiction and no global-regularity conclusion are introduced.
