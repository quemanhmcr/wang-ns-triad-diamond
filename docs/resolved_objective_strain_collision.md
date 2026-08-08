# Resolved objective strain: the pressure-Hessian near-field gap disappears

The affine/Kelvin packet is transported by the strict low-pass field

\[
V=S_{N/4}u.
\]

Therefore the objective strain in the service-or-flat theorem should be the strain of this **resolved transporter**, not an unrelated raw full-velocity strain. This changes the source calculus in an important way: pressure and SGS forcing remain filtered objects.

## 1. Exact filtered velocity-gradient identity

The resolved equation is

\[
D_t^V V
=-\nabla P-\nabla\!\cdot R+\nu\Delta V,
\]

with

\[
R=S_{N/4}(u\otimes u)-V\otimes V.
\]

Let

\[
A=\nabla V=S+\Omega.
\]

Differentiating the resolved equation,

\[
D_t^V A
=-A^2-\nabla^2P-\nabla\nabla\!\cdot R+\nu\Delta A.
\]

In the frame corotating with \(\Omega\),

\[
\boxed{
\mathring S_V
=-S^2-\Omega^2+[S,\Omega]
-\nabla^2P
-\operatorname{sym}\nabla\nabla\!\cdot R
+\nu\Delta S.
}
\]

This is the correct source identity for the strain that actually drives the affine/Kelvin packet frame.

## 2. All nonlocal sources are band limited

The strict filter gives

\[
\operatorname{supp}\widehat V\subset B_{N/4},
\qquad
\operatorname{supp}\widehat R\subset B_{N/2}.
\]

Taking divergence of the resolved equation gives

\[
-\Delta P=\partial_i\partial_j(V_iV_j+R_{ij}),
\]

hence

\[
\operatorname{supp}\widehat P\subset B_{N/2}.
\]

For a Hilbert-valued field supported in \(B_{N/2}\), unitary Fourier Hausdorff--Young followed by Hölder gives

\[
\|D^2f\|_\infty
\le C_2N^4\|f\|_{3/2},
\]

where Arb certifies

\[
\boxed{C_2<1/380.}
\]

For \(V\), L2/L∞ interpolation and the \(N/4\) support give

\[
\boxed{
\|V\|_3^2<\frac1{15}\,N\|V\|_2^2
=\frac1{15}\mu_V.
}
\]

## 3. Four physical source collisions

Write

\[
\mathfrak d_V=N^{-1}\|\nabla V\|_2^2.
\]

### Quadratic stretching/vorticity

Since \(\|Q(A)\|\le4\|A\|^2\) and

\[
\|\nabla V\|_\infty
\le\frac{N^{3/2}}{8\sqrt6\pi}\|\nabla V\|_2,
\]

one obtains

\[
\boxed{
N^{-4}\|Q(A)\|_\infty
\le\frac{\mathfrak d_V}{96\pi^2}.
}
\]

Thus source level \(\rho_Q\) forces

\[
\boxed{\mathfrak d_V\ge96\pi^2\rho_Q.}
\]

### Resolved SGS strain source

The symbol \(R\mapsto\operatorname{sym}\nabla\nabla\cdot R\) has Hilbert operator norm at most \(|\xi|^2\). Therefore

\[
\boxed{
N^{-4}\|\operatorname{sym}\nabla\nabla\cdot R\|_\infty
\le\frac1{380}\|R\|_{3/2}.
}
\]

Hence

\[
\boxed{
\|R\|_{3/2}\ge380\rho_{R,2}.
}
\]

The Germano increment inequality and Onsager/coherent collision then route this source to coherent low/base service, ancestry/entropy/cycles, or high-frequency normalized enstrophy.

### Filtered pressure Hessian

The Fourier symbol mapping the tensor \(V\otimes V+R\) directly to \(\nabla^2P\) has Hilbert operator norm at most \(|\xi|^2\). No separate Riesz-transform constant is needed. Thus

\[
\begin{aligned}
N^{-4}\|\nabla^2P\|_\infty
&\le\frac1{380}
\left(\|V\|_3^2+\|R\|_{3/2}\right)\\
&\le
\boxed{
\frac{\mu_V}{5700}
+
\frac{\|R\|_{3/2}}{380}.
}
\end{aligned}
\]

Therefore a pressure-Hessian source level \(\rho_P\) forces

\[
\boxed{
\mu_V\ge2850\rho_P
}
\]

or

\[
\boxed{
\|R\|_{3/2}\ge190\rho_P.
}
\]

This removes the old unspecified local pressure-Hessian coefficient from the **resolved transporter** source calculus. The first branch is a coherent mass/reservoir event; the second is an SGS increment/service event.

### Viscosity

Using the gradient as the L2 input,

\[
\|D^3V\|_\infty
\le C_{31}N^{7/2}\|\nabla V\|_2,
\]

and Arb certifies

\[
C_{31}<1/1500.
\]

Hence

\[
\boxed{
N^{-4}\nu\|\Delta S\|_\infty
\le\frac\nu{1500}\sqrt{\mathfrak d_V},
}
\]

so

\[
\boxed{
\mathfrak d_V
\ge
\left(\frac{1500\rho_\nu}{\nu}\right)^2.
}
\]

## 4. Source weight, not persistence

Let

\[
A_{obj}=T\int\|\mathring S_V\|dt,
\qquad
T=cN^{-2}.
\]

In scaled time \(\tau=N^2t\), one of the four source channels carries

\[
\boxed{
\Sigma_*
:=\int\rho_*d\tau
\ge\frac{A_{obj}}{4c}.
}
\]

The quadratic and SGS channels are linear in source density after collision; the viscous channel becomes stronger under temporal concentration by Cauchy. Pressure splits into coherent resolved mass or SGS service. Thus a large objective-strain variation of the actual transporter is no longer an unstructured source action.

## 5. Scope and correction

The older full-velocity identity

\[
\mathring S_u=Q(\nabla u)-\nabla^2p+\nu\Delta S_u
\]

remains mathematically correct. Its local pressure-Hessian packet coefficient was not closed. But the service-or-flat affine packet is driven by \(V\), so using the full-\(u\) strain there was unnecessarily strong and mismatched to the actual transport geometry.

For the resolved Kelvin transporter, the pressure-Hessian near-field gap is replaced by the band-limited mass/SGS alternative above. This does **not** convert critical mass or normalized dissipation into a scale-independent global event count; those currencies still require sticky causal ancestry or a scale-sensitive telescope.
