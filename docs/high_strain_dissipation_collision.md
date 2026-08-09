# High-strain lifetime is critical dissipation, not a finite reset count

The service-or-flat gate sends a packet out of the low-strain branch when

\[
K:=\int_0^T\|S(t)\|_{op}\,dt>1/30.
\]

It is tempting to call this a globally finite "strain budget". That would be wrong. The strict low-pass transporter gives a sharper and physically meaningful statement.

## 1. Exact band-limited Bernstein scaling

The transporter is

\[
V=S_{N/4}u,
\qquad
\operatorname{supp}\widehat V\subset B_{N/4}.
\]

With the unitary Fourier convention,

\[
\begin{aligned}
\|\nabla V\|_{L^\infty_x}
&\le (2\pi)^{-3/2}|B_{N/4}|^{1/2}\|\nabla V\|_2\\
&=\boxed{
\frac{N^{3/2}}{8\sqrt6\,\pi}\|\nabla V\|_2.
}
\end{aligned}
\]

Since \(\|S\|_{op}\le\|\nabla V\|_F\), on the natural lifetime \(T=cN^{-2}\),

\[
\begin{aligned}
K
&\le \frac{N^{3/2}}{8\sqrt6\pi}
\int_0^T\|\nabla V\|_2dt\\
&\le \frac{N^{3/2}\sqrt T}{8\sqrt6\pi}
\left(\int_0^T\|\nabla V\|_2^2dt\right)^{1/2}.
\end{aligned}
\]

Define the scale-critical normalized dissipation

\[
D_V=N\int_0^T\|\nabla V\|_2^2dt.
\]

All powers of \(N\) cancel:

\[
\boxed{
K\le\frac{\sqrt{cD_V}}{8\sqrt6\pi}.
}
\]

Therefore

\[
\boxed{
D_V\ge\frac{384\pi^2}{c}K^2.
}
\]

At the existing branch threshold,

\[
\boxed{
K>1/30
\Longrightarrow
D_V>\frac{32\pi^2}{75c}.
}
\]

So high strain is not an undefined action branch: it is a definite scale-critical resolved-dissipation event.

## 2. Why this does not give a global event count

This collision must **not** be inserted into the physical multi-currency master as a scale-independent additive reset budget.

If a geometric cascade has

\[
N_j=N_0q^j,
\qquad q>1,
\]

and every generation pays the same normalized dissipation \(D_{V,j}=d_0\), its actual viscous energy cost is

\[
\nu\int\|\nabla V_j\|_2^2dt
=\frac{\nu d_0}{N_j}.
\]

Hence

\[
\boxed{
\sum_{j\ge0}\frac{\nu d_0}{N_j}
=
\frac{\nu d_0}{N_0}\frac{q}{q-1}
<\infty.
}
\]

Exactly the same obstruction holds for fresh critical energy: if \(N_jE_j=\mu\),

\[
\boxed{
\sum_{j\ge0}E_j
=
\frac\mu{N_0}\frac{q}{q-1}<\infty.
}
\]

This is not a defect of the theorem. It is the Navier--Stokes scaling obstruction that motivated the sticky/branching programme in the first place.

## 3. Correct master interpretation

The physical multi-currency telescope is exact **only for resources whose event threshold is uniform in the globally bounded resource itself**. Critical mass and normalized dissipation do not satisfy that condition along a single geometric scale path.

Thus:

- `high_strain_lifetime` routes to positive critical resolved dissipation;
- it cannot be declared a finite-count additive reset merely from the energy inequality;
- the companion resolved-ancestor theorem disintegrates at least half of the actual `D_V` law onto lower-frequency shell-time atoms with fixed critical mass;
- the companion heat-increment theorem converts the same law, up to the support factor `e^(-1/32)`, into positive coherent spatial increment edges while retaining at least `e^(-1/32)/2` of the critical-ancestor mark;
- material ownership and renewed-slab attachment of those dissipation-seeded edges remain a continuum problem.

So the high-strain event is no longer an anonymous scalar recursion label, but it is still **not** a globally finite reset.  This prevents both false global counting and the equally false step of naming a shell/coherent edge as a selected transfer parent before the PDE supplies material renewal.
