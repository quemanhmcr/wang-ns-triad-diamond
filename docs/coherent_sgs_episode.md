# Source-weighted SGS to coherent ancestry without a persistence hypothesis

The coherent increment theorem is pointwise in the filter displacement/service event, while the H1 source theorem naturally supplies an **integrated source weight**.  The exponents match exactly, so no time-superlevel persistence is needed.

Let

\[
\rho_R(\tau)=N^{-4}\|S_R(\tau)\|,
\qquad
\Sigma_R=\int\rho_R(\tau)d\tau.
\]

The filtered SGS collision gives

\[
Q_{inc}(\tau)\ge c_Q\rho_R(\tau)^{3/2}.
\]

The coherent increment theorem defines

\[
Y(\tau)=\frac{(Q_{inc}/g_1)^{2/3}}{(C_{LP}C_B)^2}.
\]

Hence

\[
\boxed{Y(\tau)\ge c_Y\rho_R(\tau)}
\]

with an explicit `c_Y` built only from the already-certified SGS, Germano, LP and Bernstein constants.

After the existing large-radius versus scale-matched source-weight split, the scale-matched branch carries at least `Sigma_R/2`, so

\[
\int Y\,d\tau\ge c_Y\Sigma_R/2.
\]

At each time, either `d_high>=Y/4` or actual low coherent service is at least `Y/2`.  Pigeonholing by **Y-weight** gives

\[
\boxed{D_{high}\ge c_Y\Sigma_R/16}
\]

or

\[
\boxed{\int S_{low}\,d\tau\ge c_Y\Sigma_R/8.}
\]

If the material old pool has integrated service capacity

\[
cC_{old}\le c_Y\Sigma_R/32,
\]

then either

\[
\boxed{\int\Xi_{cell}\,d\tau\ge c_Y\Sigma_R/32}
\]

or

\[
\boxed{\int S_{new-new}\,d\tau\ge c_Y\Sigma_R/16.}
\]

A quarter-dominant new service edge yields

\[
\boxed{\int\mu_{coh,new}\,d\tau\ge c_Y\Sigma_R/128,}
\]

and therefore, on a scaled lifetime of length at most `c`,

\[
\boxed{\sup_\tau\mu_{coh,new}\ge c_Y\Sigma_R/(128c).}
\]

If the new service is fragmented instead, the service-edge probability law gives the existing `log 2` ancestry entropy or `1/4` same-ancestry cycle mass.

This is the source-weighted form of coherent reservoir synchronization: temporal concentration of the differentiated SGS source cannot evade phase-space ancestry, because the `3/2` source-to-increment power and `2/3` increment-to-service power cancel exactly.
