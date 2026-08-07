# Affine-grain polarization curvature

For `A(x)=A0+H[Lz]`, let `C_ijc=sym_ij(H_ijk L_kc)`.  Differentiated
incompressibility makes each matrix `C_c` trace free.  The certified extremal
relative-polarization tomography theorem applies to every `C_c`, and Gaussian
orthogonality gives

`E_z Q_rel(S(z)) = sum_c Q_rel(C_c) >= (1/2)||C||_F^2`.

- random affine Hessian checks: `50000`
- worst observed RMS observability ratio: `0.506683513`
- worst Monte-Carlo expectation residual: `6.124e-02`
- worst scalar third-Hermite symmetry residual for tested swirl kernels: `8.639e-16`
- minimum sampled swirl polarization signal / `||B||^2`: `7.985e-02`

Thus quadratic curvature has two distinct packet channels: full-symmetric
curvature creates third-Hermite envelope forcing, while physical symmetric
gradient variation creates transfer-distinguishable shape/polarization forcing.
The five-dimensional swirl kernel belongs to the second channel rather than the
first.  No claim is made here of an aspect-independent lower bound comparing the
second channel directly to the affine-normalized `||B||` norm.
