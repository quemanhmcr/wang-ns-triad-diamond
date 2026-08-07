# Helical spin / triad-normal transport

- exact local Berry curvature: `F_s=-s sin(theta) dtheta wedge dphi`
- exact Chern numbers: `c1(s=+1)=-2`, `c1(s=-1)=+2`
- single-triad normal gauge: coupling phase is quadrature (`+/- pi/2`) away from coupling zeros
- random 3D checks: `50000`
- worst SO(3)-covariance residual: `9.575e-16`
- worst spin-1 normal-transition residual: `9.813e-16`
- worst coupling real/absolute ratio: `2.714e-12`
- worst coupling-magnitude residual: `4.031e-15`
- worst Berry-connection finite-difference residual: `3.704e-09`
- worst transverse-strain/helicity-matrix residual: `1.415e-15`

A single triad's absolute Berry phase is not charged as a cascade defect: the
triad's moving normal supplies an SO(3)-covariant gauge in which rigid rotation is
phase-free and the Waleffe coupling phase is constant.  The physical geometric
obstruction appears only when one Fourier mode is reused by incident triads with
different normals; their transition function is the spin-1 dihedral phase
`exp(-i s psi)`.
