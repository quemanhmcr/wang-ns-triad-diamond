# Coherent transfer cells: nonlinear work without decomposition gauge

The Moyal ledger gives positive coherent energy, but the nonlinear Navier--Stokes work must also be assigned to those cells without choosing arbitrary Gaussian synthesis coefficients.  The polarized Moyal identity does exactly this.

For a normalized coherent window `g`, write `a=V_g f` and `b=V_g F`.  Then

\[
\boxed{\int a(z)\overline{b(z)}\,d\mu(z)=\langle f,F\rangle.}
\]

For any measurable phase-space partition `C_alpha`, define

\[
\boxed{\mathcal W_\alpha=2\Re\int_{C_\alpha}a\overline b\,d\mu.}
\]

Then

\[
\boxed{\sum_\alpha\mathcal W_\alpha=2\Re\langle f,F\rangle.}
\]

Thus the work of a projected nonlinear Navier--Stokes forcing, a source term, or a band-to-band transfer has an exact coherent-cell decomposition.  It is signed, as physical energy exchange must be; the positive cells provide a canonical positive service law and the negative cells measure actual backflow/cancellation rather than packet-gauge ambiguity.

## Common affine motion is not an interface

For a common affine map `M`, evolve

\[
L\mapsto ML,\qquad X\mapsto MX,\qquad k\mapsto M^{-T}k.
\]

The intrinsic coherent coordinate

\[
\zeta=(L^{-1}X/2,L^Tk)
\]

is unchanged.  The phase-space Jacobian is one.  Hence a coherent partition transported by the common affine flow has no artificial cell-boundary or relabeling term.  This is the coherent-state analogue of the zero moving-multiplier Heisenberg residual already proved for affine frequency cells.

## Relinking has a positive Moyal price

For cell energies `E_alpha>=0`, changing a selected material set from `S_-` to `S_+` obeys

\[
|E(S_+)-E(S_-)|\le E(S_+\triangle S_-).
\]

For a piecewise material selected family, let `P_+` and `P_-` be the integrated positive and negative coherent-cell work, and let `R_switch` be the sum of symmetric-difference Moyal energies at selection changes.  Exact energy balance gives

\[
\boxed{P_+\le E_{final}+P_-+R_{switch}.}
\]

Therefore

\[
\boxed{P_+>0\Longrightarrow
E_{final}\ge P_+/3\ \lor\ P_-\ge P_+/3\ \lor\ R_{switch}\ge P_+/3.}
\]

The three exits have the intended physical meanings: actual coherent reservoir energy, actual backflow/cancellation, or actual relinking mass.  Common affine motion does not contribute to `R_switch` because the cells themselves are transported covariantly.

This theorem does not yet identify the low-frequency increment reservoir selected by the SGS collision with one particular coherent cell.  That is supplied by the companion coherent increment-service theorem, where the increment itself generates a positive edge measure between nearby coherent cells.
