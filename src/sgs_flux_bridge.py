from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


def wavevectors(n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    k = np.fft.fftfreq(n, d=1.0 / n)
    kx, ky, kz = np.meshgrid(k, k, k, indexing="ij")
    k2 = kx * kx + ky * ky + kz * kz
    return kx, ky, kz, k2


def project_div_free_hat(uhat: np.ndarray) -> np.ndarray:
    if uhat.ndim != 4 or uhat.shape[0] != 3 or not (uhat.shape[1] == uhat.shape[2] == uhat.shape[3]):
        raise ValueError("expected shape (3,n,n,n)")
    n = uhat.shape[1]
    kx, ky, kz, k2 = wavevectors(n)
    ks = (kx, ky, kz)
    dot = sum(ks[i] * uhat[i] for i in range(3))
    out = uhat.copy()
    mask = k2 > 0.0
    for i in range(3):
        corr = np.zeros_like(dot)
        corr[mask] = ks[i][mask] * dot[mask] / k2[mask]
        out[i] -= corr
    return out


def divergence_free_random_field(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(3, n, n, n))
    uhat = np.fft.fftn(raw, axes=(1, 2, 3))
    # Strict band-limit before forming quadratic products.  With coordinate
    # modes <=(n-1)//6, all cubic contractions used by the energy identity
    # remain below Nyquist, so the pseudospectral regression has no aliasing.
    kx, ky, kz, _ = wavevectors(n)
    kmax = max(1, (n - 1) // 6)
    mask = (np.abs(kx) <= kmax) & (np.abs(ky) <= kmax) & (np.abs(kz) <= kmax)
    uhat *= mask[None, ...]
    uhat = project_div_free_hat(uhat)
    return np.fft.ifftn(uhat, axes=(1, 2, 3)).real


def radial_filter_multiplier(n: int, cutoff: float) -> np.ndarray:
    if cutoff <= 0.0:
        raise ValueError("positive cutoff required")
    _, _, _, k2 = wavevectors(n)
    # Velocity filter G. The resolved quadratic energy weight is |G|^2.
    return np.exp(-0.5 * (k2 / (cutoff * cutoff)) ** 2)


def filter_scalar(f: np.ndarray, ghat: np.ndarray) -> np.ndarray:
    return np.fft.ifftn(np.fft.fftn(f) * ghat).real


def spectral_derivative(f: np.ndarray, axis: int) -> np.ndarray:
    n = f.shape[0]
    ks = wavevectors(n)[:3]
    return np.fft.ifftn(1j * ks[axis] * np.fft.fftn(f)).real


def sgs_flux_identity(u: np.ndarray, cutoff: float) -> dict[str, float]:
    """Compare mean SGS transfer with resolved nonlinear spectral work.

    Pi_SGS=-grad ubar : (overline{u u}-ubar ubar). On a periodic box,
    <Pi_SGS>=<ubar . div overline{u u}>. Leray projection of the nonlinear
    term gives the same mean work because ubar is divergence free.
    """
    if u.ndim != 4 or u.shape[0] != 3 or not (u.shape[1] == u.shape[2] == u.shape[3]):
        raise ValueError("expected shape (3,n,n,n)")
    n = u.shape[1]
    g = radial_filter_multiplier(n, cutoff)
    ubar = np.stack([filter_scalar(u[i], g) for i in range(3)])

    uu_bar = np.empty((3, 3, n, n, n), dtype=float)
    for i in range(3):
        for j in range(3):
            uu_bar[i, j] = filter_scalar(u[i] * u[j], g)
    tau = uu_bar - ubar[:, None] * ubar[None, :]

    grad = np.empty((3, 3, n, n, n), dtype=float)
    for i in range(3):
        for j in range(3):
            grad[i, j] = spectral_derivative(ubar[i], j)
    pi = -np.sum(grad * tau, axis=(0, 1))
    mean_sgs = float(np.mean(pi))

    nbar = np.zeros_like(ubar)
    for i in range(3):
        for j in range(3):
            nbar[i] += spectral_derivative(uu_bar[i, j], j)
    mean_resolved = float(np.mean(np.sum(ubar * nbar, axis=0)))

    nhat = np.fft.fftn(nbar, axes=(1, 2, 3))
    pnhat = project_div_free_hat(nhat)
    pn = np.fft.ifftn(pnhat, axes=(1, 2, 3)).real
    mean_projected = float(np.mean(np.sum(ubar * pn, axis=0)))

    div_hat = sum(wavevectors(n)[j] * np.fft.fftn(ubar[j]) for j in range(3))
    div_norm = float(np.linalg.norm(div_hat) / max(np.linalg.norm(np.fft.fftn(ubar, axes=(1,2,3))), 1e-300))
    return {
        "mean_sgs_flux": mean_sgs,
        "mean_resolved_nonlinear_work": mean_resolved,
        "mean_leray_projected_work": mean_projected,
        "sgs_vs_resolved_error": mean_sgs - mean_resolved,
        "projection_work_error": mean_projected - mean_resolved,
        "relative_divergence": div_norm,
    }


def stress(samples: int = 200, n: int = 12, seed: int = 20260807) -> dict[str, float]:
    worst_sgs = 0.0
    worst_proj = 0.0
    worst_div = 0.0
    scale = 0.0
    for i in range(samples):
        u = divergence_free_random_field(n, seed + i)
        cutoff = 0.9 + 1.8 * ((i * 0.61803398875) % 1.0)
        row = sgs_flux_identity(u, cutoff)
        worst_sgs = max(worst_sgs, abs(row["sgs_vs_resolved_error"]))
        worst_proj = max(worst_proj, abs(row["projection_work_error"]))
        worst_div = max(worst_div, row["relative_divergence"])
        scale = max(scale, abs(row["mean_sgs_flux"]), abs(row["mean_resolved_nonlinear_work"]), 1e-30)
    return {
        "samples": samples,
        "grid": n,
        "worst_sgs_identity_abs_error": worst_sgs,
        "worst_projection_work_abs_error": worst_proj,
        "worst_relative_divergence": worst_div,
        "reference_flux_scale": scale,
        "worst_sgs_identity_relative_error": worst_sgs / scale,
        "worst_projection_work_relative_error": worst_proj / scale,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=200)
    ap.add_argument("--grid", type=int, default=12)
    ap.add_argument("--outdir", type=Path, default=Path("results-sgs-flux"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples, args.grid)
    (args.outdir / "sgs_flux_bridge.json").write_text(json.dumps(out, indent=2))
    md = f'''# SGS / graded spectral flux regression

The equality is analytic; this pseudospectral run validates conventions.

- random divergence-free fields: `{out['samples']}` on `{out['grid']}^3`
- worst SGS-vs-resolved relative error: `{out['worst_sgs_identity_relative_error']:.3e}`
- worst Leray-projection work relative error: `{out['worst_projection_work_relative_error']:.3e}`
- worst relative divergence diagnostic: `{out['worst_relative_divergence']:.3e}`

The tested identity is `<Pi_SGS>=<ubar . div overline(u tensor u)>`; applying
the Leray projector to the nonlinear term leaves the global work unchanged.
'''
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
