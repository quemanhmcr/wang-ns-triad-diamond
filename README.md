# Wang–NS Triad Diamond

A reproducible numerical experiment for a research programme inspired by the extremizer/rigidity methodology used in the Wang–Zahl Kakeya work.

The experiment studies a minimal Fourier interaction diamond

- `a + b -> m`
- `m + c -> d`
- `b + c -> n`
- `a + n -> d`

for helical Navier–Stokes triads. It asks whether all four forward-transfer edges can be simultaneously near-extremal while the gauge-invariant phase constraints close around the diamond.

This repository does **not** claim to solve Navier–Stokes. It is a falsifiable computational probe of one proposed rigidity mechanism.

## Reproduce

```bash
python -m pip install -r requirements.txt
pytest -q
python -m src.optimize_diamond --quick
```

The full GitHub Actions run performs a broader search and uploads JSON/Markdown artifacts.
