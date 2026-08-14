import pytest

from src.native_stock_no_mint import StockStep, stock_telescope, theorem_certificate


def test_physical_stock_chain_telescopes_without_owner_reset_currency():
    rows = (
        StockStep(energy_in=5.0, energy_out=6.0, dissipation=0.5, signed_physical_work=1.5),
        StockStep(energy_in=6.0, energy_out=4.25, dissipation=1.0, signed_physical_work=-0.75),
        StockStep(energy_in=4.25, energy_out=4.5, dissipation=0.25, signed_physical_work=0.5),
    )
    out = stock_telescope(rows)
    assert out["initial"] == pytest.approx(5.0)
    assert out["final"] == pytest.approx(4.5)
    assert out["dissipation"] == pytest.approx(1.75)
    assert out["signed_work"] == pytest.approx(1.25)
    assert out["residual"] == pytest.approx(0.0)
    assert out["owner_reset"] is False


def test_owner_change_cannot_reset_between_time_stock_endpoint():
    first = StockStep(energy_in=2.0, energy_out=2.5, dissipation=0.25, signed_physical_work=0.75)
    second = StockStep(energy_in=3.0, energy_out=3.0, dissipation=0.0, signed_physical_work=0.0)
    with pytest.raises(ValueError, match="owner change attempted to reset endpoint stock"):
        stock_telescope((first, second))


def test_stock_certificate_keeps_same_time_work_distinct_from_persistent_stock():
    cert = theorem_certificate()
    assert "mode-set continuity" in cert["primitive_below_stock"]
    assert "E_A(t1)+D_A+Phi_out" in cert["mode_stock_continuity"]
    assert "cannot reset physical stock ancestry" in cert["no_mint"]
    assert "actual signed physical work" in cert["no_mint"]
    assert "same-time donor provenance" in cert["same_time_guard"]
    assert "conservative circulation" in cert["scope"]
