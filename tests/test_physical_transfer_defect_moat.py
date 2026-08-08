import math

from src.physical_transfer_defect_moat import (
    PHYSICAL_DEFECT_MEAN_FACTOR,
    SQUARE_SCHEDULE_SUM_CLEAN,
    defect_moat_certificate,
    physical_good_core_defect_mean_upper,
    square_schedule_total_upper,
)


def test_physical_good_core_mean_constant():
    eps=1e-6
    assert math.isclose(physical_good_core_defect_mean_upper(eps),PHYSICAL_DEFECT_MEAN_FACTOR*eps)
    assert PHYSICAL_DEFECT_MEAN_FACTOR==106/25


def test_clean_square_schedule_constant():
    assert math.pi**2/6-1<SQUARE_SCHEDULE_SUM_CLEAN
    assert math.isclose(square_schedule_total_upper(.2,10,4),SQUARE_SCHEDULE_SUM_CLEAN*(.1+.1))


def test_cross_edges_are_moat_or_tail():
    verts=[('a','b','c'),('c','d','e'),('a','e','f'),('x','y','z')]
    defects=[.01,.02,.8,4.]
    weights=[4.,3.,.2,.1]
    row=defect_moat_certificate(verts,defects,weights,1.0,4)
    assert row['cross_mass']<=row['bound']+1e-13
    assert abs(row['cross_mass']-row['moat_cross_mass']-row['tail_cross_mass'])<1e-13
