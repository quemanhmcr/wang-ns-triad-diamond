import math

from src.physical_branch_compiler import CauseHit, MasterDisposition, PhysicalCause, PhysicalCurrency
from src.physical_pair_weighted_productivity import physical_log_productivity_constant
from src.recursive_physical_witness_constructor import (
    GeneratedPairEvent,
    RegenerationHit,
    compile_generated_pair_measure,
    compile_generated_pair_master_measure,
)


def test_all_good_events_continue_with_full_physical_productivity():
    events=(GeneratedPairEvent(2.0,0,True,True),GeneratedPairEvent(1.0,1,True,True))
    out=compile_generated_pair_measure(events=events,pair_cells_upper=2)
    assert math.isclose(out.continuation_fraction,1.0)
    assert math.isclose(out.conditioned_productivity,physical_log_productivity_constant(2))
    assert out.regeneration_mass==0.0


def test_exact_source_regeneration_tie_splits_only_supplied_rn_weights():
    event=GeneratedPairEvent(
        1.0,0,True,False,
        physical_hits=(CauseHit(0.3,PhysicalCause.RESOLVED_SOURCE,1.0,"source"),),
        regeneration_hits=(RegenerationHit(0.3,3.0,"regen"),),
    )
    out=compile_generated_pair_measure(events=(event,),pair_cells_upper=1)
    assert math.isclose(out.currency_mass[PhysicalCurrency.RESOLVED_SOURCE_SGS.value],0.25)
    assert math.isclose(out.regeneration_mass,0.75)
    assert out.exact_tie_events==1


def test_initial_boundary_absorbs_zero_time_regeneration_tie():
    event=GeneratedPairEvent(
        1.0,0,True,False,
        physical_hits=(CauseHit(0.0,PhysicalCause.INITIAL_BOUNDARY,0.1,"boundary"),),
        regeneration_hits=(RegenerationHit(0.0,10.0,"regen"),),
    )
    out=compile_generated_pair_measure(events=(event,),pair_cells_upper=1)
    assert out.currency_mass[PhysicalCurrency.INITIAL_BOUNDARY.value]==1.0
    assert out.regeneration_mass==0.0


def test_bad_young_event_cannot_escape_without_transfer_loss_hit():
    bad=GeneratedPairEvent(1.0,0,False,True)
    try:
        compile_generated_pair_measure(events=(bad,),pair_cells_upper=1)
    except ValueError as exc:
        assert "transfer-loss" in str(exc)
    else:
        raise AssertionError("bad Young event escaped without a physical stop")


def test_half_survivor_costs_half_productivity_factor():
    events=(
        GeneratedPairEvent(1.0,0,True,True),
        GeneratedPairEvent(1.0,0,False,True,physical_hits=(CauseHit(0.2,PhysicalCause.TRANSFER_WORK_LOSS,1.0,"deficit"),)),
    )
    out=compile_generated_pair_measure(events=events,pair_cells_upper=1)
    assert math.isclose(out.continuation_fraction,0.5)
    assert out.majority_continues
    assert math.isclose(out.conditioned_productivity,0.5*physical_log_productivity_constant(1))


def test_preferred_joint_master_ignores_exact_tie_weights():
    a=GeneratedPairEvent(1.0,0,True,False,
        physical_hits=(CauseHit(0.3,PhysicalCause.RESOLVED_SOURCE,1e-30,"source"),),
        regeneration_hits=(RegenerationHit(0.3,1e30,"regen"),))
    b=GeneratedPairEvent(1.0,0,True,False,
        physical_hits=(CauseHit(0.3,PhysicalCause.RESOLVED_SOURCE,1e30,"source"),),
        regeneration_hits=(RegenerationHit(0.3,1e-30,"regen"),))
    oa=compile_generated_pair_master_measure(events=(a,),pair_cells_upper=1)
    ob=compile_generated_pair_master_measure(events=(b,),pair_cells_upper=1)
    assert oa.master_mass==ob.master_mass
    assert oa.continuation_mass==ob.continuation_mass==0.0


def test_preferred_joint_master_does_not_validate_dummy_weights():
    # Weight is legacy fine-subledger metadata and is ignored by preferred joint master.
    event=GeneratedPairEvent(1.0,0,True,False,
        physical_hits=(CauseHit(0.2,PhysicalCause.RESOLVED_SOURCE,float("nan"),"source"),))
    out=compile_generated_pair_master_measure(events=(event,),pair_cells_upper=1)
    assert out.master_mass[MasterDisposition.RECURSE_CRITICAL.value]==1.0
