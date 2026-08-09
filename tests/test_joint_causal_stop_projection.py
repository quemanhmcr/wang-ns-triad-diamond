from src.joint_causal_stop_projection import (
    InternalHit,
    joint_stop_master_projection,
    weight_invariance_countercheck,
)
from src.physical_branch_compiler import (
    CauseHit,
    MasterDisposition,
    PhysicalCause,
    UniformResourceCertificate,
)


def test_source_strain_relink_tie_is_one_recursive_joint_stop():
    hits = (
        CauseHit(0.4, PhysicalCause.RESOLVED_SOURCE, 1e-12),
        CauseHit(0.4, PhysicalCause.HIGH_STRAIN_DISSIPATION, 1e12),
        CauseHit(0.4, PhysicalCause.MATERIAL_RELINK, 3.0),
    )
    out = joint_stop_master_projection(physical_hits=hits)
    assert out.master_disposition == MasterDisposition.RECURSE_CRITICAL.value
    assert not out.fine_rn_split_required
    assert len(out.joint_physical_causes) == 3


def test_transfer_certificate_terminates_simultaneous_recursive_stop():
    hits = (
        CauseHit(0.4, PhysicalCause.RESOLVED_SOURCE, 1e20),
        CauseHit(0.4, PhysicalCause.CAUSAL_REUSE, 1e-20),
    )
    out = joint_stop_master_projection(physical_hits=hits)
    assert out.master_disposition == MasterDisposition.TRANSFER_COST.value
    assert out.terminal_certificate_used == "any_simultaneous_fixed_transfer_cost"


def test_exact_tie_projection_is_independent_of_dummy_weights():
    causes = (PhysicalCause.RESOLVED_SOURCE, PhysicalCause.HIGH_STRAIN_DISSIPATION, PhysicalCause.CAUSAL_REUSE)
    assert weight_invariance_countercheck(causes, (1e-30, 2.0, 1e30), (1e30, 5.0, 1e-30))


def test_source_and_hh_regeneration_need_no_split():
    out = joint_stop_master_projection(
        physical_hits=(CauseHit(0.7, PhysicalCause.RESOLVED_SOURCE, 1.0),),
        internal_hits=(InternalHit(0.7),),
    )
    assert out.master_disposition == MasterDisposition.RECURSE_CRITICAL.value
    assert out.joint_internal_causes == ("earlier_high_high_regeneration",)


def test_initial_boundary_absorbs_zero_time_joint_stop():
    out = joint_stop_master_projection(
        physical_hits=(
            CauseHit(0.0, PhysicalCause.INITIAL_BOUNDARY, 1e-50),
            CauseHit(0.0, PhysicalCause.CAUSAL_REUSE, 1e50),
        ),
        internal_hits=(InternalHit(0.0),),
    )
    assert out.master_disposition == MasterDisposition.BOUNDARY.value


def test_valid_uniform_resource_can_terminate_nontransfer_joint_stop():
    cert = UniformResourceCertificate(0.2, 3.0, True, True)
    out = joint_stop_master_projection(
        physical_hits=(
            CauseHit(0.3, PhysicalCause.UNIFORM_GLOBAL_RESOURCE, 1e-20),
            CauseHit(0.3, PhysicalCause.RESOLVED_SOURCE, 1e20),
        ),
        uniform_certificates={PhysicalCause.UNIFORM_GLOBAL_RESOURCE: cert},
    )
    assert out.master_disposition == MasterDisposition.ADDITIVE_RESET.value
