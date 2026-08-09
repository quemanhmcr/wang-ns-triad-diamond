import numpy as np

from src.event_anchored_role_registration import (
    envelope_low_low_gap,
    envelope_registration_residual,
    hard_role_projectors_from_masks_and_fibers,
    pointwise_projector_lp_contraction,
    role_partition_residual,
)


def test_envelope_buffer_survives_low_strain():
    assert envelope_low_low_gap() > 0.0


def test_hard_frequency_helicity_roles_are_orthogonal():
    masks=np.array([[1,0],[0,1]],float)
    H=np.zeros((2,2,2,2),complex)
    for k in range(2):
        H[0,k]=np.diag([1,0]); H[1,k]=np.diag([0,1])
    Ps=hard_role_projectors_from_masks_and_fibers(masks,H)
    out=role_partition_residual(Ps,np.eye(4))
    assert max(out.values())<1e-12


def test_hard_role_coefficient_equals_smooth_envelope_coefficient():
    P=np.diag([1,0,0]).astype(complex)
    Q=np.diag([1,0.4,0.2]).astype(complex)
    u=np.array([1+2j,3-1j,-2j]); phi=np.array([2-1j,1,4j])
    assert abs(envelope_registration_residual(P,Q,u,phi))<1e-12


def test_pointwise_helicity_projection_contracts_frequency_lp_norms():
    vals=np.array([[1+1j,2],[3j,-1]],complex)
    H=np.zeros((2,2,2),complex)
    H[:,0,0]=1.0
    for p in (2.0,3.0,1.5):
        a,b=pointwise_projector_lp_contraction(vals,H,p)
        assert a<=b+1e-12
