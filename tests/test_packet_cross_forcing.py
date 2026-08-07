import math
from src.packet_cross_forcing import base_child_parent_remainder_loss, child_representation_loss, split_replacement_identity
from src.transfer_profile_extraction import trilinear_replacement_loss


def test_replacement_loss_splits_exactly():
    for e in [0.,1e-4,.01,.05]:
        assert math.isclose(split_replacement_identity(e),trilinear_replacement_loss(e),rel_tol=0,abs_tol=1e-15)


def test_one_percent_cross_forcing_numbers():
    assert math.isclose(base_child_parent_remainder_loss(.01),.020301,abs_tol=1e-15)
    assert math.isclose(child_representation_loss(.01),.01,abs_tol=1e-15)
