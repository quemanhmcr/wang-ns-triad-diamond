import numpy as np

from src.coherent_localization_operators import (
    cell_energy,
    cell_localization_operator,
    partition_operators,
    synthesis_piece_energy,
    trilinear_partition_sum,
    trilinear_tensor_value,
)


def test_positive_partition_resolves_identity_for_orthonormal_frame():
    F=np.eye(4,dtype=complex)
    ops=partition_operators(F,[[0,2],[1,3]])
    assert np.linalg.norm(sum(ops)-np.eye(4))<1e-14
    assert all(np.linalg.eigvalsh(A).min()>-1e-14 for A in ops)


def test_piece_energy_below_moyal_cell_energy():
    F=np.eye(3,dtype=complex); A=cell_localization_operator(F,[0,2]); f=np.array([1+1j,2.,3j])
    assert synthesis_piece_energy(A,f)<=cell_energy(A,f)+1e-14


def test_trilinear_partition_is_exact():
    F=np.eye(2,dtype=complex); ops=partition_operators(F,[[0],[1]])
    f=np.array([1.,2.]); g=np.array([-.3,1.]); h=np.array([2.,.5]); T=np.arange(8,dtype=float).reshape(2,2,2)
    assert abs(trilinear_tensor_value(T,f,g,h)-trilinear_partition_sum(T,f,g,h,ops,ops,ops))<1e-13
