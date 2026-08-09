import numpy as np

from src.heat_edge_material_ownership import (
    IntrinsicBox,
    affine_heat_edge_endpoint_residual,
    heat_edge_intrinsic_endpoints,
    old_pool_membership,
    ownership_class,
    ownership_local_capacity,
    partition_positive_edge_measure,
)


def test_positive_edge_measure_partitions_exactly_without_field_cross_terms():
    w=np.array([1.,2.,3.,4.])
    a=np.array([True,True,False,False])
    b=np.array([True,False,True,False])
    out=partition_positive_edge_measure(w,a,b)
    assert out['old_old']==1
    assert out['old_new_interface']==5
    assert out['new_new']==4
    assert out['partition_residual']==0


def test_ownership_is_unoriented():
    w=[1.,2.,3.]; a=[True,False,True]; b=[False,False,True]
    x=partition_positive_edge_measure(w,a,b)
    y=partition_positive_edge_measure(w,b,a)
    assert x['old_old']==y['old_old']
    assert x['old_new_interface']==y['old_new_interface']
    assert x['new_new']==y['new_new']


def test_both_intrinsic_heat_edge_endpoints_are_affine_invariant():
    M=np.array([[1.1,.2,0],[0,.9,.1],[.1,0,1.02]])
    L=np.array([[1.2,.1,0],[0,.8,.2],[.1,0,1.1]])
    X=np.array([.3,-.7,1.2]); k=np.array([2.,-1.,.5]); r=np.array([.2,.1,-.4])
    assert affine_heat_edge_endpoint_residual(M,L,X,k,r)<1e-12


def test_old_pool_membership_is_defined_on_intrinsic_endpoints_not_cell_names():
    L=np.eye(3); X=np.array([.4,.2,.1]); k=np.array([.3,-.2,.1]); r=np.array([.1,0,0])
    z0,z1=heat_edge_intrinsic_endpoints(L,X,k,r)
    box=IntrinsicBox(tuple([-1]*6),tuple([1]*6))
    assert old_pool_membership(z0,[box])
    assert old_pool_membership(z1,[box])
    assert ownership_class(True,True)=='old_old'
    assert ownership_class(True,False)=='old_new_interface'
    assert ownership_class(False,False)=='new_new'


def test_each_ownership_class_has_its_own_endpoint_energy_capacity():
    A0=np.array([1+1j,2-1j,-1+.5j])
    A1=np.array([.2-.1j,-2+1j,.5+2j])
    ph=np.array([1,1j,-1],complex)
    old0=np.array([True,True,False]); old1=np.array([True,False,False])
    out=ownership_local_capacity(A0,A1,ph,old0,old1,[1,2,.5])
    for name in ('old_old','old_new_interface','new_new'):
        assert out[f'{name}_margin']>=-1e-13
