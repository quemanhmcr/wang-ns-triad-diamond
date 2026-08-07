import numpy as np
from src.affine_polarization_curvature import physical_strain_gradient, rms_transfer_relevant_curvature


def test_rms_tomography_lower_bound():
    rng=np.random.default_rng(2)
    H=rng.normal(size=(3,3,3));H=.5*(H+H.swapaxes(1,2))
    for k in range(3): H[2,2,k]=H[2,k,2]=-(H[0,0,k]+H[1,1,k])
    L=np.diag([2.,1.,.7])
    Q,N=rms_transfer_relevant_curvature(H,L)
    assert Q>=.5*N-1e-12
    C=physical_strain_gradient(H,L)
    assert abs(N-np.sum(C*C))<1e-10
