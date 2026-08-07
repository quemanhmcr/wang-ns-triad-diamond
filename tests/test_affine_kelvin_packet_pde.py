import numpy as np
from src.affine_kelvin_packet_pde import transverse_frame, kelvin_amplitude_rate, affine_fourier_rhs_from_pressure, objective_frame_rate, coordinate_rate_from_3d, objective_coordinate_generator, unit


def test_pressure_elimination_gives_kelvin_generator():
    A=np.array([[.2,.3,-.1],[-.4,.1,.2],[.1,-.2,-.3]])
    A-=np.trace(A)/3*np.eye(3)
    k=np.array([1.,2.,-.5]); E=transverse_frame(unit(k)); c=np.array([1.+.2j,-.3+.7j]); a=E@c
    assert np.linalg.norm(kelvin_amplitude_rate(A,k,a,.1)-affine_fourier_rhs_from_pressure(A,k,a,.1))<1e-12


def test_objective_coordinate_equation_matches_3d_kelvin():
    A=np.array([[.1,.2,.3],[-.2,-.4,.1],[.05,.2,.3]])
    A-=np.trace(A)/3*np.eye(3)
    k=np.array([.7,-.4,1.2]); E=transverse_frame(unit(k)); c=np.array([.4+.1j,-.2+.3j]); nu=.05
    lhs=coordinate_rate_from_3d(A,k,E,c,nu)
    rhs=objective_coordinate_generator(A,E,nu*np.dot(k,k))@c
    assert np.linalg.norm(lhs-rhs)<1e-12


def test_objective_frame_preserves_constraints_to_first_order():
    A=np.array([[.2,-.1,.3],[.4,-.3,.2],[.1,.2,.1]])
    A-=np.trace(A)/3*np.eye(3); n=unit(np.array([1.,.4,-.2])); E=transverse_frame(n)
    dE=objective_frame_rate(A,n,E)
    assert np.linalg.norm(dE.T@E+E.T@dE)<1e-12
