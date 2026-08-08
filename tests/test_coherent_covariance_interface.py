import numpy as np
from src.coherent_covariance_interface import (
    clean_window_distance_upper, covariance_energy_tv_upper, discrete_energy_density_tv,
    discrete_window_change_residual, gaussian_covariance_overlap, gaussian_window_distance,
    generalized_log_covariance_distance,
)


def test_scalar_covariance_overlap():
    S=np.eye(3); T=4*np.eye(3)
    expected=(2*np.sqrt(2)/(5**0.5))**3
    assert abs(gaussian_covariance_overlap(S,T)-expected)<1e-13


def test_clean_log_distance_bound():
    S=np.eye(3); T=np.diag([1.1,.9,1.02])
    assert gaussian_window_distance(S,T)<=clean_window_distance_upper(S,T)+1e-14


def test_window_slot_moyal_discrete():
    rng=np.random.default_rng(7); n=23
    f=rng.normal(size=n)+1j*rng.normal(size=n); g=rng.normal(size=n)+1j*rng.normal(size=n); h=rng.normal(size=n)+1j*rng.normal(size=n); g/=np.linalg.norm(g); h/=np.linalg.norm(h)
    assert abs(discrete_window_change_residual(f,g,h))<1e-10
    tv,b=discrete_energy_density_tv(f,g,h); assert tv<=b+1e-10


def test_covariance_tv_formula_positive():
    S=np.eye(3); T=np.diag([1.2,.8,1.1])
    assert generalized_log_covariance_distance(S,T)>0
    assert covariance_energy_tv_upper(2.0,S,T)>0
