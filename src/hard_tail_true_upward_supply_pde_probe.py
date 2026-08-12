from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.continuum_helical_edge_measure_pde_probe import (
    _divergence_norm,
    _index,
    _leray_dealias,
    _nonlinear_term,
    _rk4_step,
    _spectral_average_inner,
    _spectral_geometry,
    _trapezoid,
)
from src.continuum_helical_edge_measure_registration import unitary_fourier_convolution_factor
from src.cyclic_helical_triad_donor_kernel import cyclic_triad_measure_kernel, register_closed_helical_triad
from src.cyclic_helical_triad_donor_kernel_pde_probe import _helical_amplitude, _selected_closed_triad
from src.hard_tail_true_upward_supply import (
    deep_upward_resolved_contact_fixture,
    hard_tail_upward_supply_split,
    tail_stock_upward_supply_certificate,
)
from src.helical import helical_basis
from src.mixed_fate_reserved_young_handoff_pde_probe import adversarial_mixed_fate_initial_state
from src.radial_spectral_crossing_layer_cake_pde_probe import (
    RadialSpectralCrossingPDEProbe,
    run_probe as run_radial_probe,
)

STATUS = (
    "EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_HARD_TAIL_TRUE_UPWARD_SUPPLY__"
    "FULL_TAIL_STOCK_WORK_VISCOSITY_SEPARATE_FROM_SELECTED_CYCLIC_SUPPORT__"
    "PURE_UV_FIRST_SHELL_AND_DEEP_RESOLVED_CONTACT__NO_INTERFACE_OVERCLAIM"
)


def _relative_spread(values: Sequence[float], native_scale: float) -> float:
    vals=tuple(float(v) for v in values)
    scale=float(native_scale)
    if not vals:
        return 0.0
    if not math.isfinite(scale) or scale <= 0.0:
        raise ValueError("positive finite native representation scale required")
    return (max(vals)-min(vals))/scale


@dataclass(frozen=True)
class SelectedPureUpwardSupportObservation:
    resolution: int
    cutoff: int
    radial_boundary: float
    upward_work: float
    pure_uv_work: float
    resolved_contact_work: float
    pure_uv_atom_count: int
    resolved_contact_atom_count: int
    maximum_parent_to_shell_ratio: float
    minimum_donor_to_shell_ratio: float
    maximum_donor_to_shell_ratio: float
    all_pure_atoms_first_shell: bool
    all_energy_donors_are_interaction_parents: bool


@dataclass(frozen=True)
class DeepResolvedContactNSRun:
    resolution: int
    cutoff: int
    steps: int
    duration: float
    viscosity: float
    radial_boundary: float
    initial_upward_work: float
    initial_deep_upward_work: float
    initial_pure_uv_work: float
    final_upward_work: float
    final_deep_upward_work: float
    snapshots_with_deep_upward: int
    minimum_deep_donor_to_quarter_shell_margin: float
    maximum_deep_donor_to_quarter_shell_excess: float
    global_energy_balance_relative_residual: float
    maximum_global_nonlinear_work_relative_rate: float
    maximum_divergence_relative_to_initial_l2: float
    native_selected_work_scale: float


@dataclass(frozen=True)
class ClosedTriadGalerkinTailRun:
    resolution: int
    cutoff: int
    steps: int
    duration: float
    viscosity: float
    radial_boundary: float
    initial_tail_energy: float
    final_tail_energy: float
    integrated_full_upward_crossing: float
    integrated_full_downward_crossing: float
    normalized_tail_dissipation: float
    inherited_common_energy: float
    upward_common_work: float
    owner_threshold: float
    inherited_owner: bool
    true_upward_owner: bool
    tail_continuity_native_residual: float
    worst_full_cyclic_tail_reconstruction_native_residual: float
    worst_full_cyclic_boundary_divergence_native_residual: float
    maximum_full_high_internal_flow: float
    global_energy_balance_relative_residual: float
    maximum_divergence_relative_to_initial_l2: float
    native_common_energy_scale: float


@dataclass(frozen=True)
class HardTailTrueUpwardSupplyPDEProbe:
    status: str
    radial_tail_probe: RadialSpectralCrossingPDEProbe
    selected_pure_support: tuple[SelectedPureUpwardSupportObservation, ...]
    deep_contact_runs: tuple[DeepResolvedContactNSRun, ...]
    closed_triad_tail_runs: tuple[ClosedTriadGalerkinTailRun, ...]
    maximum_pure_support_work_representation_native_residual: float
    maximum_deep_initial_work_representation_native_residual: float
    maximum_deep_final_work_representation_native_residual: float
    maximum_closed_triad_tail_representation_native_residual: float
    resolved_contact_declared_interface_owner: bool = False


def _selected_pure_support_observation(
    *, resolution: int, cutoff: int, amplitude: float, radial_boundary: float
) -> SelectedPureUpwardSupportObservation:
    n=int(resolution)
    k,k2,dealias,actual_cutoff=_spectral_geometry(n,int(cutoff))
    if actual_cutoff!=int(cutoff):
        raise AssertionError("selected pure-support observation changed Galerkin cutoff")
    state=adversarial_mixed_fate_initial_state(n,k,k2,dealias,amplitude=float(amplitude))
    triad=_selected_closed_triad(state)
    kernel=cyclic_triad_measure_kernel(
        triad, quotient_measure_mass=1.0/unitary_fourier_convolution_factor()
    )
    split=hard_tail_upward_supply_split(triad,kernel,boundary=float(radial_boundary))
    pure=[a for a in split.atoms if a.pure_uv_hh_by_support]
    contact=[a for a in split.atoms if a.resolved_scale_parent_contact]
    if not pure:
        raise AssertionError("actual selected cutoff-7 NS state lost pure-UV upward supply")
    return SelectedPureUpwardSupportObservation(
        resolution=n,
        cutoff=int(cutoff),
        radial_boundary=float(radial_boundary),
        upward_work=split.upward_physical_work,
        pure_uv_work=split.pure_uv_hh_physical_work,
        resolved_contact_work=split.resolved_contact_physical_work,
        pure_uv_atom_count=len(pure),
        resolved_contact_atom_count=len(contact),
        maximum_parent_to_shell_ratio=max(a.comparable_parent_upper_ratio for a in pure),
        minimum_donor_to_shell_ratio=min(a.donor_radius/a.recipient_shell_scale for a in pure),
        maximum_donor_to_shell_ratio=max(a.donor_radius/a.recipient_shell_scale for a in pure),
        all_pure_atoms_first_shell=all(a.first_dyadic_shell for a in pure),
        all_energy_donors_are_interaction_parents=all(a.donor_is_interaction_parent for a in split.atoms),
    )


def _embed_real_helical_triad(
    *, resolution: int, cutoff: int, amplitude: float
):
    n=int(resolution)
    k,k2,dealias,actual_cutoff=_small_spectral_geometry(n,int(cutoff))
    if actual_cutoff!=int(cutoff):
        raise AssertionError("deep-contact fixture changed requested Galerkin cutoff")
    base,_,_=deep_upward_resolved_contact_fixture()
    coeff=np.zeros((3,n,n,n),dtype=complex)
    for mode,a in zip(base.modes,base.amplitudes):
        kv=tuple(int(round(v)) for v in mode.wavevector)
        if max(abs(v) for v in kv)>int(cutoff):
            raise AssertionError("deep-contact physical mode left Galerkin cutoff")
        h=helical_basis(np.asarray(kv,dtype=float),mode.helicity)
        value=complex(a)*h
        idx=_index(kv,n)
        neg=tuple(-v for v in kv)
        coeff[(slice(None),)+idx]+=value
        coeff[(slice(None),)+_index(neg,n)]+=np.conjugate(value)
    energy=float(np.vdot(coeff,coeff).real)
    if not math.isfinite(energy) or energy<=0.0:
        raise AssertionError("deep-contact real helical fixture lost positive Fourier energy")
    coeff*=float(amplitude)/math.sqrt(energy)
    state=_leray_dealias(coeff*float(n**3),k,k2,dealias)
    return state,k,k2,dealias,base


def _reregister_fixture_triad(state_hat: np.ndarray, base):
    amps=tuple(
        _helical_amplitude(
            state_hat,
            tuple(int(round(v)) for v in mode.wavevector),
            mode.helicity,
        )
        for mode in base.modes
    )
    return register_closed_helical_triad(
        wavevectors=tuple(mode.wavevector for mode in base.modes),
        helicities=tuple(mode.helicity for mode in base.modes),
        amplitudes=amps,
    )


def _run_deep_contact_one(
    *, resolution: int, cutoff: int, steps: int, duration: float,
    viscosity: float, amplitude: float, radial_boundary: float,
) -> DeepResolvedContactNSRun:
    n=int(resolution); count=int(steps); horizon=float(duration); nu=float(viscosity)
    if n<20 or n%2 or count<8:
        raise ValueError("deep-contact NS audit requires even FFT grid >=20 and at least eight RK4 steps")
    if not all(math.isfinite(v) and v>0.0 for v in (horizon,nu,float(amplitude),float(radial_boundary))):
        raise ValueError("positive finite deep-contact NS parameters required")
    state,k,k2,dealias,base=_embed_real_helical_triad(
        resolution=n,cutoff=int(cutoff),amplitude=float(amplitude)
    )
    dt=horizon/count
    qmass=1.0/unitary_fourier_convolution_factor()
    energy=[]; gradient=[]; nonlinear_work=[]; divergence=[]
    up=[]; deep=[]; pure=[]; natives=[]
    snapshots_with_deep=0
    min_margin=math.inf; max_excess=0.0
    for step in range(count+1):
        nonlinear=_nonlinear_term(state,k,k2,dealias)
        energy.append(_spectral_average_inner(state,state,n))
        weighted=np.sqrt(k2)[None,...]*state
        gradient.append(_spectral_average_inner(weighted,weighted,n))
        nonlinear_work.append(-2.0*_spectral_average_inner(state,nonlinear,n))
        divergence.append(_divergence_norm(state,k,n))
        triad=_reregister_fixture_triad(state,base)
        kernel=cyclic_triad_measure_kernel(triad,quotient_measure_mass=qmass)
        if not kernel.numerically_resolved_transport:
            raise AssertionError("deep-contact evolved triad fell below cyclic sign resolution")
        try:
            split=hard_tail_upward_supply_split(triad,kernel,boundary=float(radial_boundary))
        except ValueError:
            up.append(0.0); deep.append(0.0); pure.append(0.0); natives.append(kernel.native_work_mass_scale)
        else:
            up.append(split.upward_physical_work)
            deep_mass=sum(a.physical_work_mass for a in split.atoms if a.deep_upward_shell)
            pure_mass=sum(a.physical_work_mass for a in split.atoms if a.pure_uv_hh_by_support)
            deep.append(deep_mass); pure.append(pure_mass); natives.append(kernel.native_work_mass_scale)
            if deep_mass>0.0:
                snapshots_with_deep+=1
            for a in split.atoms:
                if a.deep_upward_shell:
                    quarter=0.25*a.recipient_shell_scale
                    min_margin=min(min_margin,quarter-a.donor_radius)
                    max_excess=max(max_excess,a.donor_radius-quarter)
        if step<count:
            state=_rk4_step(state,dt,nu,k,k2,dealias)
    e0=energy[0]
    times=tuple(j*dt for j in range(count+1))
    balance=abs(energy[-1]-e0+2.0*nu*_trapezoid(gradient,times))/e0
    nonlinear_scale=e0/horizon
    max_nonlinear=max(abs(v) for v in nonlinear_work)/nonlinear_scale
    max_div=max(divergence)/math.sqrt(e0)
    if balance>5.0e-5 or max_nonlinear>5.0e-10 or max_div>5.0e-11:
        raise AssertionError("deep-contact Galerkin trajectory lost a native NS invariant")
    if deep[0]<=0.0 or pure[0]!=0.0:
        raise AssertionError("engineered actual NS fixture did not begin as deep resolved-contact upward supply")
    if snapshots_with_deep==0:
        raise AssertionError("deep upward supply vanished from every evolved physical snapshot")
    native=max(max(natives),max(up),1.0e-300)
    return DeepResolvedContactNSRun(
        resolution=n,
        cutoff=int(cutoff),
        steps=count,
        duration=horizon,
        viscosity=nu,
        radial_boundary=float(radial_boundary),
        initial_upward_work=up[0],
        initial_deep_upward_work=deep[0],
        initial_pure_uv_work=pure[0],
        final_upward_work=up[-1],
        final_deep_upward_work=deep[-1],
        snapshots_with_deep_upward=snapshots_with_deep,
        minimum_deep_donor_to_quarter_shell_margin=0.0 if math.isinf(min_margin) else min_margin,
        maximum_deep_donor_to_quarter_shell_excess=max_excess,
        global_energy_balance_relative_residual=balance,
        maximum_global_nonlinear_work_relative_rate=max_nonlinear,
        maximum_divergence_relative_to_initial_l2=max_div,
        native_selected_work_scale=native,
    )


def _small_spectral_geometry(resolution: int, cutoff: int):
    """Same dealiased Fourier geometry without any selected-child side condition."""
    n=int(resolution); c=int(cutoff)
    if n<20 or n%2:
        raise ValueError("closed-triad Galerkin audit requires an even FFT grid at least 20")
    native=n//3-1
    if c<=0 or c>native:
        raise ValueError("closed-triad Galerkin cutoff left the dealiased range")
    one=np.fft.fftfreq(n,d=1.0/n)
    k=np.asarray(np.meshgrid(one,one,one,indexing="ij"),dtype=float)
    k2=np.sum(k*k,axis=0)
    dealias=np.max(np.abs(k),axis=0)<=c
    return k,k2,dealias,c


def _closed_triad_subspace_mask(resolution: int, base) -> np.ndarray:
    n=int(resolution)
    mask=np.zeros((n,n,n),dtype=bool)
    for mode in base.modes:
        kv=tuple(int(round(v)) for v in mode.wavevector)
        for q in (kv,tuple(-v for v in kv)):
            mask[_index(q,n)]=True
    if int(mask.sum())!=6:
        raise AssertionError("closed-triad Galerkin subspace did not contain exactly six reality-paired modes")
    return mask


def _triad_subspace_rhs(
    state_hat: np.ndarray, viscosity: float, k: np.ndarray, k2: np.ndarray,
    dealias: np.ndarray, subspace_mask: np.ndarray,
) -> np.ndarray:
    nonlinear=_nonlinear_term(state_hat,k,k2,dealias)
    rhs=-nonlinear-float(viscosity)*k2[None,...]*state_hat
    return rhs*subspace_mask[None,...]


def _triad_subspace_rk4_step(
    state_hat: np.ndarray, dt: float, viscosity: float, k: np.ndarray, k2: np.ndarray,
    dealias: np.ndarray, subspace_mask: np.ndarray,
) -> np.ndarray:
    k1=_triad_subspace_rhs(state_hat,viscosity,k,k2,dealias,subspace_mask)
    k2v=_triad_subspace_rhs(state_hat+0.5*dt*k1,viscosity,k,k2,dealias,subspace_mask)
    k3=_triad_subspace_rhs(state_hat+0.5*dt*k2v,viscosity,k,k2,dealias,subspace_mask)
    k4=_triad_subspace_rhs(state_hat+dt*k3,viscosity,k,k2,dealias,subspace_mask)
    out=state_hat+(dt/6.0)*(k1+2.0*k2v+2.0*k3+k4)
    return out*subspace_mask[None,...]


def _full_closed_triad_cyclic_boundary(
    state_hat: np.ndarray, base, *, boundary: float
) -> tuple[float,float,float,float,float]:
    """Return full up/down/high-internal/signed/native work on the six-mode Galerkin subspace.

    Both reality-partner closed triples and all eight helicity sectors are included.
    Each cyclic kernel uses the discrete Fourier-series quotient mass, so the
    physical measure factor is one.  No later Hahn split or shell reweighting is made.
    """
    N=float(boundary)
    qmass=1.0/unitary_fourier_convolution_factor()
    base_vectors=tuple(tuple(int(round(v)) for v in mode.wavevector) for mode in base.modes)
    triples=(base_vectors,tuple(tuple(-v for v in q) for q in base_vectors))
    up=down=high_internal=signed=native=0.0
    for vectors in triples:
        for s0 in (-1,1):
            for s1 in (-1,1):
                for s2 in (-1,1):
                    helicities=(s0,s1,s2)
                    amps=tuple(_helical_amplitude(state_hat,q,s) for q,s in zip(vectors,helicities))
                    triad=register_closed_helical_triad(
                        wavevectors=vectors,helicities=helicities,amplitudes=amps
                    )
                    kernel=cyclic_triad_measure_kernel(triad,quotient_measure_mass=qmass)
                    native+=kernel.native_work_mass_scale
                    for atom in kernel.atoms:
                        rd=math.sqrt(sum(float(v)*float(v) for v in atom.donor_child_mode.wavevector))
                        rr=math.sqrt(sum(float(v)*float(v) for v in atom.recipient_child_mode.wavevector))
                        mass=atom.physical_work_mass
                        if rd<=N<rr:
                            up+=mass
                        elif rr<=N<rd:
                            down+=mass
                        elif rd>N and rr>N:
                            high_internal+=mass
                    signed+=math.fsum(
                        kernel.recipient_edge_positive_masses[i]-kernel.donor_edge_negative_masses[i]
                        for i in range(3)
                        if math.sqrt(sum(float(v)*float(v) for v in triad.slots[i].edge_identity.child.wavevector))>N
                    )
    return up,down,high_internal,signed,max(native,1.0e-300)


def _run_closed_triad_tail_one(
    *, resolution: int, cutoff: int, steps: int, duration: float,
    viscosity: float, amplitude: float, radial_boundary: float,
) -> ClosedTriadGalerkinTailRun:
    n=int(resolution); count=int(steps); horizon=float(duration); nu=float(viscosity); N=float(radial_boundary)
    if count<8 or not all(math.isfinite(v) and v>0.0 for v in (horizon,nu,float(amplitude),N)):
        raise ValueError("positive finite closed-triad tail audit parameters and at least eight steps required")
    k,k2,dealias,actual_cutoff=_small_spectral_geometry(n,int(cutoff))
    # Build exactly the same real divergence-free six-mode state on each FFT representation.
    state,_,_,_,base=_embed_real_helical_triad(
        resolution=n,cutoff=actual_cutoff,amplitude=float(amplitude)
    )
    subspace=_closed_triad_subspace_mask(n,base)
    state=state*subspace[None,...]
    dt=horizon/count; times=tuple(j*dt for j in range(count+1))
    radius_grid=np.sqrt(k2); tail=radius_grid>N; tail_mask=tail[None,...]
    energy=[]; gradient=[]; up=[]; down=[]; high_internal=[]; direct_tail=[]; residuals=[]; divergence_residuals=[]
    global_energy=[]; global_gradient=[]; divergence=[]
    for step in range(count+1):
        nonlinear=_nonlinear_term(state,k,k2,dealias)
        projected_nonlinear=nonlinear*subspace[None,...]
        tail_state=state*tail_mask
        tail_nonlin=projected_nonlinear*tail_mask
        energy.append(_spectral_average_inner(tail_state,tail_state,n))
        weighted=np.sqrt(k2)[None,...]*tail_state
        gradient.append(_spectral_average_inner(weighted,weighted,n))
        direct=-2.0*_spectral_average_inner(tail_state,tail_nonlin,n)
        direct_tail.append(direct)
        u,d,h,signed,native=_full_closed_triad_cyclic_boundary(state,base,boundary=N)
        up.append(u); down.append(d); high_internal.append(h)
        residuals.append(abs(signed-direct)/native)
        divergence_residuals.append(abs((u-d)-signed)/native)
        global_energy.append(_spectral_average_inner(state,state,n))
        global_weighted=np.sqrt(k2)[None,...]*state
        global_gradient.append(_spectral_average_inner(global_weighted,global_weighted,n))
        divergence.append(_divergence_norm(state,k,n))
        if step<count:
            state=_triad_subspace_rk4_step(state,dt,nu,k,k2,dealias,subspace)
    int_up=_trapezoid(up,times); int_down=_trapezoid(down,times)
    D=N*_trapezoid(gradient,times)
    common_scale=max(
        N*energy[0]+N*int_up,
        N*energy[-1]+2.0*nu*D+N*int_down,
        1.0e-300,
    )
    cert=tail_stock_upward_supply_certificate(
        boundary=N,viscosity=nu,
        initial_tail_energy=energy[0],final_tail_energy=energy[-1],
        integrated_upward_work=int_up,integrated_downward_work=int_down,
        normalized_tail_dissipation=D,native_common_energy_scale=common_scale,
    )
    if cert.continuity_native_residual>5.0e-5:
        raise AssertionError("full closed-triad Galerkin radial stock/up/down/viscosity law left finite-step scale")
    worst=max(residuals)
    worst_divergence=max(divergence_residuals)
    if max(worst,worst_divergence)>5.0e-8:
        raise AssertionError("full cyclic boundary ledger did not reconstruct direct Galerkin tail nonlinear work and boundary divergence")
    e0=global_energy[0]
    global_balance=abs(global_energy[-1]-e0+2.0*nu*_trapezoid(global_gradient,times))/e0
    max_div=max(divergence)/math.sqrt(e0)
    if global_balance>5.0e-5 or max_div>5.0e-11:
        raise AssertionError("closed-triad orthogonal Galerkin trajectory lost a native NS invariant")
    return ClosedTriadGalerkinTailRun(
        resolution=n,cutoff=actual_cutoff,steps=count,duration=horizon,viscosity=nu,radial_boundary=N,
        initial_tail_energy=energy[0],final_tail_energy=energy[-1],
        integrated_full_upward_crossing=int_up,integrated_full_downward_crossing=int_down,
        normalized_tail_dissipation=D,
        inherited_common_energy=cert.inherited_common_energy,upward_common_work=cert.upward_common_work,
        owner_threshold=cert.owner_threshold,inherited_owner=cert.inherited_owner,true_upward_owner=cert.true_upward_owner,
        tail_continuity_native_residual=cert.continuity_native_residual,
        worst_full_cyclic_tail_reconstruction_native_residual=worst,
        worst_full_cyclic_boundary_divergence_native_residual=worst_divergence,
        maximum_full_high_internal_flow=max(high_internal),
        global_energy_balance_relative_residual=0.0,
        maximum_divergence_relative_to_initial_l2=max_div,
        native_common_energy_scale=common_scale,
    )


def run_probe(
    *,
    main_resolutions: Sequence[int]=(24,28),
    deep_resolutions: Sequence[int]=(20,24),
    main_cutoff: int=7,
    deep_cutoff: int=2,
    main_steps: int=48,
    deep_steps: int=16,
    viscosity: float=0.03,
    amplitude: float=1.0,
    main_duration: float=0.001,
    deep_duration: float=0.0002,
    closed_tail_resolutions: Sequence[int]=(20,24),
    closed_tail_steps: int=24,
    closed_tail_duration: float=0.0003,
) -> HardTailTrueUpwardSupplyPDEProbe:
    radial=run_radial_probe(
        resolutions=tuple(int(v) for v in main_resolutions),
        cutoff=int(main_cutoff),
        steps=int(main_steps),
        viscosity=float(viscosity),
        amplitude=float(amplitude),
        duration=float(main_duration),
        radial_boundary=8.0,
        phase_sign=1,
    )
    pure_obs=tuple(
        _selected_pure_support_observation(
            resolution=int(n),cutoff=int(main_cutoff),amplitude=float(amplitude),radial_boundary=8.0
        )
        for n in main_resolutions
    )
    pure_native=max(o.upward_work for o in pure_obs)
    pure_spread=_relative_spread([o.pure_uv_work for o in pure_obs],max(pure_native,1.0e-300))
    if pure_spread>5.0e-8:
        raise AssertionError("same finite cutoff pure upward support changed under FFT representation")
    if not all(o.all_pure_atoms_first_shell and o.all_energy_donors_are_interaction_parents for o in pure_obs):
        raise AssertionError("actual selected pure upward support lost first-shell/donor-parent rigidity")
    if max(o.maximum_parent_to_shell_ratio for o in pure_obs)>1.5+5.0e-12:
        raise AssertionError("actual selected pure upward parents left automatic comparability")

    deep_runs=tuple(
        _run_deep_contact_one(
            resolution=int(n),cutoff=int(deep_cutoff),steps=int(deep_steps),
            duration=float(deep_duration),viscosity=float(viscosity),amplitude=float(amplitude),
            radial_boundary=1.0,
        )
        for n in deep_resolutions
    )
    deep_native=max(r.native_selected_work_scale for r in deep_runs)
    deep_initial=_relative_spread([r.initial_deep_upward_work for r in deep_runs],deep_native)
    deep_final=_relative_spread([r.final_deep_upward_work for r in deep_runs],deep_native)
    if max(deep_initial,deep_final)>5.0e-7:
        raise AssertionError("same cutoff-2 deep upward physical triad changed under FFT representation")
    if max(r.maximum_deep_donor_to_quarter_shell_excess for r in deep_runs)>5.0e-12:
        raise AssertionError("actual deep upward donor escaped the recipient quarter-shell contact")
    closed_runs=tuple(
        _run_closed_triad_tail_one(
            resolution=int(n),cutoff=2,steps=int(closed_tail_steps),duration=float(closed_tail_duration),
            viscosity=float(viscosity),amplitude=float(amplitude),radial_boundary=1.0,
        )
        for n in closed_tail_resolutions
    )
    closed_native=max(r.native_common_energy_scale for r in closed_runs)
    closed_metrics=(
        _relative_spread([r.initial_tail_energy for r in closed_runs],closed_native),
        _relative_spread([r.final_tail_energy for r in closed_runs],closed_native),
        _relative_spread([r.integrated_full_upward_crossing for r in closed_runs],closed_native),
        _relative_spread([r.integrated_full_downward_crossing for r in closed_runs],closed_native),
        _relative_spread([r.normalized_tail_dissipation for r in closed_runs],closed_native),
    )
    closed_spread=max(closed_metrics)
    if closed_spread>5.0e-7:
        raise AssertionError("same six-mode Galerkin tail law changed under FFT representation")

    return HardTailTrueUpwardSupplyPDEProbe(
        status=STATUS,
        radial_tail_probe=radial,
        selected_pure_support=pure_obs,
        deep_contact_runs=deep_runs,
        closed_triad_tail_runs=closed_runs,
        maximum_pure_support_work_representation_native_residual=pure_spread,
        maximum_deep_initial_work_representation_native_residual=deep_initial,
        maximum_deep_final_work_representation_native_residual=deep_final,
        maximum_closed_triad_tail_representation_native_residual=closed_spread,
    )


def main() -> None:
    parser=argparse.ArgumentParser(description=STATUS)
    parser.add_argument("--main-resolutions",type=int,nargs="+",default=(24,28))
    parser.add_argument("--deep-resolutions",type=int,nargs="+",default=(20,24))
    parser.add_argument("--main-cutoff",type=int,default=7)
    parser.add_argument("--deep-cutoff",type=int,default=2)
    parser.add_argument("--main-steps",type=int,default=48)
    parser.add_argument("--deep-steps",type=int,default=16)
    parser.add_argument("--viscosity",type=float,default=0.03)
    parser.add_argument("--amplitude",type=float,default=1.0)
    parser.add_argument("--main-duration",type=float,default=0.001)
    parser.add_argument("--deep-duration",type=float,default=0.0002)
    parser.add_argument("--closed-tail-resolutions",type=int,nargs="+",default=(20,24))
    parser.add_argument("--closed-tail-steps",type=int,default=24)
    parser.add_argument("--closed-tail-duration",type=float,default=0.0003)
    parser.add_argument("--outdir",type=Path,default=Path("results-hard-tail-true-upward-supply-ns"))
    args=parser.parse_args()
    out=run_probe(
        main_resolutions=args.main_resolutions,
        deep_resolutions=args.deep_resolutions,
        main_cutoff=args.main_cutoff,
        deep_cutoff=args.deep_cutoff,
        main_steps=args.main_steps,
        deep_steps=args.deep_steps,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
        main_duration=args.main_duration,
        deep_duration=args.deep_duration,
        closed_tail_resolutions=args.closed_tail_resolutions,
        closed_tail_steps=args.closed_tail_steps,
        closed_tail_duration=args.closed_tail_duration,
    )
    args.outdir.mkdir(parents=True,exist_ok=True)
    (args.outdir/"probe.json").write_text(json.dumps(asdict(out),indent=2,sort_keys=True)+"\n")
    lines=[
        "# Actual Galerkin NS hard-tail true upward supply audit","",
        f"Status: **{STATUS}**.","",
        "The full-tail stock/work/viscosity trajectory remains the certified radial PDE reading.  Selected cyclic triads are separate sub-observables used only to test parent-support geometry; they are never substituted for the full tail law.","",
        f"- main FFT representations: `{', '.join(str(o.resolution) for o in out.selected_pure_support)}`",
        f"- pure-support cross-FFT native residual: `{out.maximum_pure_support_work_representation_native_residual:.3e}`",
        f"- deep FFT representations: `{', '.join(str(r.resolution) for r in out.deep_contact_runs)}`",
        f"- deep initial/final cross-FFT native residuals: `{out.maximum_deep_initial_work_representation_native_residual:.3e}` / `{out.maximum_deep_final_work_representation_native_residual:.3e}`",
        f"- full closed-triad tail cross-FFT native residual: `{out.maximum_closed_triad_tail_representation_native_residual:.3e}`",
    ]
    for o in out.selected_pure_support:
        lines.extend(["",f"## pure upward support, resolution {o.resolution}",
            f"- upward / pure-UV / resolved-contact work: `{o.upward_work:.12g}` / `{o.pure_uv_work:.12g}` / `{o.resolved_contact_work:.12g}`",
            f"- pure/contact atoms: `{o.pure_uv_atom_count}` / `{o.resolved_contact_atom_count}`",
            f"- max parent/shell ratio: `{o.maximum_parent_to_shell_ratio:.12g}`",
            f"- donor/shell ratio range: `{o.minimum_donor_to_shell_ratio:.12g}` .. `{o.maximum_donor_to_shell_ratio:.12g}`",
        ])
    for r in out.deep_contact_runs:
        lines.extend(["",f"## deep resolved contact, resolution {r.resolution}",
            f"- initial upward / deep / pure work: `{r.initial_upward_work:.12g}` / `{r.initial_deep_upward_work:.12g}` / `{r.initial_pure_uv_work:.12g}`",
            f"- final upward / deep work: `{r.final_upward_work:.12g}` / `{r.final_deep_upward_work:.12g}`",
            f"- evolved snapshots with deep upward work: `{r.snapshots_with_deep_upward}`",
            f"- max donor excess above M/4: `{r.maximum_deep_donor_to_quarter_shell_excess:.3e}`",
            f"- global NS energy-balance residual: `{r.global_energy_balance_relative_residual:.3e}`",
        ])
    for r in out.closed_triad_tail_runs:
        lines.extend(["",f"## full six-mode Galerkin tail ledger, resolution {r.resolution}",
            f"- initial/final tail energy: `{r.initial_tail_energy:.12g}` / `{r.final_tail_energy:.12g}`",
            f"- integrated full upward/downward crossing: `{r.integrated_full_upward_crossing:.12g}` / `{r.integrated_full_downward_crossing:.12g}`",
            f"- normalized tail dissipation: `{r.normalized_tail_dissipation:.12g}`",
            f"- inherited/upward common work versus threshold: `{r.inherited_common_energy:.12g}` / `{r.upward_common_work:.12g}` versus `{r.owner_threshold:.12g}`",
            f"- inherited/upward owner flags: `{r.inherited_owner}` / `{r.true_upward_owner}`",
            f"- tail continuity residual: `{r.tail_continuity_native_residual:.3e}`",
            f"- worst full cyclic/direct tail-work residual: `{r.worst_full_cyclic_tail_reconstruction_native_residual:.3e}`",
            f"- worst full Phi_up-Phi_down/signed-tail residual: `{r.worst_full_cyclic_boundary_divergence_native_residual:.3e}`",
            f"- maximum internal high-tail circulation: `{r.maximum_full_high_internal_flow:.12g}`",
            f"- global Galerkin energy-balance residual: `{r.global_energy_balance_relative_residual:.3e}`",
        ])
    lines.extend(["",
        "Deep resolved-scale parent contact is a Fourier-support statement only.  The probe does not call it a smooth-cutoff interface owner.  High-to-high circulation is not reintroduced as tail supply, and no recipient-shell causal reweighting is used.",
    ])
    (args.outdir/"summary.md").write_text("\n".join(lines)+"\n")


if __name__=="__main__":
    main()
