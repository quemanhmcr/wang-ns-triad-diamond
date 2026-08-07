import math
import numpy as np
from src.atomic_component_entropy import *

def test_chain_exact():
    w=[.1,.2,.3,.4]; lab=['A','A','B','B']; d=collision_chain(w,lab)
    assert abs(d['q_ancestry']-d['q_atom']-d['hidden_pair_mass'])<1e-14
    assert abs(d['hidden_pair_mass']-d['hidden_pair_decomposition'])<1e-14
    assert abs(d['hidden_entropy'] + math.log(d['conditional_collision_mean']))<1e-14

def test_probabilistic_identity_direct():
    w=normalize([1,2,3,4]); lab=['A','A','B','B']; brute=0
    for i in range(len(w)):
      for j in range(len(w)):
        if i!=j and lab[i]==lab[j]: brute += w[i]*w[j]
    assert abs(brute-collision_chain(w,lab)['hidden_pair_mass'])<1e-14

def test_entropy_pair_dichotomy_random():
    rng=np.random.default_rng(3)
    for _ in range(1000):
      n=int(rng.integers(3,30)); w=rng.dirichlet(np.ones(n)); lab=rng.integers(0,5,size=n)
      assert verify_entropy_pair_dichotomy(w,lab)[0]

def test_core_markov_bound():
    rng=np.random.default_rng(4)
    for _ in range(500):
      w=rng.dirichlet(np.ones(20)*.4)
      mass,_,_=core_mass_bound(w,4)
      assert mass >= .75-1e-12

def test_attachment_cycle_gain_example():
    e1,e2=example_graphs()
    assert e1['rank_gain'] >= e1['attachment_lower_bound']
    assert e2['rank_gain'] >= e2['attachment_lower_bound']

def test_single_old_group_exact_gain_on_tree_chain():
    triads=[('a','x','m'),('b','m','n'),('c','n','d')]
    lab={p:'FRESH' for t in triads for p in t}; lab.update(a='A',b='A',c='A')
    d=ancestry_cycle_gain(triads,lab)
    assert d['rank_gain']==2
    assert d['attachment_lower_bound']==2


def test_pair_biased_multiplicity_certificate():
    rng=np.random.default_rng(9)
    for _ in range(500):
        n=int(rng.integers(5,40)); w=rng.dirichlet(np.ones(n)*.5); lab=rng.integers(0,6,size=n).tolist()
        c=pair_biased_multiplicity_certificate(w,lab,lam=2.0)
        assert c['pair_biased_good_mass'] >= .5-1e-12
        if c['good_labels']:
            assert c['minimum_actual_multiplicity'] + 1e-12 >= c['multiplicity_lower_bound']
