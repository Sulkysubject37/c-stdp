import pytest
import numpy as np
from src.cstdp.core import CausalSTDP

def test_stdp_initialization():
    stdp = CausalSTDP(A_pos=0.01, A_neg=0.02, tau_pos=10, tau_neg=20)
    assert stdp.A_pos == 0.01
    assert stdp.A_neg == 0.02
    assert stdp.tau_pos == 10
    assert stdp.tau_neg == 20

def test_stdp_update_potentiation():
    stdp = CausalSTDP(A_pos=1.0, tau_pos=10)
    # Pre at 10, Post at 15 -> dt = 5 (Causal)
    dw = stdp.stdp_update(0, 10, 15)
    expected = 1.0 * np.exp(-5 / 10)
    assert np.isclose(dw, expected)
    assert dw > 0

def test_stdp_update_depression():
    stdp = CausalSTDP(A_neg=1.0, tau_neg=10)
    # Pre at 15, Post at 10 -> dt = -5 (Anti-Causal)
    dw = stdp.stdp_update(0, 15, 10)
    # Note: logic is -A_neg * exp(dt/tau_neg) where dt is negative
    expected = -1.0 * np.exp(-5 / 10)
    assert np.isclose(dw, expected)
    assert dw < 0

def test_compute_spike_times():
    stdp = CausalSTDP()
    # Simple step function: 0, 0, 10, 10
    # Derivative will have a spike at index 2
    data = np.array([[0, 0, 10, 10]])
    time = np.array([0, 1, 2, 3])
    # Threshold < 10
    thresh = np.array([5.0])
    
    spikes = stdp.compute_spike_times(data, time, thresh)
    assert len(spikes) == 1
    # Gradient of [0, 0, 10, 10]
    # at 0: (0-0)/1 = 0
    # at 1: (10-0)/2 = 5 (central diff) -> might not cross >5
    # at 2: (10-0)/2 = 5
    # Let's check numpy gradient behavior precisely or use sharper step
    
    # Sharp step
    data = np.array([[0, 0, 20, 20]])
    thresh = np.array([5.0])
    spikes = stdp.compute_spike_times(data, time, thresh)
    
    assert len(spikes[0]) > 0 # Should detect at least one

def test_run_cstdp_causality():
    # 2 genes: Gene 0 spikes at 10, Gene 1 spikes at 15
    # Should infer 0 -> 1
    stdp = CausalSTDP(A_pos=0.1, A_neg=0.1)
    spikes = [np.array([10.0]), np.array([15.0])]
    
    weights = stdp.run_cstdp(spikes, 2)
    
    # 0->1 (causal) should be positive
    assert weights[0, 1] > 0
    # 1->0 (anti-causal) should be negative or zero (clipped)
    # The run_cstdp clips to [0, w_max]
    assert weights[1, 0] == 0 
