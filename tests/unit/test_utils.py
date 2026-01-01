import pytest
import numpy as np
from src.cstdp.utils.simulate_grn import generate_synthetic_grn, simulate_expression
from src.cstdp.utils.evaluate import calculate_metrics

def test_generate_grn():
    n = 10
    adj, delays = generate_synthetic_grn(n, connection_prob=0.5, seed=42)
    assert adj.shape == (n, n)
    assert delays.shape == (n, n)
    assert np.any(adj > 0) # Should have edges

def test_simulate_expression():
    n = 5
    adj, delays = generate_synthetic_grn(n)
    expr = simulate_expression(n, 100, adj, delays)
    assert expr.shape == (n, 100)
    assert np.all(expr >= 0) # Non-negative constraints

def test_calculate_metrics():
    true_adj = np.array([[0, 1], [0, 0]])
    pred_adj = np.array([[0, 0.8], [0.2, 0]])
    
    # Threshold 0.5: 
    # Pred bin: [[0, 1], [0, 0]]
    # Matches truth exactly
    metrics = calculate_metrics(true_adj, pred_adj, threshold=0.5)
    
    assert metrics['precision'] == 1.0
    assert metrics['recall'] == 1.0
    assert metrics['shd'] == 0
    
    # Threshold 0.1:
    # Pred bin: [[0, 1], [1, 0]]
    # FP = 1 (edge 1->0)
    metrics_loose = calculate_metrics(true_adj, pred_adj, threshold=0.1)
    assert metrics_loose['precision'] == 0.5 # 1 TP / 2 Positive Preds
