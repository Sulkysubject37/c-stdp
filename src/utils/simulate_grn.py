import numpy as np
import networkx as nx
from typing import Tuple, Dict

def generate_synthetic_grn(n_genes: int = 15, 
                           connection_prob: float = 0.2, 
                           seed: int = 42) -> Tuple[np.ndarray, Dict]:
    """
    Generates a random directed Gene Regulatory Network.
    
    Args:
        n_genes: Number of genes.
        connection_prob: Probability of an edge.
        seed: Random seed.
        
    Returns:
        adj_matrix: Weighted adjacency matrix (Ground Truth).
        delays: Dictionary or matrix of delays for each edge.
    """
    np.random.seed(seed)
    
    # Create a random graph (Erdos-Renyi variant for directed)
    # We want it to be somewhat sparse
    adj_matrix = np.zeros((n_genes, n_genes))
    delays = np.zeros((n_genes, n_genes))
    
    for i in range(n_genes):
        for j in range(n_genes):
            if i == j: continue
            if np.random.random() < connection_prob:
                # Weight in [0.5, 1.0]
                weight = np.random.uniform(0.5, 1.0)
                adj_matrix[i, j] = weight
                
                # Delay in discrete steps (e.g., 5 to 20 steps)
                # Assuming dt=1 for simplicity in simulation
                delay = np.random.randint(5, 21)
                delays[i, j] = delay
                
    return adj_matrix, delays

def simulate_expression(n_genes: int, 
                        n_timepoints: int, 
                        adj_matrix: np.ndarray, 
                        delays: np.ndarray,
                        noise_level: float = 0.1,
                        burst_prob: float = 0.05,
                        decay: float = 0.1,
                        seed: int = 42) -> np.ndarray:
    """
    Simulates time-series gene expression based on the GRN.
    Uses a simple phenomenological model:
    dX_j/dt = -decay * X_j + sum(w_ij * input_from_i) + noise
    
    Where input_from_i is a "burst" at t - delay.
    
    Args:
        n_genes: Number of genes.
        n_timepoints: Duration of simulation.
        adj_matrix: Ground truth weights.
        delays: Delays in time steps.
        noise_level: Std dev of additive noise.
        burst_prob: Probability of spontaneous activation (for source nodes).
        decay: Decay rate of expression.
        seed: Random seed.
        
    Returns:
        expression: Shape (n_genes, n_timepoints)
    """
    np.random.seed(seed)
    expression = np.zeros((n_genes, n_timepoints))
    
    # Pre-generate spontaneous bursts (inputs)
    # This represents external signals or master regulators
    spontaneous_input = np.zeros((n_genes, n_timepoints))
    for t in range(n_timepoints):
        if np.random.random() < burst_prob:
            # Random gene gets a burst
            target = np.random.randint(0, n_genes)
            spontaneous_input[target, t] = 1.0
            
    # Iterate through time
    # We need history for delays, so we simply index carefully
    
    for t in range(1, n_timepoints):
        dxdt = np.zeros(n_genes)
        
        # 1. Decay
        dxdt -= decay * expression[:, t-1]
        
        # 2. External Input
        dxdt += spontaneous_input[:, t]
        
        # 3. Network Interactions
        # Gene i regulates Gene j with delay
        for i in range(n_genes): # Source
            for j in range(n_genes): # Target
                w = adj_matrix[i, j]
                if w > 0:
                    tau = int(delays[i, j])
                    if t - tau >= 0:
                        # Simple model: if source was high, target gets input
                        # Use a threshold or sigmoid, or linear coupling?
                        # Linear coupling implies correlation. 
                        # We want bursts to propagate.
                        # Let's say input is proportional to Source(t-tau)
                        input_val = expression[i, t - tau]
                        dxdt[j] += w * input_val
                        
        # Update state (Euler integration with dt=1)
        expression[:, t] = expression[:, t-1] + dxdt
        
        # Add noise
        expression[:, t] += np.random.normal(0, noise_level, n_genes)
        
        # Ensure non-negative (gene expression count)
        expression[:, t] = np.maximum(expression[:, t], 0)
        
    return expression
