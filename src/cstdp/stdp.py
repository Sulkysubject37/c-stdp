import numpy as np
from typing import List, Tuple, Optional, Dict

class CausalSTDP:
    """
    Implements the core Causal-STDP algorithm for inferring directed graphs from time-series.
    
    Based on the mathematical foundation:
    
    Spike Encoding:
        Spike_i(t_k) = 1 if dX_i/dt(t_k) > theta_i else 0
        
    STDP Rule (Pair-based):
        Delta w_ij = 
            A_+ * exp(-(t_j - t_i)/tau_+)  if t_i < t_j (Pre before Post -> Potentiation)
            -A_- * exp(-(t_i - t_j)/tau_-) if t_i > t_j (Post before Pre -> Depression)
            
    Weights are directed, bounded [0, w_max], and normalized.
    """
    
    def __init__(self, 
                 A_pos: float = 0.01, 
                 A_neg: float = 0.01, 
                 tau_pos: float = 10.0, 
                 tau_neg: float = 10.0,
                 w_max: float = 1.0):
        """
        Initialize STDP parameters.
        
        Args:
            A_pos: Amplitude of potentiation.
            A_neg: Amplitude of depression.
            tau_pos: Time constant for potentiation (decays with delay).
            tau_neg: Time constant for depression.
            w_max: Maximum weight bound.
        """
        self.A_pos = A_pos
        self.A_neg = A_neg
        self.tau_pos = tau_pos
        self.tau_neg = tau_neg
        self.w_max = w_max

    def compute_spike_times(self, 
                          data: np.ndarray, 
                          time_points: np.ndarray, 
                          thresholds: np.ndarray) -> List[np.ndarray]:
        """
        Converts continuous gene expression data into spike trains based on the derivative.
        
        Equation:
            Spike_i(t_k) = 1 if dX_i/dt(t_k) > theta_i else 0
            
        Args:
            data: Shape (n_genes, n_time_points). Gene expression X.
            time_points: Shape (n_time_points,). Time values t.
            thresholds: Shape (n_genes,). Threshold theta_i for each gene.
            
        Returns:
            List of size n_genes, where each element is an array of spike times (t_k).
        """
        n_genes, n_points = data.shape
        spike_trains = []
        
        # Calculate derivative dX/dt using finite differences
        # We assume uniform sampling or use gradients if time_points vary
        # Using numpy gradient for potentially non-uniform time
        gradients = np.zeros_like(data)
        for i in range(n_genes):
            gradients[i, :] = np.gradient(data[i, :], time_points)
            
        for i in range(n_genes):
            # Identify indices where derivative exceeds threshold
            spike_indices = np.where(gradients[i, :] > thresholds[i])[0]
            # Convert indices to times
            spike_times = time_points[spike_indices]
            spike_trains.append(spike_times)
            
        return spike_trains

    def stdp_update(self, 
                    w_ij: float, 
                    t_pre: float, 
                    t_post: float) -> float:
        """
        Calculates the weight change delta_w_ij for a single pair of spikes.
        
        Equation:
            Delta w = A_+ * exp(-(t_post - t_pre)/tau_+)  if t_pre < t_post
            Delta w = -A_- * exp(-(t_pre - t_post)/tau_-) if t_pre > t_post
            
        Args:
            w_ij: Current weight (unused in additive STDP but kept for signature).
            t_pre: Spike time of presynaptic neuron (source gene i).
            t_post: Spike time of postsynaptic neuron (target gene j).
            
        Returns:
            Delta w_ij
        """
        dt = t_post - t_pre
        
        if dt > 0:
            # Pre before Post (Causal) -> Potentiation
            return self.A_pos * np.exp(-dt / self.tau_pos)
        elif dt < 0:
            # Post before Pre (Anti-Causal) -> Depression
            return -self.A_neg * np.exp(dt / self.tau_neg) # dt is negative here, so exponent is negative
        else:
            return 0.0

    def normalize_weights(self, weights: np.ndarray) -> np.ndarray:
        """
        Applies bounding and normalization to the weight matrix.
        
        1. Clip weights to [0, w_max].
        2. Optional: Row/Column normalization could be added here.
           Current strategy: Simple clipping as per 'Bounded' requirement.
        
        Args:
            weights: The weight matrix W.
            
        Returns:
            Normalized weight matrix.
        """
        # Bounded [0, w_max]
        return np.clip(weights, 0, self.w_max)

    def run_cstdp(self, 
                  spike_trains: List[np.ndarray], 
                  n_genes: int,
                  initial_weights: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Runs the pair-based STDP learning rule over all spike trains.
        
        Args:
            spike_trains: List of spike time arrays for each gene.
            n_genes: Number of genes.
            initial_weights: Optional starting weights. Defaults to zeros.
            
        Returns:
            Final weight matrix W where W_ij is influence of i on j.
        """
        if initial_weights is None:
            weights = np.zeros((n_genes, n_genes))
        else:
            weights = initial_weights.copy()
            
        # Iterate over all pairs of genes (i -> j)
        for i in range(n_genes):
            for j in range(n_genes):
                if i == j:
                    continue
                    
                spikes_i = spike_trains[i]
                spikes_j = spike_trains[j]
                
                # All-to-all pair interactions
                # In strict STDP, often nearest neighbor is used, but for GRN inference 
                # all pairs within a window are usually considered to capture integrated causality.
                # We implement all-to-all here.
                
                delta_w_sum = 0.0
                
                for t_i in spikes_i:
                    for t_j in spikes_j:
                        delta_w_sum += self.stdp_update(weights[i, j], t_i, t_j)
                        
                weights[i, j] += delta_w_sum
                
        return self.normalize_weights(weights)
