import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Optional

def plot_raster(spike_trains: List[np.ndarray], 
                time_range: Tuple[float, float],
                ax: Optional[plt.Axes] = None):
    """
    Plots a raster plot of spike trains.
    
    Args:
        spike_trains: List of spike times for each gene.
        time_range: (start, end) of the plot.
        ax: Matplotlib axes. If None, creates new.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
        
    n_genes = len(spike_trains)
    
    for i, spikes in enumerate(spike_trains):
        # Filter spikes within range
        valid_spikes = spikes[(spikes >= time_range[0]) & (spikes <= time_range[1])]
        ax.vlines(valid_spikes, i - 0.4, i + 0.4, color='black')
        
    ax.set_ylim(-1, n_genes)
    ax.set_xlabel("Time")
    ax.set_ylabel("Gene Index")
    ax.set_title("Spike Raster Plot")
    ax.set_yticks(range(n_genes))

def calculate_adaptive_thresholds(data: np.ndarray, 
                                  sigma: float = 2.0) -> np.ndarray:
    """
    Calculates thresholds for spike detection based on signal statistics.
    threshold = mean(derivative) + sigma * std(derivative)
    
    Args:
        data: Gene expression data (n_genes, n_time).
        sigma: Number of standard deviations.
        
    Returns:
        thresholds: Shape (n_genes,)
    """
    n_genes = data.shape[0]
    thresholds = np.zeros(n_genes)
    
    # Calculate global or per-gene derivative stats
    # Using simple finite difference along axis 1
    derivatives = np.gradient(data, axis=1)
    
    for i in range(n_genes):
        d = derivatives[i, :]
        mean_d = np.mean(d)
        std_d = np.std(d)
        thresholds[i] = mean_d + sigma * std_d
        
    return thresholds
