import numpy as np
from typing import Dict

def calculate_metrics(ground_truth_adj: np.ndarray, 
                      inferred_adj: np.ndarray, 
                      threshold: float = 0.0) -> Dict[str, float]:
    """
    Calculates classification metrics for GRN inference.
    
    Args:
        ground_truth_adj: True binary or weighted matrix.
        inferred_adj: Inferred weighted matrix.
        threshold: Threshold for binarizing inferred weights.
        
    Returns:
        Dictionary of metrics: precision, recall, f1, accuracy, shd.
    """
    # Binarize
    true_bin = (np.abs(ground_truth_adj) > 0).astype(int)
    pred_bin = (np.abs(inferred_adj) > threshold).astype(int)
    
    # Flatten
    y_true = true_bin.flatten()
    y_pred = pred_bin.flatten()
    
    tp = np.sum((y_true == 1) & (y_pred == 1))
    fp = np.sum((y_true == 0) & (y_pred == 1))
    fn = np.sum((y_true == 1) & (y_pred == 0))
    tn = np.sum((y_true == 0) & (y_pred == 0))
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
    accuracy = (tp + tn) / len(y_true)
    
    # Structural Hamming Distance
    # SHD = FP + FN (for binary graphs)
    shd = fp + fn
    
    return {
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "accuracy": accuracy,
        "shd": shd
    }

def evaluate_directionality(ground_truth_adj: np.ndarray, 
                            inferred_adj: np.ndarray,
                            threshold: float = 0.0) -> float:
    """
    Checks if correctly identified edges have correct direction.
    Only considers edges present in both Ground Truth and Prediction (TP).
    
    Args:
        ground_truth_adj: True matrix.
        inferred_adj: Inferred matrix.
        threshold: Cutoff.
        
    Returns:
        Fraction of TPs with correct direction.
        (Note: In adj matrix, (i,j) implies i->j. If both (i,j) are 1, it matches.)
        Since we compare adjacency matrices directly, standard Precision/Recall 
        already accounts for direction (i->j is different from j->i).
        
        This function calculates "Reverse Edge Rate" among False Positives? 
        Or just reiterates standard accuracy?
        
        Let's calculate: Of the edges that exist in GT (undirected), 
        how many were inferred with correct direction?
    """
    # Simply return precision as our adjacency matrices are directed.
    # A true positive in adj[i,j] means correct direction.
    # If I inferred j->i but truth is i->j, 
    # adj[j,i] is FP, adj[i,j] is FN.
    
    # Let's return a metric specific to reverse edges
    true_edges = np.argwhere(np.abs(ground_truth_adj) > 0)
    
    correct_dir = 0
    total_true_edges_found = 0
    
    for i, j in true_edges:
        # Check if we found an edge between i and j
        # Found i->j
        found_fwd = np.abs(inferred_adj[i, j]) > threshold
        # Found j->i
        found_rev = np.abs(inferred_adj[j, i]) > threshold
        
        if found_fwd:
            correct_dir += 1
            total_true_edges_found += 1
        elif found_rev:
            # Found edge but wrong direction
            total_true_edges_found += 1
            
    if total_true_edges_found == 0:
        return 0.0
        
    return correct_dir / total_true_edges_found
