import pandas as pd
import numpy as np
import os
import sys
from typing import List, Optional

def load_expression_matrix(file_path: str) -> pd.DataFrame:
    """
    Parses a Gene Expression CSV file into a DataFrame.
    
    Args:
        file_path: Path to .csv or .csv.gz file.
        
    Returns:
        DataFrame with Genes as rows and Samples as columns.
    """
    print(f"Loading {file_path}...")
    
    try:
        # Check extension for separator
        sep = ',' if 'csv' in file_path else '\t'
        
        # Read data
        df = pd.read_csv(file_path, 
                         sep=sep, 
                         index_col=0, 
                         compression='infer') # Auto-detect gzip
        
        print(f"Loaded raw data: {df.shape[0]} genes, {df.shape[1]} samples")
        return df
        
    except Exception as e:
        raise ValueError(f"Failed to load Matrix: {e}")

def normalize_expression(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalizes expression data.
    
    Steps:
    1. Log2(x+1) transform (if not already logged).
    2. Z-score normalization per gene (row-wise).
    
    Args:
        df: (Genes x Samples)
        
    Returns:
        Normalized DataFrame.
    """
    # Check max value to guess if log transform is needed
    # LogCPM is usually < 20. Raw counts > 100.
    if df.values.max() > 100:
        print(f"Max value {df.values.max()} > 100. Applying Log2(x+1)...")
        df = np.log2(df + 1)
    else:
        print(f"Max value {df.values.max()} < 100. Assuming already logged.")
        
    # Z-score normalization per gene
    print("Applying Z-score normalization per gene...")
    mean = df.mean(axis=1)
    std = df.std(axis=1)
    
    # Avoid division by zero
    std[std == 0] = 1.0
    
    df_norm = df.sub(mean, axis=0).div(std, axis=0)
    
    # Fill NaNs with 0.0 (mean) to ensure clean data for STDP
    df_norm = df_norm.fillna(0.0)
    
    return df_norm

def select_genes(df: pd.DataFrame, 
                 selection_strategy: str = 'variance', 
                 n_genes: int = 50,
                 target_genes: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Selects a subset of genes.
    """
    if selection_strategy == 'list' and target_genes is not None:
        available = [g for g in target_genes if g in df.index]
        print(f"Selected {len(available)}/{len(target_genes)} target genes found.")
        return df.loc[available]
    
    elif selection_strategy == 'variance':
        print(f"Selecting top {n_genes} genes by variance...")
        variances = df.var(axis=1)
        top_genes = variances.nlargest(n_genes).index
        return df.loc[top_genes]
        
    else:
        return df

def save_processed_data(df: pd.DataFrame, output_path: str):
    """Saves DataFrame to CSV."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path)
    print(f"Saved processed data to {output_path}")

def run_preprocessing(input_path: str, output_path: str, n_genes: int = 20):
    """
    Main pipeline function.
    """
    if not os.path.exists(input_path):
        print(f"Error: Input file {input_path} not found.")
        sys.exit(1)
        
    df = load_expression_matrix(input_path)
    
    # Normalize
    df_norm = normalize_expression(df)
    
    # Select subset 
    df_final = select_genes(df_norm, n_genes=n_genes)
    
    save_processed_data(df_final, output_path)

if __name__ == "__main__":
    # Default to GSE215865 LogCPM
    input_file = "data/raw/GSE215865/GSE215865_rnaseq_logCPM_matrix.csv.gz"
    output_file = "data/processed/GSE215865_subset.csv"
    
    run_preprocessing(input_file, output_file, n_genes=50)