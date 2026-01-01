import pandas as pd
import numpy as np
import os
import sys
from typing import Optional

def load_expression_matrix(file_path: str) -> pd.DataFrame:
    """
    Parses a Gene Expression file (TPM matrix).
    """
    print(f"Loading {file_path}...")
    try:
        # It's .txt.gz, likely tab-separated
        df = pd.read_csv(file_path, 
                         sep='\t', 
                         index_col=0, 
                         compression='infer')
        
        print(f"Loaded raw data: {df.shape[0]} genes, {df.shape[1]} samples")
        return df
    except Exception as e:
        # Try csv if tab fails?
        try:
            print("Tab separator failed, trying comma...")
            df = pd.read_csv(file_path, sep=',', index_col=0, compression='infer')
            return df
        except:
            raise ValueError(f"Failed to load Matrix: {e}")

def normalize_expression(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalizes expression data (TPM).
    """
    # Check max value. TPM can be large (e.g. 10,000+).
    if df.values.max() > 100:
        print(f"Max value {df.values.max()} > 100. Applying Log2(x+1)...")
        df = np.log2(df + 1)
    else:
        print(f"Max value {df.values.max()} < 100. Assuming already logged.")
        
    # Z-score
    print("Applying Z-score normalization per gene...")
    mean = df.mean(axis=1)
    std = df.std(axis=1)
    std[std == 0] = 1.0
    
    df_norm = df.sub(mean, axis=0).div(std, axis=0)
    df_norm = df_norm.fillna(0.0)
    
    return df_norm

def select_genes(df: pd.DataFrame, n_genes: int = 50) -> pd.DataFrame:
    print(f"Selecting top {n_genes} genes by variance...")
    variances = df.var(axis=1)
    top_genes = variances.nlargest(n_genes).index
    return df.loc[top_genes]

def run_preprocessing(input_path: str, output_path: str, n_genes: int = 50):
    if not os.path.exists(input_path):
        print(f"Error: Input file {input_path} not found.")
        sys.exit(1)
        
    df = load_expression_matrix(input_path)
    df_norm = normalize_expression(df)
    df_final = select_genes(df_norm, n_genes=n_genes)
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df_final.to_csv(output_path)
    print(f"Saved processed data to {output_path}")

if __name__ == "__main__":
    input_file = "data/raw/GSE157859/GSE157859_TPM_matrix.txt.gz"
    output_file = "data/processed/GSE157859_subset.csv"
    
    run_preprocessing(input_file, output_file, n_genes=50)
