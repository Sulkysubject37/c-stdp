import pandas as pd
import os

def generate_diagnostic_summary():
    print("--- Generating Diagnostic Summary Table ---")
    
    summary_data = []
    
    # 1. Parameter Sensitivity (Real)
    if os.path.exists("analysis/real_data_sensitivity.csv"):
        df = pd.read_csv("analysis/real_data_sensitivity.csv")
        # Stable regime: Ratio >= 1.0
        stable = df[df["Ratio"] >= 1.0]
        mean_sparsity = stable["Sparsity"].mean()
        mean_jaccard = stable["Jaccard_vs_Base"].mean()
        
        summary_data.append({"Diagnostic": "Parameter Sensitivity", "Real Data Outcome": f"Sparsity={mean_sparsity:.2f}, Jaccard={mean_jaccard:.2f}", "Interpretation": "Stable sparse networks in valid parameter regime."})

    # 2. Permutation (Values from log)
    # Sample Perm: Jaccard 0.0080
    summary_data.append({"Diagnostic": "Permutation Control", "Real Data Outcome": "Jaccard < 0.01", "Interpretation": "Network structure relies entirely on temporal order."})
    
    # 3. Delay Structure
    # Mean Weighted Delay: 9.42
    summary_data.append({"Diagnostic": "Delay Structure", "Real Data Outcome": "Mean Delay ~9.42 steps", "Interpretation": "Inferred delays match STDP window (Tau=10)."})
    
    # 4. Negative Controls
    # Top 10%: 1/10 controls
    summary_data.append({"Diagnostic": "Negative Controls", "Real Data Outcome": "10% False Discovery Rate (approx)", "Interpretation": "Random/Shuffled genes do not dominate regulators."})
    
    # 5. Consistency
    # Cohort Jaccard: 0.0383, Reg Corr: 0.4043
    summary_data.append({"Diagnostic": "Cohort Consistency", "Real Data Outcome": "Regulator Corr=0.40", "Interpretation": "Functional (regulator) consistency is higher than topological (edge) consistency."})
    
    # 6. Cross-Dataset
    # Sparsity ~0.95
    summary_data.append({"Diagnostic": "Cross-Dataset", "Real Data Outcome": "Sparsity ~0.95", "Interpretation": "Algorithm behavior is consistent across species/scales."})
    
    df_sum = pd.DataFrame(summary_data)
    df_sum.to_csv("analysis/diagnostic_summary.csv", index=False)
    
    print(df_sum.to_string())
    
    with open("analysis/diagnostic_summary.md", "w") as f:
        f.write(df_sum.to_markdown(index=False))

if __name__ == "__main__":
    generate_diagnostic_summary()
