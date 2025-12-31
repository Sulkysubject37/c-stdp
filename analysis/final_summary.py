import pandas as pd
import numpy as np
import os

def generate_final_summary():
    print("--- Generating Final Quantitative Metrics Summary ---")
    
    summary_data = []
    
    # 1. Synthetic Stability (from parameter_sensitivity.csv)
    if os.path.exists("analysis/parameter_sensitivity.csv"):
        df = pd.read_csv("analysis/parameter_sensitivity.csv")
        best_f1 = df["F1"].max()
        sparse_regime = df[df["Ratio_An_Ap"] >= 1.0]
        avg_precision_sparse = sparse_regime["Precision"].mean()
        
        summary_data.append({"Category": "Synthetic", "Metric": "Max F1 Score", "Value": f"{best_f1:.4f}"})
        summary_data.append({"Category": "Synthetic", "Metric": "Avg Precision (Sparse Regime)", "Value": f"{avg_precision_sparse:.4f}"})
        
    # 2. Causality Controls (from delay_stress_test.csv)
    if os.path.exists("analysis/delay_stress_test.csv"):
        df = pd.read_csv("analysis/delay_stress_test.csv")
        summary_data.append({"Category": "Causality", "Metric": "Direction Accuracy (Overall)", "Value": f"{df['Dir_Accuracy'].mean():.4f}"})
        summary_data.append({"Category": "Causality", "Metric": "Max Detectable Delay (Tau=10)", "Value": f"{df[df['F1'] > 0]['Delay'].max()} steps"})

    # 3. Baseline Comparison (from baseline_comparison.csv)
    if os.path.exists("analysis/baseline_comparison.csv"):
        df = pd.read_csv("analysis/baseline_comparison.csv")
        stdp_f1 = df[df["Method"] == "STDP"]["F1"].values[0]
        granger_f1 = df[df["Method"] == "Granger"]["F1"].values[0]
        summary_data.append({"Category": "Baseline", "Metric": "STDP vs Granger F1 Ratio", "Value": f"{stdp_f1 / (granger_f1 + 1e-6):.2f}x"})

    # 4. Real Data (Values extracted from logs/manual logic)
    # Cohort Overlap (Jaccard): 0.0383
    # Regulator Correlation: 0.4043
    summary_data.append({"Category": "Real Data", "Metric": "Cohort Jaccard Overlap", "Value": "0.0383"})
    summary_data.append({"Category": "Real Data", "Metric": "Regulator Consistency (Corr)", "Value": "0.4043"})
    
    # 5. Output
    summary_df = pd.DataFrame(summary_data)
    print("\n--- FINAL QUANTITATIVE SUMMARY ---")
    print(summary_df.to_string(index=False))
    
    summary_df.to_csv("analysis/final_metrics_summary.csv", index=False)
    
    # Generate Markdown table for report
    md_table = summary_df.to_markdown(index=False)
    with open("analysis/final_metrics_summary.md", "w") as f:
        f.write(md_table)

if __name__ == "__main__":
    generate_final_summary()
