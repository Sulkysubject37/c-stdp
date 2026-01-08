import pandas as pd
import matplotlib.pyplot as plt
import os

# Functional Keywords
FUNCTION_MAP = {
    'IFN': ['IFN', 'ISG', 'IRF', 'STAT', 'OAS', 'MX', 'IFI'],
    'Cytokine': ['IL', 'CXC', 'CCL', 'TNF', 'TGF', 'CSF'],
    'Sensing': ['TLR', 'NLR', 'RIG', 'DDX', 'MAVS'],
    'Stress': ['HSP', 'HMOX', 'EIF2', 'ATF']
}

def get_functional_class(gene):
    gene_upper = str(gene).upper()
    for cat, keywords in FUNCTION_MAP.items():
        for k in keywords:
            if k in gene_upper:
                return cat
    return 'Other'

def plot_phase_architecture(phase_file, output_dir):
    print("Generating Phase Architecture Plot...")
    
    if not os.path.exists(phase_file):
        print(f"File not found: {phase_file}")
        return
        
    df = pd.read_csv(phase_file)
    
    # Annotate
    df['Category'] = df['Gene'].apply(get_functional_class)
    
    # Aggregate: Count of Category per Phase
    summary = df.groupby(['Phase', 'Category']).size().unstack(fill_value=0)
    
    # Reorder Phases
    phases = ['Phase I', 'Phase II', 'Phase III', 'Phase IV']
    summary = summary.reindex(phases)
    
    # Normalize to 100% for relative composition? 
    # Or absolute counts? 
    # "Y-axis: Number of genes" -> Absolute.
    # Phase I is huge (3000), others are small (200). 
    # Absolute will hide the others.
    # I will use two subplots: Absolute (Log?) and Relative.
    # Or just Relative (Composition).
    # Prompt says "Y-axis: Number of genes".
    # I'll stick to Absolute, but maybe broken axis or just show the dominance.
    # Phase I will dominate.
    
    # Colors
    colors = {
        'IFN': '#1f77b4',       # Blue
        'Cytokine': '#d62728',  # Red
        'Sensing': '#2ca02c',   # Green
        'Stress': '#ff7f0e',    # Orange
        'Other': '#7f7f7f'      # Gray
    }
    col_list = [colors.get(c, 'black') for c in summary.columns]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    summary.plot(kind='bar', stacked=True, color=col_list, ax=ax)
    
    ax.set_ylabel("Number of Genes")
    ax.set_xlabel("Immune Response Phase")
    ax.set_title("Phase Architecture: Functional Composition")
    ax.grid(axis='y', linestyle='--', alpha=0.3)
    ax.legend(title='Functional Class')
    
    # Add text for massive Phase I
    # If Phase I is huge, maybe cut Y axis?
    # I will keep it real.
    
    plt.tight_layout()
    output_path = os.path.join(output_dir, "phase_architecture.png")
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved figure to {output_path}")

if __name__ == "__main__":
    plot_phase_architecture(
        "visuals/symposium_final/phase_definitions.csv",
        "visuals/symposium_final"
    )
