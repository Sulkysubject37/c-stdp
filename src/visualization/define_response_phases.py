import pandas as pd
import numpy as np
import os
import sys

def define_phases(rank_file, output_dir):
    print("Defining Immune Response Phases (Adaptive)...")
    
    if not os.path.exists(rank_file):
        print(f"Error: Rank file not found at {rank_file}")
        sys.exit(1)
        
    df = pd.read_csv(rank_file)
    
    # Adaptive Phase Definition due to T=0 clump (85% genes)
    # Phase I: Immediate (OnsetIndex == 0)
    # Phase II, III, IV: Tertiles of the remaining genes (OnsetIndex > 0)
    
    # Create a mask for T=0
    immediate_mask = df['OnsetIndex'] == 0
    delayed_df = df[~immediate_mask].copy()
    
    # Re-rank delayed genes for binning
    delayed_df['DelayedPercentile'] = delayed_df['Percentile'].rank(pct=True) * 100
    
    # Assign Phase I
    df['Phase'] = 'Unknown'
    df.loc[immediate_mask, 'Phase'] = 'Phase I'
    
    # Assign Phase II, III, IV based on delayed percentile
    # II: 0-33%, III: 33-66%, IV: 66-100% of delayed
    # We map back to original indices
    
    phase_ii_idx = delayed_df[delayed_df['DelayedPercentile'] <= 33.3].index
    phase_iii_idx = delayed_df[(delayed_df['DelayedPercentile'] > 33.3) & (delayed_df['DelayedPercentile'] <= 66.6)].index
    phase_iv_idx = delayed_df[delayed_df['DelayedPercentile'] > 66.6].index
    
    df.loc[phase_ii_idx, 'Phase'] = 'Phase II'
    df.loc[phase_iii_idx, 'Phase'] = 'Phase III'
    df.loc[phase_iv_idx, 'Phase'] = 'Phase IV'
    
    # Validation
    counts = df['Phase'].value_counts()
    total = len(df)
    print("\nAdaptive Phase Distribution:")
    print(counts)
    
    order = ['Phase I', 'Phase II', 'Phase III', 'Phase IV']
    for phase in order:
        count = counts.get(phase, 0)
        pct = (count / total) * 100
        print(f"  {phase}: {count} genes ({pct:.1f}%)")
        
    # Check if any (delayed) phase is empty
    if len(delayed_df) > 0 and (counts.get('Phase II', 0) == 0):
         print("WARNING: Splitting delayed genes failed.")
            
    output_path = os.path.join(output_dir, "phase_definitions.csv")
    df.to_csv(output_path, index=False)
    print(f"Phase definitions saved to {output_path}")

if __name__ == "__main__":
    define_phases(
        "visuals/symposium/GSE215865/activation_ranks.csv",
        "visuals/symposium_final"
    )
