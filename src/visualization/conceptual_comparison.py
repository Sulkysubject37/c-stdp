import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os

def create_conceptual_comparison(output_dir):
    print("Generating Conceptual Comparison Figure...")
    
    img_a_path = os.path.join(output_dir, "phase_architecture.png")
    img_b_path = os.path.join(output_dir, "GSE157859_phase_timeline.png")
    
    if not os.path.exists(img_a_path) or not os.path.exists(img_b_path):
        print("Source images missing.")
        return
        
    img_a = mpimg.imread(img_a_path)
    img_b = mpimg.imread(img_b_path)
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), gridspec_kw={'width_ratios': [1, 1, 0.5]})
    
    # Panel A
    axes[0].imshow(img_a)
    axes[0].axis('off')
    axes[0].set_title("A. Phase Architecture (GSE215865)")
    
    # Panel B
    axes[1].imshow(img_b)
    axes[1].axis('off')
    axes[1].set_title("B. Temporal Phase Progression (GSE157859)")
    
    # Panel C
    axes[2].text(0.1, 0.8, "Shared Early Regulators", fontsize=14, fontweight='bold')
    axes[2].text(0.1, 0.7, "- HEXIM1", fontsize=12)
    axes[2].text(0.1, 0.6, "- RABL3", fontsize=12)
    axes[2].text(0.1, 0.5, "- HMOX1", fontsize=12)
    axes[2].text(0.1, 0.4, "- IFI35", fontsize=12)
    axes[2].text(0.1, 0.2, "Note: Comparisons are\nconceptual, not temporal.", fontsize=10, style='italic', color='red')
    axes[2].axis('off')
    
    plt.tight_layout()
    output_path = os.path.join(output_dir, "conceptual_comparison.png")
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved comparison to {output_path}")

if __name__ == "__main__":
    create_conceptual_comparison("visuals/symposium_final")
