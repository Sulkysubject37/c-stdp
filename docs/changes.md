## Phase 1: Immune-Aware Data Preprocessing (Executed)
- **Step 1:** Detected Ensembl IDs and Symbols in datasets.
- **Step 2:** Implemented `src/utils/gene_id_mapping.py` using `mygene`.
- **Step 3:** Defined Immune Gene Universe in `src/utils/immune_gene_sets.py`. Universe size: 4352 genes.
- **Step 4:** Applied robust immune filtering in `src/utils/preprocess.py` for both datasets.
    - `GSE215865`: Retained 4326 genes.
    - `GSE157859`: Retained 3113 genes.

## Phase 2: Execute C-STDP on Immune-Filtered Data
- **Step 6:** Implemented Derivative-based Z-score Spike Encoding (`sigma=1.5`).
- **Step 7:** Executed C-STDP on top 200 spiking immune genes for both datasets.
    - Optimized with vectorized STDP runner in `scripts/inference/run_real_data_cstdp.py`.
    
### Primary Run: GSE215865 (COVID-19 Blood)
- **Top Regulators:** `ENSG00000186834`, `ENSG00000144840`, `ENSG00000102871`.
- **Adjacency Sparsity:** 0.74.

### Supplementary Run: GSE157859 (Infection Response)
- **Top Regulators:** `ENSG00000100292` (HMOX1), `ENSG00000120949` (IFI35), `ENSG00000245848` (IL1RL1).
- **Adjacency Sparsity:** 0.50.

## Phase 3: Biological Contextualization
- **Step 9:** Pathway Annotation (`src/utils/pathway_annotation.py`)
    - Regulators linked to "Viral Messenger RNA Synthesis", "Interferon Induction", and "Cytokine Signaling".
- **Step 10:** Druggability Intersection (`src/utils/drug_intersection.py`)
    - Annotations found for `HEXIM1`, `RABL3` (GSE215865) and `HMOX1`, `TNFRSF8` (GSE157859).
    - Note: Coverage sparse, primarily ChEMBL hits.

## Phase 4: Visualization
- **Step 11:** Generated Infection Cascade Plots (`src/visualization/temporal_cascade.py`)
    - Visualized top 30 strongest temporal links against gene onset time.
    - Saved to `visuals/real_data/*/temporal_cascade.png`.