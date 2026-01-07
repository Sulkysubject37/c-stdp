## Phase 1: Immune-Aware Data Preprocessing (Executed)
- **Step 1:** Detected Ensembl IDs with version suffixes in `GSE215865`.
- **Step 2:** Implemented `src/utils/gene_id_mapping.py` using `mygene` to map Ensembl to Symbol. Mapped 44778 genes.
- **Step 3:** Defined Immune Gene Universe in `src/utils/immune_gene_sets.py` sourcing from KEGG, Reactome, GO via `gseapy`. Universe size: 4352 genes.
- **Step 4:** Applied immune filtering in `src/utils/preprocess.py`. Retained 4326 genes (7.34% of original).

## Phase 2: Execute C-STDP on Immune-Filtered Data
- **Step 6:** Implemented Derivative-based Z-score Spike Encoding (`sigma=1.5`).
    - Handled NaNs (imputed 0.0).
    - Detected 475,731 total spikes (Avg ~110/gene).
- **Step 7:** Executed C-STDP on top 200 spiking immune genes.
    - Optimized with vectorized STDP runner in `scripts/inference/run_real_data_cstdp.py`.
    - Identified top regulators (e.g., `ENSG00000186834`).
    - Adjacency sparsity: 0.74.