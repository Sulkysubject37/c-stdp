## Phase 1: Immune-Aware Data Preprocessing (Executed)
- **Step 1:** Detected Ensembl IDs with version suffixes in `GSE215865`.
- **Step 2:** Implemented `src/utils/gene_id_mapping.py` using `mygene` to map Ensembl to Symbol. Mapped 44778 genes.
- **Step 3:** Defined Immune Gene Universe in `src/utils/immune_gene_sets.py` sourcing from KEGG, Reactome, GO via `gseapy`. Universe size: 4352 genes.
- **Step 4:** Applied immune filtering in `src/utils/preprocess.py`. Retained 4326 genes (7.34% of original).