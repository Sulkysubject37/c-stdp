# Project Changelog

## 2025-12-31

### Step 1: Core STDP Algorithm Implementation

**Changes:**
- Created `src/cstdp/` directory (resolving typo `cstpd` from initial structure).
- Created `src/cstdp/__init__.py`.
- Implemented `src/cstdp/stdp.py` containing the `CausalSTDP` class.

**Details:**
- Implemented `compute_spike_times` using `np.gradient` for derivative calculation.
- Implemented `stdp_update` using the exponential pair-based rule.
- Implemented `run_cstdp` with all-to-all spike pairing.
- Implemented `normalize_weights` with simple clipping [0, w_max].

**Assumptions:**
- `src/cstpd` in the initial file list was a typo; used `src/cstdp`.
- Time points are provided or can be handled by `np.gradient` (supports non-uniform).
- "Bounded" weights interpreted as hard clipping [0, w_max].

### Step 2: Synthetic GRN Simulator

**Changes:**
- Created `src/utils/simulate_grn.py`.

**Details:**
- Implemented `generate_synthetic_grn` to create random sparse directed graphs with delays.
- Implemented `simulate_expression` using a delay-differential equation approximation (Euler method).
- Model includes spontaneous bursts, exponential decay, and linear delayed coupling.

**Assumptions:**
- `dt=1` for simulation steps.
- Delays are integer multiples of `dt`.
- Expression is constrained to be non-negative.

### Step 3: Spike Encoding Utilities

**Changes:**
- Created `src/utils/spike_encoding.py`.

**Details:**
- Implemented `plot_raster` for visualizing spike timing.
- Implemented `calculate_adaptive_thresholds` using mean + sigma * std of derivatives.

**Assumptions:**
- Simple gradient `np.gradient` used for derivative estimation in threshold calculation.

### Step 4: Evaluation Utilities

**Changes:**
- Created `src/utils/evaluate.py`.

**Details:**
- Implemented `calculate_metrics` (Precision, Recall, F1, SHD).
- Implemented `evaluate_directionality` to specifically assess directed edge recovery among found edges.

**Assumptions:**
- Adjacency matrix is directed (A[i, j] means i -> j).

### Step 5: Synthetic Pipeline Test

- **Observation:** Failure: Precision 0.1786, Recall 1.0. High False Positives.
- **Diagnosis:** Strong LTD deficit. `A_pos * tau_pos` (0.5) > `A_neg * tau_neg` (0.4).
- **Fix:** Increased `A_neg` to 0.06.
- **Result (Run 2):** Precision 0.2857 (> chance 0.15), Recall 0.1143. SHD 41. Sparsity achieved (75% weights = 0).
- **Visuals:** Generated `raster_plot.png` and `comparison_matrices.png`.

### Step 6: Explanatory Notebooks

**Changes:**
- Created `notebooks/01_synthetic_demo.ipynb`.
- Created `notebooks/02_interpretability_trace.ipynb`.

**Details:**
- `01_synthetic_demo`: Replicates the synthetic test pipeline in a notebook format for demonstration.
- `02_interpretability_trace`: Traces the weight evolution of a single edge to explain the STDP mechanism.

### Step 7: Download COVID-19 Dataset

**Changes:**
- Created `data/raw/GSE215865/`.
- Downloaded `GSE215865_series_matrix.txt.gz` from GEO.

**Details:**
- Dataset: GSE215865 ("Molecular states during acute COVID-19...").
- Format: GEO Series Matrix (metadata + expression table).
- Contains longitudinal whole-blood RNA-seq from COVID-19 patients.
- **Secondary Dataset:** Downloaded `GSE157859` (GPL16791 and GPL20301). Multi-platform series.

### Step 8: Data Preprocessing

**Changes:**
- Created `src/utils/preprocess.py`.
- Processed `GSE215865` LogCPM matrix.

**Details:**
- Detected `GSE215865_series_matrix.txt.gz` was empty of data; used `suppl/GSE215865_rnaseq_logCPM_matrix.csv.gz` instead.
- Implemented `load_expression_matrix`, `normalize_expression` (Z-score), `select_genes` (Variance).
- Filled NaNs with 0.0 (mean imputation).
- Selected top 50 genes by variance for initial analysis.
- Saved to `data/processed/GSE215865_subset.csv`.
- **Secondary Dataset:** Processed `GSE157859` TPM matrix.
    - Loaded `GSE157859_TPM_matrix.txt.gz`.
    - Applied Log2(x+1) transform (max value > 100).
    - Saved to `data/processed/GSE157859_subset.csv`.
