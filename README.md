# C-STDP: Causal Spike-Timing Dependent Plasticity for GRN Inference

**C-STDP** is a computational biology pipeline designed to infer directed, causal Gene Regulatory Networks (GRNs) from transcriptomic time-series data. It adapts the neurobiological principle of Spike-Timing Dependent Plasticity (STDP) to gene expression, treating transcriptional bursts as "spikes" to identify temporal precedence patterns.

## Symposium Pivot: Mapping the COVID-19 Host Response
The project has been specialized for **Infectious Disease Analysis**, specifically mapping the temporal cascade of the host immune response to SARS-CoV-2. By focusing on the *timing* of activation events, C-STDP identifies upstream regulators that drive inflammatory and antiviral pathways.

## Key Features
-   **Sudden-Onset Detection:** Implements derivative-based Z-score encoding to identify transcriptional "bursts" (rate of change) rather than absolute amplitude.
-   **Immune-Aware Preprocessing:** Automated pipeline for gene ID harmonization (Ensembl -> Symbol) and filtering via a curated "Immune Gene Universe" (KEGG, Reactome, GO).
-   **Vectorized C-STDP:** High-performance inference engine capable of processing thousands of gene pairs in minutes.
-   **Temporal Cascade Visualization:** Generates infection timelines showing the chronological order of host factor activation and their causal dependencies.

## Repository Structure

```
C-STDP/
├── data/                   # Raw and processed transcriptomic data (GSE215865, GSE157859)
├── src/                    # Core Library Package
│   ├── core.py             # Main CausalSTDP algorithm
│   ├── utils/              # Preprocessing, mapping, and pathway annotation
│   └── visualization/      # Temporal cascade plotting
├── scripts/                # Functional Scripts
│   ├── inference/          # Spike encoding and vectorized GRN inference
│   └── visualization/      # Cascade generation
├── results/                # Output Artifacts
│   └── visuals/            # Generated cascade plots and heatmaps
├── docs/                   # Documentation
    ├── algorithm.md        # Detailed description of the infectious disease framing
    └── pipeline/           # Historical algorithmic documentation
```

## Methodology
The C-STDP algorithm treats temporal precedence as a proxy for causality. For any pair of genes (A, B), if Gene A consistently activates *before* Gene B, the link A $\to$ B is strengthened. This allows for the reconstruction of the host factor activation timeline.

Detailed algorithmic documentation is available in [docs/algorithm.md](docs/algorithm.md).

## Usage

### 1. Preprocess Immune Subset
```bash
python3 src/utils/preprocess.py
```

### 2. Infer GRN (Vectorized)
```bash
python3 scripts/inference/run_real_data_cstdp.py
```

### 3. Generate Visualizations
```bash
python3 scripts/visualization/generate_cascade.py
```

## Results
The pipeline identifies top regulators associated with:
- **Viral mRNA Synthesis** (e.g., *POLR2C*, *NUP153*)
- **Interferon Induction** (e.g., *NLRC5*, *IFI35*)
- **Stress Response** (e.g., *HMOX1*)

## Author
**MD. Arshad**
*Email:* arshad10867c@gmail.com
*Affiliation:* Jamia Millia Islamia

## License
MIT License. See [LICENSE](LICENSE) for details.