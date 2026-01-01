# C-STDP: Causal Spike-Timing Dependent Plasticity for GRN Inference

**C-STDP** is a computational biology pipeline designed to infer directed, causal Gene Regulatory Networks (GRNs) from transcriptomic time-series data. It adapts the neurobiological principle of Spike-Timing Dependent Plasticity (STDP) to gene expression, treating transcriptional bursts as "spikes" to identify temporal precedence patterns.

## Project Goal
To move beyond correlation-based association networks and identifying **driver genes** (regulators) that temporally precede downstream effects in disease states (e.g., COVID-19).

## Key Features
-   **Event-Based:** Transforms continuous expression data into discrete activation events (spikes).
-   **Causality via Timing:** Infers direction ($A \to B$) based on consistent time delays, not regression.
-   **Interpretable:** Every inferred edge can be traced back to specific, observable event pairs.
-   **Robust:** Validated against synthetic ground truth, noise controls, and permutation tests.

## Repository Structure

```
C-STDP/
├── data/                   # Raw and processed transcriptomic data (GSE215865, GSE157859)
├── src/
│   ├── cstdp/              # Core algorithm implementation (STDP rule)
│   └── utils/              # Preprocessing, simulation, and evaluation tools
├── pipelines/
│   ├── run_cstdp_pipeline.sh       # Full end-to-end workflow
│   ├── run_real_data_cstdp.py      # Inference on real data
│   └── run_synthetic_test.py       # Validation on synthetic data
├── analysis/               # Diagnostic and validation scripts
│   ├── parameter_sensitivity.py    # Stability analysis
│   ├── negative_controls.py        # Null model testing
│   └── visuals/                    # Generated plots and heatmaps
└── docs/
    └── pipeline/
        └── algorithmic_approach.md # Detailed method explanation
```

## Installation

1.  Clone the repository.
2.  Set up the Python environment:
    ```bash
    python3 -m venv casual-stdp
    source casual-stdp/bin/activate
    pip install numpy pandas matplotlib seaborn networkx scipy statsmodels tabulate
    ```

## Usage

### Run the Full Pipeline
To execute the entire workflow (Synthetic Validation -> Real Data Preprocessing -> Inference):
```bash
./pipelines/run_cstdp_pipeline.sh
```

### Run Specific Analysis
To test the algorithm's robustness:
```bash
python analysis/parameter_sensitivity.py
```

To infer a GRN from a specific dataset:
```bash
python pipelines/run_real_data_cstdp.py
```

## License
MIT License.

```