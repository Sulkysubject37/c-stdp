# C-STDP: Causal Spike-Timing Dependent Plasticity for GRN Inference

**C-STDP** is a computational biology pipeline designed to infer directed, causal Gene Regulatory Networks (GRNs) from transcriptomic time-series data. It adapts the neurobiological principle of Spike-Timing Dependent Plasticity (STDP) to gene expression, treating transcriptional bursts as "spikes" to identify temporal precedence patterns.

## Project Goal
To move beyond correlation-based association networks and identify **driver genes** (regulators) that temporally precede downstream effects in disease states (e.g., COVID-19).

## Key Features
-   **Event-Based:** Transforms continuous expression data into discrete activation events (spikes).
-   **Causality via Timing:** Infers direction ($A \to B$) based on consistent time delays, not regression.
-   **Interpretable:** Every inferred edge can be traced back to specific, observable event pairs.
-   **Robust:** Validated against synthetic ground truth, noise controls, and permutation tests.

## Repository Structure

```
C-STDP/
├── data/                   # Raw and processed transcriptomic data (GSE215865, GSE157859)
├── src/                    # Core Library Package
│   ├── core.py         # Main CausalSTDP algorithm
│   └── utils/          # Utilities (preprocessing, simulation, evaluation)
├── scripts/                # Functional Scripts
│   ├── data_prep/          # Data download and preprocessing
│   ├── validation/         # Synthetic validation & baseline comparisons
│   ├── inference/          # Real-data spike encoding & GRN inference
│   ├── diagnostics/        # Robustness & sensitivity analysis
│   └── reporting/          # Summary generation
├── results/                # Output Artifacts
│   ├── visuals/            # Generated plots and heatmaps
│   └── *.csv               # Metrics and inferred networks
├── tests/                  # Unit and Integration tests
└── docs/                   # Documentation
    └── pipeline/
        └── algorithmic_approach.md # Detailed method explanation
```

## Installation

1.  Clone the repository.
2.  Set up the Python environment:
    ```bash
    python3 -m venv causal-stdp
    source causal-stdp/bin/activate
    pip install -r requirements.txt
    ```

## Usage

### Run the Full Pipeline
To execute the entire workflow (Synthetic Validation -> Real Data Preprocessing -> Inference -> Diagnostics):
```bash
./run_pipeline.sh
```

### Run Specific Modules
To test the algorithm's robustness:
```bash
python scripts/validation/parameter_sensitivity.py
```

To infer a GRN from a specific dataset:
```bash
python scripts/inference/run_real_data_cstdp.py
```

## Author
**MD. Arshad**
*Email:* arshad10867c@gmail.com
*Affiliation:* Jamia Millia Islamia

## License
MIT License. See [LICENSE](LICENSE) for details.

```
