#!/bin/bash
set -e

# C-STDP End-to-End Pipeline
# --------------------------

echo "=========================================="
echo "   C-STDP: Causal Spike-Timing Dependent  "
echo "   Plasticity for GRN Inference           "
echo "=========================================="

# Use python3 directly (assuming environment is active or system python is sufficient)
PYTHON="python3"

echo "[1/6] Fetching Data..."
$PYTHON scripts/data_prep/fetch_data.py

echo "[2/6] Preparing Immune Context..."
$PYTHON src/utils/gene_id_mapping.py
$PYTHON src/utils/immune_gene_sets.py

echo "[3/6] Preprocessing Datasets..."
$PYTHON src/utils/preprocess.py

echo "[4/6] Inferring GRNs (Vectorized C-STDP)..."
$PYTHON scripts/inference/run_real_data_cstdp.py

echo "[5/6] Annotating Results..."
$PYTHON src/utils/pathway_annotation.py
$PYTHON src/utils/drug_intersection.py

echo "[6/6] Generating Visualizations..."
$PYTHON scripts/visualization/generate_cascade.py

echo "=========================================="
echo "✅ Pipeline Completed Successfully."
echo "   Check visuals/real_data/ for results."
echo "=========================================="