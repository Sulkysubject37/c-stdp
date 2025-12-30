#!/bin/bash
set -e

# C-STDP End-to-End Pipeline
# --------------------------

echo "=========================================="
echo "   C-STDP: Causal Spike-Timing Dependent  "
echo "   Plasticity for GRN Inference           "
echo "=========================================="

# Activate Environment (Assuming user has set it up or it's active)
# If running in this specific session, we use the venv path explicitly in python calls
# or just assume python points to the venv.
# For robustness in this environment, I'll use the venv path.
PYTHON="casual-stdp/bin/python"

if [ ! -f "$PYTHON" ]; then
    echo "Python venv not found at $PYTHON. Please run setup first."
    exit 1
fi

echo "[1/4] Running Synthetic Validation..."
$PYTHON pipelines/run_synthetic_test.py

echo "[2/4] Preprocessing Real Data..."
# Primary
$PYTHON src/utils/preprocess.py
# Secondary
$PYTHON src/utils/preprocess_secondary.py

echo "[3/4] Encoding Spikes (Real Data)..."
$PYTHON pipelines/run_real_data_spike_encoding.py

echo "[4/4] Inferring GRNs (Real Data)..."
$PYTHON pipelines/run_real_data_cstdp.py

echo "=========================================="
echo "✅ Pipeline Completed Successfully."
echo "   Check visuals/ for results."
echo "=========================================="
