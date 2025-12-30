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
