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
