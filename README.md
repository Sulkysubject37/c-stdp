# C-STDP

**Causal-STDP (C-STDP):**  
Event-Timing based causal gene regulatory network inference for COVID-19 transcriptomics.

## Structure

```
.
├── data
│   ├── raw
│   └── processed
├── docs
├── notebooks
├── pipelines
├── src
│   ├── cstpd
│   └── utils
├── tests
│   ├── unit
│   └── pipeline
├── visuals
└── README.md
```

## Setup

Install requirements:
```
pip install -r requirements.txt
```

## Modules

- `src/cstpd`: Core STDP + GRN algorithms  
- `src/utils`: Helpers (data, encoding, evaluation)  
- `pipelines`: End-to-end workflow scripts  
- `notebooks`: EDA + demos  
- `tests`: Unit + pipeline tests  
```
