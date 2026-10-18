# ED4DSE — *Experimental Design for Data Science and Engineering*

> A complete computational reproduction of **V. Roshan Joseph, *Experimental Design for Data Science and Engineering*, Chapman & Hall/CRC, 2026** — runnable R code for every section of all 11 chapters, with figures, cached results, and locally maintained companion R packages.

## What this repository contains

```
ed4dse/
├── ED4DSE.pdf            # Full compiled textbook (~56 MB)
├── codes/                # 79 section-level R scripts, named "N.M.R"
│   └── data/             # Cached results: "N.M-plot.RData" per script,
│                         #   run logs (run1–run51.txt), and 5 CSV datasets
├── figures/              # 158 figures: one PDF + one SVG per script
└── Rpackages/            # Three locally maintained/patched CRAN packages
```

## Chapter coverage (79 scripts, all 11 chapters)

| # | Chapter | Topics | Scripts |
|---|---------|--------|---------|
| 1 | Experiments | Response surfaces, simulation vs. physical experiments | 1.4–1.5 |
| 2 | Modeling Techniques | Kriging / Gaussian process regression | 2.1–2.9 |
| 3 | Model-based Designs | Prediction-based designs, maximum-entropy designs | 3.1–3.4 |
| 4 | Space-filling Designs | MaxPro, minimum-energy, Latin hypercube designs | 4.1–4.27 |
| 5 | Representative Points | Uniform designs, support points, uncertainty propagation | 5.1–5.9 |
| 6 | Screening Designs | Sensitivity analysis, Morris, MOFAT | 6.1–6.8 |
| 7 | Sequential Designs | Simulation optimization, Bayesian optimization, inverse design | 7.2–7.6 |
| 8 | Fractional Factorial Designs | Two-level and Bayesian-inspired designs | 8.1–8.5 |
| 9 | Model Calibration | Nonlinear optimal designs, robust designs | 9.1–9.4 |
| 10 | Data Subsampling | Support-point subsampling, data twins | 10.1–10.6 |
| 11 | Data Analysis | Factor selection, twin Gaussian processes | 11.1–11.3 |

## Script conventions

Every script follows the same pattern:

1. `setwd(...)` to the script directory;
2. If `data/N.M-plot.RData` does not exist, run the computation (which can be expensive — Bayesian optimization runs, GP fits, optimal-design searches) and cache the results;
3. Otherwise load the cache and skip the computation;
4. Produce the figure(s) into `figures/` as both PDF and SVG.

This makes every script **re-runnable and resumable**: delete a cache file to force recomputation of that section only.

## Locally maintained packages (`Rpackages/`)

Three companion CRAN packages are kept locally with fixes:

- **SFDesign** (v0.1.5, Shangkun Wang & Roshan Joseph) — space-filling designs: LHD, clustered, maximin, maximum projection, uniform designs (Rcpp-backed);
- **mined** (v1.0-3, Dianpeng Wang & V. Roshan Joseph) — minimum-energy designs (MinED);
- **support** (v0.1.7, Simon Mak) — support points (`sp()`, `sp_seq()`; Mak & Joseph 2018).

## Data

`codes/data/` holds five CSV datasets used across chapters: concrete compressive strength, sorption/desorption experiments, and wind speed data, plus the cached `RData` results and run logs.

## Language

Code and comments are in English. There is no formal build system — scripts are run individually in R.

## Related repositories

See also [dace](https://github.com/wenje-ma/dace) (reading notes for Santner–Williams–Notz) and [singapore](https://github.com/wenje-ma/singapore) (multi-fidelity Bayesian optimization research). Full list on the [profile page](https://github.com/wenje-ma).
