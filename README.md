# Reproduction materials — EDGE goodness-of-fit paper

Simulation code, results, and figure scripts for:

> **EDGE: A Closed-Form Directed Goodness-of-Fit Test for Sparse Logistic Regression**
> Ebrahim Khaled Ebrahim and Ahmed El-Kotory (submitted to *Statistics and Computing*, 2026).

**EDGE** (Ebrahim Directed Goodness-of-fit Evaluation) keeps the grouped, sparse-data-robust
construction of the Hosmer–Lemeshow / Ebrahim–Farrington tests, but *directs* its few degrees of
freedom onto a low-dimensional basis of calibration-curve shapes by projecting the grouped
standardized residuals onto that basis. Its null distribution is a weighted sum of χ² variables with
an exact, closed-form (Ω-projection) calibration — one eigen-decomposition, no refit, no resampling,
and no tuning.

## The test and comparators

EDGE and every partition-family comparator are implemented in the R package **`ebrahim.gof`**, on CRAN:

```r
install.packages("ebrahim.gof")
```

<https://CRAN.R-project.org/package=ebrahim.gof>

(The exported function for the EDGE test is `def.gof()`; `run.all.gof()` returns the full panel of tests.)

## Datasets

All real datasets are **public** and loaded from named R packages — none are redistributed here,
except the UMARU IMPACT Study (UIS) data, included as `results/uis_data.rds`:

| dataset | source |
|---|---|
| low birth weight (`birthwt`) | `MASS` |
| `kyphosis` | `rpart` |
| `nodal` | `boot` |
| vasoconstriction (`vaso`) | `robustbase` |
| ICU (`icu`), GLOW (`glow500`) | `aplore3` |
| UIS (UMARU IMPACT Study), n = 575 | `results/uis_data.rds` |
| flour-beetle mortality (Bliss 1935) | entered in script |

## Contents

### `R/` — scripts

| script | role |
|---|---|
| `_dgp_library.R` | data-generating processes for all scenarios |
| `_harness.R` | parallel simulation scaffolding |
| `_ek_theme.R` | shared figure theme |
| `_proj_test.R` | reimplementation of the Escanciano–Liu covariate-space projection test |
| `grid_*.R` | power / size / boundary simulation drivers |
| `bench_compute.R`, `bench_slow_timing.R` | timing benchmarks |
| `null_calibration_checks.R` | KS/AD null-uniformity + Imhof-vs-Satterthwaite checks |
| `make_headline_recount.R` | the 19-of-22 cell-by-cell census |
| `uis_reproduce.R` | UIS real-data vignette |
| `_fig1..8_*.R` | figures 1–8 |
| `run_all.R` | master driver |

### `results/` — result files

| file | feeds |
|---|---|
| `sim_power_broad.csv` | power table + Figs 2/4/6 + the 19-of-22 headline |
| `sim_null.csv`, `sim_null_pvalues.csv.gz` | size table + min-p multiplicity |
| `null_ks_table.csv`, `null_bootstrap_agreement.csv`, `imhof_satterthwaite.csv` | null-calibration checks |
| `sim_g_sensitivity.csv` | grouping-sensitivity sweep |
| `sim_edge_loses.csv` | honest-boundary scenario |
| `bench_*_summary.csv`, `bench_compute_scaling_fits.csv`, `bench_stukel_failure.csv` | timings + Stukel separation |
| `uis_gof_reproducible.csv`, `uis_realdata_gof.csv` | UIS vignette |
| `headline_recount.csv` | 19-of-22 enumeration |

### `figures/` — the paper's figures (PDF)

## Reproducing

With the R packages above installed, from the repo root:

```r
source("R/run_all.R")
```

Randomness uses fixed `L'Ecuyer-CMRG` seeds, so every table and figure regenerates exactly.
Benchmarks were run on an AMD Ryzen 9 3900X (12 cores), 32 GB RAM, Windows 11, R 4.4.1.

**Note:** `results/sim_null_pvalues.csv.gz` is gzip-compressed (≈128 MB uncompressed). `read.csv()`
reads it directly through a `gzfile()` connection — e.g. `read.csv("results/sim_null_pvalues.csv.gz")`.

## License

MIT — see [LICENSE](LICENSE).
