# Reproduction materials — EDGE goodness-of-fit paper

[![DOI](https://zenodo.org/badge/1292650009.svg)](https://doi.org/10.5281/zenodo.21247541)

Simulation code, results, and figure scripts for:

> **A directed goodness-of-fit test for detecting smooth miscalibration in logistic risk-prediction models**
> Ebrahim Khaled Ebrahim and Ahmed El-Kotory (2026). Preprint: arXiv:2608.20511.

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
| `_pstar_giviti_harness.R` | shared harness for the 2026-09 runs (HL verified vs `ResourceSelection`, Stukel refit, GiViTI) |
| `run_A_pstar_giviti.R`, `run_C_ef_gate.R` | covariate-count sweep + GiViTI head-to-head; the HL_F-vs-HL gate |
| `run_D_probit_G_giviti.R`, `run_E_probit_highG.R` | probit vs n; G sweep to 60 and to 80 |
| `run_F_grule.R`, `run_G_validate_rule.R`, `analyse_paired.R` | the m = n/G sweep, the hold-out of the G = n/25 rule (GMAX = 400), McNemar on paired rejections |
| `run_H_overconfidence.R` | external-validation mode (frozen predictions, Ω = I) |
| `run_I_bigdata.R` | Diabetes-130 demonstration (UCI, doi 10.24432/C5230J; data not redistributed here) |
| `run_J_weights_tinym.R` | Theorem 2 check: plug-in weights and size down to 2 observations per group |
| `run_K_basis_score.R`, `run_K2_ao_basis.R` | score-form EDGE and the Aranda–Ordaz link-family basis |
| `probe_probit.R` | population decomposition of the probit signal |
| `emit_tables_A.R`, `emit_tables_FGHI.R`, `emit_tables_K2.R` | the paper's Tables 4–5, 7, 8, 10 (S1) and 9, generated from the CSVs |
| `_fig9_grule.R`, `_fig10_diabetes.R` | figures 9 and 10 (Figures 7 and 8 of the SiM numbering) |

### `results/` — result files

| file | feeds |
|---|---|
| `sim_power_broad.csv` | power table + Figs 2/4/6 + the 19-of-22 headline |
| `sim_null.csv`, `sim_null_pvalues.csv.gz` | size table + min-p multiplicity |
| `null_ks_table.csv`, `null_bootstrap_agreement.csv`, `imhof_satterthwaite.csv` | null-calibration checks |
| `sim_g_sensitivity.csv` | grouping-sensitivity sweep |
| `sim_edge_loses.csv` | honest-boundary scenario |
| `bench_*_summary.csv`, `bench_compute_scaling_fits.csv`, `bench_stukel_failure.csv` | timings + Stukel separation |
| `proj_power_grid.csv`, `proj_power_grid_pvalues.csv`, `proj_timing.txt` | EDGE-vs-projection power head-to-head (24 cells, §5) + per-rep p-values + projection timing |
| `uis_gof_reproducible.csv`, `uis_realdata_gof.csv` | UIS vignette |
| `headline_recount.csv` | 19-of-22 enumeration |
| `runA_pstar_giviti_*.csv`, `runC_ef_gate_*.csv` | covariate-count sweep, GiViTI, the HL_F gate |
| `runD_probit_G_giviti.csv`, `runE_probit_highG_*.csv` | probit and G sweeps |
| `runF_grule_*.csv`, `runG_validate_*.csv`, `paired_edge_vs_rivals.csv` | the G rule (Table 7, Fig 7), its hold-out, McNemar. `runG_validate_*_GMAX200.csv` is the first hold-out run, capped at G = 200, kept for the record |
| `runH_overconfidence_summary.csv` | external mode (Table S1) |
| `runI_bigdata_{dev,val,calcurve}.csv` | Diabetes-130 (Table 10, Fig 8) |
| `runJ_weights_tinym.csv` | Theorem 2 check |
| `runK_basis_score_*.csv`, `runK2_ao_basis_*.csv` | basis experiments (Table 9) |

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

## Version

**v1.1.0** (2026-09-13) — the *Statistics in Medicine* version: adds the 2026-09 runs A–K2 with their
per-replicate p-values, the table/figure emitters, and Figures 9–10. Every table and figure of the
manuscript and its Supporting Information regenerates from these files. Previous: **v1.0.2** — final projection head-to-head grid complete (24/24 cells). Adds
`results/proj_power_grid.csv` (EDGE-poly3 vs the Escanciano–Liu projection test vs EF and HL
across all pre-declared scenarios at n = 500 and 1000), its per-replicate p-values
(`proj_power_grid_pvalues.csv`), and the projection timing (`proj_timing.txt`), feeding the
new head-to-head subsection of §5. Previous: v1.0.1 (benchmark environment freeze + raw timings),
v1.0.0 (initial release).

## License

MIT — see [LICENSE](LICENSE).

## Environment

All benchmarks and simulation grids were run on R 4.4.1 (Windows 11, AMD Ryzen 9 3900X); the exact package versions are frozen in `results/sessionInfo.txt`. Per-grid seed bases are set at the top of each `grid_*.R` script (L'Ecuyer-CMRG streams, invariant to core count).
