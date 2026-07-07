## _harness.R — shared 24-core parallel scaffolding + incremental CSV for the EDGE grids.
## Sources _dgp_library.R inside every worker. Per-rep set.seed(seed_base + rep) makes results
## INVARIANT to worker count (reproducible on the 24-core box). Honor REPS_SCALE for smoke runs.
suppressMessages(library(parallel))
SIMDIR <- "C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/simulations"
DGP    <- file.path(SIMDIR, "_dgp_library.R")
## EK_NCORES defaults to (logical cores - 1); a self-test / co-running agent may cap it
## by setting env EK_NCORES_CAP (e.g. "2") to avoid oversubscribing the box.
EK_NCORES <- max(1L, detectCores(logical = TRUE) - 1L)
{ .cap <- suppressWarnings(as.integer(Sys.getenv("EK_NCORES_CAP", ""))); if (!is.na(.cap) && .cap >= 1L) EK_NCORES <- min(EK_NCORES, .cap) }

## scale B down for a smoke pass, e.g. REPS_SCALE=0.02; always >= 2 reps
ek_reps <- function(B) max(2L, as.integer(round(B * as.numeric(Sys.getenv("REPS_SCALE", "1")))))

## cluster with ebrahim.gof + the DGP library loaded in every worker, RNG stream set
ek_cluster <- function(seed) {
  cl <- makeCluster(EK_NCORES)
  clusterSetRNGStream(cl, seed)
  clusterExport(cl, "DGP", envir = environment())
  invisible(clusterEvalQ(cl, { suppressMessages(library(ebrahim.gof)); source(DGP); NULL }))
  cl
}

## incremental append: writes the header only when the file does not yet exist
ek_append <- function(df, path) {
  ex <- file.exists(path)
  write.table(df, path, sep = ",", row.names = FALSE, col.names = !ex, append = ex, qmethod = "double")
}

ek_mcse <- function(p, B) sqrt(pmax(0, p * (1 - p)) / B)

## Run ONE cell: B reps in parallel, each seeded set.seed(seed_base + rep). one_rep(rep, cell)
## must return a NAMED numeric vector with identical names on every call. Returns a data.frame
## of B rows x named columns. Give each cell a DISTINCT seed_base (e.g. grid_seed + id*1e5).
ek_run_cell <- function(cl, B, seed_base, one_rep, cell) {
  clusterExport(cl, c("one_rep", "cell", "seed_base"), envir = environment())
  M <- parSapply(cl, seq_len(B), function(rep) { set.seed(seed_base + rep); one_rep(rep, cell) })
  if (is.null(dim(M))) M <- matrix(M, nrow = length(one_rep(1, cell)))
  as.data.frame(t(M))
}
