# Project paths, resolved so scripts work whether the working directory is the
# repo root or src/ (RStudio tends to give one, `Rscript src/foo.R` the other).
#
# Every script starts with the same two-candidate bootstrap:
#   source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
# after which src_source() and the *_path() helpers are location-independent.

.find_root <- function() {
  for (p in c(".", "..", "../..")) {
    if (dir.exists(file.path(p, "data")) && dir.exists(file.path(p, "src")))
      return(normalizePath(p, winslash = "/"))
  }
  stop("cannot locate project root: expected sibling data/ and src/ directories")
}

PROJ_ROOT <- .find_root()

data_path <- function(...) file.path(PROJ_ROOT, "data", ...)
src_path  <- function(...) file.path(PROJ_ROOT, "src", ...)
out_path  <- function(...) {
  d <- file.path(PROJ_ROOT, "output")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  file.path(d, ...)
}

# Source a sibling script by bare name, from anywhere.
src_source <- function(file) source(src_path(file))

## ---- shared fitting cluster ----------------------------------------------------

# PSOCK workers with the whole fitting stack loaded, created on first use and
# reused for the rest of the process. Fits are independent across benchmarks
# (and model families), so the fitting loops hand them to these workers; R and
# its bundled BLAS are single-threaded, so without this a many-core machine
# fits one model at a time. FIT_WORKERS=1 in the environment forces the serial
# path -- worth using when a fit needs debugging, because an error inside a
# PSOCK worker comes back as a one-line summary, not a traceback.
#
# The state lives in an environment guarded by exists(): every script sources
# this file several times over (each library file's bootstrap re-sources it),
# and an unguarded assignment would orphan a running cluster on each pass.
FIT_WORKERS <- {
  w <- suppressWarnings(as.integer(Sys.getenv("FIT_WORKERS", "")))
  if (is.na(w)) max(1L, min(8L, parallel::detectCores(logical = FALSE) - 2L))
  else max(1L, w)
}

if (!exists(".fit_cluster_env", inherits = FALSE))
  .fit_cluster_env <- new.env(parent = emptyenv())

# The cluster, sized min(FIT_WORKERS, n_tasks) at FIRST use and reused as-is
# thereafter (every current caller has one task per benchmark, so the sizes
# agree; a later call wanting more workers than the first simply queues).
# Returns NULL when parallelism is off or pointless, and callers fall back to
# lapply. Workers get FIT_WORKERS=1 so a fit that itself calls a *_by helper
# cannot spawn clusters of clusters, and they die with this process -- the
# finalizer just closes them promptly on a normal exit.
fit_cluster <- function(n_tasks = FIT_WORKERS) {
  n <- min(FIT_WORKERS, n_tasks)
  if (n < 2L) return(NULL)
  if (is.null(.fit_cluster_env$cl)) {
    cl <- parallel::makePSOCKcluster(n)
    parallel::clusterCall(cl, function(root) {
      Sys.setenv(FIT_WORKERS = "1")
      setwd(root)
      source("src/boxcox_frontier.R")   # sources the entire fitting stack
      invisible(NULL)
    }, PROJ_ROOT)
    .fit_cluster_env$cl <- cl
    reg.finalizer(.fit_cluster_env, function(e) {
      if (!is.null(e$cl)) try(parallel::stopCluster(e$cl), silent = TRUE)
    }, onexit = TRUE)
  }
  .fit_cluster_env$cl
}
