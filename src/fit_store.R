# One computation per heavy fit, shared by every consumer in the process.
#
# The output scripts no longer fit models themselves: they ask this store,
# which computes a CANONICAL recipe on first request and hands the same fitted
# objects to every later requester in the process. run_all.R sources the
# scripts into a single process for exactly this reason, so the Box-Cox
# profiles and SFA fits -- minutes each across eleven benchmarks -- are
# computed once instead of once per script. Run alone, a script still gets
# everything it asks for, computed fresh, so every script remains standalone.
#
# Nothing is written to disk: the store lives and dies with the process, so a
# cached fit can never go stale against edited code or data. The price is that
# a fresh process always pays for its own fits, which is the deliberate
# trade-off against silently reusing yesterday's estimates.
#
# The recipes are canonical on purpose. Figure and table code used to fit the
# same model with slightly different seeding -- regression_tables.R seeded the
# envelope and Pareto-grid Box-Cox profiles from S's profiled lambdas while
# the plot scripts profiled unseeded -- so on a multimodal profile the figure
# and the table beside it could disagree about the lambdas. Here there is one
# recipe per fit: every family's profile is seeded from S's, the cheap one,
# matching what the table code did.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("boxcox_frontier.R")   # accuracy-direction stack, fit_specs, viz
src_source("cost_frontier.R")     # cost-direction stack

# Guarded like .fit_cluster_env (paths.R): library files are re-sourced many
# times per process, and an unguarded assignment would empty the store on
# each pass -- which under run_all.R would mean once per script, recreating
# the very duplication the store exists to remove.
if (!exists(".fit_store", inherits = FALSE))
  .fit_store <- new.env(parent = emptyenv())

fits_once <- function(name, compute) {
  if (is.null(.fit_store[[name]])) {
    t0 <- Sys.time()
    .fit_store[[name]] <- compute()
    cat(sprintf("[fit store] %s: fitted in %.1fs\n", name,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  .fit_store[[name]]
}

# The store fits on its OWN copy of the data rather than any script's, so a
# fit cannot depend on which script asked first. load_runs() is deterministic,
# so this equals every script's local `d`.
store_runs <- function() fits_once("runs", load_runs)

# S/A/B, linear and quadratic: the full specification grid of fit_specs.R,
# keyed "<family>_<time>" with $fits inside, exactly as fit_all_specs returns.
store_specs <- function() fits_once("specs", function() {
  fit_all_specs(store_runs())
})

# paretologit / paretologitenv: list(lin = , quad = ) of per-benchmark fit
# lists.
store_grid <- function(key) fits_once(paste0(key, "_grid"), function() {
  d <- store_runs()
  bs <- bench_levels(d$benchmark)
  fitter <- if (key == "paretologitenv") fit_pareto_logit_env else
    fit_pareto_logit
  lapply(TIME_FORMS, function(form)
    setNames(lapply(bs, function(b)
      fitter(d[d$benchmark == b, ], formula = form)), bs))
})

# Box-Cox profiles, accuracy direction, per-benchmark. S first (its glm inner
# loop makes its profile the cheap one), and its profiled lambdas seed every
# slower family's search -- a starting point, not a constraint.
store_bc <- function(key) fits_once(paste0("bc_", key), function() {
  seeds <- if (key == "S") list() else
    lapply(store_bc("S"), function(f) unname(attr(f, "bc_lambda"))[1:2])
  fit_bc_by(key, store_runs(), seeds)
})

# Cost-direction duals: list(lin = , quad = ) of per-benchmark fit lists,
# through fit_cost_model (cost_frontier.R), the same recipe for figure and
# table.
store_cost <- function(key) fits_once(paste0(key, "_grid"), function() {
  d <- store_runs()
  bs <- bench_levels(d$benchmark)
  lapply(COST_FORMS, function(form)
    setNames(lapply(bs, function(b)
      fit_cost_model(key, d[d$benchmark == b, ], form)), bs))
})

# Box-Cox profiles, cost direction, on the worker cluster via fit_cost_bc_by.
store_cost_bc <- function(key) fits_once(paste0("bc_", key), function() {
  fit_cost_bc_by(key, store_runs())
})
