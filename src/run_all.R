# Run every output-producing script in the analysis.
#
#     Rscript src/run_all.R                     everything, animations included
#     Rscript src/run_all.R --skip-animations   everything except the mp4s
#
# The scripts run in THIS process, each sourced into its own private
# environment. One process, because the scripts draw their fits from the
# shared in-memory store (fit_store.R): the profiled Box-Cox fits and SFA fits
# -- minutes each across the benchmarks -- are computed once and every script
# reuses the same objects, where separate Rscript children each paid for their
# own. (Fits are never cached to disk: the store dies with this process, so a
# stale fit cannot survive a code or data edit.) The private environments keep
# one script's globals from leaking into the next, and each script still runs
# standalone by hand -- it just computes what it asks the store for.
#
# All scripts run even if one fails; failures are summarised at the end and
# the exit code is nonzero if there were any, so a wrapper (or a person
# skimming the tail of the log) cannot mistake a partial rebuild for a full
# one. The isolation is per-script tryCatch rather than per-script process,
# so a crash of R itself (native code, out-of-memory) now aborts the whole
# run -- the price of sharing the fits.
#
# The runnable scripts, in order. Scripts whose fits overlap share them
# through the store, so the order only affects which script's banner the
# fitting time is logged under -- no script consumes another's OUTPUT FILES.
# The remaining .R files in src/ (paths.R, prepare_data.R, frontier_viz.R,
# fit_specs.R, panel_frontier.R, fractional_frontier.R, envelope_frontier.R,
# boxcox_frontier.R, cost_frontier.R, fit_store.R) are libraries the scripts
# below source; they produce nothing and are not run directly. "Model
# inference cost.do" is Stata and out of scope here.
SCRIPTS <- c(
  "pareto_frontiers.R",    # nonparametric staircase figure
  "record_timelines.R",    # cost-record timelines table (HTML + CSV)
  "plot_frontiers.R",      # SFA + S figures, both views, + fit diagnostics
  "plot_paretologit.R",    # Pareto-frontier logit figures, both views
  "plot_paretologitenv.R", # envelope-constrained variant, both views
  "plot_cost_frontier.R",  # cost-direction dual figures, both views
  "plot_surfaces_3d.R",    # interactive 3-D surfaces for the Pareto-grid pair
  "cost_frontier_report.R",# cost-direction dual fits: rate comparison (console)
  "regression_tables.R",   # HTML + RTF tables in output/tables/
  "animate_frontiers.R"    # semiannual-accumulation movies (slow; skippable)
)

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")

if ("--skip-animations" %in% commandArgs(trailingOnly = TRUE))
  SCRIPTS <- setdiff(SCRIPTS, "animate_frontiers.R")

# Each script's opening bootstrap (`source(if (file.exists("src/paths.R"))
# ... else "paths.R")`) only resolves from the root or from src/. Pin the
# root so run_all works no matter where it was launched from.
setwd(PROJ_ROOT)

status <- integer(0)
t_all <- Sys.time()
for (s in SCRIPTS) {
  cat(sprintf("\n===== %s =====\n", s))
  t0 <- Sys.time()
  # A private environment per script: top-level assignments (d, fits, p, ...)
  # land there and cannot collide across scripts, while the libraries the
  # script sources land in the global environment as always -- which is where
  # fit_store.R's store lives, shared by every script in this loop.
  rc <- tryCatch({
    source(file.path(PROJ_ROOT, "src", s),
           local = new.env(parent = globalenv()))
    0L
  }, error = function(e) {
    cat(sprintf("ERROR in %s: %s\n", s, conditionMessage(e)))
    1L
  })
  status[s] <- rc
  cat(sprintf("----- %s: %s in %.1f min\n", s,
              if (rc == 0) "ok" else "FAILED",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

cat(sprintf("\n===== summary (%.1f min total) =====\n",
            as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
for (s in names(status))
  cat(sprintf("  %-24s %s\n", s, if (status[s] == 0) "ok" else "FAILED"))

if (any(status != 0)) {
  cat("\nsome scripts FAILED -- outputs are a partial rebuild\n")
  quit(status = 1)
}
cat("\nall scripts completed\n")
