# Run every output-producing script in the analysis.
#
#     Rscript src/run_all.R                     everything, animations included
#     Rscript src/run_all.R --skip-animations   everything except the mp4s
#
# Each script runs in its OWN fresh Rscript process, not source()d into this
# one. That is deliberate: the scripts are written to be standalone (each does
# its own sourcing and fitting), a fresh process is exactly the environment
# they are tested in by hand, and a crash in one cannot leave half-built
# globals lying around to corrupt the next. The cost is that shared fits are
# recomputed in several scripts, which was already true of running them by
# hand and keeps every script's output reproducible in isolation.
#
# All scripts run even if one fails; failures are summarised at the end and
# the exit code is nonzero if there were any, so a wrapper (or a person
# skimming the tail of the log) cannot mistake a partial rebuild for a full
# one.
#
# The runnable scripts, in order. The order only matters for readability of
# the log -- no script consumes another's OUTPUT FILES; everything downstream
# of the data is refitted from the source CSVs each time via prepare_data.R.
# The remaining .R files in src/ (paths.R, prepare_data.R, frontier_viz.R,
# fit_specs.R, panel_frontier.R, fractional_frontier.R, envelope_frontier.R,
# boxcox_frontier.R) are libraries the scripts below source; they produce
# nothing and are not run directly. "Model inference cost.do" is Stata and out
# of scope here.
SCRIPTS <- c(
  "pareto_frontiers.R",    # nonparametric staircase figure
  "plot_frontiers.R",      # SFA + S figures, both views, + fit diagnostics
  "plot_envelope.R",       # envelope figures, both views, + touch diagnostics
  "plot_paretologit.R",    # Pareto-frontier logit figures, both views
  "regression_tables.R",   # HTML + RTF tables in output/tables/
  "animate_frontiers.R"    # quarterly-accumulation movies (slow; skippable)
)

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")

if ("--skip-animations" %in% commandArgs(trailingOnly = TRUE))
  SCRIPTS <- setdiff(SCRIPTS, "animate_frontiers.R")

# The Rscript that is running THIS file, not whatever happens to be on PATH --
# so the children are guaranteed the same R version and library paths.
rscript <- file.path(R.home("bin"), "Rscript")

# Children inherit this process's working directory, and each script's opening
# bootstrap (`source(if (file.exists("src/paths.R")) ... else "paths.R")`) only
# resolves from the root or from src/. Pin the root so run_all works no matter
# where it was launched from.
setwd(PROJ_ROOT)

status <- integer(0)
t_all <- Sys.time()
for (s in SCRIPTS) {
  cat(sprintf("\n===== %s =====\n", s))
  t0 <- Sys.time()
  rc <- system2(rscript, shQuote(file.path(PROJ_ROOT, "src", s)),
                stdout = "", stderr = "")
  status[s] <- rc
  cat(sprintf("----- %s: %s in %.1f min\n", s,
              if (rc == 0) "ok" else sprintf("FAILED (exit %d)", rc),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

cat(sprintf("\n===== summary (%.1f min total) =====\n",
            as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
for (s in names(status))
  cat(sprintf("  %-24s %s\n", s, if (status[s] == 0) "ok" else
    sprintf("FAILED (exit %d)", status[s])))

if (any(status != 0)) {
  cat("\nsome scripts FAILED -- outputs are a partial rebuild\n")
  quit(status = 1)
}
cat("\nall scripts completed\n")
