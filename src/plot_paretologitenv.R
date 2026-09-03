# The envelope-constrained frontier logit (fit_pareto_logit_env,
# envelope_frontier.R): the Pareto-grid logit's objective minimised subject to
# the envelope's constraints -- at or above every run, monotone in cost and
# date. Drawn exactly as its unconstrained parent is, so the two can be read
# side by side in the viewer:
#
#   plot_paretologit.R  the surface THROUGH the staircase, nothing imposed
#   this                THROUGH the staircase, as closely as ABOVE-the-runs
#                       allows
#
# Where the two agree, the constraints are slack; where they diverge, the
# constraints are doing the work, and the binding report below says which one.
#
# Three specifications: linear and quadratic below, and the Box-Cox alternative
# (fit_bc with key "paretologitenv") in its own section at the end.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("fit_store.R")   # brings the whole fitting stack with it

LEVELS <- seq(0.10, 0.90, by = 0.20)

d <- load_runs()
tbar    <- bench_tbar(d)
dates   <- bench_dates(d)
benches <- bench_levels(d$benchmark)

fits_by_spec <- store_grid("paretologitenv")

# The staircase the fit is tracking, and its iso-accuracy counterpart, drawn
# underneath exactly as the two parent scripts draw them.
steps     <- pareto_curves(d, dates)
iso_steps <- iso_pareto_curves(d, LEVELS)

axis_ranges <- do.call(rbind, lapply(benches, function(b) {
  s <- d[d$benchmark == b, ]
  data.frame(benchmark = b, cost = range(s$cost), value = c(0, 1))
}))
iso_ranges <- do.call(rbind, lapply(benches, function(b) {
  s <- d[d$benchmark == b, ]
  data.frame(benchmark = b, date = range(s$releasedate), cost = range(s$cost))
}))

NOTES_FRONTIER <- c(
  paste("Solid: the Pareto-grid logit refitted subject to the envelope's",
        "constraints -- the surface lies at or above every run and cannot",
        "slope backwards in cost or date."),
  PARETO_STEP_NOTE,
  paste("Where the constraints are slack this matches the unconstrained",
        "Pareto-grid logit; where they bind, the surface is held above the",
        "staircase it is trying to run through."))

NOTES_ISO <- c(
  paste("Contours are accuracy targets from 10% to 90%; a falling contour",
        "means the same performance costs less over time. Dots are observed",
        "runs."),
  ISO_PARETO_NOTE,
  paste("Read off the Pareto-grid logit refitted subject to the envelope's",
        "constraints: at or above every run, monotone in cost and date."),
  paste("Cut at the observed cost range; blanked only where both earlier than",
        "a level's first record and dearer than its dearest, so each contour",
        "reaches its record's start or cost ceiling, whichever is more",
        "generous. Never-achieved levels show no contour."),
  ISO_BRANCH_NOTE)

## ---- the pooled pseudo-benchmark, sixth panel ----------------------------------------
# Exactly as in plot_paretologit.R (see the comment there); this script's
# pooled fits are the envelope-constrained ones.
pd <- pooled_acc_display(d, grid_dates(min(d$releasedate), max(d$releasedate)))
axis_ranges_p <- rbind(axis_ranges, data.frame(
  benchmark = "pooled", cost = range(pd$sa$cost), value = range(pd$spx$acc)))
iso_ranges_p <- rbind(iso_ranges, data.frame(
  benchmark = "pooled", date = range(pd$sa$releasedate),
  cost = range(pd$sa$cost)))
PC <- c("benchmark", "cost", "acc", "year")
IC <- c("benchmark", "releasedate", "cost", "acc")
pts_frontier <- rbind(d[, PC], pd$spx[, PC])
pts_iso      <- rbind(d[, IC], pd$sa[, IC])
POOL_NOTE <- paste(
  "Sixth panel: the five primaries pooled on the anchored ECI capability",
  "scale (2PL: C = logit(a)/alpha_b + D_b, Claude 3.5 Sonnet = 130) with",
  "benchmark fixed effects; its value axis is in ECI points and its curve is",
  "the shared capability surface at the fixed effects' anchored mean.")
POOL_ISO_NOTE <- paste(
  "Sixth panel: the pooled primaries; contours are ECI capability levels,",
  "labelled in ECI points and shaded by position within the pooled",
  "capability range on the shared ramp; dashes their record staircases; dots",
  "keep each run's own-benchmark accuracy colour.")
pool_iso_layers <- function(pool_iso)
  pooled_iso_layers(pool_iso, crng = range(pd$spx$acc), steps = pd$iso_steps)

for (tt in names(TIME_FORMS)) {
  fits <- fits_by_spec[[tt]]
  pf <- store_pooled_acc("paretologitenv")[[tt]]

  curves <- rbind(frontier_curves(fits, d, dates, tbar),
                  pooled_acc_frontier_curves(pf, pd$sa, pd$dates))
  p <- frontier_plot(
    curves, pts_frontier, ranges = axis_ranges_p,
    labels = LABELS_POOLED, free_value = TRUE,
    ylab = "Fitted frontier accuracy",
    notes = c(NOTES_FRONTIER, POOL_NOTE)) +
    pareto_step_layer(rbind(steps, pd$steps), labels = LABELS_POOLED)

  f <- sprintf("paretologitenv_%s.png", tt)
  ggsave(out_path(f), p, width = 10, height = fig_height(length(benches)), dpi = 200,
         device = ragg::agg_png)
  cat("wrote", f, "\n")

  iso <- iso_acc_curves(fits, d, tbar, levels = LEVELS, cost_cap = iso_steps)
  pi <- iso_acc_plot(
    iso, pts_iso, ranges = iso_ranges_p, labels = LABELS_POOLED,
    notes = c(NOTES_ISO, POOL_ISO_NOTE)) +
    iso_pareto_layer(iso_steps, labels = LABELS_POOLED) +
    pool_iso_layers(pooled_acc_iso_curves(pf, pd$sa, pd$levels,
                                          cost_cap = pd$iso_steps))

  fi <- sprintf("isoaccuracy_paretologitenv_%s.png", tt)
  ggsave(out_path(fi), pi, width = 10, height = fig_height(length(benches)), dpi = 200,
         device = ragg::agg_png)
  cat("wrote", fi, "\n")
}

## ---- the Box-Cox specification ------------------------------------------------------
# fit_bc with key "paretologitenv": the Box-Tidwell family, its lambdas
# profiled against this model's own objective -- the probability-scale
# quasi-likelihood on the grid, under the envelope's constraints. The contours
# have one branch (phi is monotone), so ISO_BRANCH_NOTE is replaced by the BC
# note, exactly as in the parent scripts.
NOTES_BC <- paste("Box-Tidwell: logit accuracy is linear in phi(cost),",
                  "phi(years since mid-2020) and their product, both lambdas",
                  "profiled against the constrained probability-scale",
                  "quasi-likelihood. Only the REGRESSORS are transformed --",
                  "the response keeps the plain logit link.")
fits_bc <- store_bc("paretologitenv")
pf_bc <- store_pooled_acc_bc("paretologitenv")

curves <- rbind(frontier_curves(fits_bc, d, dates, tbar),
                pooled_acc_frontier_curves(pf_bc, pd$sa, pd$dates))
p <- frontier_plot(
  curves, pts_frontier, ranges = axis_ranges_p,
  labels = LABELS_POOLED, free_value = TRUE,
  ylab = "Fitted frontier accuracy",
  notes = c(NOTES_FRONTIER, NOTES_BC, POOL_NOTE)) +
  pareto_step_layer(rbind(steps, pd$steps), labels = LABELS_POOLED)
ggsave(out_path("paretologitenv_bc.png"), p, width = 10,
       height = fig_height(length(benches)), dpi = 200, device = ragg::agg_png)
cat("wrote paretologitenv_bc.png\n")

iso <- iso_acc_curves(fits_bc, d, tbar, levels = LEVELS, cost_cap = iso_steps)
pi <- iso_acc_plot(
  iso, pts_iso, ranges = iso_ranges_p, labels = LABELS_POOLED,
  notes = c(head(NOTES_ISO, -1), NOTES_BC, POOL_ISO_NOTE)) +
  iso_pareto_layer(iso_steps, labels = LABELS_POOLED) +
  pool_iso_layers(pooled_acc_iso_curves(pf_bc, pd$sa, pd$levels,
                                        cost_cap = pd$iso_steps))
ggsave(out_path("isoaccuracy_paretologitenv_bc.png"), pi, width = 10,
       height = fig_height(length(benches)), dpi = 200, device = ragg::agg_png)
cat("wrote isoaccuracy_paretologitenv_bc.png\n")

## ---- which constraint is doing the work --------------------------------------------
# slack_envelope ~ 0: the surface touches a run -- the ABOVE-the-runs side
# binds. slack_mono ~ 0: the unconstrained fit wanted to slope backwards
# (simpleqa's time slope, for instance) and monotonicity is what holds it flat
# -- read the fitted rate there as the constraint's, not the data's.
for (tt in names(TIME_FORMS)) {
  cat(sprintf("\n== %s: which constraints bind ==\n", tt))
  cat(sprintf("%-10s %8s %12s %12s\n", "bench", "n touch", "slack env",
              "slack mono"))
  for (b in benches) {
    f <- fits_by_spec[[tt]][[b]]
    cat(sprintf("%-10s %8d %12.2e %12.2e\n", b,
                sum(f$env_slack < 0.01), f$slack_envelope, f$slack_mono))
  }
}
