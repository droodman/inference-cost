# Logistic surface fitted to the empirical Pareto frontier, sampled on a grid.
#
# The Pareto frontier at date t is
#
#     P_t(c) = max { a_i : c_i <= c, t_i <= t }
#
# fit_pareto_logit() (envelope_frontier.R) evaluates it at every node of the
# SAME fixed (log cost, tc) grid the envelope scores its objective on, and fits
# model S's fractional logit to those sampled values. An earlier version fitted
# the logit to the frontier-defining runs themselves; that summed the objective
# over wherever the Pareto points cluster -- half of them inside 21-30% of the
# cost range, almost none in the dearest deciles -- while the envelope weighted
# the whole rectangle uniformly, so the two models disagreed about weighting
# before they ever disagreed about anything interesting. On a common grid the
# remaining difference is the one worth having:
#
#   envelope  the lowest monotone logistic surface ABOVE the frontier
#   this      the logistic surface THROUGH the frontier, nothing imposed
#
# Consequences worth keeping in view --
#
#   * grid nodes are not observations, so there are NO standard errors here,
#     exactly as for the envelope; point estimates only;
#   * free disposal is not imposed, so the cost slope can come out negative and
#     the iso-accuracy contours blank out where it does;
#   * P is a running maximum: a node's value is the best run at-or-before it in
#     both coordinates, so the fitted surface summarises the frontier, not the
#     cloud of runs beneath it.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("fit_store.R")          # brings the whole fitting stack with it

LEVELS <- seq(0.10, 0.90, by = 0.20)

d <- load_runs()
tbar    <- bench_tbar(d)
dates   <- bench_dates(d)
benches <- bench_levels(d$benchmark)

## ---- the staircase being sampled ------------------------------------------------

# The frontier's jump points -- what the grid values resample. Selected on RAW
# accuracy over ALL runs, which is P_t(c) as defined above. The
# envelope-constrained variant feeds pareto_binding() something slightly
# different (clipped logit, zero scores dropped, because logit(0) = -Inf
# cannot be a constraint); the counts of both sets are reported below rather
# than assumed away.
frontier_rows <- function(s) s[pareto_binding(s$lncost, s$tc, s$acc), , drop = FALSE]

cat("staircase corners, and the grid that resamples them\n")
cat(sprintf("%-6s %8s %9s %11s %10s\n", "bench", "all runs", "corners",
            "constraints", "grid nodes"))
fits_by_spec <- store_grid("paretologit")

for (b in benches) {
  s <- d[d$benchmark == b, ]
  L <- qlogis(clip_acc(s$acc, s$n_samples))
  pos <- which(is.finite(L))
  n_env <- length(pareto_binding(s$lncost[pos], s$tc[pos], L[pos]))
  f <- fits_by_spec$lin[[b]]
  cat(sprintf("%-6s %8d %9d %11d %10d\n", b, nrow(s),
              attr(f, "n_corners"), n_env, attr(f, "n_grid")))
}

# The grid manufactures rows, not information: 4000 nodes resampling 27 corners
# are still 27 data points. Warn where the corner count is thin relative to the
# quadratic's six coefficients, since the glm itself will no longer complain.
for (b in benches) {
  nc <- attr(fits_by_spec$lin[[b]], "n_corners")
  if (nc < 8)
    warning(sprintf("%s staircase has only %d corners; the quadratic fit there rests on them however many grid nodes resample them",
                    b, nc))
}

## ---- figures ----------------------------------------------------------------------

NOTES_FRONTIER <- c(
  paste("Solid: model S's logit fitted to P_t(c) = max{a_i : c_i <= c, t_i <= t}",
        "sampled on the same fixed (cost, date) grid the envelope uses, so the",
        "two models weight the rectangle identically."),
  PARETO_STEP_NOTE,
  paste("Grid nodes are not observations, so this model reports no standard",
        "errors; the surface runs THROUGH the frontier rather than above it, and",
        "monotonicity is not imposed."))

NOTES_ISO <- c(
  paste("Contours are accuracy targets from 10% to 90%; a falling contour means the",
        "same performance costs less over time. Dots are observed runs."),
  ISO_PARETO_NOTE,
  paste("Read off a logit fitted to the empirical Pareto frontier sampled on a",
        "uniform (cost, date) grid -- the frontier, not the full sample."),
  paste("Cut at the observed cost range; blanked only where both earlier than",
        "a level's first record and dearer than its dearest, so each contour",
        "reaches its record's start or cost ceiling, whichever is more",
        "generous. Never-achieved levels show no contour."),
  paste("Free disposal is not imposed, so a contour also blanks out at dates where",
        "the fitted cost slope is too flat to invert."),
  ISO_BRANCH_NOTE)

iso_ranges <- do.call(rbind, lapply(benches, function(b) {
  s <- d[d$benchmark == b, ]
  data.frame(benchmark = b, date = range(s$releasedate), cost = range(s$cost))
}))

# Pinned the same way plot_paretologitenv.R pins its panels, so the two figures
# can be laid side by side and read panel for panel.
axis_ranges <- do.call(rbind, lapply(benches, function(b) {
  s <- d[d$benchmark == b, ]
  data.frame(benchmark = b, cost = range(s$cost), value = c(0, 1))
}))

# The staircase the fit is tracking, drawn underneath it exactly as
# plot_paretologitenv.R draws it -- and its iso-accuracy counterpart at the
# contour levels, for the swapped-axes figure.
steps <- pareto_curves(d, dates)
iso_steps <- iso_pareto_curves(d, LEVELS)

## ---- the pooled pseudo-benchmark, sixth panel ----------------------------------------
# The five primaries pooled on the anchored ECI scale with benchmark fixed
# effects (pooled_acc_runs / fit_pooled_acc, envelope_frontier.R), drawn as a
# sixth panel: its value axis is in ECI points, its dots and staircases are
# the runs' own capabilities C = logit(a)/alpha_b + D_b. The iso view's
# pooled contours are ECI levels, shaded by position within the pooled
# capability range (pooled_iso_layers, frontier_viz.R) and labelled in place.
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
  "the shared capability surface at its highest benchmark copy (the maximum",
  "fixed effect), the model's counterpart of the capability record.")
POOL_ISO_NOTE <- paste(
  "Sixth panel: the pooled primaries; contours are ECI capability levels,",
  "labelled in ECI points and shaded by position within the pooled",
  "capability range on the shared ramp; dashes their record staircases; dots",
  "keep each run's own-benchmark accuracy colour.")
pool_iso_layers <- function(pool_iso)
  pooled_iso_layers(pool_iso, crng = range(pd$spx$acc), steps = pd$iso_steps)

for (tt in names(TIME_FORMS)) {
  fits <- fits_by_spec[[tt]]
  pf <- store_pooled_acc("paretologit")[[tt]]

  curves <- rbind(frontier_curves(fits, d, dates, tbar),
                  pooled_acc_frontier_curves(pf, pd$sa, pd$dates))
  p <- frontier_plot(
    curves, pts_frontier, ranges = axis_ranges_p,
    labels = LABELS_POOLED, free_value = TRUE,
    ylab = "Fitted frontier accuracy",
    notes = c(NOTES_FRONTIER, POOL_NOTE)) +
    pareto_step_layer(rbind(steps, pd$steps), labels = LABELS_POOLED)

  f <- sprintf("paretologit_%s.png", tt)
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

  fi <- sprintf("isoaccuracy_paretologit_%s.png", tt)
  ggsave(out_path(fi), pi, width = 10, height = fig_height(length(benches)), dpi = 200,
         device = ragg::agg_png)
  cat("wrote", fi, "\n")
}

## ---- the Box-Cox specification ------------------------------------------------------
# fit_bc with key "paretologit" profiles the transform parameters against the
# grid deviance, on the same fixed grid as the other two specifications. The
# contours have one branch (phi is monotone), so ISO_BRANCH_NOTE is replaced by
# the BC note rather than joined by it -- the caption budget is 6 lines.
NOTES_BC <- paste("Box-Tidwell: logit accuracy is linear in phi(cost),",
                  "phi(years since October 2020) and their product, both lambdas",
                  "profiled against the probability-scale quasi-likelihood.",
                  "Only the REGRESSORS are transformed -- the response keeps",
                  "the plain logit link.")
fits_bc <- store_bc("paretologit")
pf_bc <- store_pooled_acc_bc("paretologit")

curves <- rbind(frontier_curves(fits_bc, d, dates, tbar),
                pooled_acc_frontier_curves(pf_bc, pd$sa, pd$dates))
p <- frontier_plot(
  curves, pts_frontier, ranges = axis_ranges_p,
  labels = LABELS_POOLED, free_value = TRUE,
  ylab = "Fitted frontier accuracy",
  notes = c(NOTES_FRONTIER, NOTES_BC, POOL_NOTE)) +
  pareto_step_layer(rbind(steps, pd$steps), labels = LABELS_POOLED)
ggsave(out_path("paretologit_bc.png"), p, width = 10, height = fig_height(length(benches)), dpi = 200,
       device = ragg::agg_png)
cat("wrote paretologit_bc.png\n")

iso <- iso_acc_curves(fits_bc, d, tbar, levels = LEVELS, cost_cap = iso_steps)
pi <- iso_acc_plot(
  iso, pts_iso, ranges = iso_ranges_p, labels = LABELS_POOLED,
  notes = c(head(NOTES_ISO, -1), NOTES_BC, POOL_ISO_NOTE)) +
  iso_pareto_layer(iso_steps, labels = LABELS_POOLED) +
  pool_iso_layers(pooled_acc_iso_curves(pf_bc, pd$sa, pd$levels,
                                        cost_cap = pd$iso_steps))
ggsave(out_path("isoaccuracy_paretologit_bc.png"), pi, width = 10, height = fig_height(length(benches)),
       dpi = 200, device = ragg::agg_png)
cat("wrote isoaccuracy_paretologit_bc.png\n")

## ---- what does moving from the cloud to the frontier change? ------------------------
# Model S on all runs vs the same functional form fitted to the sampled frontier.
# beta_x is the cost slope in logits, beta_t the improvement rate per year at the
# benchmark's reference date.

cat("\ncoefficients: all runs (S) vs Pareto frontier on the grid, linear in time\n")
cat(sprintf("%-6s %9s %9s %9s %9s\n", "bench", "b_x all", "b_x front",
            "b_t all", "b_t front"))
all_fits <- store_specs()$S_lin$fits
for (b in benches) {
  ca <- frontier_coefs(all_fits[[b]])
  cf <- frontier_coefs(fits_by_spec$lin[[b]])
  cat(sprintf("%-6s %9.3f %9.3f %9.3f %9.3f\n", b,
              ca[["bx"]], cf[["bx"]], ca[["bt"]], cf[["bt"]]))
}

# Neither monotonicity condition is imposed by a glm, so check both rather than
# assume them. With the full quadratic each derivative varies over both cost and
# date, so the worst case is a corner of the observed rectangle, not an endpoint
# of one range -- frontier_slope_bounds() does exactly what
# report_quadratic_pathologies() does for the SFA families, on the same helpers.
cat("\nmonotonicity over the observed (cost, date) rectangle\n")
cat(sprintf("%-6s %-5s %12s %12s  %s\n", "bench", "spec", "min dz/dlnc",
            "min dz/dtc", "status"))
for (tt in names(TIME_FORMS)) {
  for (b in benches) {
    co <- frontier_coefs(fits_by_spec[[tt]][[b]])
    s  <- d[d$benchmark == b, ]
    bd <- frontier_slope_bounds(co, s$lncost, s$tc)
    cat(sprintf("%-6s %-5s %12.4f %12.4f  %s\n", b, tt,
                bd[["dcost"]], bd[["dtime"]],
                if (all(bd >= 0)) "ok" else "NON-MONOTONE"))
  }
}
