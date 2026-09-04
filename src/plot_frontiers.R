# Parametric-specification figures: BOTH views from ONE set of fits.
#
#   frontier_progression_<fam>_<spec>.png  fitted frontier by cost, one curve
#                                          per half-year, colour = elapsed time
#   isoaccuracy_<fam>_<spec>.png           the same fits with the axes swapped:
#                                          what a fixed accuracy costs over time,
#                                          one contour per accuracy target
#
# Nine fit sets: three inefficiency families crossed with linear / quadratic /
# Box-Cox specifications. See fit_specs.R for the family grid and
# boxcox_frontier.R for the third specification. S (no inefficiency term) is
# the reference case: with no u its "frontier" IS the conditional mean, so
# roughly half the runs sit above it by construction. A and B are worth having
# only insofar as they beat that.
#
# This script absorbed plot_isoaccuracy.R (Aug 2026). The two views always used
# identical fits, and each script refitted them -- tolerable at seconds per SFA
# fit, not at the ~10 minutes the profiled Box-Cox specification added per
# script. Everything is fitted ONCE, via the shared store (fit_store.R), and
# the fitted objects are passed to the two figure builders; under run_all.R
# the same objects also serve the regression tables.
#
# On the iso view's reading: each contour holds ACCURACY fixed and lets cost
# vary, so these are isoquants, not isocosts (they were called isocost figures
# until Aug 2026; that name says the opposite of what is drawn). Colour carries
# accuracy there rather than date -- still a magnitude, so the same sequential
# ramp is right -- and the observed runs are plotted at their (date, cost)
# coloured by the accuracy they actually reached, so a run's shade can be read
# against the contour it lies on.
#
# marginaleffects is not used -- it dispatches on classes with predict() methods
# and cannot handle a raw maxLik fit -- but the prediction is one line.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("fit_store.R")   # brings the whole fitting stack with it

d <- load_runs()
tbar  <- bench_tbar(d)     # from t - tc, so it matches whatever the fit used
dates <- bench_dates(d)
benches <- bench_levels(d$benchmark)

LEVELS <- seq(0.10, 0.90, by = 0.20)   # iso-accuracy contour targets

## ---- every fit, once -- from the shared store -----------------------------------------
# fit_store.R computes each of these on first request; under run_all.R the
# table script then reuses the same objects instead of refitting.

specs <- store_specs()
bc_fits <- list(S = store_bc("S"), A = store_bc("A"), B = store_bc("B"))

## ---- semiannual frontier figures ------------------------------------------------------

NOTES_BASE <- paste("All models retained. Data prepared once in prepare_data.R,",
                    "which also drops the duplicate rows, so every specification",
                    "sees the same sample.")
# Describes the test rather than asserting its outcome. The previous wording
# ("the quadratic time term is not significant on any benchmark") both went stale
# -- the LR table it cited reports p < 1e-7 on aime and gpqa -- and named a
# single term the test never isolated: quad adds all three second-order terms at
# once, so the LR is a joint test of the block.
NOTES_QUAD <- paste("Quadratic adds all three second-order terms (cost^2, time^2,",
                    "cost x time) jointly; see this script's LR table and",
                    "monotonicity check for whether the block earns its keep.")
NOTES_BC <- paste("Box-Cox: phi(odds) -- profiled for S, the logit for the",
                  "SFA pair -- is linear in phi(cost), phi(years since",
                  "mid-2020) and their product, with the transform parameters",
                  "profiled per benchmark; monotone in cost at every date.")

## ---- the pooled pseudo-benchmark, S family's sixth panel ------------------------------
# The five primaries pooled on the anchored ECI scale with benchmark fixed
# effects (fit_pooled_acc, envelope_frontier.R), drawn as a sixth panel on
# the S figures only -- the SFA families are not pooled (see the pooled
# section of envelope_frontier.R for why).
pd <- pooled_acc_display(d, grid_dates(min(d$releasedate), max(d$releasedate)))
axis_ranges_all <- do.call(rbind, lapply(benches, function(b) {
  s <- d[d$benchmark == b, ]
  data.frame(benchmark = b, cost = range(s$cost), value = c(0, 1))
}))
axis_ranges_p <- rbind(axis_ranges_all, data.frame(
  benchmark = "pooled", cost = range(pd$sa$cost), value = range(pd$spx$acc)))
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
  "capability range on the shared ramp; dots keep each run's own-benchmark",
  "accuracy colour.")
pool_iso_layers <- function(pool_iso)
  pooled_iso_layers(pool_iso, crng = range(pd$spx$acc))

frontier_fig <- function(fits, fam, fname, extra_notes, pooled = NULL) {
  curves <- frontier_curves(fits, d, dates, tbar)
  if (!is.null(pooled))
    curves <- rbind(curves,
                    pooled_acc_frontier_curves(pooled, pd$sa, pd$dates))
  p <- frontier_plot(
    curves, if (is.null(pooled)) d else pts_frontier,
    labels = if (is.null(pooled)) LABELS else LABELS_POOLED,
    free_value = !is.null(pooled),
    ranges = if (is.null(pooled)) NULL else axis_ranges_p,
    ylab = if (fam == "S") "Fitted accuracy" else "Frontier accuracy",
    notes = c(NOTES_BASE, extra_notes, if (!is.null(pooled)) POOL_NOTE))
  ggsave(out_path(fname), p, width = 10, height = fig_height(length(benches)), dpi = 200,
         device = ragg::agg_png)
  cat("wrote", fname, "\n")
}

for (k in names(specs)) {
  sp <- specs[[k]]
  frontier_fig(sp$fits, sp$family, sprintf("frontier_progression_%s.png", k),
               if (sp$time == "quad") NOTES_QUAD,
               pooled = if (sp$family == "S")
                 store_pooled_acc("S")[[sp$time]])
}
for (fam in names(bc_fits))
  frontier_fig(bc_fits[[fam]], fam,
               sprintf("frontier_progression_%s_bc.png", fam), NOTES_BC,
               pooled = if (fam == "S") store_pooled_acc_bc("S"))

## ---- iso-accuracy figures -------------------------------------------------------------

# Pin the panels to the observed date and cost ranges, so blanked contour
# segments cannot pull the axes around.
iso_ranges <- do.call(rbind, lapply(benches, function(b) {
  s <- d[d$benchmark == b, ]
  data.frame(benchmark = b, date = range(s$releasedate), cost = range(s$cost))
}))

ISO_NOTES_BASE <- c(
  paste("Contours are accuracy targets from 10% to 90%; a falling contour means",
        "the same performance costs less over time. Dots are observed runs,",
        "coloured by the accuracy they reached."),
  paste("Contours are cut where they leave the benchmark's observed cost range,",
        "so a missing contour means that target lies outside the data rather",
        "than that it costs nothing."))
ISO_NOTES_QUAD <- paste("Quadratic adds all three second-order terms (cost^2,",
                        "time^2, cost x time) jointly, so curvature here is the",
                        "whole block's doing, not any one term's.")
# The BC contours have ONE branch (phi is monotone, so the inversion has one
# root), so ISO_BRANCH_NOTE does not apply and NOTES_BC takes its place.

iso_ranges_p <- rbind(iso_ranges, data.frame(
  benchmark = "pooled", date = range(pd$sa$releasedate),
  cost = range(pd$sa$cost)))

iso_fig <- function(fits, fname, extra_notes, pooled = NULL) {
  curves <- iso_acc_curves(fits, d, tbar, levels = LEVELS)
  p <- iso_acc_plot(curves, if (is.null(pooled)) d else pts_iso,
                    labels = if (is.null(pooled)) LABELS else LABELS_POOLED,
                    notes = c(ISO_NOTES_BASE, extra_notes,
                              if (!is.null(pooled)) POOL_ISO_NOTE),
                    ranges = if (is.null(pooled)) iso_ranges else iso_ranges_p)
  if (!is.null(pooled))
    p <- p + pool_iso_layers(pooled_acc_iso_curves(pooled, pd$sa, pd$levels))
  ggsave(out_path(fname), p, width = 10, height = fig_height(length(benches)), dpi = 200,
         device = ragg::agg_png)
  cat("wrote", fname, "\n")
}

for (k in names(specs)) {
  sp <- specs[[k]]
  iso_fig(sp$fits, sprintf("isoaccuracy_%s.png", k),
          if (sp$time == "quad") c(ISO_NOTES_QUAD, ISO_BRANCH_NOTE),
          pooled = if (sp$family == "S") store_pooled_acc("S")[[sp$time]])
}
for (fam in names(bc_fits))
  iso_fig(bc_fits[[fam]], sprintf("isoaccuracy_%s_bc.png", fam), NOTES_BC,
          pooled = if (fam == "S") store_pooled_acc_bc("S"))

## ---- does each fitted curve actually envelope the data? -------------------------------
# For S this is a sanity check, not a criticism: a conditional mean SHOULD have
# about half the runs above it. For A and B a high share means the inefficiency
# term is not buying a frontier.

cat("\nshare of observed runs lying ABOVE the fitted curve\n")
cat(sprintf("%-6s %s\n", "bench",
            paste(sprintf("%10s", names(specs)), collapse = "")))
for (b in benches) {
  s <- d[d$benchmark == b, ]
  sh <- vapply(specs, function(sp) {
    co <- frontier_coefs(sp$fits[[b]])
    mean(s$acc > plogis(frontier_index(co, s$lncost, s$tc)))
  }, numeric(1))
  cat(sprintf("%-6s %s\n", b, paste(sprintf("%9.1f%%", 100 * sh), collapse = "")))
}

report_quadratic_pathologies(specs, d)
report_scale_lr(specs, d)
report_quad_lr(specs)

## ---- how much of each contour is inside the data? ------------------------------------
# A contour that is mostly blank is the model extrapolating, not a cost decline.

# Rising branch only. iso_acc_curves() returns both branches, so averaging
# over all its rows would halve every figure here for reasons that have nothing
# to do with coverage.
cat("\nshare of each contour inside the observed cost range (%), model A linear\n")
cv <- iso_acc_curves(specs$A_lin$fits, d, tbar, levels = LEVELS)
cvr <- cv[cv$branch == "rising", ]
print(round(100 * tapply(!is.na(cvr$cost), list(cvr$benchmark, cvr$acc), mean)))

## ---- implied cost decline at 50% accuracy --------------------------------------------
# Cost at the first vs last observed date. Under the linear specification the
# ratio is the same at every accuracy level (the contours are parallel in logs);
# the quadratic breaks that, so the two columns can differ.

cat("\ncost of 50% accuracy, first vs last observed date\n")
cat(sprintf("%-6s %-6s %12s %12s %8s %10s\n",
            "spec", "bench", "at first", "at last", "ratio", "%/yr"))
for (k in names(specs)) {
  for (b in benches) {
    co <- frontier_coefs(specs[[k]]$fits[[b]])
    s  <- d[d$benchmark == b, ]
    ends <- range(s$releasedate)
    tc <- as_t(ends) - tbar[[b]]
    cst <- exp((qlogis(0.5) - co[["b0"]] - co[["bt"]] * tc - co[["btt"]] * tc^2) /
                 co[["bx"]])
    yrs <- diff(as_t(ends))
    cat(sprintf("%-6s %-6s %12.4f %12.4f %8.3f %9.1f%%\n", k, b, cst[1], cst[2],
                cst[2] / cst[1], 100 * (1 - (cst[2] / cst[1])^(1 / yrs))))
  }
}
