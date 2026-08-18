# Logistic regression fitted ONLY to the frontier-defining runs.
#
# The Pareto frontier at date t is
#
#     P_t(c) = max { a_i : c_i <= c, t_i <= t }
#
# and the runs that DEFINE it are its jump points: the runs no other run
# dominates, where j dominates i when j is no dearer, no later, and scored no
# worse. That is the same set envelope_frontier.R isolates, for the same reason
# -- with the surface non-decreasing in cost and date, a dominated run's
# constraint is implied by its dominator's -- so pareto_binding() is reused here
# rather than reimplemented, and the two models are guaranteed to be looking at
# the same points.
#
# What differs is what is asked of those points. The envelope wants the lowest
# S-curve passing ABOVE all of them. This wants the S-curve that runs THROUGH
# them: a conventional quasibinomial logit, nothing imposed. Consequences worth
# keeping in view --
#
#   * the fit need not lie above the points, so it is not a frontier in the
#     envelope's sense; it is a summary of how the frontier-defining runs trade
#     accuracy against cost and date;
#   * free disposal is NOT imposed, so the cost slope can come out negative and
#     the iso-accuracy contours blank out where it does;
#   * the estimator is exactly model S's (fit_family("S", ...) is called
#     literally), so S vs this model isolates ONE thing: which runs are fitted.
#     S uses all of them, this uses the frontier.
#
# Standard errors are quasi-likelihood and, more importantly, condition on a
# SELECTED sample -- the points were chosen for being maxima -- so they are
# descriptive here, not inferential.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("fit_specs.R")          # TIME_FORMS, TIME_LABEL, fit_family, viz
src_source("envelope_frontier.R")  # pareto_binding()

LEVELS <- seq(0.05, 0.95, by = 0.10)

d <- load_runs(drop_gpt4o_chess = FALSE)
tbar    <- bench_tbar(d)
dates   <- bench_dates(d)
benches <- sort(unique(d$benchmark))

## ---- the frontier-defining runs -------------------------------------------------

# Selected on RAW accuracy over ALL runs, which is P_t(c) as defined above.
# fit_envelope() feeds pareto_binding() something slightly different -- clipped
# logit, with zero-scoring runs dropped first, because logit(0) = -Inf cannot be
# a constraint. Those two differences change the set only at the margins (a
# cheapest-and-earliest zero run, and ties among perfect scores that clipping
# breaks by sample size); the count of each is reported below rather than
# assumed away.
frontier_rows <- function(s) s[pareto_binding(s$lncost, s$tc, s$acc), , drop = FALSE]

front <- do.call(rbind, lapply(benches, function(b) frontier_rows(d[d$benchmark == b, ])))

cat("frontier-defining runs, by benchmark\n")
cat(sprintf("%-6s %8s %8s %8s\n", "bench", "all", "frontier", "share"))
for (b in benches) {
  n_all <- sum(d$benchmark == b)
  n_fr  <- sum(front$benchmark == b)
  cat(sprintf("%-6s %8d %8d %7.1f%%\n", b, n_all, n_fr, 100 * n_fr / n_all))
}

# How far does this set sit from the one the envelope actually constrains on?
cat("\noverlap with the envelope's binding-candidate set\n")
cat(sprintf("%-6s %9s %9s %9s\n", "bench", "P_t set", "envelope", "common"))
for (b in benches) {
  s <- d[d$benchmark == b, ]
  a <- rownames(frontier_rows(s))
  L <- qlogis(clip_acc(s$acc, s$n_samples))
  pos <- which(is.finite(L))
  e <- rownames(s)[pos[pareto_binding(s$lncost[pos], s$tc[pos], L[pos])]]
  cat(sprintf("%-6s %9d %9d %9d\n", b, length(a), length(e),
              length(intersect(a, e))))
}

# A fit needs more points than parameters; say so loudly rather than reading a
# rank-deficient glm's coefficients as a frontier.
for (b in benches) {
  n_fr <- sum(front$benchmark == b)
  if (n_fr < 8)
    warning(sprintf("%s has only %d frontier-defining runs; the quadratic fit there is not credible",
                    b, n_fr))
}

## ---- figures ----------------------------------------------------------------------

FAMILY_SUB <- "logit fitted to frontier-defining runs only"

NOTES_FRONTIER <- c(
  paste("Solid: a conventional logistic regression fitted ONLY to the runs that",
        "define the staircase -- its jump points, not the full sample."),
  PARETO_STEP_NOTE,
  paste("The same runs that can bind the envelope: with the surface non-decreasing",
        "in cost and date, a dominated run's constraint is implied by its dominator's."),
  paste("The fit runs THROUGH those points rather than above them, and",
        "monotonicity in cost or time is not imposed."))

NOTES_ISO <- c(
  paste("Contours are accuracy targets from 5% to 95%; a falling contour means the",
        "same performance costs less over time. Dots are observed runs."),
  paste("Read off a logit fitted only to the frontier-defining runs -- the jump",
        "points of P_t(c) -- not to the full sample."),
  paste("Cut where they leave the benchmark's observed cost range, so a missing",
        "contour means that target lies outside the data, not that it is free."),
  paste("Free disposal is not imposed, so a contour also blanks out at dates where",
        "the fitted cost slope is too flat to invert."),
  ISO_BRANCH_NOTE)

iso_ranges <- do.call(rbind, lapply(benches, function(b) {
  s <- d[d$benchmark == b, ]
  data.frame(benchmark = b, date = range(s$releasedate), cost = range(s$cost))
}))

# Pinned the same way plot_envelope.R pins its panels, so the two figures can be
# laid side by side and read panel for panel.
axis_ranges <- do.call(rbind, lapply(benches, function(b) {
  s <- d[d$benchmark == b, ]
  data.frame(benchmark = b, cost = range(s$cost), value = c(0, 1))
}))

# The empirical staircase, drawn underneath the fit exactly as plot_envelope.R
# draws it. Showing the frontier itself rather than highlighting the fitted runs
# keeps the two figures readable side by side, and loses nothing: the staircase's
# jump points ARE the fitted runs.
steps <- pareto_curves(d, dates)

for (tt in names(TIME_FORMS)) {
  fits <- fit_family("S", TIME_FORMS[[tt]], front)   # same estimator as model S

  curves <- frontier_curves(fits, d, dates, tbar)
  p <- frontier_plot(
    curves, d, ranges = axis_ranges,
    title = sprintf("Accuracy frontier fitted to frontier-defining runs (%s)", tt),
    subtitle = sprintf("%s; %s", FAMILY_SUB, TIME_LABEL[[tt]]),
    ylab = "Fitted frontier accuracy",
    notes = NOTES_FRONTIER) +
    pareto_step_layer(steps)

  f <- sprintf("paretologit_%s.png", tt)
  ggsave(out_path(f), p, width = 10, height = 7.5, dpi = 200,
         device = ragg::agg_png)
  cat("wrote", f, "\n")

  iso <- iso_acc_curves(fits, d, tbar, levels = LEVELS)
  pi <- iso_acc_plot(
    iso, d, ranges = iso_ranges,
    title = sprintf("Cost of a fixed accuracy level -- frontier-run logit (%s)", tt),
    subtitle = sprintf("%s; %s", FAMILY_SUB, TIME_LABEL[[tt]]),
    notes = NOTES_ISO)

  fi <- sprintf("isoaccuracy_paretologit_%s.png", tt)
  ggsave(out_path(fi), pi, width = 10, height = 7.5, dpi = 200,
         device = ragg::agg_png)
  cat("wrote", fi, "\n")
}

## ---- what does restricting the sample actually change? ------------------------------
# The same estimator on the same functional form, all runs vs frontier runs only.
# beta_x is the cost slope in logits, beta_t the improvement rate per year at the
# benchmark's reference date.

cat("\ncoefficients: all runs (S) vs frontier-defining runs only, linear in time\n")
cat(sprintf("%-6s %9s %9s %9s %9s\n", "bench", "b_x all", "b_x front",
            "b_t all", "b_t front"))
all_fits <- fit_family("S", TIME_FORMS$lin, d)
frt_fits <- fit_family("S", TIME_FORMS$lin, front)
for (b in benches) {
  ca <- frontier_coefs(all_fits[[b]])
  cf <- frontier_coefs(frt_fits[[b]])
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
  ff <- fit_family("S", TIME_FORMS[[tt]], front)
  for (b in benches) {
    co <- frontier_coefs(ff[[b]])
    s  <- d[d$benchmark == b, ]
    bd <- frontier_slope_bounds(co, s$lncost, s$tc)
    cat(sprintf("%-6s %-5s %12.4f %12.4f  %s\n", b, tt,
                bd[["dcost"]], bd[["dtime"]],
                if (all(bd >= 0)) "ok" else "NON-MONOTONE"))
  }
}
