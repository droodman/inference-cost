# Figures for the cost-direction fits (cost_frontier.R): the mirror of the
# accuracy-direction figure scripts, with the two views' roles swapped.
#
# For the accuracy-direction models the quarterly-frontier quartets are the
# fits' NATIVE direction and the iso-accuracy figures are inversions -- which
# is why iso_acc_curves() carries the two-root fold bookkeeping. Here it is
# the other way round: the iso-accuracy view reads contours straight off the
# fitted cost surface (no inversion at all, whatever the specification), and
# the frontier view is produced by SWEEPING the surface along the observed
# accuracy range -- exact for every specification, with points blanked where
# the surface bends back in accuracy rather than drawn as a fold.
#
# Five models, the duals of the five in the plot viewer, each drawn in three
# specifications (linear, full quadratic, Box-Cox -- COST_FORMS and
# fit_cost_bc), on the same panels, dates, staircases and axis pins as their
# accuracy-direction counterparts so the viewer can flip between directions
# panel for panel:
#
#   costols       lm of ln cost on the accuracy terms and date over all
#                 positive-accuracy runs -- the reverse regression of model
#                 S: the TYPICAL cost of a run scoring a, not a frontier
#   costsfa       stochastic cost frontier, constant inefficiency scale
#                 (family A's dual)
#   costsfab      the same with log sigma_u linear in date (B's dual)
#   costgridols   least squares to the record cost ln C_a(t) on the fixed
#                 (logit accuracy, date) grid -- the Pareto-grid logit's dual
#   costenvelope  the highest surface under every run's log cost -- the
#                 strict envelope's dual

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("fit_store.R")   # fits come from the shared store

LEVELS <- seq(0.10, 0.90, by = 0.20)

d <- load_runs()
benches <- bench_levels(d$benchmark)
tbar    <- bench_tbar(d)
dates   <- bench_dates(d)

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

# One caption line saying what kind of surface each model's is; the rest of
# each caption is shared. The fitting recipes live in fit_cost_model
# (cost_frontier.R) and the fits come from the shared store, so the tables use
# the same objects. The SFA line states the cloud-versus-record caveat
# documented at length in cost_frontier.R -- on this data its time slope
# tracks the dense cheap edge of the model-effort cells, which grows dearer
# as reasoning configurations arrive, so read it against the dashed record.
MODELS <- list(
  costols = list(
    note = paste("Least squares of log cost over all positive-accuracy runs:",
                 "the TYPICAL cost of a run scoring a -- model S's reverse",
                 "regression, not a frontier.")),
  costsfa = list(
    note = paste("Stochastic cost frontier (half-normal inefficiency per",
                 "model x effort). Its time slope follows the dense cheap",
                 "edge of the cells, which grows DEARER as reasoning models",
                 "arrive -- not the record; compare the dashed staircase.")),
  costsfab = list(
    note = paste("Stochastic cost frontier with log sigma_u linear in date.",
                 "Like the constant-scale variant, its time slope follows",
                 "the dense cheap edge of the model-effort cells, not the",
                 "record; compare the dashed staircase.")),
  costgridols = list(
    note = paste("Least squares to the record cost ln C_a(t) sampled on a",
                 "uniform (logit accuracy, date) grid. Grid nodes are not",
                 "observations, so no standard errors -- point estimates",
                 "only, as for the Pareto-grid logit.")),
  costenvelope = list(
    note = paste("The highest surface in (logit accuracy, date) lying at or",
                 "below every run's log cost, monotone by constraint --",
                 "pinned by a few extreme runs, like its accuracy-direction",
                 "dual, and reporting no standard errors.")))

# Specification lines, one per caption. The frontier-view line also covers
# the zero-exclusion and (for the curved specs) the rising-branch blanking;
# the iso-view line states each specification's known limitation.
SPEC_NOTE <- c(
  lin = paste("ln cost is modeled linearly in logit accuracy and date, on",
              "runs scoring above zero (logit 0 is unusable as a",
              "coordinate)."),
  quad = paste("ln cost is a full quadratic in logit accuracy and date",
               "(zeros excluded); where the surface bends back in accuracy",
               "the curve stops rather than fold."),
  bc = paste("Box-Cox both sides: phi(cost) is linear in phi(odds), phi(years",
             "since mid-2020) and their product (zeros excluded), all three",
             "transform parameters profiled on lambda-invariant objectives."))
ISO_SPEC_NOTE <- c(
  lin = paste("The linear surface cannot bend, so its parallel contours miss",
              "the records' sharp rise toward each date's best-achieved",
              "level: functional-form misfit, priced in dollars here and so",
              "visible."),
  quad = paste("Full quadratic in logit accuracy and date. Curvature helps",
               "the fit, but a cliff whose location rides the advancing",
               "ceiling still cannot be tracked by one fixed surface."),
  bc = paste("Box-Cox both sides: phi(cost) linear in phi(odds), phi(years",
             "since mid-2020) and their product, all three lambdas profiled;",
             "a lambda on the search-box edge means the profile ran to the",
             "wall, and contours blank where the fitted index leaves phi's",
             "range."))

for (key in names(MODELS)) {
  m <- MODELS[[key]]
  for (tt in c("lin", "quad", "bc")) {
    fits <- if (tt == "bc") store_cost_bc(key) else store_cost(key)[[tt]]

    # frontier view: the surface swept along accuracy at each drawn date
    curves <- cost_frontier_curves(fits, d, dates, tbar)
    p <- frontier_plot(
      curves, d, ranges = axis_ranges,
      ylab = "Fitted accuracy (cost fit inverted)",
      notes = c(
        paste("Solid: the fitted cost surface traced over the observed",
              "accuracy range at each drawn date -- the inverted view; the",
              "iso-accuracy figure is this model's native one."),
        PARETO_STEP_NOTE, SPEC_NOTE[[tt]], m$note)) +
      pareto_step_layer(steps)
    f <- sprintf("%s_%s.png", key, tt)
    ggsave(out_path(f), p, width = 10, height = fig_height(length(benches)), dpi = 200,
           device = ragg::agg_png)
    cat("wrote", f, "\n")

    # iso-accuracy view: contours straight off the fitted surface
    iso <- cost_iso_curves(fits, d, tbar, levels = LEVELS,
                           cost_cap = iso_steps)
    p_iso <- iso_acc_plot(
      iso, d, ranges = iso_ranges,
      notes = c(
        paste("Contours are accuracy targets from 10% to 90% read DIRECTLY",
              "off the fitted cost surface -- the native view of a model of",
              "ln cost; no inversion is involved."),
        ISO_PARETO_NOTE,
        paste("Cut at the observed cost range; blanked only where both",
              "earlier than a level's first record and dearer than its",
              "dearest, so each contour reaches its record's start or cost",
              "ceiling, whichever is more generous. Never-achieved levels",
              "show no contour."),
        ISO_SPEC_NOTE[[tt]], m$note)) +
      iso_pareto_layer(iso_steps)
    fi <- sprintf("isoaccuracy_%s_%s.png", key, tt)
    ggsave(out_path(fi), p_iso, width = 10, height = fig_height(length(benches)), dpi = 200,
           device = ragg::agg_png)
    cat("wrote", fi, "\n")

    if (tt == "lin") {
      cat(sprintf("  %-12s per-quarter cost change at fixed accuracy:", key))
      for (b in benches)
        cat(sprintf("  %s %+.1f%%", b, -cost_decline_qtr(fits[[b]])))
      cat("\n")
    }
    if (tt == "quad") {
      # slopes over the observed rectangle, at its corners: dacc < 0 says
      # more accuracy is fitted as cheaper, dtime > 0 that the record rises.
      # The envelope imposes both signs; the others only get checked.
      cat(sprintf("  %-12s quad slopes at rectangle corners:", key))
      for (b in benches) {
        s <- iso_runs(d[d$benchmark == b, ])
        co <- cost_coefs(fits[[b]])
        cg <- expand.grid(la = range(s$la), tc = range(s$tc))
        cat(sprintf("  %s [%s]", b,
                    if (min(cost_dacc(co, cg$la, cg$tc)) >= -1e-6 &&
                        max(cost_dtime(co, cg$la, cg$tc)) <= 1e-6) "mono"
                    else "NON-MONO"))
      }
      cat("\n")
    }
    if (tt == "bc") {
      cat(sprintf("  %-12s BC lambdas (cost, odds, time):", key))
      for (b in benches) {
        lam <- attr(fits[[b]], "bc_lambda")
        cat(sprintf("  %s (%.2f, %.2f, %.2f)", b, lam[["lambda_cost"]],
                    lam[["lambda_odds"]], lam[["lambda_time"]]))
      }
      cat("\n")
    }
  }
}

## ---- the SFA duals' own diagnostics -------------------------------------------------

cat("\nstochastic cost frontier (constant scale, linear): robust summaries\n")
for (b in benches) {
  f <- store_cost("costsfa")$lin[[b]]
  cat(sprintf("\n== %s: code %d, logLik %.2f, %d groups, %d obs, sigma_u %.2f, sigma_v %.2f\n",
              b, f$code, as.numeric(logLik(f)), attr(f, "n_groups"),
              attr(f, "n_obs"), sigma_u_hat(f), exp(coef(f)[["logsig_v"]])))
  print(summary_robust(f))
}
