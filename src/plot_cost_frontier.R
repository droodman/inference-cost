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
#   costols         lm of ln cost on the accuracy terms and date over all
#                   positive-accuracy runs -- the reverse regression of model
#                   S: the TYPICAL cost of a run scoring a, not a frontier
#   costsfa         stochastic cost frontier, constant inefficiency scale
#                   (family A's dual)
#   costsfab        the same with log sigma_u linear in date (B's dual)
#   costgridols     least squares to the record cost ln C_a(t) on the fixed
#                   (logit accuracy, date) grid -- the Pareto-grid logit's dual
#   costgridolsenv  the grid OLS objective under the cost envelope's
#                   constraints -- the envelope-constrained logit's dual

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

## ---- the pooled pseudo-benchmark, sixth panel ----------------------------------------
#
# The five primaries pooled on the anchored ECI scale with benchmark fixed
# effects (pooled_cost_runs, cost_frontier.R), drawn as a sixth panel. The
# staircase machinery is generic in (cost, acc), so handing it la as `acc`
# yields ECI-unit staircases; the untouched-acc copy keeps the iso view's
# dots coloured by each run's own-benchmark accuracy on the shared 0-1 scale.
# The iso view's CONTOURS cannot ride that shared accuracy colour scale --
# their levels are ECI scores -- so the pooled panel draws them in a fixed
# neutral colour, labelled in place with their ECI values.
sp  <- pooled_cost_runs(d)
spx <- sp
spx$acc <- spx$la
LB <- append(LABELS, c(pooled = "Pooled primaries (ECI scale)"),
             after = match("mystery", names(LABELS)))
pool_dates <- {
  g <- grid_dates(min(d$releasedate), max(d$releasedate))
  g[g >= min(sp$releasedate)]
}
pool_levels <- {
  pr <- pretty(range(sp$la), 5)
  pr[pr > min(sp$la) & pr < max(sp$la)]
}
pool_steps     <- pareto_curves(spx, setNames(list(pool_dates), "pooled"))
pool_iso_steps <- iso_pareto_curves(spx, pool_levels)
PC <- c("benchmark", "cost", "acc", "year")           # frontier-view point columns
IC <- c("benchmark", "releasedate", "cost", "acc")    # iso-view point columns
pts_frontier <- rbind(d[, PC], spx[, PC])
pts_iso      <- rbind(d[, IC], sp[, IC])
axis_ranges_p <- rbind(axis_ranges, data.frame(
  benchmark = "pooled", cost = range(sp$cost), value = range(sp$la)))
iso_ranges_p <- rbind(iso_ranges, data.frame(
  benchmark = "pooled", date = range(sp$releasedate), cost = range(sp$cost)))
POOL_NOTE <- paste(
  "Sixth panel: the five primaries pooled on the anchored ECI capability",
  "scale (2PL: C = logit(a)/alpha_b + D_b, Claude 3.5 Sonnet = 130) with",
  "benchmark fixed effects; its value axis is in ECI points and its curves",
  "trace the fitted surface at its cheapest benchmark copy (the minimum",
  "fixed effect), the model's counterpart of the cost record.")
POOL_ISO_NOTE <- paste(
  "Sixth panel: the pooled primaries; contours are ECI capability levels",
  "from the pooled fit at its cheapest benchmark copy, labelled in ECI",
  "points and shaded by position within the pooled capability range on the",
  "shared ramp; dashes their record staircases; dots keep each run's",
  "own-benchmark accuracy colour.")

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
  costgridolsenv = list(
    note = paste("Least squares to the record cost ln C_a(t) on the grid,",
                 "constrained to lie at or below every run's log cost and",
                 "stay monotone -- the cost envelope's feasible set with the",
                 "grid OLS's loss; no standard errors.")))

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
  bc = paste("Box-Cox: phi(cost) -- profiled for the least-squares fits,",
             "fixed at ln cost for the SFA duals, whose inefficiency term the",
             "response transform would confound -- is linear in phi(odds) and",
             "phi(years since October 2020) (zeros excluded)."))
ISO_SPEC_NOTE <- c(
  lin = paste("The linear surface cannot bend, so its parallel contours miss",
              "the records' sharp rise toward each date's best-achieved",
              "level: functional-form misfit, priced in dollars here and so",
              "visible."),
  quad = paste("Full quadratic in logit accuracy and date. Curvature helps",
               "the fit, but a cliff whose location rides the advancing",
               "ceiling still cannot be tracked by one fixed surface."),
  bc = paste("Box-Cox: phi(cost) -- profiled for the least-squares fits,",
             "fixed at ln cost for the SFA duals -- linear in phi(odds) and",
             "phi(years since October 2020); a lambda on the search-box edge",
             "means the profile ran to the wall, and contours blank where",
             "the fitted index leaves phi's range."))

for (key in names(MODELS)) {
  m <- MODELS[[key]]
  for (tt in c("lin", "quad", "bc")) {
    fits <- if (tt == "bc") store_cost_bc(key) else store_cost(key)[[tt]]
    # the pooled companion fit for the sixth panel: lin and quad for every
    # key, Box-Cox for the least-squares keys only (fit_pooled_cost_bc); the
    # SFA duals' BC figures keep their five panels
    pf <- if (tt == "bc") {
      if (key %in% c("costols", "costgridols", "costgridolsenv"))
        store_pooled_cost_bc(key) else NULL
    } else store_pooled_cost(key)[[tt]]
    lbs <- if (is.null(pf)) LABELS else LB

    # frontier view: the surface swept along accuracy at each drawn date
    curves <- cost_frontier_curves(fits, d, dates, tbar)
    if (!is.null(pf))
      curves <- rbind(curves, pooled_frontier_curves(pf, sp, pool_dates))
    p <- frontier_plot(
      curves, if (is.null(pf)) d else pts_frontier,
      ranges = if (is.null(pf)) axis_ranges else axis_ranges_p,
      labels = lbs, free_value = !is.null(pf),
      ylab = "Fitted accuracy (cost fit inverted)",
      notes = c(
        paste("Solid: the fitted cost surface traced over the observed",
              "accuracy range at each drawn date -- the inverted view; the",
              "iso-accuracy figure is this model's native one."),
        PARETO_STEP_NOTE, SPEC_NOTE[[tt]], m$note,
        if (!is.null(pf)) POOL_NOTE)) +
      pareto_step_layer(if (is.null(pf)) steps else rbind(steps, pool_steps),
                        labels = lbs)
    f <- sprintf("%s_%s.png", key, tt)
    ggsave(out_path(f), p, width = 10, height = fig_height(length(benches)), dpi = 200,
           device = ragg::agg_png)
    cat("wrote", f, "\n")

    # iso-accuracy view: contours straight off the fitted surface
    iso <- cost_iso_curves(fits, d, tbar, levels = LEVELS,
                           cost_cap = iso_steps)
    p_iso <- iso_acc_plot(
      iso, if (is.null(pf)) d else pts_iso,
      ranges = if (is.null(pf)) iso_ranges else iso_ranges_p,
      labels = lbs,
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
        ISO_SPEC_NOTE[[tt]], m$note,
        if (!is.null(pf)) POOL_ISO_NOTE)) +
      iso_pareto_layer(iso_steps, labels = lbs)
    if (!is.null(pf)) {
      # contours at ECI capability levels, shaded by position within the
      # pooled capability range on the shared ramp (pooled_iso_layers,
      # frontier_viz.R), each labelled in place with its ECI value
      p_iso <- p_iso +
        pooled_iso_layers(pooled_iso_curves(pf, sp, pool_levels,
                                            cost_cap = pool_iso_steps),
                          crng = range(sp$la), steps = pool_iso_steps)
    }
    fi <- sprintf("isoaccuracy_%s_%s.png", key, tt)
    ggsave(out_path(fi), p_iso, width = 10, height = fig_height(length(benches)), dpi = 200,
           device = ragg::agg_png)
    cat("wrote", fi, "\n")

    if (tt == "lin") {
      cat(sprintf("  %-12s per-quarter cost change at fixed accuracy:", key))
      for (b in benches)
        cat(sprintf("  %s %+.1f%%", b, -cost_decline_qtr(fits[[b]])))
      if (!is.null(pf)) cat(sprintf("  pooled %+.1f%%", -cost_decline_qtr(pf)))
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
      if (!is.null(pf)) {
        si <- iso_runs(sp)
        co <- cost_coefs(pf)
        cg <- expand.grid(la = range(si$la), tc = range(si$tc))
        cat(sprintf("  pooled [%s]",
                    if (min(cost_dacc(co, cg$la, cg$tc)) >= -1e-6 &&
                        max(cost_dtime(co, cg$la, cg$tc)) <= 1e-6) "mono"
                    else "NON-MONO"))
      }
      cat("\n")
    }
    if (tt == "bc") {
      cat(sprintf("  %-12s BC lambdas (odds, time):", key))
      for (b in benches) {
        lam <- attr(fits[[b]], "bc_lambda")
        cat(sprintf("  %s (%.2f, %.2f)", b, lam[["lambda_odds"]],
                    lam[["lambda_time"]]))
      }
      if (!is.null(pf)) {
        lam <- attr(pf, "bc_lambda")
        cat(sprintf("  pooled (%.2f, %.2f)", lam[["lambda_odds"]],
                    lam[["lambda_time"]]))
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
  # Same guard as fit_specs.R's LR table and regression_tables.R's est_se: a
  # fit pinned at the sigma_u = 0 boundary has a flat direction, the Hessian
  # is singular, and the sandwich cannot be built. Say so and carry on --
  # sigma_u printed above is what diagnoses it.
  m <- tryCatch(summary_robust(f), error = function(e) NULL)
  if (is.null(m))
    cat("  robust summary unavailable: singular Hessian (flat direction)\n")
  else print(m)
}
