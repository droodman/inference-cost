# Visualising the deterministic envelope frontier against the data it must cover.
#
# Two specifications:
#   lin   z = b0 + b_x*ln c + b_t*tc
#   quad  z = b0 + b_x*ln c + b_xx*(ln c)^2 + b_t*tc      (bends in cost)
#
# Each figure draws, at a handful of annual dates so the panels stay legible:
#   solid  the fitted envelope
#   dashed the empirical Pareto staircase at the same date
# so the parametric surface can be read directly against the nonparametric
# object it is trying to approximate. Where they diverge is where the functional
# form is doing the work rather than the data.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("envelope_frontier.R")

d <- load_runs()
benches <- sort(unique(d$benchmark))

# The same two forms every other model uses (TIME_FORMS in fit_specs.R): `quad`
# is the full second-order surface, not the cost-only quadratic it used to be
# here. The envelope now bends in time and lets its cost slope drift with date,
# which is what makes it comparable with the SFA and logit fits at last.
SPECS <- list(
  lin  = list(form = acc ~ lncost + tc),
  quad = list(form = acc ~ lncost * tc + I(lncost^2) + I(tc^2))
)

# The same quarterly grid as every other figure, via bench_dates(), so the
# envelope can be compared with the SFA and Pareto panels curve for curve.
dates <- bench_dates(d)

# The empirical staircase to draw underneath the fitted envelope. Same dates as
# the fitted curves, and the same pareto_curves() the standalone Pareto figure
# and the frontier-run logit use -- see frontier_viz.R.
steps <- pareto_curves(d, dates)

# Its iso-accuracy counterpart, at the same levels as the fitted contours:
# the minimum cost of achieving each accuracy level by each date.
LEVELS <- seq(0.05, 0.95, by = 0.10)
iso_steps <- iso_pareto_curves(d, LEVELS)

axis_ranges <- do.call(rbind, lapply(benches, function(b) {
  s <- d[d$benchmark == b, ]
  data.frame(benchmark = b, cost = range(s$cost), value = c(0, 1))
}))

for (k in names(SPECS)) {
  sp <- SPECS[[k]]
  curves <- list(); fits <- list()
  for (b in benches) {
    s <- d[d$benchmark == b, ]
    f <- fit_envelope(s, formula = sp$form)
    fits[[b]] <- f
    co <- frontier_coefs_envelope(f)
    dts <- dates[[b]]
    cost <- exp(seq(min(s$lncost), max(s$lncost), length.out = 200))
    for (q in dts) {
      tc <- as_t(q) - unique(s$t - s$tc)[1]
      # frontier_index(), not a hand-written sum: this line used to spell out
      # b0 + bx*lnc + bt*tc + bxx*lnc^2 and would have silently dropped the tc^2
      # and lncost:tc terms the full quadratic added, plotting a surface that was
      # not the one fitted.
      z <- frontier_index(co, log(cost), tc)
      curves[[length(curves) + 1]] <- data.frame(
        cost = cost, value = plogis(z), qdate = q, benchmark = b,
        year = 2023 + as_t(q))
    }
  }
  curves <- do.call(rbind, curves)

  p <- frontier_plot(
    curves, d, ranges = axis_ranges,
    ylab = "Frontier accuracy",
    notes = c(
      # Pre-empts the obvious misreading: the curves float above the staircases
      # almost everywhere, and that is the fit working, not failing.
      paste("Solid: the lowest logistic surface above every run. ONE surface",
            "spans all dates, pinned at only a few runs, so a curve need not",
            "meet its staircase."),
      PARETO_STEP_NOTE,
      paste("Monotonicity in cost and time is imposed, so the surface cannot",
            "slope backwards in either."),
      paste("Tightness is scored on a fixed grid, so it does not shift with where",
            "runs happen to cluster; the data enters only as constraints."))) +
    pareto_step_layer(steps)

  f <- sprintf("envelope_%s.png", k)
  ggsave(out_path(f), p, width = 10, height = 7.5, dpi = 200,
         device = ragg::agg_png)
  cat("wrote", f, "\n")

  ## ---- iso-accuracy view of the same fits ----------------------------------
  # Axes swapped: what a fixed accuracy target costs over time. Because the
  # envelope is monotone by construction, these contours cannot double back
  # the way the SFA ones did.
  iso <- iso_acc_curves(fits, d, bench_tbar(d), levels = LEVELS,
                        cost_cap = iso_steps)
  iso_ranges <- do.call(rbind, lapply(benches, function(b) {
    s <- d[d$benchmark == b, ]
    data.frame(benchmark = b, date = range(s$releasedate), cost = range(s$cost))
  }))
  pi <- iso_acc_plot(
    iso, d, ranges = iso_ranges,
    notes = c(
      paste("Contours are accuracy targets from 5% to 95% read off the",
            "deterministic envelope; a falling contour means the same",
            "performance costs less over time."),
      ISO_PARETO_NOTE,
      paste("Cut at the observed cost range and above each level's dashed",
            "record, so a level never achieved by any run shows no contour at",
            "all."),
      paste("The envelope must clear every run, so these are upper bounds on",
            "what frontier performance costs, not central estimates."),
      ISO_BRANCH_NOTE)) +
    iso_pareto_layer(iso_steps)
  fi <- sprintf("isoaccuracy_envelope_%s.png", k)
  ggsave(out_path(fi), pi, width = 10, height = 7.5, dpi = 200,
         device = ragg::agg_png)
  cat("wrote", fi, "\n")

  # Monotonicity is IMPOSED here, unlike the SFA and logit fits, so this is a
  # check that the constraints did their job rather than a diagnostic of the
  # data. It matters because the constraint rows are built from the formula: a
  # term whose derivative the constraint code failed to account for would leave
  # the surface free to slope backwards in exactly the way the envelope promises
  # it cannot.
  # WHERE the surface touches the data, and whether that is visible on the
  # figure. It generally is not: contact happens at each binding run's OWN
  # release date, the curves are drawn on a quarterly grid, and the surface
  # rises with date -- so a curve drawn after a run it touches sits above it by
  # the time-growth in between. On fm13 that is the whole apparent gap between
  # the solid curves and the dashed staircases; the fit is tight, the slices just
  # are not taken where it is tight. Reported so "the envelope looks loose" can
  # be checked rather than eyeballed.
  cat(sprintf("\n== %s: where the envelope touches the data ==\n", k))
  cat(sprintf("%-6s %8s %10s %9s %-12s %s\n", "bench", "n touch", "worst slack",
              "at acc", "run date", "a drawn quarter?"))
  for (b in benches) {
    f <- fits[[b]]
    s <- d[d$benchmark == b, ]
    for (i in f$bind[f$env_slack < 0.01]) {
      cat(sprintf("%-6s %8d %10.2e %9.3f %-12s %s\n", b,
                  sum(f$env_slack < 0.01), f$slack_envelope, s$acc[i],
                  format(s$releasedate[i], "%Y-%m-%d"),
                  if (s$releasedate[i] %in% dates[[b]]) "yes" else "no"))
    }
  }

  cat(sprintf("\n== %s: monotonicity actually achieved (imposed, so must hold) ==\n", k))
  cat(sprintf("%-6s %12s %12s %10s  %s\n", "bench", "min dz/dlnc", "min dz/dtc",
              "slack", "status"))
  for (b in benches) {
    co <- frontier_coefs_envelope(fits[[b]])
    s  <- d[d$benchmark == b, ]
    bd <- frontier_slope_bounds(co, s$lncost, s$tc)
    cat(sprintf("%-6s %12.4f %12.4f %10.2e  %s\n", b, bd[["dcost"]], bd[["dtime"]],
                fits[[b]]$worst_slack,
                if (all(bd >= -1e-6)) "ok" else "CONSTRAINT FAILED"))
  }

  # Rising branch only -- iso_acc_curves() returns both, and averaging over both
  # would report a coverage drop that is just the second branch's NAs.
  cat(sprintf("\n== %s: share of each contour inside the observed cost range ==\n", k))
  isor <- iso[iso$branch == "rising", ]
  print(round(100 * tapply(!is.na(isor$cost), list(isor$benchmark, isor$acc), mean)))

  ## ---- in-range reporting only ------------------------------------------------
  # The cost implied for an accuracy target is meaningful only when it lands
  # inside the observed cost range at BOTH endpoints. Outside, the surface is
  # extrapolating and the "decline" is an artefact of the functional form.
  cat(sprintf("\n== %s: cost of a target accuracy, first vs last date ==\n", k))
  cat(sprintf("%-6s %6s %12s %12s %8s %9s  %s\n", "bench", "target",
              "at first", "at last", "ratio", "%/yr", "status"))
  for (b in benches) {
    co <- frontier_coefs_envelope(fits[[b]])
    s  <- d[d$benchmark == b, ]
    rng <- range(s$cost); tcr <- range(s$tc); yrs <- diff(tcr)
    for (tgt in c(0.25, 0.50, 0.75)) {
      # Solve z(u, tc) = logit(tgt) for u = ln cost, with the SAME coefficients
      # iso_acc_curves() uses -- previously this spelled out only b0, b_t, b_x
      # and b_xx, so the full quadratic's tc^2 and lncost:tc terms would have
      # been dropped here while the plotted contours kept them, and the table
      # would have disagreed with the figure beside it.
      #   aa*u^2 + bb*u + cc = 0,  aa = b_xx,  bb = b_x + b_xt*tc,
      #   cc = b0 + b_t*tc + b_tt*tc^2 - logit(tgt)
      cst <- vapply(tcr, function(tc) {
        aa <- co[["bxx"]]
        bb <- co[["bx"]] + co[["bxt"]] * tc
        cc <- co[["b0"]] + co[["bt"]] * tc + co[["btt"]] * tc^2 - qlogis(tgt)
        if (abs(aa) < 1e-10) return(exp(-cc / bb))
        disc <- bb^2 - 4 * aa * cc
        if (disc < 0) return(NA_real_)
        r <- (-bb + c(-1, 1) * sqrt(disc)) / (2 * aa)
        r <- exp(r[is.finite(r)])
        if (!length(r)) NA_real_ else min(r)
      }, numeric(1))
      ok <- all(is.finite(cst)) && all(cst >= rng[1] & cst <= rng[2])
      cat(sprintf("%-6s %6.0f%% %12.5f %12.5f %8.3f %8.1f%%  %s\n",
                  b, 100 * tgt, cst[1], cst[2], cst[2] / cst[1],
                  100 * (1 - (cst[2] / cst[1])^(1 / yrs)),
                  if (isTRUE(ok)) "in range" else "OUT OF RANGE -- ignore"))
    }
  }
}
