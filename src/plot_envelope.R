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

SPECS <- list(
  lin  = list(form = acc ~ lncost + tc,
              sub  = "envelope linear in log cost"),
  quad = list(form = acc ~ lncost + I(lncost^2) + tc,
              sub  = "envelope quadratic in log cost")
)

# The same quarterly grid as every other figure, via bench_dates(), so the
# envelope can be compared with the SFA and Pareto panels curve for curve.
dates <- bench_dates(d)

pareto_steps <- function(sub, q) {
  s <- sub[sub$releasedate <= q, ]
  if (!nrow(s)) return(NULL)
  s <- s[order(s$cost), ]
  m <- cummax(s$acc)
  keep <- c(TRUE, diff(m) > 0)
  data.frame(cost = c(s$cost[keep], max(sub$cost)),
             value = c(m[keep], m[length(m)]))
}

axis_ranges <- do.call(rbind, lapply(benches, function(b) {
  s <- d[d$benchmark == b, ]
  data.frame(benchmark = b, cost = range(s$cost), value = c(0, 1))
}))

for (k in names(SPECS)) {
  sp <- SPECS[[k]]
  curves <- list(); steps <- list(); fits <- list()
  for (b in benches) {
    s <- d[d$benchmark == b, ]
    f <- fit_envelope(s, formula = sp$form)
    fits[[b]] <- f
    co <- frontier_coefs_envelope(f)
    dts <- dates[[b]]
    cost <- exp(seq(min(s$lncost), max(s$lncost), length.out = 200))
    for (q in dts) {
      tc <- as_t(q) - unique(s$t - s$tc)[1]
      z <- co[["b0"]] + co[["bx"]] * log(cost) + co[["bt"]] * tc +
        co[["bxx"]] * log(cost)^2
      curves[[length(curves) + 1]] <- data.frame(
        cost = cost, value = plogis(z), qdate = q, benchmark = b,
        year = 2023 + as_t(q))
      st <- pareto_steps(s, q)
      if (!is.null(st)) {
        st$qdate <- q; st$benchmark <- b; st$year <- 2023 + as_t(q)
        steps[[length(steps) + 1]] <- st
      }
    }
  }
  curves <- do.call(rbind, curves)
  steps  <- do.call(rbind, steps)
  steps$benchmark <- factor(LABELS[steps$benchmark], levels = LABELS)

  p <- frontier_plot(
    curves, d, ranges = axis_ranges,
    title = sprintf("Deterministic envelope frontier -- %s", k),
    subtitle = sp$sub,
    ylab = "Frontier accuracy",
    notes = c(
      paste("Solid: fitted envelope, the lowest logistic surface passing above",
            "every run. Dashed: the empirical Pareto staircase at the same date."),
      paste("Monotonicity in cost and time is imposed, so the surface cannot",
            "slope backwards in either."),
      paste("Tightness is scored on a fixed grid, so it does not shift with where",
            "runs happen to cluster; the data enters only as constraints."))) +
    geom_step(data = steps, aes(cost, value, group = qdate, colour = year),
              direction = "hv", linewidth = 0.4, linetype = "22",
              inherit.aes = FALSE)

  f <- sprintf("envelope_%s.png", k)
  ggsave(out_path(f), p, width = 10, height = 7.5, dpi = 200,
         device = ragg::agg_png)
  cat("wrote", f, "\n")

  ## ---- iso-cost view of the same fits -------------------------------------
  # Axes swapped: what a fixed accuracy target costs over time. Because the
  # envelope is monotone by construction, these contours cannot double back
  # the way the SFA ones did.
  iso <- iso_cost_curves(fits, d, bench_tbar(d), levels = seq(0.05, 0.95, 0.10))
  iso_ranges <- do.call(rbind, lapply(benches, function(b) {
    s <- d[d$benchmark == b, ]
    data.frame(benchmark = b, date = range(s$releasedate), cost = range(s$cost))
  }))
  pi <- iso_cost_plot(
    iso, d, ranges = iso_ranges,
    title = sprintf("Cost of a fixed accuracy level -- envelope frontier (%s)", k),
    subtitle = sp$sub,
    notes = c(
      paste("Contours are accuracy targets from 5% to 95% read off the",
            "deterministic envelope; a falling contour means the same",
            "performance costs less over time."),
      paste("Cut where they leave the observed cost range, so a missing contour",
            "means that target lies outside the data rather than that it is",
            "free."),
      paste("The envelope must clear every run, so these are upper bounds on",
            "what frontier performance costs, not central estimates.")))
  fi <- sprintf("isocost_envelope_%s.png", k)
  ggsave(out_path(fi), pi, width = 10, height = 7.5, dpi = 200,
         device = ragg::agg_png)
  cat("wrote", fi, "\n")

  cat(sprintf("\n== %s: share of each contour inside the observed cost range ==\n", k))
  print(round(100 * tapply(!is.na(iso$cost), list(iso$benchmark, iso$acc), mean)))

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
      # solve b_xx*u^2 + b_x*u + (b0 + b_t*tc - logit(tgt)) = 0 for u = ln cost
      cst <- vapply(tcr, function(tc) {
        cc <- co[["b0"]] + co[["bt"]] * tc - qlogis(tgt)
        if (abs(co[["bxx"]]) < 1e-10) return(exp(-cc / co[["bx"]]))
        disc <- co[["bx"]]^2 - 4 * co[["bxx"]] * cc
        if (disc < 0) return(NA_real_)
        r <- (-co[["bx"]] + c(-1, 1) * sqrt(disc)) / (2 * co[["bxx"]])
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
