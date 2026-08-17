# The parametric specification grid, shared by plot_frontiers.R and
# plot_isocost.R so both figure families are built from identical fits.
#
# Two dimensions:
#   inefficiency  A = SFA with mu = 0        (half-normal)
#                 B = SFA with mu = delta_0  (truncated normal)
#                 S = none                   (plain fractional logit)
#   time          lin  = tc
#                 quad = tc + tc^2
#
# tc is time demeaned within benchmark, so beta_tc is the improvement rate at
# each benchmark's reference date and tc is decorrelated from its square.
# frontier_coefs() returns 0 for an absent I(tc^2), so linear and quadratic fits
# flow through the same prediction code untouched.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("panel_frontier.R")
src_source("frontier_viz.R")

# `quad` is the flexible frontier: curvature in time AND a cost-slope that may
# drift with time. The interaction is what lets the cost-accuracy relationship
# steepen as inference scaling improves -- but it also lets the fitted slope go
# NEGATIVE at some dates, which violates free disposal. report_slope_monotonicity()
# checks whether that happens inside the observed span or only outside it.
TIME_FORMS <- list(lin  = acc ~ lncost + tc,
                   quad = acc ~ lncost * tc + I(tc^2))

TIME_LABEL <- c(lin  = "linear in time",
                quad = "quadratic in time, cost slope varying with time")

# Both SFA families now fix mu = 0 (half-normal); they differ in whether the
# SCALE of inefficiency may move with time:
#
#   A  log sigma_u = g0            constant spread
#   B  log sigma_u = g0 + g1*tc    spread shrinking (or growing) over time
#
# B replaces the earlier truncated-normal variant. That one let delta_0 set the
# typical gap independently of the spread, but delta_0 traded off against beta_0
# -- the likelihood sees only beta_0 - E[u] -- so the frontier LEVEL was not
# identified for gpqa (12 starts spanned delta_0 from -2.1 to +5.8 within 0.12
# log-likelihood units). Moving the time story into the scale is better posed:
# with mu = 0, E[u] = 0.8*sigma_u, so g1 shrinks the mean and the spread
# together, and the spread change is a signature beta_t cannot imitate.
FAMILY_LABEL <- c(
  A = "SFA, mu = 0, constant scale (half-normal)",
  B = "SFA, mu = 0, log sigma_u linear in time",
  S = "plain fractional logit -- no inefficiency term")

SIGMA_FORM <- list(A = ~ 1, B = ~ tc)

U_GROUP <- c("model", "effort")

fit_family <- function(family, form, d) {
  benches <- sort(unique(d$benchmark))
  if (family == "S") {
    setNames(lapply(benches, function(b) {
      glm(form, data = d[d$benchmark == b, ], family = quasibinomial(link = "logit"))
    }), benches)
  } else {
    fit_panel_frontier_by(
      form, ~ 1, data = d, u_group = U_GROUP, by = "benchmark", dedup = FALSE,
      formula_sigma = SIGMA_FORM[[family]],
      fixed = "delta_(Intercept)")     # mu = 0 in both
  }
}

# Returns a named list; each element carries the fits plus the labels a figure
# needs, so callers loop once and never reconstruct the naming themselves.
fit_all_specs <- function(d) {
  out <- list()
  for (tt in names(TIME_FORMS)) {
    for (fam in names(FAMILY_LABEL)) {
      key <- paste0(fam, "_", tt)
      out[[key]] <- list(
        family = fam, time = tt, form = TIME_FORMS[[tt]],
        subtitle = sprintf("%s; %s", FAMILY_LABEL[[fam]], TIME_LABEL[[tt]]),
        fits = fit_family(fam, TIME_FORMS[[tt]], d))
    }
  }
  out
}

# LR test of the quadratic against the linear fit, per benchmark. Available for
# the SFA families only: S is a quasi-likelihood, so glm reports no logLik and
# the comparison has no likelihood to test.
# B nests A (gamma_1 = 0), so the LR tests whether the inefficiency SCALE moves
# with time. Reported with the minimum eigenvalue of -H and the correlation
# between beta_tc and gamma_1, the two checks that condemned the earlier
# truncated-normal variant: delta_t correlated +0.955 with beta_t and several
# benchmarks landed on saddles.
report_scale_lr <- function(specs, d) {
  cat("\nLR test, time-varying vs constant inefficiency scale (B vs A)\n")
  cat(sprintf("%-6s %-6s %8s %7s %11s %9s %9s %11s\n", "time", "bench",
              "LR", "p", "logsig_tc", "se", "min_eig", "corr(bt,g1)"))
  for (tt in names(TIME_FORMS)) {
    fa <- specs[[paste0("A_", tt)]]$fits
    fb <- specs[[paste0("B_", tt)]]$fits
    for (b in names(fa)) {
      f <- fb[[b]]
      ap <- activePar(f); cf <- coef(f); nm <- names(cf)[ap]
      V <- vcov_robust(f)[ap, ap, drop = FALSE]
      i1 <- match("beta_tc", nm); i2 <- match("logsig_tc", nm)
      ev <- min(eigen(-hessian(f)[ap, ap, drop = FALSE], symmetric = TRUE,
                      only.values = TRUE)$values)
      lr <- 2 * (as.numeric(logLik(f)) - as.numeric(logLik(fa[[b]])))
      cat(sprintf("%-6s %-6s %8.2f %7.4f %11.4f %9.4f %9.1e %+11.3f%s\n",
                  tt, b, lr, pchisq(max(lr, 0), 1, lower.tail = FALSE),
                  cf[["logsig_tc"]], sqrt(V[i2, i2]), ev,
                  V[i1, i2] / sqrt(V[i1, i1] * V[i2, i2]),
                  if (ev > 0) "" else "  SADDLE"))
    }
  }
  cat("\nimplied sigma_u at the first vs last observed date\n")
  cat(sprintf("%-6s %-6s %11s %11s %8s %9s\n", "time", "bench",
              "sigma first", "sigma last", "ratio", "%/yr"))
  for (tt in names(TIME_FORMS)) {
    for (b in names(specs[[paste0("B_", tt)]]$fits)) {
      cf <- coef(specs[[paste0("B_", tt)]]$fits[[b]])
      tcr <- range(d$tc[d$benchmark == b])
      sg <- exp(cf[["logsig_(Intercept)"]] + cf[["logsig_tc"]] * tcr)
      cat(sprintf("%-6s %-6s %11.4f %11.4f %8.3f %8.1f%%\n", tt, b, sg[1], sg[2],
                  sg[2] / sg[1], 100 * (1 - exp(cf[["logsig_tc"]]))))
    }
  }
}

# Free disposal: the frontier's slope in log cost, b_x + b_xt*tc, must stay
# positive or the fit says more spend buys less accuracy. The interaction makes
# that slope drift with date, so the question is not merely whether it ever turns
# negative but WHERE: a crossing outside the observed span is harmless
# extrapolation, one inside it contaminates the fitted region. Reported as the
# crossing date, whether it falls in range, and the share of actual runs sitting
# at a date where the slope is negative.
report_slope_monotonicity <- function(specs, d) {
  cat("\nfree disposal: where does the fitted cost slope cross zero?\n")
  cat(sprintf("%-8s %-6s %9s %9s %11s %-12s %10s %9s\n", "spec", "bench",
              "b_lncost", "b_x:tc", "slope@first", "crossing", "in range?",
              "runs < 0"))
  for (k in names(specs)) {
    if (specs[[k]]$time != "quad") next          # lin has no interaction
    for (b in names(specs[[k]]$fits)) {
      co <- frontier_coefs(specs[[k]]$fits[[b]])
      s  <- d[d$benchmark == b, ]
      rng <- range(s$tc)
      sl  <- cost_slope(co, rng)
      cross <- if (co[["bxt"]] == 0) NA_real_ else -co[["bx"]] / co[["bxt"]]
      inrng <- !is.na(cross) && cross >= rng[1] && cross <= rng[2]
      share <- mean(cost_slope(co, s$tc) <= 0)
      cdate <- if (is.na(cross)) "--" else
        format(EPOCH + (cross + mean(s$t - s$tc)) * 365.25, "%Y-%m")
      cat(sprintf("%-8s %-6s %9.4f %9.4f %11.4f %-12s %10s %8.1f%%\n",
                  k, b, co[["bx"]], co[["bxt"]], sl[1], cdate,
                  if (inrng) "INSIDE" else "outside", 100 * share))
    }
  }
}

report_quad_lr <- function(specs) {
  cat("\nLR test, quadratic vs linear time (SFA families only)\n")
  cat(sprintf("%-6s %-8s %10s %10s %8s %8s\n",
              "family", "bench", "logLik lin", "logLik quad", "LR", "p"))
  for (fam in c("A", "B")) {
    lin <- specs[[paste0(fam, "_lin")]]$fits
    qd  <- specs[[paste0(fam, "_quad")]]$fits
    for (b in names(lin)) {
      l0 <- as.numeric(logLik(lin[[b]]))
      l1 <- as.numeric(logLik(qd[[b]]))
      lr <- 2 * (l1 - l0)
      cat(sprintf("%-6s %-8s %10.3f %10.3f %8.2f %8.4f\n", fam, b, l0, l1, lr,
                  pchisq(max(lr, 0), 1, lower.tail = FALSE)))
    }
  }
}
