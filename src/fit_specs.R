# The parametric specification grid, used by plot_frontiers.R (which draws both
# figure families from one set of fits) and regression_tables.R.
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

# `quad` is the FULL second-order surface in (log cost, time): both squares and
# the interaction, six coefficients including the intercept. Nothing else.
#
# It used to mean something different here (time squared plus the interaction)
# than it did for the envelope (cost squared alone), which made "the quadratic
# fits" a comparison between models that were quadratic in different arguments.
# One complete second-order surface everywhere is the neutral starting point:
# whatever curvature the data wants is available to every model, and the
# specification stops being a place where the models differ.
#
# The flexibility is not free, and the pathologies are the point of the exercise:
#   * the cost slope b_x + 2*b_xx*ln c + b_xt*tc now moves with BOTH coordinates,
#     so free disposal can fail in a corner of the (cost, date) rectangle rather
#     than after a single crossing date;
#   * (ln c)^2 is near-collinear with ln c, which has already stalled the
#     envelope's solver once (see envelope_frontier.R);
#   * three extra parameters on benchmarks as small as fm13 (266 runs).
# report_quadratic_pathologies() checks the first, and the LR table the last.
TIME_FORMS <- list(lin  = acc ~ lncost + tc,
                   quad = acc ~ lncost * tc + I(lncost^2) + I(tc^2))

TIME_LABEL <- c(lin  = "linear in log cost and time",
                quad = "full quadratic in log cost and time")

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
  benches <- bench_levels(d$benchmark)
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
      # A thin benchmark can leave the fit with a flat direction -- family B's
      # time-varying scale is the usual culprit -- and the Hessian is then
      # singular, so the sandwich cannot be built. The LR and the curvature
      # check need no inverse, so report them with EMPTY se and correlation
      # cells rather than halting the whole run: this table is a diagnostic,
      # and the rows it CAN compute are exactly the ones that diagnose the
      # problem. regression_tables.R's est_se() takes the same line.
      V <- tryCatch(vcov_robust(f)[ap, ap, drop = FALSE],
                    error = function(e) NULL)
      i1 <- match("beta_tc", nm); i2 <- match("logsig_tc", nm)
      ev <- min(eigen(-hessian(f)[ap, ap, drop = FALSE], symmetric = TRUE,
                      only.values = TRUE)$values)
      lr <- 2 * (as.numeric(logLik(f)) - as.numeric(logLik(fa[[b]])))
      cat(sprintf("%-6s %-6s %8.2f %7.4f %11.4f %9.4f %9.1e %+11.3f%s\n",
                  tt, b, lr, pchisq(max(lr, 0), 1, lower.tail = FALSE),
                  cf[["logsig_tc"]],
                  if (is.null(V)) NA_real_ else sqrt(V[i2, i2]), ev,
                  if (is.null(V)) NA_real_ else
                    V[i1, i2] / sqrt(V[i1, i1] * V[i2, i2]),
                  paste0(if (ev > 0) "" else "  SADDLE",
                         if (is.null(V)) "  SINGULAR" else "")))
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
# With the full quadratic BOTH derivatives vary over BOTH coordinates, so the
# old "where does the cost slope cross zero" framing no longer describes the
# failure: there is no single crossing date, there is a region of the (cost,
# date) rectangle where the surface slopes the wrong way. Each derivative is
# linear in (ln c, tc) jointly, so its minimum over the rectangle is attained at
# a corner -- that minimum is what is reported, alongside the share of ACTUAL
# runs sitting where the sign is wrong, which is what decides whether the
# pathology contaminates the fitted region or lives in an empty corner.
# These fits describe the WHOLE distribution of runs, not just its upper edge,
# and the two derivatives do not deserve the same prior away from that edge:
#
#   d z / d ln c  free disposal. A strong prior anywhere: spending more should
#                 not buy less. Flagged when it fails.
#   d z / d tc    NOT expected to be positive in the interior. As effort moves to
#                 expensive reasoning models, cheap runs can genuinely get worse
#                 over time, so a negative time slope at low spend is a pattern
#                 in the data rather than a defect in the fit. Reported, not
#                 judged.
#
# So the status column speaks only to free disposal, and the time column is left
# as a number to read. Only the ENVELOPE imposes both, because only it claims to
# be a frontier; see plot_paretologitenv.R, whose binding report shows where
# its constraints hold and bind.
report_quadratic_pathologies <- function(specs, d) {
  cat("\nslopes of the fitted surface over the observed (cost, date) rectangle\n")
  cat("  dz/dlnc < 0 breaks free disposal; dz/dtc < 0 at low spend is expected\n")
  cat(sprintf("%-8s %-6s %11s %10s %11s %10s  %s\n", "spec", "bench",
              "min dz/dlnc", "runs dc<0", "min dz/dtc", "runs dt<0",
              "free disposal"))
  for (k in names(specs)) {
    for (b in names(specs[[k]]$fits)) {
      co <- frontier_coefs(specs[[k]]$fits[[b]])
      s  <- d[d$benchmark == b, ]
      bounds <- frontier_slope_bounds(co, s$lncost, s$tc)
      # at the runs themselves, not just the corners: a violation in a corner
      # where nothing was ever run is extrapolation, not a finding
      shr_c <- mean(frontier_dcost(co, s$lncost, s$tc) < 0)
      shr_t <- mean(frontier_dtime(co, s$lncost, s$tc) < 0)
      cat(sprintf("%-8s %-6s %11.4f %9.1f%% %11.4f %9.1f%%  %s\n",
                  k, b, bounds[["dcost"]], 100 * shr_c,
                  bounds[["dtime"]], 100 * shr_t,
                  if (bounds[["dcost"]] >= 0) "ok"
                  else if (shr_c < 0.01) "violated, but <1% of runs"
                  else "VIOLATED at observed runs"))
    }
  }
}

report_quad_lr <- function(specs) {
  cat("\nLR test, quadratic vs linear time (SFA families only)\n")
  cat(sprintf("%-6s %-8s %10s %10s %8s %5s %8s\n",
              "family", "bench", "logLik lin", "logLik quad", "LR", "df", "p"))
  for (fam in c("A", "B")) {
    lin <- specs[[paste0(fam, "_lin")]]$fits
    qd  <- specs[[paste0(fam, "_quad")]]$fits
    for (b in names(lin)) {
      l0 <- as.numeric(logLik(lin[[b]]))
      l1 <- as.numeric(logLik(qd[[b]]))
      lr <- 2 * (l1 - l0)
      # df counted from the fits, not hardcoded: `quad` adds three beta terms
      # now (both squares and the interaction) where it added two before, and
      # this test was left testing on 1 df through that change.
      df <- sum(activePar(qd[[b]])) - sum(activePar(lin[[b]]))
      cat(sprintf("%-6s %-8s %10.3f %10.3f %8.2f %5d %8.4f\n", fam, b, l0, l1,
                  lr, df, pchisq(max(lr, 0), df, lower.tail = FALSE)))
    }
  }
}
