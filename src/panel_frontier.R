# Panel version of the fractional-logit stochastic frontier.
#
# Difference from fractional_frontier.R: the inefficiency term u is a GROUP-level
# effect (here: model x effort, within a benchmark), constant across all budget
# levels of that group's cost curve, rather than an independent draw per row.
#
#   ln y*_ig = x_ig'beta - u_g,    u_g ~ N+(mu_g, sigma_u^2),   mu_g = w_g'delta
#
# Motivation: with u drawn per observation, sigma_u is not identified -- the
# fractional-logit kernel is a broad bump in the index (curvature bounded by 1/4,
# width ~3), so a sigma_u of 0.005-0.5 changes it only at second order while
# E[u] is absorbed by the intercept. The profile likelihood is monotone in
# sigma_u and the MLE sits on the sigma_u = 0 boundary (see fractional_frontier.R
# results). Sharing u across a group's ~11-12 curve points means u is inferred
# against the PRODUCT of those kernels, whose curvature is their sum, so the
# identification comes from real replication rather than an assumed noise scale.
#
# Side benefit: the likelihood units become groups, so the strong within-curve
# dependence sits inside a unit instead of violating independence across units,
# and the BHHH/sandwich scores are automatically clustered at the group level.
#
# Caveats worth keeping in mind when reading the output:
#   - u_g is now ENTIRELY inefficiency: there is no group-level noise term, so
#     this is a deterministic frontier at the group level.
#   - Benchmarks with few groups (fm13 has 8 models) cannot support a
#     cross-group variance and will stay degenerate.
#   - 20-46% of acc values are exactly 0, many at budgets too small to answer at
#     all; u will partly absorb that floor. Filter `data` before calling if you
#     want a minimum-budget restriction -- none is imposed here.
#
# Why formula_mu has no time term. A rising frontier (beta_t) and a shrinking
# inefficiency spread (delta_t) both make later models score better, and they are
# separated only by the shape of the truncated normal -- a thin channel that
# supports ONE mu parameter, not two:
#
#   mu = d0 + dt*t : chess/fm13/gpqa land on a SADDLE (min eigenvalue of -H is
#     -3.3e-5, -1.3e-6, -3.9e-3, confirmed by Richardson extrapolation, so not a
#     finite-difference artifact) despite reporting code 0. aime is a proper
#     maximum but multimodal -- modes at beta_t = -1.94, +1.25, +3.76 spanning
#     3.1 log-likelihood units, all inside a 95% LR region.
#   mu = dt*t (d0 fixed at 0) : proper maximum on all four, so delta_t IS
#     estimable alone -- and it is insignificant everywhere (LR vs mu = 0:
#     p = 0.099, 0.91, 0.75, 1.00 for aime/gpqa/fm13/chess).
#   mu = d0 : the better use of the single degree of freedom -- LR p = 0.0001
#     for aime, and 7.7 log-likelihood units better than mu = dt*t on equal df.
#
# So delta_t is dropped on evidence, not just on numerical grounds. Use mu = d0
# where it earns its keep (aime), otherwise mu = 0; beta_t is stable across both.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("fractional_frontier.R")   # reuse log1pexp(), gauss_legendre_01()

## ---- likelihood ---------------------------------------------------------------

make_panel_loglik <- function(y, Xb, Xmu_g, Xs_g, grp, quad) {
  K <- length(quad$nodes)
  G <- nrow(Xmu_g)
  logw <- log(quad$weights)
  log1m_p <- log1p(-quad$nodes)
  skeleton <- list(beta = numeric(ncol(Xb)), delta = numeric(ncol(Xmu_g)),
                   logsig = numeric(ncol(Xs_g)))

  function(theta) {
    p <- relist(theta, skeleton)
    beta    <- p$beta
    delta   <- p$delta
    # sigma_u is a G-vector, not a scalar: log sigma_u = s_g'gamma, the same
    # linear-predictor construction already used for mu. With Xs_g = ~1 this is
    # the constant-scale model exactly as before.
    sigma_u <- exp(as.vector(Xs_g %*% p$logsig))   # G

    eta <- as.vector(Xb %*% beta)        # n:  x_ig'beta
    mu  <- as.vector(Xmu_g %*% delta)    # G:  one truncation centre per group

    # identical PIT quadrature to the cross-sectional version, but at group level
    log_1mPhi_c <- pnorm(mu, sd = sigma_u, log.p = TRUE)
    log_1mx <- outer(log_1mPhi_c, log1m_p, "+")        # G x K
    Q <- -qnorm(log_1mx, log.p = TRUE)                 # G x K
    U <- mu + sigma_u * Q                              # G x K

    Z <- eta - U[grp, , drop = FALSE]                  # n x K
    H <- y * Z - log1pexp(Z)                           # n x K
    S <- rowsum(H, grp)                                # G x K: sum_{i in g}
    M <- S + matrix(logw, G, K, byrow = TRUE)
    ll <- matrixStats::rowLogSumExps(M)                # G: per-group log L_g

    # Score. dH/dZ = y - G(Z); dZ/ddelta and dZ/dlog_sigma_u depend only on
    # (g,k), so they multiply the group-summed residual E_gk, while dZ/dbeta = x_i
    # needs the per-observation weighting c_i = sum_k q_gk * e_ik.
    q  <- exp(M - ll)                                  # G x K
    e  <- y - plogis(Z)                                # n x K
    E  <- rowsum(e, grp)                               # G x K
    qE <- q * E
    ci <- rowSums(q[grp, , drop = FALSE] * e)          # n
    logR <- outer(-0.5 * (mu / sigma_u)^2, log1m_p, "+") + 0.5 * Q^2
    Rm1 <- expm1(logR)                                 # G x K: R - 1

    # dl_g/d(log sigma_u,g) as before; the chain rule to gamma is just a
    # multiplication by s_g, exactly as delta's is by w_g.
    dtau <- -rowSums(qE * (U + mu * Rm1))              # G
    attr(ll, "gradient") <- cbind(
      rowsum(Xb * ci, grp),                            # dl/dbeta
      Xmu_g * rowSums(qE * Rm1),                       # dl/ddelta
      Xs_g * dtau                                      # dl/dgamma
    )
    ll
  }
}

## ---- fitting wrapper ------------------------------------------------------------

fit_panel_frontier <- function(formula_beta, formula_mu, data, u_group,
                               formula_sigma = ~ 1,
                               n_quad = 64, dedup = TRUE, start = NULL,
                               method = "BFGS", ...) {
  vars <- unique(c(all.vars(formula_beta), all.vars(formula_mu),
                   all.vars(formula_sigma), u_group))
  data <- data[complete.cases(data[vars]), , drop = FALSE]

  key <- interaction(data[u_group], drop = TRUE, lex.order = TRUE)
  grp <- as.integer(key)

  mf <- model.frame(formula_beta, data)
  y  <- model.response(mf)
  if (any(y < 0 | y > 1)) stop("response must lie in [0, 1]")
  Xb  <- model.matrix(formula_beta, data)
  Xmu <- model.matrix(formula_mu, data)
  Xs  <- model.matrix(formula_sigma, data)

  # Drop rows that repeat an earlier row of the same group exactly in both design
  # and outcome. These arise where the effort cap binds and extra budget buys
  # nothing, so the row simply repeats; keeping them would count one observation
  # several times and spuriously sharpen the very kernel u is identified against.
  if (dedup) {
    dup <- duplicated(data.frame(grp, Xb, y))
    if (any(dup)) {
      message(sprintf("dropping %d duplicate row(s) of %d (%.1f%%)",
                      sum(dup), length(dup), 100 * mean(dup)))
      keep <- !dup
      data <- data[keep, , drop = FALSE]; y <- y[keep]
      Xb <- Xb[keep, , drop = FALSE]; Xmu <- Xmu[keep, , drop = FALSE]
      Xs <- Xs[keep, , drop = FALSE]
      key <- droplevels(key[keep]); grp <- as.integer(key)
    }
  }

  # mu and log sigma_u both vary at the group level, so their regressors must be
  # group-constant. Compare values only: model.matrix() carries an "assign"
  # attribute that matrix subsetting drops, so comparing the matrices themselves
  # fails on attributes.
  first <- match(seq_len(nlevels(key)), grp)
  Xmu_g <- Xmu[first, , drop = FALSE]
  Xs_g  <- Xs[first, , drop = FALSE]
  if (!isTRUE(all.equal(as.vector(Xmu), as.vector(Xmu_g[grp, , drop = FALSE])))) {
    stop("formula_mu regressors vary within a u_group; they must be constant")
  }
  if (!isTRUE(all.equal(as.vector(Xs), as.vector(Xs_g[grp, , drop = FALSE])))) {
    stop("formula_sigma regressors vary within a u_group; they must be constant")
  }

  quad <- gauss_legendre_01(n_quad)
  ll <- make_panel_loglik(y, Xb, Xmu_g, Xs_g, grp, quad)

  if (is.null(start)) {
    b0 <- coef(glm(formula_beta, data = data, family = quasibinomial(link = "logit")))
    d0 <- numeric(ncol(Xmu_g))
    g0 <- c(log(0.3), numeric(ncol(Xs_g) - 1))   # scale at the reference date
    names(b0) <- paste0("beta_", colnames(Xb))
    names(d0) <- paste0("delta_", colnames(Xmu_g))
    names(g0) <- paste0("logsig_", colnames(Xs_g))
    start <- c(b0, d0, g0)
  }

  res <- maxLik(ll, start = start, method = method, ...)
  attr(res, "loglik_fn") <- ll
  attr(res, "n_groups")  <- nlevels(key)
  attr(res, "n_obs")     <- length(y)
  res
}

# One fit per level of `by`, independent by construction, so they run on the
# shared worker cluster when one is available (fit_cluster() in paths.R -- R
# and its bundled BLAS are single-threaded, so run serially the benchmarks fit
# one after another). The dots are captured with list() and re-spliced by
# do.call because a closure shipped to a PSOCK worker cannot carry `...`
# through serialisation.
fit_panel_frontier_by <- function(formula_beta, formula_mu, data, u_group, by, ...) {
  # Force every argument before the closure ships: an unforced promise crossing
  # the serialisation boundary re-evaluates its expression on the worker, where
  # the caller's variables do not exist. (data and by are forced by the next
  # two lines; list(...) forces the dots.)
  force(formula_beta)
  force(formula_mu)
  force(u_group)
  split_var <- data[[by]]
  levels <- sort(unique(split_var[!is.na(split_var)]))
  dots <- list(...)
  one <- function(lev) do.call(fit_panel_frontier, c(
    list(formula_beta, formula_mu, data[split_var == lev, , drop = FALSE],
         u_group), dots))
  cl <- fit_cluster(length(levels))
  fits <- if (is.null(cl)) lapply(levels, one) else
    parallel::parLapply(cl, levels, one)
  setNames(fits, as.character(levels))
}

## ---- cost of fixed frontier performance over time --------------------------------
#
# Holding target accuracy y fixed and inverting the frontier for cost:
#   ln c(t) = [logit(y) - b0 - b_t*t] / (b_x + b_xt*t)
# With no lncost:t interaction b_xt = 0, the denominator is constant and
# d ln(c)/dt = -b_t/b_x -- one rate at every accuracy level and every date.
# WITH the interaction that fails: the rate depends on both the date and the
# accuracy target, and the denominator crosses zero at t = -b_x/b_xt (Sep 2023
# for aime, Dec 2023 for gpqa), where the implied cost diverges and before which
# it carries the wrong sign. So at_date is required once an interaction is
# present, and dates near that singularity are refused rather than reported.
# Delta method (numeric Jacobian) off the robust vcov.

cost_trend <- function(fit, accuracy = 0.5, at_date = NULL, level = 0.95) {
  ap <- activePar(fit)
  cf <- coef(fit)[ap]
  V  <- vcov_robust(fit)[ap, ap, drop = FALSE]
  has_ix <- "beta_lncost:t" %in% names(cf)
  if (has_ix && is.null(at_date))
    stop("this fit has a lncost:t interaction, so d ln(c)/dt varies with date and ",
         "accuracy target: supply at_date (and accuracy)")

  tt <- if (is.null(at_date)) 0 else
    as.numeric(as.Date(at_date) - as.Date("2023-01-01")) / 365.25
  L <- qlogis(accuracy)

  lnc <- function(th, t.) {
    bx <- th[["beta_lncost"]] + if (has_ix) th[["beta_lncost:t"]] * t. else 0
    (L - th[["beta_(Intercept)"]] - th[["beta_t"]] * t.) / bx
  }
  denom <- cf[["beta_lncost"]] + if (has_ix) cf[["beta_lncost:t"]] * tt else 0
  if (abs(denom) < 0.05)
    stop(sprintf("cost slope is %.4f here -- too close to the singularity for the ",
                 denom), "cost inversion to be meaningful")

  J  <- numDeriv::grad(function(x) lnc(setNames(x, names(cf)), tt), cf)
  g  <- numDeriv::grad(function(x) lnc(cf, x), tt)      # d ln(c)/dt at this date
  se <- sqrt(drop(t(J) %*% V %*% J))                    # SE of ln c at this date

  c(lncost = lnc(cf, tt), cost = exp(lnc(cf, tt)), se_lncost = se,
    dlncost_dt = g, pct_per_yr = 100 * (1 - exp(g)),
    halving_months = 12 * log(2) / -g)
}

## ---- example usage ----------------------------------------------------------------

if (sys.nframe() == 0) {
  src_source("frontier_viz.R")   # load_runs() -- same sample as the figures

  # Specification matches plot_frontiers.R exactly: quadratic in time, demeaned
  # within benchmark, all models retained. Keep the two in step -- reporting
  # coefficients here that differ from the plotted curves is how they drifted
  # apart before. Deduplication happens once, in prepare_data.R.
  d <- load_runs()
  FORM <- acc ~ lncost + tc + I(tc^2)

  # delta_t is constrained out of both (formula_mu = ~1); they differ only in
  # whether delta_0 is free. S is the no-inefficiency reference.
  fitA <- fit_panel_frontier_by(FORM, ~ 1, data = d,
    u_group = c("model", "effort"),   # u is per model-effort cell, per benchmark
    by = "benchmark", dedup = FALSE,
    fixed = "delta_(Intercept)")      # half-normal: mu = 0
  fitB <- fit_panel_frontier_by(FORM, ~ 1, data = d,
    u_group = c("model", "effort"),
    by = "benchmark", dedup = FALSE)  # truncated normal: mu = delta_0, free

  for (bench in names(fitA)) {
    cat(sprintf("\n======== benchmark: %s  (%d groups, %d obs) ========\n",
                bench, attr(fitA[[bench]], "n_groups"), attr(fitA[[bench]], "n_obs")))
    sub <- d[d$benchmark == bench, ]
    fitS <- glm(FORM, data = sub, family = quasibinomial(link = "logit"))

    for (nm in c("A: mu = 0", "B: mu = delta_0")) {
      f <- if (startsWith(nm, "A")) fitA[[bench]] else fitB[[bench]]
      cat(sprintf("\n-- %-16s logLik = %10.3f   sigma_u = %.4f\n",
                  nm, as.numeric(logLik(f)), sigma_u_hat(f)))
      print(summary_robust(f))
    }
    cat("\n-- S: plain fractional logit (no inefficiency term)\n")
    print(lmtest::coeftest(fitS, vcov = sandwich::sandwich))

    lr <- 2 * (as.numeric(logLik(fitB[[bench]])) - as.numeric(logLik(fitA[[bench]])))
    cat(sprintf("\n   LR (B vs A, delta_0 = 0): %.2f   p = %.4f\n",
                lr, pchisq(lr, 1, lower.tail = FALSE)))
  }
}
