# Cost-direction duals of the frontier-per-se models: the same frontier
# object, fitted with the penalty in the cost direction.
#
# Every "cost drop, %/qtr" row targets d ln c / dt at fixed accuracy -- a
# horizontal quantity -- but the models in envelope_frontier.R price misfit in
# the ACCURACY direction and back the rate out as -b_t/b_x. Where the surface
# is not logit-linear that back-out is a ratio of grid-average slopes standing
# in for an average of local ratios, and on this data it overstates the
# decline 1.5-3x: the fitted b_t matches the staircase's actual rise almost
# exactly, while the fitted b_x is diluted by the staircase's wide flat steps
# and understates the slope where records actually move (compare the tables
# against pareto_decline_qtr's output in pareto_frontiers.R). The mirror fits
# ln C_a(t) -- the record cost of performance a at date t, the same object
# iso_pareto_steps() traces -- directly. The decline per quarter is then a
# single coefficient rather than a ratio, and misfit is priced in log
# dollars, the units the estimand lives in.
#
# For the FRONTIER these are duals, not different models: with the linear
# family, logit acc = b0 + b_x ln c + b_t tc inverts exactly to
# ln c = g0 + g_a logit(acc) + g_t tc with g_t = -b_t/b_x. What the flip
# changes is only the loss direction and the grid weighting -- which is the
# point.
#
# Two fits, mirroring envelope_frontier.R's pair:
#   fit_lncost_grid    ln C_a(t) sampled on a fixed (logit acc, date) grid,
#                      OLS through the samples -- the dual of the Pareto-grid
#                      logit, and the regression version of the nonparametric
#                      check pareto_decline_qtr()
#   fit_cost_envelope  the HIGHEST surface in (logit acc, date) lying at or
#                      below every run's log cost -- the dual of the strict
#                      envelope: the same Pareto corners bind, and only the
#                      objective differs (mean fitted log cost over the grid,
#                      maximized, where fit_envelope() minimizes mean fitted
#                      accuracy over its grid)
#
# Three specifications throughout, mirroring the accuracy direction: linear
# and full quadratic in (logit accuracy, time) (COST_FORMS), and the Box-Cox
# alternative with profiled transform parameters (fit_cost_bc, end of file).
#
# What the flip costs, stated rather than hidden:
#   * runs scoring 0 leave the sample: logit(0) = -Inf is unusable as a
#     coordinate, and "the cost of achieving at least nothing" is just the
#     cheapest run, which says nothing about the frontier (the same reasoning
#     that drops a = 0 nodes from pareto_decline_qtr);
#   * runs scoring n/n are clipped by their own sample size via clip_acc(),
#     exactly as the accuracy envelope clips its constraints;
#   * accuracy becomes a REGRESSOR, so its binomial noise is measurement
#     error rather than response noise. For these two fits that takes the
#     familiar form of the record's one-lucky-run fragility, already
#     documented for the staircase itself.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("envelope_frontier.R")   # clip_acc, pareto_binding, the accuracy-direction fits
src_source("fractional_frontier.R") # maxLik, vcov_robust/summary_robust, for the SFA dual
src_source("boxcox_frontier.R")     # BC_BOX_*, bc_lt_free, for the BC cost specs

# The cost-direction specification pair, mirroring TIME_FORMS (fit_specs.R):
# `quad` is the FULL second-order surface in (logit accuracy, time) -- both
# squares and the interaction -- for the same reason the accuracy direction
# standardized on it: whatever curvature the data wants is available to every
# model, and the specification stops being a place where models differ. The
# same warnings transfer in mirror image: the accuracy slope g_a + 2*g_aa*la
# + g_at*tc now moves with both coordinates, so "costlier to do better" can
# fail in a corner of the rectangle (checked, or imposed by the envelope);
# and no single cost-drop rate describes a quadratic column.
COST_FORMS <- list(lin  = lncost ~ la + tc,
                   quad = lncost ~ la * tc + I(la^2) + I(tc^2))

# Runs usable in the cost direction: logit of clipped accuracy as a
# coordinate, zeros dropped. Every function below starts from this.
#
#   data  DATA FRAME of one benchmark's runs; needs acc, n_samples, lncost, tc
#
# Returns the rows with finite logit accuracy, sorted by lncost, with the
# column `la` added.
iso_runs <- function(data) {
  s <- data
  s$la <- qlogis(clip_acc(s$acc, s$n_samples))
  s <- s[is.finite(s$la), , drop = FALSE]
  s[order(s$lncost), ]
}

# The fixed grid the cost-direction objectives use: the full cross of logit
# accuracy levels and dates, uniform over each observed range -- the exact
# mirror of objective_grid(). Uniform in LOGIT accuracy, so levels are
# weighted evenly on the scale the family is linear in; pareto_decline_qtr()
# instead weights levels by the log-cost width of their staircase steps (its
# grid is uniform in log cost), so its average and these fits answer the same
# question under slightly different level weightings.
#
# An alternative that sampled the staircase's own levels at objective_grid's
# (log cost, date) nodes -- matching the check's step-width weighting exactly
# -- was tried and REVERTED: it broke the mirror symmetry with the
# accuracy-direction grids without fixing what motivated it. The contours'
# visible misfit to the records is the linear family's inability to bend
# around the frontier's shelf-and-cliff shape (see plot_cost_frontier.R and
# the moving-ceiling analysis), not a weighting artifact. Each direction's
# uniform grid carries the mirrored hazard, kept knowingly: this one weights
# the clipped top tail heavily (clip_acc stretches n/n scores to la ~ 7,
# where the record is a single best run), as the accuracy-direction grid
# weights the wide flat cost-steps that dilute its -b_t/b_x.
#
#   data     DATA FRAME as for iso_runs(); only the ranges of la and tc are used
#   n_level  SCALAR integer, nodes in logit accuracy
#   n_date   SCALAR integer, nodes in tc
#
# Returns a DATA FRAME of n_level * n_date rows, columns la and tc.
iso_grid <- function(data, n_level = 100, n_date = 40) {
  s <- iso_runs(data)
  expand.grid(la = seq(min(s$la), max(s$la), length.out = n_level),
              tc = seq(min(s$tc), max(s$tc), length.out = n_date))
}

# ln C_a(t) = min{ ln c_i : t_i <= t, la_i >= la } at every grid node: the
# iso_pareto_steps() record (frontier_viz.R) evaluated in bulk, the mirror of
# pareto_grid_response(). Over cost-sorted runs cummax(la) is the running
# best, so the first run whose running best reaches the level is the record.
# Nodes where no run had yet achieved the level are dropped, not imputed.
#
# Returns a DATA FRAME with columns la, tc, lnC -- the defined nodes only.
iso_grid_response <- function(data, n_level = 100, n_date = 40) {
  s <- iso_runs(data)
  gr <- iso_grid(data, n_level, n_date)
  gr$lnC <- NA_real_
  for (tk in unique(gr$tc)) {
    el <- s[s$tc <= tk, ]
    if (!nrow(el)) next
    m <- cummax(el$la)
    rows <- which(gr$tc == tk)
    gr$lnC[rows] <- vapply(gr$la[rows], function(a) {
      j <- which(m >= a)[1]
      if (is.na(j)) NA_real_ else el$lncost[j]
    }, numeric(1))
  }
  gr[!is.na(gr$lnC), , drop = FALSE]
}

# OLS of the sampled record on logit accuracy and date. The coefficient on tc
# IS the annual log cost decline at fixed performance -- no ratio, no delta
# method -- and misfit costs squared log dollars, so the wide flat regions
# that diluted the accuracy-direction b_x are priced at their true horizontal
# size. No standard errors are meaningful: grid nodes are not observations
# (lm will happily print SEs; the class tag is what table code must key on to
# refuse them, as est_se() does for pareto_grid_logit).
#
# Returns an lm with class "lncost_grid_ols" prepended, carrying n_grid and
# n_corners attributes with the same meanings as fit_pareto_logit()'s.
# `formula` is one of COST_FORMS (or COST_BC_FORM), written with the runs'
# response name (lncost); the grid's sampled record is aliased to it here.
# `grid_augment` adds columns derived from the grid's (la, tc) coordinates --
# how the BC terms become evaluable at grid nodes, as in fit_pareto_logit().
fit_lncost_grid <- function(data, formula = COST_FORMS$lin,
                            n_level = 100, n_date = 40, grid_augment = NULL) {
  gr <- iso_grid_response(data, n_level, n_date)
  gr$lncost <- gr$lnC
  if (!is.null(grid_augment)) gr <- grid_augment(gr)
  fit <- lm(formula, data = gr)
  s <- iso_runs(data)
  attr(fit, "n_grid")    <- nrow(gr)
  attr(fit, "n_corners") <- length(pareto_binding(s$lncost, s$tc, s$la))
  class(fit) <- c("lncost_grid_ols", class(fit))
  fit
}

# The dual strict envelope: the highest surface ln C(la, tc) lying at or
# below every run's log cost -- achieving what run i achieved, by run i's
# date, can cost no more than run i paid. Monotone by constraint:
# d lnC/d la >= 0 (better performance cannot be cheaper) and d lnC/d tc <= 0
# (a released model stays available, so the record cannot rise). With the
# quadratic each derivative is affine in (la, tc) jointly, so requiring it at
# the FOUR CORNERS of the grid rectangle enforces it throughout -- the same
# vertex argument fit_envelope() uses; for the linear specification the
# corners collapse to one row per derivative. The run constraints are the
# exact duals of fit_envelope()'s and the same Pareto reduction applies:
# run j implies run i's constraint when it is no dearer, no later, and scored
# no worse -- the identical pareto_binding() call. Everything is linear in
# the coefficients, so this is a small LP; SLSQP dispatches it exactly, and
# the best FEASIBLE candidate is kept. As in fit_envelope(), the quadratic
# also starts from the linear solution, so the richer model can never
# legitimately do worse than the one it nests.
#
# Returns an object of class "cost_envelope_frontier": a LIST with
#   coefficients   named numeric VECTOR (Intercept, la, tc)
#   value          SCALAR, mean fitted log cost over the grid (the objective)
#   worst_slack    SCALAR, minimum slack over all constraints; negative past
#                  -1e-8 means infeasible
#   n_binding      SCALAR integer, candidates surviving the Pareto reduction
#   slack_envelope SCALAR, minimum slack over the RUN constraints, in log
#                  dollars; 0 means the plane touches the data
#   slack_mono     SCALAR, minimum slack over the two monotonicity rows
#   tightest_row   SCALAR integer, row of `data` (post zero-drop) nearest the
#                  plane
#   env_slack      numeric VECTOR of run slacks, parallel to bind
#   bind           integer VECTOR of candidate row indices
fit_cost_envelope <- function(data, formula = COST_FORMS$lin,
                              n_level = 100, n_date = 40, margin = 0.05,
                              grid_augment = NULL) {
  s <- iso_runs(data)
  bind <- pareto_binding(s$lncost, s$tc, s$la)
  X  <- model.matrix(formula, s)
  Xb <- X[bind, , drop = FALSE]
  cb <- s$lncost[bind]
  gr <- iso_grid(data, n_level, n_date)
  if (!is.null(grid_augment)) gr <- grid_augment(gr)
  Xg <- model.matrix(delete.response(terms(formula)), gr)

  # monotonicity rows requiring row %*% g >= 0. e() returns a zero vector for
  # an absent term, so each block covers every specification it applies to.
  nmv <- colnames(X)
  e <- function(nm) {
    v <- numeric(ncol(X))
    if (nm %in% nmv) v[match(nm, nmv)] <- 1
    v
  }
  if ("phia" %in% nmv) {
    # Box-Cox terms (fit_cost_bc): lnC = g0 + ga*phia + gt*phit +
    # gat*phia*phit, and phi is increasing in its argument whatever lambda
    # is, so monotonicity in accuracy and date IS monotonicity in phia and
    # phit. Each derivative is linear in ONE other coordinate, so the ends of
    # that coordinate's grid range bound it -- fit_envelope()'s phic branch,
    # mirrored:
    #    d lnC/d phia =  ga + gat*phit  >= 0
    #   -d lnC/d phit = -(gt + gat*phia) >= 0
    mono <- unique(rbind(
      t(vapply(range(gr$phit), function(p) e("phia") + p * e("phiat"),
               numeric(ncol(X)))),
      t(vapply(range(gr$phia), function(p) -(e("phit") + p * e("phiat")),
               numeric(ncol(X))))))
  } else {
    # quadratic (or linear): each derivative is affine in (la, tc) JOINTLY,
    # so the four corners of the rectangle bound it -- fit_envelope()'s
    # corner blocks, mirrored:
    #    d lnC/d la  =  g_a + 2*g_aa*la + g_at*tc  >= 0
    #   -d lnC/d tc  = -(g_t + 2*g_tt*tc + g_at*la) >= 0
    # Without the second-order terms the corners collapse and unique() keeps
    # one row per derivative.
    corners <- expand.grid(la0 = range(gr$la), tc0 = range(gr$tc))
    mono <- unique(rbind(
      t(mapply(function(la0, tc0)
        e("la") + 2 * la0 * e("I(la^2)") + tc0 * e("la:tc"),
        corners$la0, corners$tc0)),
      t(mapply(function(la0, tc0)
        -(e("tc") + 2 * tc0 * e("I(tc^2)") + la0 * e("la:tc")),
        corners$la0, corners$tc0))))
  }

  # Feasible start: OLS through the runs, curvature and interactions zeroed
  # (they enter both derivatives, so leaving lm's values can hand SLSQP an
  # infeasible start), slopes forced to the monotone signs, intercept dropped
  # until the surface clears under every run.
  b0 <- coef(lm(formula, data = s))
  for (nm in c("I(la^2)", "I(tc^2)", "la:tc", "phiat"))
    if (nm %in% names(b0)) b0[[nm]] <- 0
  for (nm in c("la", "phia"))
    if (nm %in% names(b0)) b0[[nm]] <- max(b0[[nm]], 0.05)
  for (nm in c("tc", "phit"))
    if (nm %in% names(b0)) b0[[nm]] <- min(b0[[nm]], 0)
  b0[[1]] <- b0[[1]] - max(0, max(drop(Xb %*% b0) - cb)) - margin
  starts <- list(unname(b0))
  tl <- attr(terms(formula), "term.labels")
  if (any(grepl("^I\\(|:", tl))) {
    inner <- fit_cost_envelope(data, COST_FORMS$lin, n_level, n_date, margin)
    bs <- setNames(numeric(ncol(X)), nmv)
    bs[names(coef(inner))] <- coef(inner)
    starts <- c(starts, list(unname(bs)))
  }

  grad <- -colMeans(Xg)   # maximize mean height = minimize its negative
  slack <- function(g) min(c(cb - drop(Xb %*% g), drop(mono %*% g)))
  # Coefficient box: with near-collinear columns (a Box-Cox lambda that makes
  # phi(tau) nearly constant) the LP can have an unbounded improving ray, and
  # SLSQP will ride it to 1e28. The box keeps the solver finite;
  # fit_cost_bc() then treats a solution AT the box as a failed fit, so the
  # lambda profile steers away from the degenerate region.
  bound <- 1e5
  run1 <- function(x0) {
    r <- nloptr::nloptr(
      x0 = x0,
      eval_f = function(g) list(objective = -mean(Xg %*% g), gradient = grad),
      eval_g_ineq = function(g) list(
        constraints = c(drop(Xb %*% g) - cb, -drop(mono %*% g)),
        jacobian = rbind(Xb, -mono)),
      lb = rep(-bound, ncol(X)), ub = rep(bound, ncol(X)),
      opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-10,
                  maxeval = 5000, print_level = 0))
    list(g = r$solution, obj = -mean(Xg %*% r$solution),
         slack = slack(r$solution))
  }
  cand <- c(lapply(starts, run1),
            lapply(starts, function(x)
              list(g = x, obj = -mean(Xg %*% x), slack = slack(x))))
  ok <- vapply(cand, function(z) z$slack >= -1e-8, logical(1))
  if (!any(ok)) stop("no feasible cost envelope found")
  best <- cand[ok][[which.min(vapply(cand[ok], `[[`, numeric(1), "obj"))]]

  g <- setNames(best$g, nmv)
  env_slack <- cb - drop(Xb %*% g)
  structure(list(coefficients = g, value = -best$obj,
                 worst_slack = best$slack, n_binding = length(bind),
                 slack_envelope = min(env_slack),
                 slack_mono = min(drop(mono %*% g)),
                 tightest_row = bind[which.min(env_slack)],
                 env_slack = env_slack, bind = bind, formula = formula),
            class = "cost_envelope_frontier")
}

coef.cost_envelope_frontier <- function(object, ...) object$coefficients

## ---- the SFA dual: a stochastic cost frontier -----------------------------------------

# The cost-direction counterpart of the panel SFA (panel_frontier.R):
#
#   ln c_ig = x_ig'gamma + u_g + v_ig,  u_g ~ N+(0, sigma_u,g^2),  v ~ N(0, sigma_v^2)
#
# with x = (1, logit acc, tc) and u at the same model x effort group level. The
# sign flips with the direction: inefficiency now ADDS cost -- a group's u_g is
# how much it overpays, in log dollars, for the accuracy it delivers -- where
# the accuracy-direction u subtracted logits. Unlike the fractional-logit case,
# the Gaussian response gives the group likelihood in CLOSED FORM (the panel
# analogue of Aigner-Lovell-Schmidt / Pitt-Lee): with e_i = y_i - x_i'gamma,
# S1 = sum e_i, S2 = sum e_i^2 within the group, and A = n/sv2 + 1/su2,
#
#   ln L_g = ln 2 - (n/2) ln(2 pi sv2) - (1/2) ln su2 - S2/(2 sv2)
#            + S1^2/(2 sv2^2 A) - (1/2) ln A + ln Phi(S1 / (sv2 sqrt A))
#
# (completing the square in u and integrating over u >= 0; at n = 1 this
# reduces exactly to the classical (2/sigma) phi(e/sigma) Phi(e lambda/sigma)
# cost-frontier density, which the accompanying synthetic test exploits).
# No quadrature, no PIT machinery. The analytic per-GROUP score is attached so
# vcov_robust() (fractional_frontier.R) clusters the sandwich at the group
# level exactly as it does for the accuracy-direction fits; it is validated
# against a numerical Jacobian in the test suite rather than trusted blind.
#
# Same caveats as the accuracy-direction panel fit, mirrored: u_g is entirely
# inefficiency (no group-level noise), and benchmarks with few groups support
# the cross-group variance weakly. New caveat from the flip: logit accuracy is
# a REGRESSOR here, so its binomial sampling noise is measurement error, which
# attenuates gamma_la; the time coefficient is the estimand and is less
# exposed, but not immune where sampling noise correlates with date.
#
# The scale can move with time, mirroring family B: log sigma_u,g = s_g'gs
# with group-constant regressors (Xs_g), the same construction as
# panel_frontier.R's formula_sigma.
#
# READ THE TIME COEFFICIENT WITH CARE. On this data beta_tc comes out
# POSITIVE (cost rising at fixed accuracy) on three of four benchmarks,
# under BOTH scale families -- the time-varying scale was tried precisely as
# a rescue and does not change the sign (its LR against the constant scale
# is weak everywhere). That is a finding about what this model measures, not
# a defect in the fit: tc is group-constant (release date is a model
# property), so beta_tc is identified from cross-group variation that u_g
# also absorbs, and the half-normal's mode at zero places the frontier at
# the DENSE cheap edge of the model x effort distribution. That edge really
# does drift dearer over time -- ever more of the cells are expensive
# reasoning configurations, at every accuracy level -- while the record
# (the minimum, not the dense edge) collapses. So this fit answers "what
# does a typically-efficient model-effort cell pay for accuracy a", the
# cost-direction dual of the S/A/B story about the whole distribution; for
# the record's decline, the grid OLS and envelope duals above are the
# frontier-per-se mirrors.
make_cost_sfa_loglik <- function(y, Xb, Xs_g, grp) {
  G <- max(grp)
  n_g <- tabulate(grp, G)
  k <- ncol(Xb)
  function(theta) {
    gam <- theta[seq_len(k)]
    sv2 <- exp(2 * theta[[k + 1]])
    su2 <- exp(2 * as.vector(Xs_g %*% theta[k + 1 + seq_len(ncol(Xs_g))]))  # G
    e  <- y - as.vector(Xb %*% gam)
    S1 <- rowsum(e, grp)[, 1]
    S2 <- rowsum(e^2, grp)[, 1]
    A  <- n_g / sv2 + 1 / su2
    h  <- S1 / (sv2 * sqrt(A))
    logPhi <- pnorm(h, log.p = TRUE)
    P  <- S1^2 / (sv2^2 * A)
    # base::pi, defensively: the figure scripts' `pi <- iso_acc_plot(...)`
    # idiom (plot_paretologit.R) leaves a ggplot object named pi in the
    # global environment, which is where source()d code resolves names
    ll <- log(2) - (n_g / 2) * log(2 * base::pi * sv2) - 0.5 * log(su2) -
      S2 / (2 * sv2) + P / 2 - 0.5 * log(A) + logPhi

    # score, per group. Mills ratio in logs: both phi and Phi underflow
    # separately for very negative h. kap + rho = 1 by construction (the two
    # variance components' shares of A), a cheap identity the test checks.
    Mh   <- exp(dnorm(h, log = TRUE) - logPhi)
    uhat <- S1 / (sv2 * A) + Mh / sqrt(A)   # posterior mean of u_g
    kap  <- n_g / (sv2 * A)
    rho  <- 1 / (su2 * A)
    ci   <- (e - uhat[grp]) / sv2           # dl/dgamma weights, per observation
    dsu  <- -1 + rho * (P + 1 + Mh * h)     # dl/d(log sigma_u,g), per group
    attr(ll, "gradient") <- cbind(
      rowsum(Xb * ci, grp),
      logsig_v = -n_g + S2 / sv2 - 2 * P + kap * P + kap + Mh * h * (kap - 2),
      Xs_g * dsu)
    ll
  }
}

# Fit for ONE benchmark. Runs with zero accuracy are dropped by iso_runs();
# u_group members must be columns of `data`. formula_sigma follows
# panel_frontier.R's convention: ~ 1 for the constant scale (the A family's
# dual, and the default -- the time-varying scale earns no LR keep here),
# ~ tc for B's dual; its regressors must be group-constant.
#
# Returns a maxLik fit whose coefficients are beta_(Intercept), beta_la,
# beta_tc, logsig_v, then logsig_* for the u scale -- the beta_ prefix so
# cost_coefs() and the summary helpers treat it like the other maxLik fits,
# and the u-scale names matching the accuracy-direction fits' so
# sigma_u_hat() reads the reference-date scale. Attributes n_groups, n_obs.
fit_cost_sfa <- function(data, formula_beta = COST_FORMS$lin,
                         u_group = c("model", "effort"),
                         formula_sigma = ~ 1,
                         start = NULL, method = "BFGS", ...) {
  s <- iso_runs(data)
  key <- interaction(s[u_group], drop = TRUE, lex.order = TRUE)
  grp <- as.integer(key)
  y  <- model.response(model.frame(formula_beta, s))
  Xb <- model.matrix(formula_beta, s)
  Xs <- model.matrix(formula_sigma, s)
  first <- match(seq_len(nlevels(key)), grp)
  Xs_g <- Xs[first, , drop = FALSE]
  if (!isTRUE(all.equal(as.vector(Xs), as.vector(Xs_g[grp, , drop = FALSE]))))
    stop("formula_sigma regressors vary within a u_group; they must be constant")
  ll <- make_cost_sfa_loglik(y, Xb, Xs_g, grp)

  if (is.null(start)) {
    # OLS through the runs; the group means of its residuals proxy u, whose
    # spread seeds sigma_u and whose mean E[u] = sigma_u*sqrt(2/pi) is what
    # separates the frontier's intercept from OLS's.
    f0 <- lm(formula_beta, data = s)
    r  <- residuals(f0)
    gm <- tapply(r, grp, mean)
    su0 <- max(sd(gm), 0.1)
    sv0 <- max(sd(r - gm[grp]), 0.05)
    b0 <- unname(coef(f0))
    start <- c(b0[1] - su0 * sqrt(2 / base::pi), b0[-1],
               log(sv0), log(su0), numeric(ncol(Xs_g) - 1))
  }
  names(start) <- c(paste0("beta_", colnames(Xb)), "logsig_v",
                    paste0("logsig_", colnames(Xs_g)))
  res <- maxLik(ll, start = start, method = method, ...)
  attr(res, "loglik_fn") <- ll
  attr(res, "n_groups")  <- nlevels(key)
  attr(res, "n_obs")     <- length(y)
  res
}

## ---- shared accessors and curve builders ----------------------------------------------

# One accessor for every quadratic-family cost fit: all six possible
# coefficients, absent terms as 0, whatever the fit's naming convention (lm
# and the envelope use bare names, maxLik the beta_ prefix). The mirror of
# frontier_coefs(); the same warning applies -- a term present in the FIT but
# absent here would be silently dropped from every plotted surface.
cost_coefs <- function(fit) {
  cf <- coef(fit)
  names(cf) <- sub("^beta_", "", names(cf))
  get1 <- function(nm) if (nm %in% names(cf)) unname(cf[[nm]]) else 0
  c(g0 = get1("(Intercept)"), ga = get1("la"), gt = get1("tc"),
    gaa = get1("I(la^2)"), gtt = get1("I(tc^2)"), gat = get1("la:tc"))
}

# The fitted surface and its two derivatives, mirrors of frontier_index(),
# frontier_dcost() and frontier_dtime().
cost_index <- function(co, la, tc) {
  co[["g0"]] + co[["ga"]] * la + co[["gt"]] * tc + co[["gaa"]] * la^2 +
    co[["gtt"]] * tc^2 + co[["gat"]] * la * tc
}
cost_dacc  <- function(co, la, tc) co[["ga"]] + 2 * co[["gaa"]] * la + co[["gat"]] * tc
cost_dtime <- function(co, la, tc) co[["gt"]] + 2 * co[["gtt"]] * tc + co[["gat"]] * la

# A BC cost fit announces itself by the lambda attribute fit_cost_bc() stamps
# on it, exactly as is_bc_fit() detects the accuracy-direction BC fits.
is_cost_bc <- function(fit) !is.null(attr(fit, "bc_lambda"))

# One evaluation interface over all three specifications: `f(la, tc)` is the
# fitted ln cost and `dacc(la, tc)` carries the SIGN of d lnC / d la (for the
# BC fit the positive factor d phia/d la is dropped -- only the sign is
# consumed, by the rising-branch filter below). `sub` is the benchmark's
# runs, needed only to recover the BC time offset.
cost_surface <- function(fit, sub) {
  if (is_cost_bc(fit)) {
    cf <- coef(fit)
    names(cf) <- sub("^beta_", "", names(cf))
    lam <- attr(fit, "bc_lambda")
    lo <- unname(lam[[1]]); lt <- unname(lam[[2]])
    off <- (sub$year - sub$tc)[1] - BC_T0
    g0 <- cf[["(Intercept)"]]; ga <- cf[["phia"]]
    gt <- cf[["phit"]]; gat <- cf[["phiat"]]
    list(f = function(la, tc) {
           pa <- bc_tf(exp(la), lo)          # odds = exp(logit)
           pt <- bc_tf(tc + off, lt)
           g0 + ga * pa + gt * pt + gat * pa * pt
         },
         dacc = function(la, tc) ga + gat * bc_tf(tc + off, lt))
  } else {
    co <- cost_coefs(fit)
    list(f    = function(la, tc) cost_index(co, la, tc),
         dacc = function(la, tc) cost_dacc(co, la, tc))
  }
}

# The quarterly decline a LINEAR cost fit implies: gt is d ln C / dt in log
# dollars per year (negative when cost falls), so the quarterly percentage
# drop is the direct transform -- one coefficient, where the
# accuracy-direction fits need the -b_t/b_x ratio. Linear specification only:
# with curvature or transforms the time slope moves with (accuracy, date) and
# no single rate describes the fit.
cost_decline_qtr <- function(fit) 100 * (1 - exp(cost_coefs(fit)[["gt"]] / 4))

# Accuracy-versus-cost curves at each drawn date, for frontier_plot(). NOT an
# inversion: the surface is SWEPT parametrically along the observed logit
# accuracy range and plotted as (exp(lnC), plogis(la)) -- exact for every
# specification, no root-finding, no fold bookkeeping. Points where the
# surface falls in accuracy are blanked (the sweep would double back in cost
# and geom_line cannot draw a fold); the rising region is an interval in la
# for every family here, and the cost-range clip then trims its ends, so each
# date's kept points form one contiguous monotone piece.
#
# Same argument conventions as frontier_curves() (frontier_viz.R).
cost_frontier_curves <- function(fitset, data, dates_by_bench, tbar,
                                 n_la = 200) {
  do.call(rbind, lapply(names(fitset), function(b) {
    sub <- data[data$benchmark == b, ]
    srf <- cost_surface(fitset[[b]], sub)
    su <- iso_runs(sub)
    urng <- range(sub$lncost)
    g <- expand.grid(la = seq(min(su$la), max(su$la), length.out = n_la),
                     qdate = dates_by_bench[[b]])
    tc <- as_t(g$qdate) - tbar[[b]]
    u <- srf$f(g$la, tc)
    u[srf$dacc(g$la, tc) <= 0] <- NA_real_
    u[u < urng[1] | u > urng[2]] <- NA_real_
    out <- data.frame(cost = exp(u), value = plogis(g$la), qdate = g$qdate,
                      benchmark = b, year = 2023 + as_t(g$qdate))
    out[!is.na(out$cost), , drop = FALSE]
  }))
}

# Iso-accuracy contours read DIRECTLY off the cost surface -- this direction's
# native view, no inversion for ANY specification: the surface is evaluated
# at each level's logit. Clipped exactly as iso_acc_curves() clips
# (frontier_viz.R): to the observed cost range, tightened per level to the
# empirical record where `cost_cap` (in practice iso_pareto_curves()'s
# staircases) supplies one, so a level no run ever achieved gets no contour.
# Routed through iso_segments() so clip gaps break the drawn line rather than
# being bridged by a chord.
cost_iso_curves <- function(fitset, data, tbar,
                            levels = seq(0.10, 0.90, by = 0.20), n_date = 300,
                            cost_cap = NULL) {
  do.call(rbind, lapply(names(fitset), function(b) {
    sub <- data[data$benchmark == b, ]
    srf <- cost_surface(fitset[[b]], sub)
    urng <- range(sub$lncost)
    dts <- seq(min(sub$releasedate), max(sub$releasedate), length.out = n_date)
    g <- expand.grid(date = dts, acc = levels)
    u <- srf$f(qlogis(g$acc), as_t(g$date) - tbar[[b]])

    # per-level upper clip, same matching rule as iso_acc_curves(): exact on
    # the shared `levels` values, -Inf (no contour) for a level the cap data
    # never covers
    umax <- rep(urng[2], nrow(g))
    if (!is.null(cost_cap)) {
      cb <- cost_cap[cost_cap$benchmark == b, ]
      lv <- unique(cb$acc)
      caps <- vapply(lv, function(a) max(cb$cost[cb$acc == a]), numeric(1))
      cap_u <- log(caps)[match(g$acc, lv)]
      umax <- ifelse(is.na(cap_u), -Inf, pmin(umax, cap_u))
    }
    u[u < urng[1] | u > umax] <- NA_real_

    g$cost <- exp(u)
    g$branch <- "rising"        # single-valued: monotone in accuracy
    g$disc <- NA_real_
    g$benchmark <- b
    iso_segments(g)
  }))
}

## ---- the Box-Cox cost specification ----------------------------------------------------

# The third member of the cost-direction specification family, mirroring
# boxcox_frontier.R:
#
#   ln cost = g0 + ga*phi(odds; lambda_odds) + gt*phi(tau; lambda_time)
#                + gat*phi(odds)*phi(tau)
#
# with phi the Box-Cox transform, odds = a/(1-a) = exp(logit a) -- accuracy
# re-expressed as a positive quantity phi can act on -- and tau = years since
# BC_T0, exactly as in the accuracy direction. At lambda_odds = 0 phi(odds)
# IS logit accuracy, so the specification nests the linear cost model
# (lambda_odds = 0, lambda_time = 1, gat = 0), three restrictions, the same
# count as in the accuracy direction; it does not nest the quadratic. It also
# spans the raw-vs-logit regressor question empirically: the profiled
# lambda_odds lets the data place the level axis between logit-like and
# level-like shapes instead of the analyst picking one.
#
# Estimation is by PROFILE, as for the accuracy BC: at fixed lambdas the
# model is linear in its coefficients, so each family's own fitter runs on
# pre-transformed columns and the outer optimiser moves the lambdas against
# that family's own objective -- SSR for the two least-squares fits, the
# closed-form likelihood for the SFA duals, mean fitted height (maximized)
# for the envelope. Coefficient standard errors are conditional on the
# profiled lambdas, which are reported without standard errors. lambda_time
# is fixed at 1 where the observed tau span cannot identify it (bc_lt_free;
# fm13), exactly as in the accuracy direction.

COST_BC_FORM <- lncost ~ phia + phit + phiat

# lambda_odds search box: odds spans roughly e^-7..e^7 here, so the same
# considerations as BC_BOX_C apply and the same box serves.
COST_BC_BOX_O <- BC_BOX_C

# The transformed columns, on run data carrying la (i.e. AFTER iso_runs)...
cost_bc_augment <- function(s, lo, lt, off) {
  s$phia  <- bc_tf(exp(s$la), lo)
  s$phit  <- bc_tf(s$tc + off, lt)
  s$phiat <- s$phia * s$phit
  s
}

# ... and on the (la, tc) grid, for the grid OLS response and the envelope
# objective; `off` is the benchmark's tc -> tau shift, as in bc_grid_augment.
cost_bc_grid_augment <- function(lo, lt, off) {
  function(gr) {
    gr$phia  <- bc_tf(exp(gr$la), lo)
    gr$phit  <- bc_tf(gr$tc + off, lt)
    gr$phiat <- gr$phia * gr$phit
    gr
  }
}

# Profile fit for one cost model x benchmark. Returns the refitted inner fit
# at the profiled optimum, carrying the same bc_lambda / bc_lambda_free
# attributes fit_bc() stamps, so est_se_bc() and the curve builders treat
# both directions' BC fits identically.
fit_cost_bc <- function(key, data, lambda_start = c(0, 1)) {
  s <- iso_runs(data)
  off <- (s$year - s$tc)[1] - BC_T0
  lt_free <- bc_lt_free(s$year - BC_T0)
  ws <- new.env(parent = emptyenv())

  fit_at <- function(lo, lt) {
    sa <- cost_bc_augment(s, lo, lt, off)
    ga <- cost_bc_grid_augment(lo, lt, off)
    if (key == "costols") {
      f <- lm(COST_BC_FORM, data = sa)
      list(fit = f, obj = -sum(residuals(f)^2))
    } else if (key == "costgridols") {
      f <- fit_lncost_grid(sa, COST_BC_FORM, grid_augment = ga)
      list(fit = f, obj = -sum(residuals(f)^2))
    } else if (key == "costenvelope") {
      f <- fit_cost_envelope(sa, COST_BC_FORM, grid_augment = ga)
      # a solution with huge coefficients means the LP degenerated (see the
      # coefficient box in fit_cost_envelope): score these lambdas as
      # unfittable rather than letting a numerically meaningless objective
      # win the profile
      if (max(abs(coef(f))) > 1e4) stop("degenerate envelope LP")
      list(fit = f, obj = f$value)   # the HIGHEST surface under the runs
    } else {
      f <- fit_cost_sfa(sa, COST_BC_FORM,
                        formula_sigma = if (key == "costsfab") ~ tc else ~ 1,
                        start = ws$start)
      ws$start <- coef(f)   # warm-start the next profile evaluation
      list(fit = f, obj = as.numeric(logLik(f)))
    }
  }

  # A failed inner fit scores -1e6: bad, not fatal, so the outer search
  # steers away from lambdas where the model cannot be fitted.
  neg <- function(l) {
    lo <- l[1]; lt <- if (lt_free) l[2] else 1
    if (lo < COST_BC_BOX_O[1] || lo > COST_BC_BOX_O[2] ||
        lt < BC_BOX_T[1] || lt > BC_BOX_T[2]) return(1e6)
    r <- tryCatch(fit_at(lo, lt), error = function(e) NULL)
    if (is.null(r) || !is.finite(r$obj)) return(1e6)
    -r$obj
  }
  if (lt_free) {
    opt <- optim(lambda_start, neg, method = "Nelder-Mead",
                 control = list(reltol = 1e-6, maxit = 300))
    lam <- opt$par
  } else {
    opt <- optimize(function(x) neg(c(x, 1)), interval = COST_BC_BOX_O,
                    tol = 1e-4)
    lam <- c(opt$minimum, 1)
  }

  fit <- fit_at(lam[1], lam[2])$fit
  attr(fit, "bc_lambda") <- c(lambda_odds = lam[1], lambda_time = lam[2])
  attr(fit, "bc_lambda_free") <- c(lambda_odds = TRUE, lambda_time = lt_free)
  fit
}
