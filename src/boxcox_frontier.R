# Box-Cox frontier specification: the third member of the specification family,
# alongside linear and quadratic.
#
#   logit acc = b0 + bx*phi(cost; lambda_c) + bt*phi(tau; lambda_t)
#                  + bxt*phi(cost)*phi(tau)
#
# with phi(x; l) = (x^l - 1)/l, the Box-Cox transform (log at l = 0), applied to
# LEVEL cost per task and to tau = years since BC_T0 (GPT-3's release, mid-2020
# -- level time needs an origin, and log cost does not, which is why the linear
# and quadratic specifications never faced this question).
#
# Why this family: the quadratic's (ln c)^2 term lets the fitted surface bend
# back toward -- sometimes into -- the data. Here curvature comes from the
# transform instead, and phi is strictly increasing whatever lambda is, so the
# surface is monotone in cost at EVERY date; what survives of the quadratic's
# pathology is a single possible sign change over time via the product term,
# exactly the behaviour the plain lncost:tc interaction already had. The
# specification nests the linear one (lambda_c = 0, lambda_t = 1, bxt = 0) --
# three restrictions, the same count the quadratic adds -- but not the quadratic.
#
# Estimation is by PROFILE: at fixed (lambda_c, lambda_t) the model is linear in
# its coefficients, so each family's existing fitter runs on pre-transformed
# columns, and the outer optimisation moves the lambdas against that family's
# own objective -- deviance for S, the likelihood for A/B, mean fitted height
# for the envelope, the grid deviance for the Pareto logit. Coefficient standard
# errors are therefore CONDITIONAL on the profiled lambdas; the lambdas
# themselves are reported without standard errors rather than with invented
# ones.
#
# How the lambdas are searched differs by family: S/A/B use a derivative-free
# Nelder-Mead over the profiled objective, while the two frontier-per-se models
# get gradient-based searches (see the bc_lambda_* block below) that reach the
# same optima 10-80x faster -- groundwork for bootstrapping the whole pipeline.
#
# Identification of lambda_t is weak by construction: Box-Cox shape is read off
# deviation-from-linearity over the observed span, roughly (span)/(8*midpoint),
# which is ~10% for benchmarks observed 2.7-6.1 years after the origin and 0.7%
# for fm13's five months of 2026 data. So fm13's lambda_t is FIXED at 1 (any
# value fits identically there) and everyone else's should be read as a shape
# the data tolerates rather than demands.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("fit_specs.R")          # fit_panel_frontier via panel_frontier.R, U_GROUP, SIGMA_FORM
src_source("envelope_frontier.R")  # fit_envelope, fit_pareto_logit

# BC_T0, bc_tf() and bc_inv() live in frontier_viz.R (sourced via fit_specs.R):
# the model-agnostic curve builders there must evaluate and invert BC fits, and
# defining the primitives here would close a source cycle.

BC_FORM <- acc ~ phic + phit + phixt

# Search boxes for the profile. Cost spans ~5 decades, so lambda_c much below 0
# turns the cheapest runs into transform values of ~1e9 and the design matrix
# into rubble; the box stops the optimiser before the arithmetic does. tau spans
# only [2.7, 6.1], where any lambda is numerically tame.
BC_BOX_C <- c(-1, 2)
BC_BOX_T <- c(-2, 3)

# lambda_t is estimable only when the observed tau span is wide relative to its
# distance from the origin; a ratio near 1 (fm13: 6.06/5.71 = 1.06) makes phi
# affine in tau to within a fraction of a percent for every lambda.
bc_lt_free <- function(tau) max(tau) / min(tau) > 1.2

# d phi(x; l) / d l, the transform's own lambda-derivative, vectorised over x.
# The closed form is 0/0 at l = 0; the series through l^2 keeps the switch at
# |l| = 1e-4 exact to machine precision. Companion to bc_tf (frontier_viz.R),
# but only the gradient-based lambda searches below need it, so it lives here.
bc_dtf <- function(x, l) {
  lx <- log(x)
  if (abs(l) < 1e-4) lx^2 / 2 + l * lx^3 / 3 + l^2 * lx^4 / 8
  else (x^l * (l * lx - 1) + 1) / l^2
}

# The transformed columns, on run data (which carries year) ...
bc_augment_runs <- function(s, lc, lt) {
  s$phic  <- bc_tf(s$cost, lc)
  s$phit  <- bc_tf(s$year - BC_T0, lt)
  s$phixt <- s$phic * s$phit
  s
}

# ... and on the objective grid, whose coordinates are (lncost, tc); `off` is
# the benchmark's tc -> tau shift, recovered from the data as (year - tc) - BC_T0
# in the same spirit as bench_tbar().
bc_grid_augment <- function(lc, lt, off) {
  function(gr) {
    gr$phic  <- bc_tf(exp(gr$lncost), lc)
    gr$phit  <- bc_tf(gr$tc + off, lt)
    gr$phixt <- gr$phic * gr$phit
    gr
  }
}

# One fit at fixed lambdas. Returns list(fit, obj) with obj oriented so that
# LARGER is better for every family, whatever its native objective's sense.
bc_fit_at <- function(key, s, off, lc, lt, start = NULL) {
  sa <- bc_augment_runs(s, lc, lt)
  if (key == "S") {
    f <- glm(BC_FORM, data = sa, family = quasibinomial(link = "logit"))
    list(fit = f, obj = -deviance(f) / 2)
  } else if (key %in% c("A", "B")) {
    f <- fit_panel_frontier(BC_FORM, ~ 1, data = sa, u_group = U_GROUP,
                            formula_sigma = SIGMA_FORM[[key]], dedup = FALSE,
                            fixed = "delta_(Intercept)", start = start)
    list(fit = f, obj = as.numeric(logLik(f)))
  } else if (key == "envelope") {
    f <- fit_envelope(sa, BC_FORM, grid_augment = bc_grid_augment(lc, lt, off))
    list(fit = f, obj = -f$value)   # the lowest-sitting surface, now over lambda too
  } else {
    f <- fit_pareto_logit(sa, BC_FORM, grid_augment = bc_grid_augment(lc, lt, off))
    list(fit = f, obj = -deviance(f) / 2)
  }
}

## ---- dedicated lambda searches for the frontier-per-se models ---------------------
#
# The two models fitted to the Pareto frontier rather than to the runs get
# gradient-based lambda searches in place of the generic Nelder-Mead profile.
# Head-to-head comparisons (experiment_bc_paretologit.R and
# experiment_bc_envelope.R) picked a DIFFERENT winner for each -- the
# asymmetry is structural, not taste:
#
#   paretologit  profile kept, outer L-BFGS-B over the lambdas with the exact
#                gradient via the ENVELOPE THEOREM: at the inner optimum the
#                beta-score is zero, so d(profiled ll)/d lambda is the partial
#                derivative, read off the inner glm's own residuals. The inner
#                problem is a convex canonical-link glm IRLS cannot fail, and
#                the profiled surface is smooth. (Joint L-BFGS-B over
#                (beta, lambda) silently misconverges at default tolerances --
#                the b_x/lambda_c ravine -- and needs ~50x more evaluations
#                once tightened enough to be trustworthy.)
#
#   envelope     JOINT SLSQP over (beta, lambda) at once: analytic objective
#                gradient and constraint Jacobian, run constraints bilinear in
#                (beta, lambda), monotonicity still enforced at the phi range
#                endpoints because phi is increasing whatever lambda is. Two
#                starts, best feasible kept. (The profile-gradient route fails
#                HERE: the value function is only piecewise smooth in lambda --
#                kinks where the active constraint set changes -- and recovering
#                the KKT multipliers adds fragility, not robustness.)
#
# Both searches only pick the lambdas; fit_bc() then refits through bc_fit_at()
# at the optimum, so the returned object is the canonical fit_pareto_logit /
# fit_envelope result whatever route found the lambdas.

bc_lambda_paretologit <- function(s, off, lt_free, lambda_start) {
  # The grid staircase does not depend on the lambdas, so compute it once; only
  # the phi columns move inside the search.
  gr <- pareto_grid_response(s)
  y <- gr$acc
  cost <- exp(gr$lncost)
  tau <- gr$tc + off
  # One-point cache: L-BFGS-B asks for fn and gr at the same lambdas, and one
  # inner fit serves both.
  cache <- new.env(parent = emptyenv())
  at <- function(lc, lt) {
    key <- c(lc, lt)
    if (identical(cache$key, key)) return(cache$val)
    phic <- bc_tf(cost, lc)
    phit <- bc_tf(tau, lt)
    b <- unname(coef(glm(y ~ phic + phit + I(phic * phit),
                         family = quasibinomial(link = "logit"))))
    eta <- b[1] + b[2] * phic + b[3] * phit + b[4] * phic * phit
    r <- y - plogis(eta)
    cache$key <- key
    cache$val <- list(
      ll = sum(y * eta - log1pexp(eta)),
      g = c(sum(r * (b[2] + b[4] * phit) * bc_dtf(cost, lc)),
            sum(r * (b[3] + b[4] * phic) * bc_dtf(tau, lt))))
    cache$val
  }
  l0 <- pmin(pmax(lambda_start, c(BC_BOX_C[1], BC_BOX_T[1])),
             c(BC_BOX_C[2], BC_BOX_T[2]))
  if (lt_free) {
    o <- optim(l0, function(l) -at(l[1], l[2])$ll,
               function(l) -at(l[1], l[2])$g,
               method = "L-BFGS-B",
               lower = c(BC_BOX_C[1], BC_BOX_T[1]),
               upper = c(BC_BOX_C[2], BC_BOX_T[2]))
  } else {
    o <- optim(l0[1], function(l) -at(l, 1)$ll,
               function(l) -at(l, 1)$g[1],
               method = "L-BFGS-B",
               lower = BC_BOX_C[1], upper = BC_BOX_C[2])
  }
  if (o$convergence != 0)
    stop("lambda L-BFGS-B did not converge: ", o$message)
  if (lt_free) o$par else c(o$par, 1)
}

bc_lambda_envelope <- function(s, off, lt_free, lambda_start) {
  # The lambda-invariant setup, hoisted: the logit clipping, the Pareto
  # reduction to binding candidates, and the grid coordinates. Mirrors the
  # front half of fit_envelope(); the phi columns are rebuilt per evaluation.
  L <- qlogis(clip_acc(s$acc, s$n_samples))
  pos <- which(is.finite(L))
  bind <- pos[pareto_binding(s$lncost[pos], s$tc[pos], L[pos])]
  Lb <- L[bind]
  cost_b <- s$cost[bind]
  tau_b <- s$year[bind] - BC_T0
  gr <- objective_grid(s)
  cost_g <- exp(gr$lncost)
  tau_g <- gr$tc + off
  te <- range(tau_g)
  ce <- range(cost_g)

  # Everything the joint solver needs at theta = (beta, lc[, lt]): objective
  # (mean surface height on the grid) and its gradient, constraints
  # (ci - ui(lambda) beta <= 0: run rows then 4 monotonicity rows at the phi
  # range endpoints) and their Jacobian. All closed-form via bc_dtf.
  pieces <- function(b, lc, lt) {
    phicb <- bc_tf(cost_b, lc)
    phitb <- bc_tf(tau_b, lt)
    ui <- rbind(cbind(1, phicb, phitb, phicb * phitb),
                c(0, 1, 0, bc_tf(te[1], lt)), c(0, 1, 0, bc_tf(te[2], lt)),
                c(0, 0, 1, bc_tf(ce[1], lc)), c(0, 0, 1, bc_tf(ce[2], lc)))
    ci <- c(Lb, rep(0, 4))
    phicg <- bc_tf(cost_g, lc)
    phitg <- bc_tf(tau_g, lt)
    Xg <- cbind(1, phicg, phitg, phicg * phitg)
    p <- plogis(drop(Xg %*% b))
    w <- p * (1 - p)
    dg_lc <- c(-(b[2] + b[4] * phitb) * bc_dtf(cost_b, lc), 0, 0,
               -b[4] * bc_dtf(ce[1], lc), -b[4] * bc_dtf(ce[2], lc))
    dg_lt <- c(-(b[3] + b[4] * phicb) * bc_dtf(tau_b, lt),
               -b[4] * bc_dtf(te[1], lt), -b[4] * bc_dtf(te[2], lt), 0, 0)
    list(obj = mean(p),
         g_obj = c(drop(crossprod(Xg, w)) / nrow(Xg),
                   mean(w * (b[2] + b[4] * phitg) * bc_dtf(cost_g, lc)),
                   mean(w * (b[3] + b[4] * phicg) * bc_dtf(tau_g, lt))),
         g_con = as.vector(ci - ui %*% b),
         jac = cbind(-ui, dg_lc, dg_lt),
         ui = ui, ci = ci)
  }

  # Feasible start at given lambdas, replicating fit_envelope()'s b0: glm on
  # all runs, curvature zeroed, slopes forced positive (so the monotonicity
  # rows hold), intercept lifted to clear every run constraint.
  cold <- function(lc, lt) {
    phic <- bc_tf(s$cost, lc)
    phit <- bc_tf(s$year - BC_T0, lt)
    b0 <- unname(coef(glm(s$acc ~ phic + phit + I(phic * phit),
                          family = quasibinomial(link = "logit"))))
    b0[4] <- 0
    b0[2] <- max(b0[2], 0.05)
    b0[3] <- max(b0[3], 0.05)
    pp <- pieces(b0, lc, lt)
    b0[1] <- b0[1] + max(0, max(pp$g_con)) + 0.05
    b0
  }

  free_i <- if (lt_free) 1:6 else 1:5
  run1 <- function(lam0) {
    th0 <- c(cold(lam0[1], lam0[2]), if (lt_free) lam0 else lam0[1])
    unpack <- function(th) list(b = th[1:4], lc = th[5],
                                lt = if (lt_free) th[6] else 1)
    r <- nloptr::nloptr(
      x0 = th0,
      eval_f = function(th) {
        u <- unpack(th)
        pp <- pieces(u$b, u$lc, u$lt)
        list(objective = pp$obj, gradient = pp$g_obj[free_i])
      },
      eval_g_ineq = function(th) {
        u <- unpack(th)
        pp <- pieces(u$b, u$lc, u$lt)
        list(constraints = pp$g_con, jacobian = pp$jac[, free_i, drop = FALSE])
      },
      lb = c(rep(-Inf, 4), BC_BOX_C[1], if (lt_free) BC_BOX_T[1]),
      ub = c(rep( Inf, 4), BC_BOX_C[2], if (lt_free) BC_BOX_T[2]),
      opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-10,
                  maxeval = 5000, print_level = 0))
    u <- unpack(r$solution)
    pp <- pieces(u$b, u$lc, u$lt)
    list(lam = c(u$lc, u$lt), obj = pp$obj, slack = -max(pp$g_con))
  }

  # Two starts, best feasible kept: agreement between them is the working
  # check against the silent-misconvergence failure mode the joint approach
  # showed on the (unconstrained) Pareto logit.
  l0 <- pmin(pmax(lambda_start, c(BC_BOX_C[1], BC_BOX_T[1])),
             c(BC_BOX_C[2], BC_BOX_T[2]))
  starts <- unique(list(l0, c(1, 2)))
  cand <- lapply(starts, run1)
  ok <- vapply(cand, function(z) z$slack >= -1e-8, logical(1))
  if (!any(ok)) stop("no feasible joint envelope solution from either start")
  cand[ok][[which.min(vapply(cand[ok], `[[`, numeric(1), "obj"))]]$lam
}

# Profile fit for one family x benchmark. `s` is that benchmark's rows of the
# analysis data; `lambda_start` seeds the outer optimiser (the caller passes S's
# profiled lambdas to the slower families -- a starting point, not a constraint).
#
# Returns the refitted inner fit at the profiled optimum, carrying attributes
#   bc_lambda       named numeric length 2, the (lambda_cost, lambda_time) used
#   bc_lambda_free  named logical length 2, FALSE where fixed rather than profiled
# so the table code can print the lambdas and count LR degrees of freedom.
fit_bc <- function(key, s, lambda_start = c(0, 1)) {
  off <- (s$year - s$tc)[1] - BC_T0
  lt_free <- bc_lt_free(s$year - BC_T0)
  ws <- new.env(parent = emptyenv())

  # The frontier-per-se models take their dedicated gradient-based searches
  # (bc_lambda_* above); a failure there falls through to the generic profile
  # below rather than aborting -- the slow path is always a correct path.
  lam <- if (key %in% c("paretologit", "envelope")) {
    tryCatch(
      if (key == "paretologit")
        bc_lambda_paretologit(s, off, lt_free, lambda_start)
      else bc_lambda_envelope(s, off, lt_free, lambda_start),
      error = function(e) {
        message("bc_lambda_", key, " failed (", conditionMessage(e),
                "); falling back to the Nelder-Mead profile")
        NULL
      })
  } else NULL

  if (is.null(lam)) {
    # Warm starts: consecutive outer evaluations sit at nearby lambdas, so the
    # previous coefficient vector is close to the next optimum and cuts the SFA
    # fits from dozens of BFGS iterations to a few. A failed inner fit scores
    # 1e6 - bad, not fatal - so the outer search simply steers away from
    # lambdas where the model cannot be fitted.
    neg <- function(l) {
      lc <- l[1]; lt <- if (lt_free) l[2] else 1
      if (lc < BC_BOX_C[1] || lc > BC_BOX_C[2] ||
          lt < BC_BOX_T[1] || lt > BC_BOX_T[2]) return(1e6)
      r <- tryCatch(bc_fit_at(key, s, off, lc, lt, start = ws$start),
                    error = function(e) NULL)
      if (is.null(r) || !is.finite(r$obj)) return(1e6)
      if (key %in% c("A", "B")) ws$start <- coef(r$fit)
      -r$obj
    }

    if (lt_free) {
      opt <- optim(lambda_start, neg, method = "Nelder-Mead",
                   control = list(reltol = 1e-6, maxit = 300))
      lam <- c(opt$par[1], opt$par[2])
    } else {
      # optimize(), not 1-D Nelder-Mead, which optim() itself warns is
      # unreliable
      opt <- optimize(function(x) neg(c(x, 1)), interval = BC_BOX_C, tol = 1e-4)
      lam <- c(opt$minimum, 1)
    }
  }

  r <- bc_fit_at(key, s, off, lam[1], lam[2], start = ws$start)
  fit <- r$fit
  attr(fit, "bc_lambda") <- c(lambda_cost = lam[1], lambda_time = lam[2])
  attr(fit, "bc_lambda_free") <- c(lambda_cost = TRUE, lambda_time = lt_free)
  fit
}

# All benchmarks of one family, on the shared worker cluster when one is
# available (fit_cluster() in paths.R) -- the profiles are independent across
# benchmarks, so the wall clock is the slowest one rather than the sum.
# `lambda_starts` is a list benchmark -> c(lambda_c, lambda_t) seeding the
# outer optimiser; missing entries start at (0, 1), the linear specification.
fit_bc_by <- function(key, data, lambda_starts = list()) {
  # Force the arguments BEFORE the closure ships: an unforced promise crossing
  # the serialisation boundary re-evaluates its expression on the worker, where
  # the caller's variables do not exist. as.list() both forces and lets a
  # caller pass an environment (regression_tables.R's seed cache) directly.
  force(key)
  lambda_starts <- as.list(lambda_starts)
  bs <- sort(unique(data$benchmark))
  one <- function(b) fit_bc(key, data[data$benchmark == b, ],
                            lambda_start = lambda_starts[[b]] %||% c(0, 1))
  cl <- fit_cluster(length(bs))
  fits <- if (is.null(cl)) lapply(bs, one) else parallel::parLapply(cl, bs, one)
  setNames(fits, bs)
}
