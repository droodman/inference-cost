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

  # Warm starts: consecutive outer evaluations sit at nearby lambdas, so the
  # previous coefficient vector is close to the next optimum and cuts the SFA
  # fits from dozens of BFGS iterations to a few. A failed inner fit scores
  # 1e6 - bad, not fatal - so the outer search simply steers away from
  # lambdas where the model cannot be fitted.
  ws <- new.env(parent = emptyenv())
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
    # optimize(), not 1-D Nelder-Mead, which optim() itself warns is unreliable
    opt <- optimize(function(x) neg(c(x, 1)), interval = BC_BOX_C, tol = 1e-4)
    lam <- c(opt$minimum, 1)
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
