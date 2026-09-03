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
# own objective -- deviance for S, the likelihood for A/B, the grid deviance
# for the Pareto logit, constrained or not. Coefficient standard
# errors are therefore CONDITIONAL on the profiled lambdas; the lambdas
# themselves are reported without standard errors rather than with invented
# ones.
#
# How the lambdas are searched differs by family: S/A/B use a derivative-free
# Nelder-Mead over the profiled objective, while the frontier-per-se models
# seed theirs from a gradient-based search (bc_lambda_paretologit below) that
# reaches the neighbourhood 10-80x faster.
#
# Identification of lambda_t is weak by construction: Box-Cox shape is read off
# deviation-from-linearity over the observed span, roughly (span)/(8*midpoint),
# which is ~10% for benchmarks observed 2.7-6.1 years after the origin and 0.7%
# for fm13's five months of 2026 data. So lambda_t is FIXED at 1 where the span
# cannot identify it at all (bc_lt_free: fm13), profiled elsewhere -- and when
# the profiled value rides a BC_BOX_T edge, non-identification in practice, the
# fit falls back to lambda_t = 1 too (bc_lt_stuck, in bc_lambda_search).

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("fit_specs.R")          # fit_panel_frontier via panel_frontier.R, U_GROUP, SIGMA_FORM
src_source("envelope_frontier.R")  # fit_pareto_logit, fit_pareto_logit_env

# BC_T0, bc_tf() and bc_inv() live in frontier_viz.R (sourced via fit_specs.R):
# the model-agnostic curve builders there must evaluate and invert BC fits, and
# defining the primitives here would close a source cycle.

BC_FORM <- acc ~ phic + phit + phixt

# Search boxes for the profile. Cost spans ~5 decades, so lambda_c much below 0
# turns the cheapest runs into transform values of ~1e9 and the design matrix
# into rubble; the box stops the optimiser before the arithmetic does. tau spans
# only [2.7, 6.1], where any lambda is numerically tame -- but not statistically:
# past roughly |lambda_t| = 8, phi(tau) is nearly constant over the span (its
# relative variation ~1e-4 by -8), a near-degenerate column. The old box of
# [-2, 3] amputated genuine interior optima in the cost-direction grid fits
# (aime/gpqa/fm13 peak near -4..-3, chess near +4, with R2 falling away steeply
# on both sides), which is what made those fits ride its edges.
BC_BOX_C <- c(-1, 2)
BC_BOX_T <- c(-8, 8)

# These are BOX-TIDWELL fits: phi() acts on the REGRESSORS only. There is no
# response-side lambda in either direction -- the accuracy models use the
# plain logit link, the cost models model ln cost -- so the search box above
# applies to covariate transforms alone. The reasoning, and what was removed,
# is in frontier_viz.R.

# The gradient seed search (bc_lambda_paretologit) is IDENTICAL for the two
# frontier-per-se keys -- same benchmark, same S-seeded lambda start -- and
# store_bc() profiles the two families back-to-back in the same processes
# (the shared PSOCK workers persist across fit_bc_by calls), so memoize it
# per benchmark rather than paying its ~50 inner glms twice. Guarded like
# the other re-source-safe environments.
if (!exists(".bc_seed_memo", inherits = FALSE))
  .bc_seed_memo <- new.env(parent = emptyenv())

# lambda_t is estimable only when the observed tau span is wide relative to its
# distance from the origin; a ratio near 1 (fm13: 6.06/5.71 = 1.06) makes phi
# affine in tau to within a fraction of a percent for every lambda.
bc_lt_free <- function(tau) max(tau) / min(tau) > 1.2

# Even where the span test passes, the profiled lambda_t often lands on a
# BC_BOX_T edge -- the objective monotone across the whole box, a bound-riding
# shape the data tolerates rather than demands. The box edges are arbitrary
# numerical guards, so an edge value is not an estimate; bc_lambda_search()
# treats one as non-identification and refits with lambda_t locked at 1. The
# tolerance is loose because the penalty wall (1e6 outside the box) can stall
# the simplex slightly short of the edge itself.
bc_lt_stuck <- function(lt, tol = 0.05)
  lt - BC_BOX_T[1] < tol || BC_BOX_T[2] - lt < tol

# The one lambda search every profile below runs: 2-D Nelder-Mead over
# (lambda_c, lambda_t) when lambda_t is free, falling back to the 1-D
# optimize() over `box1` alone -- not 1-D Nelder-Mead, which optim() itself
# warns against -- when lambda_t is fixed a priori (bc_lt_free) OR comes back
# riding a BC_BOX_T edge (bc_lt_stuck). `neg` is the family's profiled
# objective over l = c(lambda_c, lambda_t), smaller-better; the fallback
# passes c(x, 1) explicitly so neg's own lt_free gate is irrelevant there.
# Returns list(lam, lt_free): lt_free FALSE whenever the fallback decided.
bc_lambda_search <- function(neg, start, lt_free, box1 = BC_BOX_C) {
  if (lt_free) {
    opt <- optim(start[1:2], neg, method = "Nelder-Mead",
                 control = list(reltol = 1e-6, maxit = 300))
    if (!bc_lt_stuck(opt$par[2]))
      return(list(lam = c(opt$par[1], opt$par[2]), lt_free = TRUE))
  }
  opt <- optimize(function(x) neg(c(x, 1)), interval = box1, tol = 1e-4)
  list(lam = c(opt$minimum, 1), lt_free = FALSE)
}

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

# bc_qll(), the probability-scale quasi-log-likelihood these fits profile on,
# lives in frontier_viz.R: fit_pareto_logit_env() (envelope_frontier.R) needs
# it too, and this file sources that one.
#
# There was also a fit_pareto_bclink() here: the Pareto-grid fit under a
# parametric link, mu = phi^-1(eta; lambda_odds), hand-rolled on BFGS because
# glm could not take a link with a free shape parameter in its inner loop. It
# is gone with the rest of the response-transform machinery (see
# frontier_viz.R): at lambda_odds = 0 -- the only value now used, and the
# only one without a pole -- it was bit-identical to the plain glm path this
# file already took, so deleting it removes a second implementation of the
# same fit rather than a capability.

# One fit at fixed lambdas. Returns list(fit, obj) with obj oriented so that
# LARGER is better for every family, whatever its native objective's sense.
# Both lambdas are REGRESSOR-side: lc on cost, lt on time. The response is
# never transformed, so every family here uses the plain logit.
# `inv`, when given (the frontier-per-se profiles), is the list of
# lambda-invariant precomputations fit_bc() hoists once per benchmark:
#   gr0   the staircase pareto_grid_response(s)
#   bind  the envelope constraint candidates (row indices into s)
#   ncor  the Pareto corner count for the n_corners attribute
# `start` doubles as the SFA warm start (A/B, as before) and the constrained
# grid fit's SLSQP warm start (paretologitenv).
bc_fit_at <- function(key, s, off, lc, lt, start = NULL,
                      final = TRUE, inv = NULL) {
  if (is.null(inv)) inv <- list()
  sa <- bc_augment_runs(s, lc, lt)
  if (key == "S") {
    # an ordinary quasibinomial glm: IRLS does the estimation, deviance() the
    # profile objective, sandwich/coeftest the robust standard errors. The
    # quasibinomial deviance differs from the probability-scale quasi-ll only
    # by the saturated term, constant in data alone, so -deviance/2 is
    # comparable across both lambdas -- and no Jacobian is needed, the
    # response being untransformed. An IRLS that ran out its iterations is
    # unfittable at these lambdas, not evidence.
    f <- glm(BC_FORM, data = sa, family = quasibinomial(link = "logit"),
             control = glm.control(maxit = 100))
    if (!f$converged) stop("IRLS did not converge at these lambdas")
    list(fit = f, obj = -deviance(f) / 2)
  } else if (key %in% c("A", "B")) {
    # `final = FALSE` inside the lambda profile: maxLik's finite-difference
    # Hessian at each inner optimum costs more than the warm-started
    # optimisation itself, and the profile reads only the objective. The
    # optimiser's path is untouched, so the profiled lambdas are identical;
    # only the refit at the optimum pays for the Hessian the tables need.
    f <- fit_panel_frontier(BC_FORM, ~ 1, data = sa, u_group = U_GROUP,
                            formula_sigma = SIGMA_FORM[[key]], dedup = FALSE,
                            fixed = "delta_(Intercept)", start = start,
                            finalHessian = final)
    list(fit = f, obj = as.numeric(logLik(f)))
  } else if (key == "paretologitenv") {
    # the constrained grid fit's own objective: `value` is the mean negative
    # quasi-ll per node, and n_grid does not move with the lambdas, so -value
    # profiles identically to the sum bc_qll the unconstrained fit reports
    f <- fit_pareto_logit_env(sa, BC_FORM,
                              grid_augment = bc_grid_augment(lc, lt, off),
                              gr0 = inv$gr0,
                              bind = inv$bind, n_corners = inv$ncor,
                              start = start)
    list(fit = f, obj = -f$value)
  } else {
    f <- fit_pareto_logit(sa, BC_FORM,
                          grid_augment = bc_grid_augment(lc, lt, off),
                          gr0 = inv$gr0, n_corners = inv$ncor)
    list(fit = f, obj = bc_qll(f$y, fitted(f)))
  }
}

## ---- dedicated lambda search for the frontier-per-se models -----------------------
#
# The models fitted to the Pareto frontier rather than to the runs seed their
# lambda profiles from a gradient-based search on the UNCONSTRAINED grid
# logit: outer L-BFGS-B over the lambdas with the exact gradient via the
# ENVELOPE THEOREM -- at the inner optimum the beta-score is zero, so
# d(profiled ll)/d lambda is the partial derivative, read off the inner glm's
# own residuals. The inner problem is a convex canonical-link glm IRLS cannot
# fail, and the profiled surface is smooth. (Joint L-BFGS-B over
# (beta, lambda) silently misconverges at default tolerances -- the
# b_x/lambda_c ravine -- and needs ~50x more evaluations once tightened
# enough to be trustworthy; see experiment_bc_paretologit.R.)
#
# The search only picks seed lambdas; fit_bc() then profiles the full triple
# by Nelder-Mead and refits through bc_fit_at() at the optimum, so the
# returned object is the canonical inner fit whatever route found the lambdas.

bc_lambda_paretologit <- function(s, off, lt_free, lambda_start, gr0 = NULL) {
  # The grid staircase does not depend on the lambdas, so compute it once (or
  # take the caller's); only the phi columns move inside the search.
  gr <- if (is.null(gr0)) pareto_grid_response(s) else gr0
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

  # Every family profiles the same two REGRESSOR lambdas, (lambda_cost,
  # lambda_time), over the plain logit link. The frontier-per-se pair still
  # gets its own branch below, but only for machinery: the dedicated gradient
  # search above seeds it, and Nelder-Mead refines from there.
  if (key %in% c("paretologit", "paretologitenv")) {
    # the lambda-invariant pieces, hoisted out of the profile loop and handed
    # to every inner fit (see bc_fit_at's `inv`): the staircase, the envelope
    # constraint candidates, and the corner count -- together they were about
    # half of each profile's time when recomputed per evaluation
    L0 <- qlogis(clip_acc(s$acc, s$n_samples))
    pos <- which(is.finite(L0) & s$acc > 0)
    inv <- list(gr0  = pareto_grid_response(s),
                bind = pos[pareto_binding(s$lncost[pos], s$tc[pos], L0[pos])],
                ncor = length(pareto_binding(s$lncost, s$tc, s$acc)))

    # the envelope-constrained grid fit seeds from the UNCONSTRAINED profile
    # search too: same objective family, and a seed need only land near --
    # which is also why the two keys can share one memoized search
    mk <- paste(s$benchmark[1], lt_free,
                paste(signif(lambda_start[1:2], 12), collapse = ","))
    seed <- .bc_seed_memo[[mk]]
    if (is.null(seed)) {
      seed <- tryCatch(
        bc_lambda_paretologit(s, off, lt_free, lambda_start[1:2],
                              gr0 = inv$gr0),
        error = function(e) {
          message("bc_lambda_paretologit seed for ", key, " failed (",
                  conditionMessage(e),
                  "); seeding the profile from lambda_start")
          c(lambda_start[1], if (lt_free) lambda_start[2] else 1)
        })
      .bc_seed_memo[[mk]] <- seed
    }
    neg <- function(l) {
      lc <- l[1]; lt <- if (lt_free) l[2] else 1
      if (lc < BC_BOX_C[1] || lc > BC_BOX_C[2] ||
          lt < BC_BOX_T[1] || lt > BC_BOX_T[2]) return(1e6)
      r <- tryCatch(bc_fit_at(key, s, off, lc, lt, start = ws$start,
                              inv = inv),
                    error = function(e) NULL)
      if (is.null(r) || !is.finite(r$obj)) return(1e6)
      # warm-start the next constrained solve from this one, as the SFA
      # profiles do; the unconstrained keys' inner fitters self-start
      if (key == "paretologitenv") ws$start <- coef(r$fit)
      -r$obj
    }
    sr <- bc_lambda_search(neg, seed, lt_free)
    lam <- sr$lam; lt_free <- sr$lt_free
    # the canonical refit at the optimum: cold-started, so the returned
    # object is exactly what a standalone fit at these lambdas produces
    r <- bc_fit_at(key, s, off, lam[1], lam[2], inv = inv)
    fit <- r$fit
    attr(fit, "bc_lambda") <- c(lambda_cost = lam[1], lambda_time = lam[2])
    attr(fit, "bc_lambda_free") <- c(lambda_cost = TRUE,
                                     lambda_time = lt_free)
    return(fit)
  }

  if (key == "S") {
    # S profiles the same two regressor lambdas as every other family, so the
    # generic block below would serve -- it keeps its own only because it
    # takes neither a gradient seed nor a warm start (it IS the cheap family
    # the others are seeded from).
    neg <- function(l) {
      lc <- l[1]; lt <- if (lt_free) l[2] else 1
      if (lc < BC_BOX_C[1] || lc > BC_BOX_C[2] ||
          lt < BC_BOX_T[1] || lt > BC_BOX_T[2]) return(1e6)
      r <- tryCatch(bc_fit_at("S", s, off, lc, lt),
                    error = function(e) NULL)
      if (is.null(r) || !is.finite(r$obj)) return(1e6)
      -r$obj
    }
    sr <- bc_lambda_search(neg, lambda_start, lt_free)
    lam <- sr$lam; lt_free <- sr$lt_free
    r <- bc_fit_at("S", s, off, lam[1], lam[2])
    fit <- r$fit
    attr(fit, "bc_lambda") <- c(lambda_cost = lam[1], lambda_time = lam[2])
    attr(fit, "bc_lambda_free") <- c(lambda_cost = TRUE,
                                     lambda_time = lt_free)
    return(fit)
  }

  {
    # Warm starts: consecutive outer evaluations sit at nearby lambdas, so the
    # previous coefficient vector is close to the next optimum and cuts the SFA
    # fits from dozens of BFGS iterations to a few. A failed inner fit scores
    # 1e6 - bad, not fatal - so the outer search simply steers away from
    # lambdas where the model cannot be fitted.
    neg <- function(l) {
      lc <- l[1]; lt <- if (lt_free) l[2] else 1
      if (lc < BC_BOX_C[1] || lc > BC_BOX_C[2] ||
          lt < BC_BOX_T[1] || lt > BC_BOX_T[2]) return(1e6)
      r <- tryCatch(bc_fit_at(key, s, off, lc, lt, start = ws$start,
                              final = FALSE),
                    error = function(e) NULL)
      if (is.null(r) || !is.finite(r$obj)) return(1e6)
      if (key %in% c("A", "B")) ws$start <- coef(r$fit)
      -r$obj
    }

    sr <- bc_lambda_search(neg, lambda_start, lt_free)
    lam <- sr$lam; lt_free <- sr$lt_free
  }

  r <- bc_fit_at(key, s, off, lam[1], lam[2], start = ws$start)
  fit <- r$fit
  attr(fit, "bc_lambda") <- c(lambda_cost = lam[1], lambda_time = lam[2])
  attr(fit, "bc_lambda_free") <- c(lambda_cost = TRUE, lambda_time = lt_free)
  fit
}

# The pooled (ECI-units, benchmark-fixed-effects) accuracy Box-Cox profile,
# for S and the frontier-per-se pair -- the accuracy mirror of
# fit_pooled_cost_bc (cost_frontier.R). phi acts on each run's own cost and
# tau exactly as per benchmark; the transformed terms are then alpha_b-scaled
# (xphic = alpha * phic, ...) so the index is alpha_b * C(phi_c, phi_t) plus
# the fixed effects, and envelope_constraints()'s pooled branch supplies the
# monotonicity rows. No SFA keys, as for the other pooled accuracy fits.
fit_pooled_acc_bc <- function(key, d, lambda_start = c(0, 1)) {
  stopifnot(key %in% c("S", "paretologit", "paretologitenv"))
  sa <- pooled_acc_runs(d)
  off <- (sa$year - sa$tc)[1] - BC_T0
  lt_free <- bc_lt_free(sa$year - BC_T0)
  form <- acc ~ xphic + xphit + xphixt + bench
  gr0 <- if (key == "S") NULL else pooled_acc_grid(sa)
  bind <- if (key == "S") NULL else pooled_acc_binding(sa)
  ws <- new.env(parent = emptyenv())
  aug_runs <- function(x, lc, lt) {
    x <- bc_augment_runs(x, lc, lt)
    x$xphic  <- x$alpha * x$phic
    x$xphit  <- x$alpha * x$phit
    x$xphixt <- x$alpha * x$phixt
    x
  }
  aug_grid <- function(gr, lc, lt) {
    gr$phic   <- bc_tf(exp(gr$lncost), lc)
    gr$phit   <- bc_tf(gr$tc + off, lt)
    gr$phixt  <- gr$phic * gr$phit
    gr$xphic  <- gr$alpha * gr$phic
    gr$xphit  <- gr$alpha * gr$phit
    gr$xphixt <- gr$alpha * gr$phixt
    gr
  }
  fit_at <- function(lc, lt) {
    sr <- aug_runs(sa, lc, lt)
    if (key == "S") {
      f <- suppressWarnings(glm(form, data = sr,
                                family = quasibinomial(link = "logit")))
      if (!f$converged) stop("IRLS did not converge")
      return(list(fit = f, obj = -deviance(f) / 2))
    }
    ga <- aug_grid(gr0, lc, lt)
    if (key == "paretologit") {
      f <- fit_pareto_logit(sr, form, gr0 = ga, n_corners = length(bind))
      list(fit = f, obj = -deviance(f) / 2)
    } else {
      f <- fit_pareto_logit_env(sr, form, gr0 = ga, bind = bind,
                                n_corners = length(bind), start = ws$start)
      ws$start <- coef(f)   # warm-start the next profile evaluation
      # value is the mean negative quasi-ll over the grid; larger obj better
      list(fit = f, obj = -f$value * attr(f, "n_grid"))
    }
  }
  neg <- function(l) {
    lc <- l[1]; lt <- if (lt_free) l[2] else 1
    if (lc < BC_BOX_C[1] || lc > BC_BOX_C[2] ||
        lt < BC_BOX_T[1] || lt > BC_BOX_T[2]) return(1e6)
    r <- tryCatch(fit_at(lc, lt), error = function(e) NULL)
    if (is.null(r) || !is.finite(r$obj)) return(1e6)
    -r$obj
  }
  sr <- bc_lambda_search(neg, lambda_start, lt_free)
  lam <- sr$lam; lt_free <- sr$lt_free
  # canonical refit, cold-started, with the same degenerate-refit fallback to
  # lambda_time = 1 as fit_cost_bc
  ws$start <- NULL
  fit <- tryCatch(fit_at(lam[1], lam[2])$fit, error = function(e) e)
  if (inherits(fit, "error") && lam[2] != 1) {
    o <- optimize(function(x) neg(c(x, 1)), interval = BC_BOX_C, tol = 1e-4)
    lam <- c(o$minimum, 1); lt_free <- FALSE
    ws$start <- NULL
    fit <- fit_at(lam[1], lam[2])$fit
  } else if (inherits(fit, "error")) stop(fit)
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
  bs <- bench_levels(data$benchmark)
  one <- function(b) fit_bc(key, data[data$benchmark == b, ],
                            lambda_start = lambda_starts[[b]] %||% c(0, 1))
  cl <- fit_cluster(length(bs))
  fits <- if (is.null(cl)) lapply(bs, one) else parallel::parLapply(cl, bs, one)
  setNames(fits, bs)
}
