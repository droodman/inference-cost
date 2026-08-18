# Deterministic envelope frontier: the lowest-sitting logistic surface that
# passes above every observed run.
#
#   minimise   sum over a FIXED grid of G(z(c,t))        (accuracy units)
#   subject to z_i >= logit(y_i)  for every run i         (the envelope)
#              dz/d ln c >= 0                             (free disposal)
#              dz/dt     >= 0                             (models stay available)
#
# Why this rather than the SFA families:
#
#   * The data enters ONLY through the constraints. Tightness is scored on a
#     fixed grid, so the standard cannot drift as the population of runs
#     changes -- and it did change hugely: gpqa's cheapest tercile went from 66
#     runs with 14% zeros to 638 runs with 56% zeros, which is what dragged the
#     mean-based fits' cheap end the wrong way.
#   * Monotonicity is imposed, not checked afterwards. The frontier cannot slope
#     backwards in cost or time, which is the pathology the lncost:tc
#     interaction produced.
#   * No inefficiency distribution, so nothing hinges on sigma_u, which never
#     identified cleanly in any specification we tried.
#
# The cost of all that: no standard errors, and the fit is pinned by a handful
# of binding runs, so a single lucky run tilts the whole surface. Fit the
# best-per-model variant alongside to see whether the answer rests on flukes.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("frontier_viz.R")

## ---- constraint set ----------------------------------------------------------------

# logit is infinite at 0 and 1, and aime really does contain runs scoring exactly
# 1.000, which would make the problem infeasible. Clip using each run's OWN
# sample size: k/n = n/n is consistent with a true rate near 1 - 1/(2n), not 1.
# Clip ONLY at the top. logit(1) = +Inf makes the problem infeasible, so a run
# scoring n/n is treated as 1 - 1/(2n). But logit(0) = -Inf is harmless: "the
# frontier must exceed 0%" is no constraint at all, and such runs are dropped
# below. Clipping the bottom as well INVENTS a floor -- with n = 40 a zero
# becomes 0.0125, logit -4.37 -- and that spurious floor was what pinned fm13's
# envelope at its cheapest costs and pushed the whole surface far above the data.
#
#   y  numeric VECTOR, one accuracy per run, a proportion in [0, 1]
#   n  numeric VECTOR parallel to y (a scalar recycles), the run's n_samples;
#      NA or less than 1 is treated as 1, so a missing sample size clips to 0.5
#      rather than propagating NA into the constraint matrix
#
# Returns a numeric vector the length of y.
clip_acc <- function(y, n) {
  n <- ifelse(is.na(n) | n < 1, 1, n)
  pmin(y, 1 - 1 / (2 * n))
}

# Only 2-D Pareto-efficient runs can bind. With the frontier non-decreasing in
# both cost and date, a run that is cheaper AND earlier AND scored at least as
# well implies the constraint of any run it dominates, so the dominated ones can
# be dropped -- typically thousands down to a few dozen. This is the same
# efficient set the staircase in pareto_frontier.png traces, in two dimensions.
# Run i is redundant when some j is no dearer, no later, and scored no worse:
# then z_i >= z_j >= L_j >= L_i, so i's constraint is implied. A cost-only sweep
# is NOT enough -- a dearer but earlier run is not dominated -- so this checks
# both coordinates. Exact ties are broken by index so identical runs do not
# eliminate each other.
#
# All three arguments are numeric VECTORS of the same length, element-parallel:
# element i of each describes the same run.
#
#   cost  the cost coordinate. Only its ORDER is used, so any monotone transform
#         serves; fit_envelope() passes log cost, plot_paretologit.R the same.
#   t     the date coordinate, on any numeric scale increasing with time (tc,
#         centred years, in every caller here)
#   L     the level being dominated -- logit of clipped accuracy for the
#         envelope, raw accuracy where the caller wants the P_t(c) definition.
#         Must be finite: -Inf entries would dominate nothing and be dominated by
#         everything, so callers drop them first.
#
# Returns an INTEGER VECTOR of positions into those vectors, ascending: the runs
# whose constraints are not implied by another's. Length is typically a few dozen
# out of thousands.
pareto_binding <- function(cost, t, L) {
  n <- length(cost)
  keep <- logical(n)
  for (i in seq_len(n)) {
    dom    <- cost <= cost[i] & t <= t[i] & L >= L[i]
    strict <- cost < cost[i] | t < t[i] | L > L[i]
    idx <- seq_len(n) < i
    keep[i] <- !any(dom & (strict | idx))
  }
  which(keep)
}

## ---- the objective grid --------------------------------------------------------------

# The fixed grid the envelope scores its objective on, and which the companion
# Pareto-frontier logit below samples its response on. ONE function, used by
# both, so "the two models use the same grid" is true by construction rather
# than by two seq() calls staying in sync.
#
#   data    DATA FRAME with numeric columns lncost and tc; only their RANGES are
#           used, so any superset of columns serves
#   n_cost  SCALAR integer, number of nodes in log cost
#   n_date  SCALAR integer, number of nodes in tc
#
# Returns a DATA FRAME of n_cost * n_date rows, numeric columns lncost and tc:
# the full cross of two uniform sequences over the observed ranges, lncost
# varying fastest.
objective_grid <- function(data, n_cost = 100, n_date = 40) {
  expand.grid(lncost = seq(min(data$lncost), max(data$lncost), length.out = n_cost),
              tc     = seq(min(data$tc),     max(data$tc),     length.out = n_date))
}

## ---- companion model: logit fitted to the Pareto frontier on that grid ----------------

# The empirical Pareto frontier
#
#     P_t(c) = max { acc_i : c_i <= c, t_i <= t }
#
# evaluated at every node of objective_grid(), and a fractional logit fitted to
# those sampled values as if they were data. This replaces fitting the logit to
# the frontier-defining runs themselves: summing the objective over the runs
# weights it by where Pareto points happen to cluster (half of them sit inside
# 21-30% of the cost range), while the envelope's objective is uniform over the
# (log cost, tc) rectangle. Sampling P on the SAME grid gives both models the
# same weighting, so what remains between them is only ABOVE (envelope) versus
# THROUGH (this).
#
# The price is inference: a grid node is not an observation, so the glm's
# standard errors, dispersion and any test read off it are fictions. Callers
# must report the point estimates alone, exactly as they do for the envelope.
#
#   data     DATA FRAME, one row per run, a single benchmark. Needs acc (raw,
#            in [0,1] -- P_t(c) is defined on raw accuracy, so zero and perfect
#            scores participate, unlike the envelope's clipped constraints),
#            lncost and tc, plus anything else `formula` names.
#   formula  two-sided FORMULA with response acc, terms as in fit_envelope()
#   n_cost   SCALAR integer, forwarded to objective_grid()
#   n_date   SCALAR integer, forwarded to objective_grid()
#
# Returns a glm object with class "pareto_grid_logit" prepended, carrying two
# attributes:
#   n_grid     SCALAR integer, grid nodes actually fitted -- nodes cheaper and
#              earlier than every run have P undefined (a max over the empty
#              set) and are dropped, not imputed
#   n_corners  SCALAR integer, undominated runs defining the staircase, i.e.
#              how many real data points the n_grid fitted values resample
fit_pareto_logit <- function(data, formula = acc ~ lncost + tc,
                             n_cost = 100, n_date = 40) {
  gr <- objective_grid(data, n_cost, n_date)
  s  <- data[order(data$lncost), ]
  gr$acc <- NA_real_
  # one staircase per grid date: runs released by then, best accuracy at each
  # cost, looked up at the grid's cost nodes by findInterval (rightmost run with
  # lncost <= node, whose cummax is the running best)
  for (tk in unique(gr$tc)) {
    el <- s[s$tc <= tk, ]
    if (!nrow(el)) next
    rows <- which(gr$tc == tk)
    m <- cummax(el$acc)
    idx <- findInterval(gr$lncost[rows], el$lncost)
    gr$acc[rows] <- ifelse(idx >= 1, m[pmax(idx, 1)], NA_real_)
  }
  gr <- gr[!is.na(gr$acc), , drop = FALSE]
  fit <- glm(formula, data = gr, family = quasibinomial(link = "logit"))
  attr(fit, "n_grid")    <- nrow(gr)
  attr(fit, "n_corners") <- length(pareto_binding(s$lncost, s$tc, s$acc))
  class(fit) <- c("pareto_grid_logit", class(fit))
  fit
}

## ---- the fit -------------------------------------------------------------------------

# Fit the envelope for ONE benchmark. Callers loop over benchmarks themselves;
# nothing here splits by group.
#
#   data     DATA FRAME, one row per run, all of a single benchmark. Must carry
#            the response and every term named in `formula`, plus `lncost` and
#            `tc` (used directly for the Pareto reduction and the objective grid)
#            and `n_samples` (used by clip_acc). Extra columns are ignored.
#   formula  two-sided FORMULA, response ~ terms. Terms may be any subset of
#            lncost, I(lncost^2), tc, I(tc^2) and lncost:tc; the monotonicity
#            block and frontier_coefs_envelope() both key off exactly those
#            names, so a term spelled differently would be silently unconstrained
#            and silently dropped from the plotted surface.
#   n_cost   SCALAR integer, grid points in log cost for the objective
#   n_date   SCALAR integer, grid points in tc. The grid is the objective's
#            weighting only -- the data enters through the constraints.
#   margin   SCALAR numeric, logit units by which the starting surface is lifted
#            clear of the highest constraint, so the solver begins strictly
#            inside the feasible set rather than exactly on its boundary.
#
# Returns an object of class "envelope_frontier": a LIST with
#
#   coefficients    named numeric VECTOR, one per column of model.matrix(formula)
#   value           SCALAR, the objective: mean fitted accuracy over the grid
#   worst_slack     SCALAR, minimum slack across ALL constraints, run and
#                   monotonicity together. Negative beyond -1e-8 means infeasible.
#   n_binding       SCALAR integer, how many run constraints survived the Pareto
#                   reduction (candidates, not necessarily active)
#   slack_envelope  SCALAR, minimum slack over the RUN constraints only, in logit
#                   units; 0 means the surface touches the data
#   slack_mono      SCALAR, minimum slack over the MONOTONICITY constraints, or
#                   NA when the specification has none
#   tightest_row    SCALAR integer, row index into `data` of the closest run
#   env_slack       numeric VECTOR, one entry per candidate, in logit units
#   bind            integer VECTOR of row indices into `data`, parallel to
#                   env_slack, naming the candidates
#   formula         the formula as supplied, so the fit carries its own spec
fit_envelope <- function(data, formula = acc ~ lncost + tc,
                         n_cost = 100, n_date = 40, margin = 0.05) {
  mf <- model.frame(formula, data)
  y  <- model.response(mf)
  X  <- model.matrix(formula, data)
  L  <- qlogis(clip_acc(y, data$n_samples))

  # runs scoring zero impose nothing (logit -Inf), so drop them before the
  # Pareto reduction rather than letting -Inf into the constraint matrix
  pos <- which(is.finite(L))
  bind <- pos[pareto_binding(data$lncost[pos], data$tc[pos], L[pos])]
  Xb <- X[bind, , drop = FALSE]
  Lb <- L[bind]

  # fixed grid: the objective's weighting, independent of where runs cluster;
  # shared with fit_pareto_logit() via objective_grid()
  gr <- objective_grid(data, n_cost, n_date)
  Xg <- model.matrix(delete.response(terms(formula)), gr)

  # monotonicity, evaluated at the corners (both derivatives are linear in beta)
  # Both derivatives are linear in beta, so requiring them non-negative at the
  # ENDS of each range enforces it throughout (a linear function of one variable
  # is non-negative on an interval iff it is at both ends).
  nmv <- colnames(X)
  # e(nm): nm a SCALAR string naming a model-matrix column. Returns a numeric
  # VECTOR of length ncol(X), the indicator for that column, or all zeros if the
  # term is absent -- which is what lets one constraint block cover every
  # specification without branching on which terms exist.
  e <- function(nm) {
    v <- numeric(ncol(X))
    if (nm %in% nmv) v[match(nm, nmv)] <- 1
    v
  }
  ui <- Xb                                            # the envelope itself
  ci <- Lb
  # add(row): row a numeric VECTOR of length ncol(X), a constraint row requiring
  # row %*% beta >= 0. Appends to ui/ci in the enclosing frame; returns nothing
  # useful and is called for that effect.
  add <- function(row) { ui <<- rbind(ui, row); ci <<- c(ci, 0) }

  # free disposal:  dz/d ln c = b_x + 2*b_xx*ln c + b_xt*tc   >= 0
  # stays available: dz/dtc   = b_t + 2*b_tt*tc  + b_xt*ln c  >= 0
  #
  # Each is linear in (ln c, tc) JOINTLY, so requiring it at the four corners of
  # the grid rectangle enforces it throughout. The corners, not the ends of one
  # range: with an lncost:tc term the cost slope depends on the date too, so
  # sweeping cost at a single date leaves the surface free to slope backwards at
  # another. e() returns a zero vector for an absent term, so this one block
  # covers every specification -- linear, cost-quadratic, time-quadratic, full.
  corners <- expand.grid(lc0 = range(gr$lncost), tc0 = range(gr$tc))
  mono <- rbind(
    t(mapply(function(lc0, tc0)
      e("lncost") + 2 * lc0 * e("I(lncost^2)") + tc0 * e("lncost:tc"),
      corners$lc0, corners$tc0)),
    t(mapply(function(lc0, tc0)
      e("tc") + 2 * tc0 * e("I(tc^2)") + lc0 * e("lncost:tc"),
      corners$lc0, corners$tc0)))
  # Without the cross terms all four corners collapse to the same row; drop the
  # duplicates rather than hand the solver the same constraint four times.
  mono <- unique(mono)
  for (i in seq_len(nrow(mono))) add(mono[i, ])

  # mean fitted accuracy over the grid -- the surface's average height.
  # fn(b):    b a numeric VECTOR of length ncol(X). Returns a SCALAR.
  # gr_fn(b): the same argument; returns the gradient, a numeric VECTOR of
  #           length ncol(X). Supplied analytically because SLSQP asks for it on
  #           every iteration and a finite-difference gradient over ~200 grid
  #           points would dominate the run time.
  fn <- function(b) mean(plogis(drop(Xg %*% b)))
  gr_fn <- function(b) {
    p <- plogis(drop(Xg %*% b))
    drop(crossprod(Xg, p * (1 - p))) / nrow(Xg)
  }

  # Feasible start: a quasibinomial fit with curvature zeroed and linear slopes
  # forced positive (so monotonicity holds whatever glm returned), then lifted
  # until it clears every envelope constraint.
  b0 <- coef(glm(formula, data = data, family = quasibinomial(link = "logit")))
  # Zero the interaction too, not just the squares: it also enters both
  # derivatives, so leaving glm's value in place can hand SLSQP a start that
  # violates monotonicity, and an infeasible start is how "no feasible envelope
  # found" happens on a problem that is perfectly feasible.
  for (nm in c("I(lncost^2)", "I(tc^2)", "lncost:tc"))
    if (nm %in% names(b0)) b0[[nm]] <- 0
  b0[["lncost"]] <- max(b0[["lncost"]], 0.05)
  b0[["tc"]]     <- max(b0[["tc"]], 0.05)
  b0[1] <- b0[1] + max(0, max(Lb - drop(Xb %*% b0))) + margin
  starts <- list(b0)

  # A richer model nests a simpler one, so it can never legitimately do worse --
  # yet both constrOptim and SLSQP stalled well inside the feasible set when the
  # near-collinear (ln c)^2 column was present, returning objectives ABOVE the
  # linear fit. Starting the richer model from the simpler model's solution makes
  # the nesting property hold by construction, whatever the solver does next.
  curv <- grep("^I\\(", colnames(X), value = TRUE)
  if (length(curv)) {
    tl <- attr(terms(formula), "term.labels")
    simpler <- reformulate(tl[!grepl("^I\\(", tl)],
                           response = all.vars(formula)[1])
    inner <- fit_envelope(data, simpler, n_cost, n_date, margin)
    bs <- setNames(numeric(ncol(X)), colnames(X))
    bs[names(coef(inner))] <- coef(inner)
    starts <- c(starts, list(bs))
  }

  # SLSQP, not constrOptim. constrOptim uses a logarithmic barrier that blows up
  # near the boundary; with ~100-200 constraint rows it dominated the objective
  # and the optimiser halted well INSIDE the feasible set -- slack of +0.07 to
  # +0.17 where a tight envelope must sit at 0, and objectives worse than the
  # nested linear fit, which is impossible for a correct solver. SLSQP handles
  # linear inequality constraints directly and can terminate exactly on them.
  # run(x0): x0 a numeric VECTOR of length ncol(X), the starting coefficients.
  # Returns a LIST of b (named numeric vector, the solution), obj (SCALAR
  # objective at b) and slack (SCALAR, minimum over all constraint rows).
  run <- function(x0) {
    r <- nloptr::nloptr(
      x0 = x0, eval_f = function(b) list(objective = fn(b), gradient = gr_fn(b)),
      eval_g_ineq = function(b) list(constraints = as.vector(ci - ui %*% b),
                                     jacobian = -ui),
      opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-10,
                  maxeval = 5000, print_level = 0))
    b <- setNames(r$solution, colnames(X))
    list(b = b, obj = fn(b), slack = min(drop(ui %*% b) - ci))
  }

  # keep the best FEASIBLE candidate; a lower objective reached by leaving the
  # feasible set is not an envelope
  cand <- c(lapply(starts, run),
            lapply(starts, function(x) list(b = x, obj = fn(x),
                                            slack = min(drop(ui %*% x) - ci))))
  ok <- vapply(cand, function(z) z$slack >= -1e-8, logical(1))
  if (!any(ok)) stop("no feasible envelope found")
  best <- cand[ok][[which.min(vapply(cand[ok], `[[`, numeric(1), "obj"))]]

  # worst_slack pools two quite different constraints, and which of them is tight
  # is the whole story about whether the surface touches the data. Report them
  # apart:
  #   slack_envelope  how far above the NEAREST RUN the surface sits, in logit
  #                   units. Zero means it touches; large means it floats.
  #   slack_mono      how much room is left in the monotonicity constraints.
  # If slack_mono is 0 while slack_envelope is not, monotonicity is what is
  # holding the surface up and the data is not binding anywhere -- the fit is
  # then the lowest MONOTONE logistic surface above the runs, which is a strictly
  # stronger requirement than the lowest logistic surface above them, and it can
  # sit well clear of the Pareto staircase.
  env_slack <- drop(Xb %*% best$b) - Lb
  structure(list(coefficients = best$b, value = best$obj,
                 worst_slack = best$slack, n_binding = length(bind),
                 slack_envelope = min(env_slack),
                 slack_mono = if (nrow(mono)) min(drop(mono %*% best$b)) else NA_real_,
                 tightest_row = bind[which.min(env_slack)],
                 env_slack = env_slack, bind = bind,
                 formula = formula),
            class = "envelope_frontier")
}

# coef() method, so the fit answers to the same accessor as glm and maxLik fits.
#
#   object  an "envelope_frontier" as returned by fit_envelope()
#   ...     ignored; present only to match the generic's signature
#
# Returns a named numeric VECTOR, one entry per model-matrix column.
coef.envelope_frontier <- function(object, ...) object$coefficients

# so frontier_coefs()/frontier_index() in frontier_viz.R work unchanged
#
#   fit  an "envelope_frontier" as returned by fit_envelope()
#
# Returns a named numeric VECTOR of length 6 in the FIXED order b0, bx, bt, btt,
# bxt, bxx, whatever subset of terms the fit actually has -- a term the formula
# omitted comes back as 0, so linear, cost-quadratic and full-quadratic fits all
# evaluate through one prediction path. The order and names must match
# frontier_coefs() in frontier_viz.R, which is what frontier_index() consumes.
frontier_coefs_envelope <- function(fit) {
  cf <- coef(fit)
  # get1(nm): nm a SCALAR string naming a coefficient. Returns a SCALAR -- the
  # estimate, or 0 when the fit does not carry that term.
  get1 <- function(nm) if (nm %in% names(cf)) unname(cf[[nm]]) else 0
  c(b0 = get1("(Intercept)"), bx = get1("lncost"), bt = get1("tc"),
    btt = get1("I(tc^2)"), bxt = get1("lncost:tc"),
    bxx = get1("I(lncost^2)"))
}
