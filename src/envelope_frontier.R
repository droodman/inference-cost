# Frontier-per-se fits in the accuracy direction: the empirical Pareto
# staircase P_t(c) sampled on a fixed (log cost, date) grid, model S's
# fractional logit fitted THROUGH those samples (fit_pareto_logit), and the
# same objective minimised subject to envelope constraints
# (fit_pareto_logit_env):
#
#   z_i >= logit(y_i)  for every run i     the surface must clear everything
#                                          ever observed
#   dz/d ln c >= 0                         free disposal
#   dz/dt     >= 0                         models stay available
#
# Why constraints rather than the SFA families: monotonicity is imposed, not
# checked afterwards -- the frontier cannot slope backwards in cost or time,
# the pathology the lncost:tc interaction produced -- and there is no
# inefficiency distribution, so nothing hinges on sigma_u, which never
# identified cleanly in any specification we tried.
#
# The constrained fit replaced an earlier strict envelope (the lowest-sitting
# surface above every run, its objective never looking at the staircase).
# The two agreed closely wherever the constraints bound, but the envelope was
# pinned by a handful of binding runs, so a single lucky run tilted the whole
# surface; under the deviance objective every Pareto-efficient run pulls on
# the fit in proportion to the grid area its staircase tread covers.
#
# The cost that remains: no standard errors -- grid nodes are not
# observations, and the constrained optimum has no likelihood behind it.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("frontier_viz.R")

## ---- constraint set ----------------------------------------------------------------

# logit is infinite at 0 and 1, and aime really does contain runs scoring exactly
# 1.000, which would make the problem infeasible. Clip using each run's OWN
# sample size: k/n = n/n is consistent with a true rate near 1 - 1/(2n), not 1.
# Clip ONLY at the top. logit(1) = +Inf makes the problem infeasible, so a run
# scoring n/n is treated as 1 - 1/(2n). But logit(0) = -Inf is harmless: "the
# frontier must exceed 0%" is no constraint at all, and such runs are dropped
# below. Clipping the bottom as well INVENTS a floor -- with n = 100 a zero
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
#         serves; every caller here passes log cost.
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
objective_grid <- function(data, n_cost = 100, n_date = 100) {
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
#   formula  two-sided FORMULA with response acc, terms as in
#            fit_pareto_logit_env()
#   n_cost   SCALAR integer, forwarded to objective_grid()
#   n_date   SCALAR integer, forwarded to objective_grid()
#   grid_augment  optional FUNCTION applied to the grid data frame after
#            construction, adding columns derived from its (lncost, tc)
#            coordinates -- how the Box-Cox terms (boxcox_frontier.R) become
#            evaluable at grid nodes. The staircase itself is untouched.
#
# Returns a glm object with class "pareto_grid_logit" prepended, carrying two
# attributes:
#   n_grid     SCALAR integer, grid nodes actually fitted -- nodes cheaper and
#              earlier than every run have P undefined (a max over the empty
#              set) and are dropped, not imputed
#   n_corners  SCALAR integer, undominated runs defining the staircase, i.e.
#              how many real data points the n_grid fitted values resample
# The staircase resampled on the objective grid: P_t(c) at every node, nodes
# where P is undefined dropped. Split out of fit_pareto_logit() so the Box-Cox
# lambda search (boxcox_frontier.R) can compute it ONCE per benchmark -- the
# staircase does not depend on the lambdas, only the regressor columns do.
pareto_grid_response <- function(data, n_cost = 100, n_date = 100) {
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
  gr[!is.na(gr$acc), , drop = FALSE]
}

# The node set the nonparametric rate check below and its smoothed twin
# (surface_decline_qtr, cost_frontier.R) SHARE, built once here so the two
# cannot drift apart: at each of n_date dates with a full dt-year horizon
# still ahead of it, n_level accuracy levels placed uniformly WITHIN the
# state of the art at that date --
#
#     a_j(t) = SOTA(t) * (j - 1/2) / n_level,   j = 1..n_level
#
# the midpoint lattice on (0, SOTA(t)]. Shrinking the lattice to the achieved
# range, rather than clipping a fixed one, keeps every date's node count
# equal and makes every node well-defined by construction: no undefined
# frontiers, no a = 0 rows (the midpoint offset keeps the lattice off both
# endpoints, so C_a exists at the base date and, a fortiori, dt later), and
# no exact-1 levels the fitted surfaces cannot be asked about. Dates before
# a benchmark's first positive score contribute nothing.
#
# Returns a DATA FRAME of columns tc and a, or NULL when no date has a full
# horizon.
decline_nodes <- function(data, n_level = 100, n_date = 100, dt = 1) {
  tgrid <- seq(min(data$tc), max(data$tc), length.out = n_date)
  tgrid <- tgrid[tgrid + dt <= max(data$tc)]
  if (!length(tgrid)) return(NULL)
  do.call(rbind, lapply(tgrid, function(tk) {
    sota <- max(data$acc[data$tc <= tk])
    if (sota <= 0) return(NULL)
    data.frame(tc = tk, a = sota * (seq_len(n_level) - 0.5) / n_level)
  }))
}

# Nonparametric sense check on the fitted models' "cost drop, %/qtr" summary
# (regression_tables.R): the same quantity read straight off the Pareto
# staircase, with no functional form anywhere.
#
# At each decline_nodes() node (a, t), ask what a costs dt years on:
# C_a(t + dt) = min{c_i : t_i <= t + dt, acc_i >= a}, the record cost
# iso_pareto_steps() (frontier_viz.R) traces, against the record C_a(t) now.
# Record at BOTH ends, never the cost of a run that merely achieves a: the
# staircase is flat between jumps, so a run's own cost typically overpays for
# a, and measuring from it would count that slack as decline even over a
# stretch where no new model appeared. Each node's change is the horizontal
# shift of the staircase at that node's level -- exactly what -b_t/b_x
# measures on a fitted surface.
#
# The weighting is the node set's: levels uniform in ACCURACY within each
# date's achieved range, dates uniform. An earlier version scored the check
# on the fits' own (log cost, date) objective grid, which weighted each level
# by the log-cost WIDTH of its staircase step -- concentrating the average on
# long-plateau levels pinned by single lucky cheap runs, and (it turned out)
# on exactly the regions where fitted surfaces' misfit drifts over the
# horizon. The price of the change: a fixed node index j is a RISING absolute
# level as the SOTA climbs, so the aggregate averages fixed-level changes
# under a date-varying level measure; each node's dln is still a fixed-level
# quantity.
#
# Measured over a YEAR, expressed per QUARTER: at a fixed level the record is
# a step function, and over one quarter it usually has not moved, so the
# average is a few real drops diluted by zeros. A year of change dilutes
# less; the geometric mean log change is then compounded down (divide by 4dt)
# into the quarterly rate the tables print. A benchmark observed for less
# than dt has no measurable horizon and returns NULL -- an honest gap until
# its data span grows.
#
#   data     DATA FRAME, one benchmark's runs; needs acc, lncost, tc
#   n_level  SCALAR integer, accuracy levels per date (decline_nodes)
#   n_date   SCALAR integer, dates on the lattice (decline_nodes)
#   dt       SCALAR, the measurement horizon in years of tc
#
# Returns NULL if no node survives, else a LIST:
#   pct_qtr      SCALAR, 100*(1 - exp(mean(dln) / 4dt)): the geometric-mean
#                decline as a compound quarterly rate, comparable to the
#                tables' 100*(1 - exp(-b_t/b_x/4))
#   share_moved  SCALAR in [0, 1], fraction of nodes whose record moved at all
#                over the horizon
#   n_nodes      SCALAR integer, nodes entering the average
pareto_decline_qtr <- function(data, n_level = 100, n_date = 100, dt = 1) {
  nd <- decline_nodes(data, n_level, n_date, dt)
  if (is.null(nd) || !nrow(nd)) return(NULL)
  s <- data[order(data$lncost), ]
  # log record cost of each level at date t: cheapest run no later than t
  # scoring at least the level. Over cost-sorted runs cummax(acc) is the
  # running best, so the first index where it reaches the level is the record.
  # Every level handed in is within the SOTA at its base date, so a match
  # exists at both dates and which() cannot come back empty.
  rec <- function(t, levels) {
    el <- s[s$tc <= t, ]
    m <- cummax(el$acc)
    vapply(levels, function(a) el$lncost[which(m >= a)[1]], numeric(1))
  }
  dln <- unlist(lapply(unique(nd$tc), function(tk) {
    a <- nd$a[nd$tc == tk]
    rec(tk + dt, a) - rec(tk, a)
  }))
  list(pct_qtr     = 100 * (1 - exp(mean(dln) / (4 * dt))),
       share_moved = mean(dln < 0),
       n_nodes     = length(dln))
}

# gr0 and n_corners let a caller that fits MANY times on one benchmark -- the
# Box-Cox lambda profile -- pass the staircase and corner count in once
# instead of recomputing them per fit: neither depends on the lambdas (the
# staircase is raw accuracy; the Pareto set is invariant to any monotone
# transform), and together they were a measured third of the profile's time.
# Left NULL (every one-shot caller), behavior is unchanged.
fit_pareto_logit <- function(data, formula = acc ~ lncost + tc,
                             n_cost = 100, n_date = 100, grid_augment = NULL,
                             gr0 = NULL, n_corners = NULL) {
  gr <- if (is.null(gr0)) pareto_grid_response(data, n_cost, n_date) else gr0
  if (!is.null(grid_augment)) gr <- grid_augment(gr)
  fit <- glm(formula, data = gr, family = quasibinomial(link = "logit"))
  attr(fit, "n_grid")    <- nrow(gr)
  attr(fit, "n_corners") <- if (is.null(n_corners)) {
    s <- data[order(data$lncost), ]
    length(pareto_binding(s$lncost, s$tc, s$acc))
  } else n_corners
  class(fit) <- c("pareto_grid_logit", class(fit))
  fit
}

## ---- the fit -------------------------------------------------------------------------

# The constraint system of fit_pareto_logit_env(), stacked into one
# ui %*% beta >= ci:  every surviving run's level from below, then the
# monotonicity rows.
#
# Runs scoring zero impose nothing (logit -Inf), so they are dropped before the
# Pareto reduction rather than letting -Inf into the constraint matrix. The
# explicit y > 0 matters when lambda_odds > 0: phi(0 odds) is then FINITE
# (-1/lambda), and keeping it would smuggle the zero runs back in as
# constraints the logit version deliberately excludes.
#
#   data     DATA FRAME of runs, as in fit_pareto_logit_env()
#   formula  as in fit_pareto_logit_env(); the monotonicity block keys off its
#            term names
#   gr       DATA FRAME of grid nodes carrying the formula's coordinate columns.
#            Only the RANGES of its coordinates are used here, to place the
#            corner rows that make monotonicity hold over the whole rectangle.
#   lambda_odds  as in fit_pareto_logit_env()
#
# Returns a LIST:
#   ui, ci  the stacked system: ui %*% beta >= ci
#   Xb, Lb  the run rows alone -- model-matrix rows and clipped levels
#   bind    integer VECTOR of row indices into `data`, the runs behind Xb
#   mono    the monotonicity rows alone (0 rows when the spec has none)
envelope_constraints <- function(data, formula, gr, lambda_odds = 0,
                                 bind = NULL) {
  mf <- model.frame(formula, data)
  y  <- model.response(mf)
  X  <- model.matrix(formula, data)
  pc <- clip_acc(y, data$n_samples)
  L  <- if (abs(lambda_odds) < 1e-8) qlogis(pc) else
    bc_tf(pc / (1 - pc), lambda_odds)
  # `bind` may be passed in by the Box-Cox profile, which calls this hundreds
  # of times per benchmark: the O(n^2) Pareto reduction depends on L only
  # through order comparisons, and every lambda_odds transforms L
  # monotonically, so the binding set never moves across the profile
  if (is.null(bind)) {
    pos <- which(is.finite(L) & y > 0)
    bind <- pos[pareto_binding(data$lncost[pos], data$tc[pos], L[pos])]
  }
  Xb <- X[bind, , drop = FALSE]
  Lb <- L[bind]

  # monotonicity, evaluated at the corners: both derivatives are linear in beta,
  # so requiring them non-negative at the ENDS of each range enforces it
  # throughout (a linear function of one variable is non-negative on an interval
  # iff it is at both ends).
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

  # free disposal:  dz/d ln c = b_x + 2*b_xx*ln c + b_xt*tc   >= 0
  # stays available: dz/dtc   = b_t + 2*b_tt*tc  + b_xt*ln c  >= 0
  #
  # Each is linear in (ln c, tc) JOINTLY, so requiring it at the four corners of
  # the grid rectangle enforces it throughout. The corners, not the ends of one
  # range: with an lncost:tc term the cost slope depends on the date too, so
  # sweeping cost at a single date leaves the surface free to slope backwards at
  # another. e() returns a zero vector for an absent term, so this one block
  # covers every specification -- linear, cost-quadratic, time-quadratic, full.
  if ("phic" %in% nmv) {
    # Box-Cox terms (boxcox_frontier.R): z = b0 + bx*phic + bt*phit +
    # bxt*phic*phit, and phi is strictly increasing in its argument whatever
    # lambda is, so monotonicity in cost and date IS monotonicity in phic and
    # phit:
    #   dz/dphic = b_phic + b_phixt*phit >= 0
    #   dz/dphit = b_phit + b_phixt*phic >= 0
    # Each is linear in ONE other coordinate, so the ends of that coordinate's
    # grid range enforce it throughout.
    mono <- rbind(
      t(vapply(range(gr$phit), function(p) e("phic") + p * e("phixt"),
               numeric(ncol(X)))),
      t(vapply(range(gr$phic), function(p) e("phit") + p * e("phixt"),
               numeric(ncol(X)))))
  } else {
    corners <- expand.grid(lc0 = range(gr$lncost), tc0 = range(gr$tc))
    mono <- rbind(
      t(mapply(function(lc0, tc0)
        e("lncost") + 2 * lc0 * e("I(lncost^2)") + tc0 * e("lncost:tc"),
        corners$lc0, corners$tc0)),
      t(mapply(function(lc0, tc0)
        e("tc") + 2 * tc0 * e("I(tc^2)") + lc0 * e("lncost:tc"),
        corners$lc0, corners$tc0)))
  }
  # Without the cross terms all four corners collapse to the same row; drop the
  # duplicates rather than hand the solver the same constraint four times.
  mono <- unique(mono)
  list(ui = rbind(Xb, mono), ci = c(Lb, numeric(nrow(mono))),
       Xb = Xb, Lb = Lb, bind = bind, mono = mono)
}

# A starting point strictly inside the feasible set, made from any glm's
# coefficients. Zero the curvature and interaction terms -- they enter both
# monotonicity derivatives, so glm's values can hand the solver an infeasible
# start, and an infeasible start is how "no feasible fit found" happens on a
# problem that is perfectly feasible. Force the linear slopes positive (so
# monotonicity holds whatever glm returned), then lift the intercept until
# the surface clears every run constraint with `margin` to spare.
#
#   b0      named numeric VECTOR, glm coefficients on `formula`'s model matrix
#   Xb, Lb  the run constraints, as envelope_constraints() returns them
#   margin  SCALAR, logit units of clearance above the tightest run
feasible_start <- function(b0, Xb, Lb, margin) {
  b0[is.na(b0)] <- 0   # an aliased column would otherwise poison the lift
  for (nm in c("I(lncost^2)", "I(tc^2)", "lncost:tc", "phixt"))
    if (nm %in% names(b0)) b0[[nm]] <- 0
  for (nm in c("lncost", "tc", "phic", "phit"))
    if (nm %in% names(b0)) b0[[nm]] <- max(b0[[nm]], 0.05)
  b0[1] <- b0[1] + max(0, max(Lb - drop(Xb %*% b0))) + margin
  b0
}

# coef() method, so the fit answers to the same accessor as glm and maxLik fits.
#
#   object  an "envelope_frontier" as returned by fit_pareto_logit_env()
#   ...     ignored; present only to match the generic's signature
#
# Returns a named numeric VECTOR, one entry per model-matrix column.
coef.envelope_frontier <- function(object, ...) object$coefficients

## ---- the fit: the companion's loss under the envelope's constraints -------------------

# fit_pareto_logit()'s estimator inside the envelope's feasible set: the
# quasibinomial deviance against the staircase P_t(c) sampled on the grid,
# minimised subject to the surface lying at or above every run and to
# monotonicity in cost and date. Every Pareto-efficient run pulls on the
# surface in proportion to the grid area its staircase tread covers, at the
# price that the fit depends on the staircase's interior, which sits below
# the true frontier wherever the best available run underperformed it.
#
# No multi-start machinery, deliberately: at lambda_odds = 0 the deviance is
# convex in beta (canonical logit link) and the constraints are linear, so one
# SLSQP run from a feasible start finds the global minimum. (The strict
# envelope this fit replaced needed candidate-racing -- its mean-of-sigmoids
# objective was not convex and its solver stalled more than once. With
# lambda_odds != 0 convexity is not guaranteed here either, but the
# unconstrained analogue fit_pareto_bclink runs single-start BFGS on the same
# smooth surface, and the same trust extends here.) SLSQP rather than
# constrOptim throughout this file's history: constrOptim's logarithmic
# barrier blew up near the boundary and halted well inside the feasible set,
# where SLSQP handles linear inequality constraints directly and can
# terminate exactly on them.
#
# Zeros and ones: the staircase response keeps both, as in fit_pareto_logit()
# -- bc_qll() clamps mu away from 0/1, so P = 0 and P = 1 terms stay finite.
# The constraints exclude zeros (logit -Inf is unusable, and "the frontier
# must exceed 0%" is no constraint) and clip ones by the run's own sample
# size via clip_acc().
#
#   data     DATA FRAME, one row per run, all of a single benchmark. Must carry
#            the response and every term named in `formula`, plus `lncost` and
#            `tc` (used directly for the Pareto reduction and the grid) and
#            `n_samples` (used by clip_acc). Extra columns are ignored.
#   formula  two-sided FORMULA, response ~ terms. Terms may be any subset of
#            lncost, I(lncost^2), tc, I(tc^2) and lncost:tc, OR the Box-Cox trio
#            phic, phit, phixt (boxcox_frontier.R); the monotonicity block and
#            frontier_coefs() both key off exactly those names, so a term
#            spelled differently would be silently unconstrained and silently
#            dropped from the plotted surface.
#   n_cost   SCALAR integer, grid points in log cost
#   n_date   SCALAR integer, grid points in tc
#   margin   SCALAR numeric, logit units by which the starting surface is
#            lifted clear of the highest constraint, so the solver begins
#            strictly inside the feasible set rather than exactly on its
#            boundary.
#   grid_augment  optional FUNCTION applied to the grid data frame after
#            construction, adding columns derived from its (lncost, tc)
#            coordinates so terms like the Box-Cox trio are evaluable at grid
#            nodes. The run data must already carry the same columns.
#   lambda_odds  SCALAR, the response-side Box-Cox parameter of the doubly-
#            transformed family (fit_bc): the surface's index is phi(odds;
#            lambda_odds) rather than the logit, which is the special case
#            lambda_odds = 0 (the default -- every non-BC caller). The run
#            constraints become surface >= phi(odds_i) -- the same set of
#            accuracy surfaces, reparameterized -- and the objective stays on
#            the PROBABILITY scale via the inverse link (bc_mu), so it is
#            comparable across lambda_odds values by construction -- what
#            makes the Box-Cox profile well-posed.
#
# Returns an object of class c("pareto_logit_env", "envelope_frontier"),
# carrying fit_pareto_logit()'s n_grid / n_corners attributes: a LIST with
#
#   coefficients    named numeric VECTOR, one per column of model.matrix(formula)
#   value           SCALAR, the objective: mean deviance per grid node
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
#
# slack_envelope and slack_mono split worst_slack because which of the two is
# tight is the whole story about what shapes the fit where it leaves the
# staircase: slack_envelope = 0 means the surface touches a run;
# slack_mono = 0 means the unconstrained fit wanted to slope backwards and
# monotonicity is what holds it flat, so the fitted rate there is the
# constraint's, not the data's.
# gr0, bind and n_corners are the lambda-invariant precomputations a Box-Cox
# profile passes in (see fit_pareto_logit); `start` warm-starts SLSQP from the
# previous profile iteration's solution -- consecutive evaluations sit at
# nearby lambdas, the same argument the SFA profiles' warm starts make. A
# warm-started solve that ends infeasible falls back to the cold start rather
# than failing.
fit_pareto_logit_env <- function(data, formula = acc ~ lncost + tc,
                                 n_cost = 100, n_date = 100, margin = 0.05,
                                 grid_augment = NULL, lambda_odds = 0,
                                 gr0 = NULL, bind = NULL, n_corners = NULL,
                                 start = NULL) {
  gr <- if (is.null(gr0)) pareto_grid_response(data, n_cost, n_date) else gr0
  if (!is.null(grid_augment)) gr <- grid_augment(gr)
  Xg <- model.matrix(delete.response(terms(formula)), gr)
  P  <- gr$acc
  cs <- envelope_constraints(data, formula, gr, lambda_odds, bind = bind)

  # negative mean Bernoulli quasi-log-likelihood of the staircase values, on
  # the probability scale whatever lambda_odds is. The score has the closed
  # form (P - mu)/(1 + lo*eta) per node (the canonical (P - mu) at lo = 0),
  # exactly as in fit_pareto_bclink (boxcox_frontier.R); 0 outside the
  # clamped link's domain.
  fn <- function(b) -bc_qll(P, bc_mu(drop(Xg %*% b), lambda_odds)) / nrow(Xg)
  gr_fn <- function(b) {
    eta <- drop(Xg %*% b)
    resid <- P - bc_mu(eta, lambda_odds)
    if (abs(lambda_odds) >= 1e-8) {
      base <- 1 + lambda_odds * eta
      w <- numeric(length(eta))
      i <- which(base > 0)
      w[i] <- 1 / base[i]
      resid <- resid * w
    }
    -drop(crossprod(Xg, resid)) / nrow(Xg)
  }

  solve1 <- function(x0) nloptr::nloptr(
    x0 = x0, eval_f = function(b) list(objective = fn(b), gradient = gr_fn(b)),
    eval_g_ineq = function(b) list(constraints = as.vector(cs$ci - cs$ui %*% b),
                                   jacobian = -cs$ui),
    opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-10,
                maxeval = 5000, print_level = 0))
  r <- NULL
  if (!is.null(start) && length(start) == ncol(Xg)) {
    r <- solve1(unname(start))
    if (min(drop(cs$ui %*% r$solution) - cs$ci) < -1e-8) r <- NULL
  }
  if (is.null(r)) {
    b0 <- feasible_start(coef(glm(formula, data = gr,
                                  family = quasibinomial(link = "logit"))),
                         cs$Xb, cs$Lb, margin)
    r <- solve1(b0)
  }
  b <- setNames(r$solution, colnames(Xg))
  env_slack <- drop(cs$Xb %*% b) - cs$Lb
  fit <- structure(
    list(coefficients = b, value = fn(b),
         worst_slack = min(drop(cs$ui %*% b) - cs$ci),
         n_binding = length(cs$bind),
         slack_envelope = min(env_slack),
         slack_mono = if (nrow(cs$mono)) min(drop(cs$mono %*% b)) else NA_real_,
         tightest_row = cs$bind[which.min(env_slack)],
         env_slack = env_slack, bind = cs$bind, formula = formula),
    class = c("pareto_logit_env", "envelope_frontier"))
  if (fit$worst_slack < -1e-8) stop("no feasible constrained Pareto logit found")
  attr(fit, "n_grid")    <- nrow(gr)
  attr(fit, "n_corners") <- if (is.null(n_corners))
    length(pareto_binding(data$lncost, data$tc, data$acc)) else n_corners
  fit
}
