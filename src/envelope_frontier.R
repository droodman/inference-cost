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
pareto_binding <- function(cost, t, L) {
  n <- length(cost)
  keep <- logical(n)
  for (i in seq_len(n)) {
    dom <- cost <= cost[i] & t <= t[i] & L >= L[i]
    strict <- cost < cost[i] | t < t[i] | L > L[i]
    idx <- seq_len(n) < i
    keep[i] <- !any(dom & (strict | idx))
  }
  which(keep)
}

## ---- the fit -------------------------------------------------------------------------

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

  # fixed grid: the objective's weighting, independent of where runs cluster
  gcost <- exp(seq(min(data$lncost), max(data$lncost), length.out = n_cost))
  gtc   <- seq(min(data$tc), max(data$tc), length.out = n_date)
  gr <- expand.grid(lncost = log(gcost), tc = gtc)
  Xg <- model.matrix(delete.response(terms(formula)), gr)

  # monotonicity, evaluated at the corners (both derivatives are linear in beta)
  # Both derivatives are linear in beta, so requiring them non-negative at the
  # ENDS of each range enforces it throughout (a linear function of one variable
  # is non-negative on an interval iff it is at both ends).
  nmv <- colnames(X)
  e <- function(nm) {
    v <- numeric(ncol(X))
    if (nm %in% nmv) v[match(nm, nmv)] <- 1
    v
  }
  ui <- Xb                                            # the envelope itself
  ci <- Lb
  add <- function(row) { ui <<- rbind(ui, row); ci <<- c(ci, 0) }

  # free disposal: dz/d ln c = b_x + 2*b_xx*ln c >= 0
  if ("I(lncost^2)" %in% nmv) {
    for (lc0 in range(log(gcost))) {
      add(e("lncost") + 2 * lc0 * e("I(lncost^2)"))
    }
  } else {
    add(e("lncost"))
  }
  # frontier non-decreasing in time: dz/dtc = b_t + 2*b_tt*tc >= 0
  if ("I(tc^2)" %in% nmv) {
    for (tc0 in range(gtc)) add(e("tc") + 2 * tc0 * e("I(tc^2)"))
  } else {
    add(e("tc"))
  }

  # mean fitted accuracy over the grid -- the surface's average height
  fn <- function(b) mean(plogis(drop(Xg %*% b)))
  gr_fn <- function(b) {
    p <- plogis(drop(Xg %*% b))
    drop(crossprod(Xg, p * (1 - p))) / nrow(Xg)
  }

  # Feasible start: a quasibinomial fit with curvature zeroed and linear slopes
  # forced positive (so monotonicity holds whatever glm returned), then lifted
  # until it clears every envelope constraint.
  b0 <- coef(glm(formula, data = data, family = quasibinomial(link = "logit")))
  for (nm in c("I(lncost^2)", "I(tc^2)")) if (nm %in% names(b0)) b0[[nm]] <- 0
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

  structure(list(coefficients = best$b, value = best$obj,
                 worst_slack = best$slack, n_binding = length(bind),
                 formula = formula),
            class = "envelope_frontier")
}

coef.envelope_frontier <- function(object, ...) object$coefficients

# so frontier_coefs()/frontier_index() in frontier_viz.R work unchanged
frontier_coefs_envelope <- function(fit) {
  cf <- coef(fit)
  get1 <- function(nm) if (nm %in% names(cf)) unname(cf[[nm]]) else 0
  c(b0 = get1("(Intercept)"), bx = get1("lncost"), bt = get1("tc"),
    btt = get1("I(tc^2)"), bxt = get1("lncost:tc"),
    bxx = get1("I(lncost^2)"))
}
