# Regression tables: one per statistical model PER SPECIFICATION (linear,
# full quadratic, Box-Cox), in RTF and HTML -- regression_<key>_<spec>.*.
#
# Columns are one per benchmark, estimates with standard errors in
# parentheses beneath, plus a final Pooled pair combining the benchmarks on
# the ECI capability scale in the linear and quadratic tables (see the
# pooling section below; the Box-Cox specification is excluded from pooling
# -- its per-benchmark transforms leave no common scale).
# Model order and naming follow the plot viewer, which is the version of these
# names that has actually been read by a human:
#
#   S               Logistic, all tests
#   A               Stochastic frontier
#   B               Stochastic frontier (time-dependent inefficiency spread)
#   paretologit     Logistic, Pareto points
#   paretologitenv  Logistic, Pareto points, envelope-constrained
#                   (fit_pareto_logit_env: the frontier logit's objective under
#                   the envelope's constraints)
#
# Three things differ across models and the tables must not paper over them:
#
#   * Standard errors. A and B are maxLik fits with robust (sandwich) errors --
#     the Papke-Wooldridge point that the objective is a QUASI-likelihood, so
#     inverse-Hessian errors are wrong. S is a quasibinomial glm, given HC1
#     errors for the same reason. The envelope and the Pareto-frontier logit
#     have NO standard errors at all: the envelope is a constrained optimisation
#     with no likelihood behind it, and the frontier logit is fitted to P_t(c)
#     sampled on a fixed grid, whose nodes are not observations. Inventing
#     errors for either would be worse than leaving the cells empty.
#
#   * The quadratic-vs-linear test. A and B have real likelihoods, so this is a
#     likelihood ratio test on 3 df (quad adds lncost^2, tc^2 and lncost:tc).
#     S's quasibinomial glm reports no logLik, so it gets the Wald analogue on
#     the same three terms, LABELLED as Wald rather than passed off as LR. The
#     envelope and the frontier logit get neither, for the same reason they get
#     no standard errors.
#
#   * Sample. The Pareto-frontier logit's N is the number of grid nodes at
#     which P_t(c) is defined -- nodes cheaper and earlier than every run have
#     no frontier value. Those nodes resample a few dozen staircase corners, so
#     N there measures the sampling resolution, not independent information.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("fit_store.R")   # every fit comes from the shared store

d <- load_runs()
benches <- bench_levels(d$benchmark)

## ---- benchmark discriminations for pooling --------------------------------------
#
# ALPHA (prepare_data.R) holds each benchmark's 2PL discrimination in logits
# per ECI point; it is defined beside TASK_LABEL, its join key, because
# cost_frontier_report.R pools with it too.
#
# Only the PRIMARY benchmarks are pooled, so only they must carry an ECI
# discrimination; a new secondary benchmark missing from the ECI file is
# reported unpooled rather than halting the build.
stopifnot(!anyNA(ALPHA[intersect(PRIMARY_BENCHES, benches)]))

## ---- row layout ---------------------------------------------------------------

# term -> (plain label for RTF, HTML label). Order is the table's row order; a
# term absent from a model is simply skipped, so one layout serves all five.
TERMS <- list(
  list(t = "(Intercept)",      p = "Intercept",          h = "Intercept"),
  list(t = "lncost",           p = "ln cost",            h = "ln cost"),
  # the cost-direction tables' regressor: clipped logit accuracy
  list(t = "la",               p = "logit accuracy",     h = "logit accuracy"),
  list(t = "tc",               p = "time",               h = "time"),
  list(t = "I(lncost^2)",      p = "ln cost^2",          h = "ln cost<sup>2</sup>"),
  list(t = "I(la^2)",          p = "logit accuracy^2",   h = "logit accuracy<sup>2</sup>"),
  list(t = "I(tc^2)",          p = "time^2",             h = "time<sup>2</sup>"),
  list(t = "lncost:tc",        p = "ln cost x time",     h = "ln cost &times; time"),
  list(t = "la:tc",            p = "logit accuracy x time", h = "logit accuracy &times; time"),
  # the Box-Cox specification's rows (boxcox_frontier.R): transformed cost and
  # time, their product, and the profiled transform parameters themselves
  list(t = "phic",             p = "BC cost",            h = "BC cost"),
  # the cost-direction BC terms (fit_cost_bc): transformed odds a/(1-a)
  list(t = "phia",             p = "BC accuracy odds",            h = "BC accuracy odds"),
  list(t = "phit",             p = "BC time",            h = "BC time"),
  list(t = "phixt",            p = "BC cost x BC time",  h = "BC cost &times; BC time"),
  list(t = "phiat",            p = "BC accuracy odds x BC time",  h = "BC accuracy odds &times; BC time"),
  list(t = "lambda_cost",      p = "lambda_cost",        h = "&lambda;<sub>cost</sub>"),
  list(t = "lambda_odds",      p = "lambda_odds",        h = "&lambda;<sub>accuracy</sub>"),
  list(t = "lambda_time",      p = "lambda_time",        h = "&lambda;<sub>time</sub>"),
  # the cost-direction SFA's noise scale; the u rows below are shared with the
  # accuracy-direction SFA tables
  list(t = "logsig_v",         p = "log sigma_v",        h = "log &sigma;<sub>v</sub>"),
  list(t = "logsig_(Intercept)", p = "log sigma_u",      h = "log &sigma;<sub>u</sub>"),
  list(t = "logsig_tc",        p = "log sigma_u x time", h = "log &sigma;<sub>u</sub> &times; time")
)

MODELS <- list(
  list(key = "S",           label = "Logistic, all tests"),
  list(key = "A",           label = "Stochastic frontier"),
  list(key = "B",           label = "Stochastic frontier (time-dependent inefficiency spread)"),
  list(key = "paretologit", label = "Logistic, Pareto points"),
  list(key = "paretologitenv",
       label = "Logistic, Pareto points, envelope-constrained"),
  # the cost-direction duals (cost_frontier.R), in the same order as their
  # accuracy-direction counterparts above
  list(key = "costols",      label = "Least squares on log cost, all tests"),
  list(key = "costsfa",      label = "Stochastic cost frontier"),
  list(key = "costsfab",     label = "Stochastic cost frontier (time-dependent inefficiency spread)"),
  list(key = "costgridols",  label = "Least squares on log cost, Pareto grid"),
  list(key = "costgridolsenv",
       label = "Least squares on log cost, Pareto grid, envelope-constrained")
)

COST_KEYS <- c("costols", "costsfa", "costsfab", "costgridols",
               "costgridolsenv")

## ---- fitting and extraction ---------------------------------------------------

# One fit per specification x benchmark, from the shared store (fit_store.R):
# under run_all.R these are the very objects the figure scripts drew, so table
# and figure cannot drift apart and nothing heavy is fitted twice.
fit_grid <- function(key) {
  if (key %in% c("paretologit", "paretologitenv"))
    return(store_grid(key))
  lapply(setNames(nm = names(TIME_FORMS)), function(tt)
    store_specs()[[paste0(key, "_", tt)]]$fits)
}

# estimate/SE table with one row per term. maxLik keeps a beta_/logsig_ prefix;
# strip beta_ only, so the sigma terms stay distinguishable from the frontier's.
est_se <- function(fit) {
  # This branch must precede the glm/lm ones: the frontier logit IS a glm and
  # the grid OLS IS an lm, but their rows are grid nodes, so any SE that
  # machinery reports is a fiction. The two envelopes have no likelihood at
  # all.
  if (inherits(fit, c("envelope_frontier", "pareto_grid_logit",
                      "lncost_grid_ols", "cost_envelope_frontier"))) {
    cf <- coef(fit)
    return(data.frame(term = names(cf), est = unname(cf),
                      se = NA_real_, stringsAsFactors = FALSE))
  }
  # glm before lm: a glm inherits "lm" too
  if (inherits(fit, "glm") || inherits(fit, "lm")) {
    m <- lmtest::coeftest(fit, vcov = sandwich::vcovHC(fit, type = "HC1"))
    return(data.frame(term = rownames(m), est = m[, 1], se = m[, 2],
                      stringsAsFactors = FALSE))
  }
  # A thin benchmark can leave a maxLik fit with a flat direction (e.g. the
  # time-varying inefficiency scale where the dates barely vary): the Hessian
  # is then singular and the sandwich cannot be built. Report the point
  # estimates with EMPTY standard-error cells rather than crash or invent.
  m <- tryCatch(summary_robust(fit), error = function(e) NULL)
  if (is.null(m)) {
    cf <- coef(fit)[activePar(fit)]
    return(data.frame(term = sub("^beta_", "", names(cf)), est = unname(cf),
                      se = NA_real_, stringsAsFactors = FALSE))
  }
  se <- m[, 2]
  se[!is.finite(se)] <- NA_real_   # near-singular: negative variances
  data.frame(term = sub("^beta_", "", rownames(m)), est = m[, 1], se = se,
             stringsAsFactors = FALSE)
}

# For the frontier logit the "sample" is the grid nodes carrying a defined
# frontier value, read off the fit itself; for everything else it is the runs.
n_obs <- function(key, b, fit = NULL) {
  if (key %in% c("paretologit", "paretologitenv")) attr(fit, "n_grid")
  else sum(d$benchmark == b)
}

## ---- rate of cost decline ------------------------------------------------------
#
# Holding accuracy fixed in z = b0 + b_x*ln c + b_t*tc:
#
#     0 = b_x d(ln c) + b_t d(tc)   =>   d ln c / d tc = -b_t / b_x
#
# so the summary is minus the TIME coefficient over the LN COST one, in log
# dollars per year. Units settle the direction: b_t/b_x is (logits/yr) over
# (logits/log-$) = log-$/yr, while b_x/b_t would be years per log-dollar, which
# is a cost-time tradeoff rather than a rate of decline.
#
# Reported PER QUARTER: tc is in years, so the annual log decline -b_t/b_x is
# divided by 4 before the transform x -> 100*(1 - exp(x)), which turns log
# points into a percentage drop and flips the sign so a genuine decline reads
# positive. (Compounding, not division, relates it to the annual figure
# plot_frontiers.R prints: 1 - drop_yr = (1 - drop_qtr)^4.)
#
# The transformation is applied INSIDE the function differentiated, not to the
# ratio afterwards: the delta method is not invariant to reparametrisation, and
# an SE computed on -b_t/b_x is not the SE of 1 - exp(-b_t/b_x). (Composing by
# hand would give 100*exp(x)*SE_x, which the accompanying test confirms.)
#
# Linear specifications only. With the quadratic's lncost:tc interaction the cost
# slope is b_x + b_xt*tc, so the ratio is no longer a single number for the
# column and reporting one would be a fiction.
#
# The standard error is the delta method. car::deltaMethod() is the usual package
# for this and msm::deltamethod() the other; neither is installed here, so the
# gradient comes from numDeriv (which is) and is checked against the closed form
#     d/db_x (-b_t/b_x) =  b_t/b_x^2
#     d/db_t (-b_t/b_x) = -1/b_x
# in the accompanying test rather than trusted blind.
DECLINE <- function(p) 100 * (1 - exp(-p[["tc"]] / p[["lncost"]] / 4))

# Coefficients (beta_ prefix stripped) and the covariance that goes with them,
# or V = NULL where none exists: the envelope has no covariance matrix, and the
# frontier logit's would be a fiction built on grid nodes.
coef_vcov <- function(fit) {
  b <- coef(fit)
  names(b) <- sub("^beta_", "", names(b))
  if (inherits(fit, c("envelope_frontier", "pareto_grid_logit",
                      "lncost_grid_ols", "cost_envelope_frontier")))
    return(list(b = b, V = NULL))
  V <- if (inherits(fit, "glm") || inherits(fit, "lm"))
    sandwich::vcovHC(fit, type = "HC1") else
    tryCatch(vcov_robust(fit), error = function(e) NULL)   # singular Hessian
  if (is.null(V)) return(list(b = b, V = NULL))
  if (is.null(rownames(V))) dimnames(V) <- list(names(coef(fit)), names(coef(fit)))
  dimnames(V) <- list(sub("^beta_", "", rownames(V)),
                      sub("^beta_", "", colnames(V)))
  list(b = b, V = V)
}

# The delta-method core, on coefficients rather than a fit, so the pooled
# column's slopes flow through the identical calculation. V = NULL means a
# point estimate standing alone.
decline_from <- function(b, V) {
  need <- c("lncost", "tc")
  if (!all(need %in% names(b))) return(NULL)
  est <- DECLINE(b)
  if (is.null(V)) return(list(est = est, se = NA_real_))
  g <- function(p) { names(p) <- need; DECLINE(p) }
  J <- numDeriv::grad(g, b[need])
  list(est = est, se = sqrt(drop(t(J) %*% V[need, need] %*% J)))
}

cost_decline <- function(fit) {
  cv <- coef_vcov(fit)
  decline_from(cv$b, cv$V)
}

# The cost-direction version: the response is already ln cost, so the decline
# is the TIME COEFFICIENT's direct transform, 100*(1 - exp(b_tc/4)) -- no
# ratio, no b_x in the denominator. The SE is the delta method on that
# transform, |d/db| = exp(b_tc/4)/4 in closed form (checked against numDeriv
# in the accompanying test); NA where the fit carries no covariance.
DECLINE_DUAL <- function(g) 100 * (1 - exp(g / 4))

decline_dual_from <- function(g, se_g) {
  list(est = DECLINE_DUAL(g),
       se = if (is.na(se_g)) NA_real_ else 100 * exp(g / 4) * se_g / 4)
}

cost_decline_dual <- function(fit) {
  cv <- coef_vcov(fit)
  if (!"tc" %in% names(cv$b)) return(NULL)
  se_g <- if (is.null(cv$V) || !"tc" %in% rownames(cv$V)) NA_real_ else
    sqrt(cv$V["tc", "tc"])
  decline_dual_from(cv$b[["tc"]], se_g)
}

# percentage points, so one decimal rather than the coefficients' three
decl_cell <- function(cl, which = "est") {
  if (is.null(cl$decline)) return("")
  if (which == "est") fmt(cl$decline$est, 1)
  else if (is.na(cl$decline$se)) "" else sprintf("(%s)", fmt(cl$decline$se, 1))
}

# The at/past-SOTA decline rows (cost tables; post_sota_decline_qtr,
# cost_frontier.R): one label per entry of the qtrs vector, cells in the
# decline row's percentage format
SOTA_ROW_LABELS <- c("cost drop, %/qtr, at SOTA",
                     sprintf("%d quarter%s past SOTA", 1:8,
                             ifelse(1:8 > 1, "s", "")))
sota_cell <- function(cl, i) {
  if (is.null(cl$sota) || is.na(cl$sota[i])) "" else fmt(cl$sota[i], 1)
}

# Quadratic block test. LR where a likelihood exists, Wald where it does not,
# nothing for the envelope. Returns label + statistic + df + p, plus `row`
# naming which test row of the table the result belongs to.
quad_test <- function(key, b, fit_lin, fit_quad) {
  extra <- c("I(lncost^2)", "I(tc^2)", "lncost:tc")
  # no test for the frontier logits: no likelihood, and no sampling
  # distribution for a Wald statistic built on grid nodes
  if (key %in% c("paretologit", "paretologitenv")) return(NULL)
  if (key %in% c("A", "B")) {
    df <- sum(activePar(fit_quad)) - sum(activePar(fit_lin))
    stat <- 2 * (as.numeric(logLik(fit_quad)) - as.numeric(logLik(fit_lin)))
    return(list(kind = "LR", stat = max(stat, 0), df = df,
                p = pchisq(max(stat, 0), df, lower.tail = FALSE), row = "quad"))
  }
  V <- sandwich::vcovHC(fit_quad, type = "HC1")
  bq <- coef(fit_quad)
  i <- intersect(extra, names(bq))
  i <- i[!is.na(bq[i])]
  if (!length(i)) return(NULL)
  stat <- drop(t(bq[i]) %*% solve(V[i, i, drop = FALSE]) %*% bq[i])
  list(kind = "Wald", stat = stat, df = length(i),
       p = pchisq(stat, length(i), lower.tail = FALSE), row = "quad")
}

# Box-Cox vs linear. The BC specification nests the linear one at lambda_cost =
# 0, lambda_time = 1, b_phixt = 0, so for the likelihood families the profile
# LR has one df per free lambda plus one for the product term -- 3, or 2 where
# lambda_time is fixed (fm13). S gets nothing: its lambdas are profiled on the
# quasi-deviance, which supports no LR, and a Wald statistic has no covariance
# for the lambdas to draw on.
bc_test <- function(key, fit_lin, fit_bc) {
  if (!key %in% c("A", "B")) return(NULL)
  df <- sum(activePar(fit_bc)) - sum(activePar(fit_lin)) +
    sum(attr(fit_bc, "bc_lambda_free"))
  stat <- 2 * (as.numeric(logLik(fit_bc)) - as.numeric(logLik(fit_lin)))
  list(kind = "LR", stat = max(stat, 0), df = df,
       p = pchisq(max(stat, 0), df, lower.tail = FALSE), row = "bc")
}

# The BC fit's estimate/SE rows, plus the profiled lambdas as rows of their own:
# point values only (the profile provides no covariance for them), and only the
# FREE ones -- printing fm13's fixed lambda_time as if estimated would be a lie
# the notes would then have to walk back.
est_se_bc <- function(fit) {
  lam  <- attr(fit, "bc_lambda")
  free <- attr(fit, "bc_lambda_free")
  rbind(est_se(fit),
        data.frame(term = names(lam)[free], est = unname(lam[free]),
                   se = NA_real_, stringsAsFactors = FALSE))
}

# At most THREE SIGNIFICANT DIGITS, under two constraints. Digits left of
# the decimal point are never dropped -- no scientific notation, so a
# coefficient like -107297.415 prints as -107297, all of its integer part
# and nothing after. And `digits` still caps the DECIMALS as before, so a
# cell type keeps its accustomed budget (rates one decimal, coefficients
# three) and small values are not padded out to three significant digits
# they never displayed. The 2 is 3 significant digits minus the one digit
# floor(log10) already counts.
fmt <- function(x, digits = 3) {
  vapply(x, function(v) {
    if (is.na(v)) return("")
    d <- if (v == 0) digits else
      max(0, min(digits, 2 - floor(log10(abs(v)))))
    formatC(v, format = "f", digits = d)
  }, character(1))
}
fmt_p <- function(p) {
  if (is.na(p)) "" else if (p < 0.0001) "<0.0001" else formatC(p, format = "f", digits = 4)
}

## ---- the averages column ----------------------------------------------------------
#
# The tables' one summary column: the SIMPLE AVERAGE over the primary
# benchmarks of the quarterly decline rates, taken in log space (the mean of
# the implied annual log-cost changes, converted back to a quarterly rate) --
# the same summary the rate-comparison table leads with. It replaced two
# pooled columns, an ECI-weighted (alpha^2 x history x p(1-p)) pooling of
# per-benchmark coefficients and a pooled-data regression with benchmark
# fixed effects: a delete-one-model jackknife put per-benchmark sampling
# errors well above the gaps between the candidate weightings, so the
# transparent, calibration-free average was elevated over both. It carries
# NO coefficient estimates: slopes sit on benchmark-specific scales, and
# this column deliberately stays out of that business.
avg_rate <- function(r) {
  ok <- is.finite(r)
  if (!any(ok)) return(NA_real_)
  100 * (1 - exp(mean(4 * log(1 - r[ok] / 100)) / 4))
}
avg_col <- function(cols) {
  pri <- Filter(function(cl) cl$bench %in% PRIMARY_BENCHES, cols)
  dec <- avg_rate(vapply(pri, function(cl)
    if (is.null(cl$decline)) NA_real_ else cl$decline$est, numeric(1)))
  sota <- NULL
  if (any(vapply(pri, function(cl) !is.null(cl$sota), logical(1))))
    sota <- vapply(seq_along(SOTA_ROW_LABELS), function(i)
      avg_rate(vapply(pri, function(cl)
        if (is.null(cl$sota)) NA_real_ else cl$sota[i], numeric(1))),
      numeric(1))
  list(bench = "avg", spec = NA_character_, head = "Average, primary",
       es = data.frame(term = character(0), est = numeric(0),
                       se = numeric(0), stringsAsFactors = FALSE),
       n = NA_integer_,
       decline = if (is.na(dec)) NULL else list(est = dec, se = NA_real_),
       sota = sota, test = NULL)
}

# Quadratic-vs-linear for the cost models: LR where a likelihood exists (the
# SFA duals), the Wald analogue on the HC1 covariance for the all-runs OLS
# (mirroring model S), nothing for the two grid-based fits.
cost_quad_test <- function(key, fit_lin, fit_quad) {
  if (key %in% c("costgridols", "costgridolsenv")) return(NULL)
  if (key %in% c("costsfa", "costsfab")) {
    df <- sum(activePar(fit_quad)) - sum(activePar(fit_lin))
    stat <- 2 * (as.numeric(logLik(fit_quad)) - as.numeric(logLik(fit_lin)))
    return(list(kind = "LR", stat = max(stat, 0), df = df,
                p = pchisq(max(stat, 0), df, lower.tail = FALSE), row = "quad"))
  }
  V <- sandwich::vcovHC(fit_quad, type = "HC1")
  bq <- coef(fit_quad)
  i <- intersect(c("I(la^2)", "I(tc^2)", "la:tc"), names(bq))
  i <- i[!is.na(bq[i])]
  if (!length(i)) return(NULL)
  stat <- drop(t(bq[i]) %*% solve(V[i, i, drop = FALSE]) %*% bq[i])
  list(kind = "Wald", stat = stat, df = length(i),
       p = pchisq(stat, length(i), lower.tail = FALSE), row = "quad")
}

# Box-Cox vs linear, SFA duals only: the profile LR on the product term plus
# one df per free lambda, exactly as bc_test() counts for A and B. The
# least-squares fits get nothing -- their lambdas are profiled on SSR, which
# supports no LR without a distributional commitment the accuracy-direction
# S model also declines to make.
cost_bc_test <- function(key, fit_lin, fit_bc) {
  if (!key %in% c("costsfa", "costsfab")) return(NULL)
  df <- sum(activePar(fit_bc)) - sum(activePar(fit_lin)) +
    sum(attr(fit_bc, "bc_lambda_free"))
  stat <- 2 * (as.numeric(logLik(fit_bc)) - as.numeric(logLik(fit_lin)))
  list(kind = "LR", stat = max(stat, 0), df = df,
       p = pchisq(max(stat, 0), df, lower.tail = FALSE), row = "bc")
}

# The cost-direction mirror of build_model(): one table per specification.
# The pooled column exists for the LINEAR table only -- with curvature or
# transforms the time slope moves with (accuracy, date), so no single-number
# pooling is honest.
build_model_cost <- function(key) {
  fits <- store_cost(key)
  fits$bc <- store_cost_bc(key)

  out <- list()
  for (tt in c(names(COST_FORMS), "bc")) {
    cols <- lapply(benches, function(b) {
      f <- fits[[tt]][[b]]
      tst <- if (tt == "quad")
        cost_quad_test(key, fits$lin[[b]], f)
      else if (tt == "bc") cost_bc_test(key, fits$lin[[b]], f)
      else NULL
      # the grid fit's "sample" is defined grid nodes, resampling a few
      # dozen staircase corners; everything else fits the positive-accuracy
      # runs (zeros carry no logit coordinate)
      list(bench = b, spec = tt, head = LABELS[[b]],
           es = if (tt == "bc") est_se_bc(f) else est_se(f),
           n = if (key %in% c("costgridols", "costgridolsenv"))
             attr(f, "n_grid") else nrow(iso_runs(d[d$benchmark == b, ])),
           decline = if (tt == "lin") cost_decline_dual(f) else NULL,
           # the at/past-SOTA decline rates (post_sota_decline_qtr,
           # cost_frontier.R): every specification gets them -- on the linear
           # fit they repeat the cost drop row by construction, which is
           # itself the point of comparison with the curved columns. The
           # pooled column stays NULL: achievement dates are benchmark-
           # specific, so no single SOTA line serves the stack.
           sota = post_sota_decline_qtr(f, d[d$benchmark == b, ]),
           test = tst)
    })
    # The one summary column: the simple average of the primaries' decline
    # rates (avg_col above), every specification -- the linear tables average
    # the cost drop row, and all three average the at/past-SOTA rows. The
    # pooled-runs and coefficient-pooling columns it replaced live on in the
    # git history and in fit_pooled_cost (still drawn by the figures).
    cols <- c(cols, list(avg_col(cols)))
    out[[tt]] <- cols
  }
  out
}

# Build the columns for one model, one table PER SPECIFICATION: with eleven
# benchmarks a combined linear/quadratic/Box-Cox table no longer fits any
# page, so each specification gets its own table of one column per benchmark
# -- the quadratic table's columns carrying the quadratic-vs-linear test, the
# BC table's the BC-vs-linear test -- plus a pooled column for the linear and
# quadratic tables (the BC transforms leave no common scale to pool on).
# Returns a LIST keyed lin/quad/bc, each element a column list the renderers
# consume.
build_model <- function(key) {
  if (key %in% COST_KEYS) return(build_model_cost(key))
  grid <- fit_grid(key)
  # store_bc seeds every slower family's profile from S's lambdas, the same
  # seeding this script used to do locally with its BC_LSTART cache
  bc_fits <- store_bc(key)
  out <- list()
  for (tt in c(names(TIME_FORMS), "bc")) {
    cols <- lapply(benches, function(b) {
      f <- if (tt == "bc") bc_fits[[b]] else grid[[tt]][[b]]
      tst <- if (tt == "quad") quad_test(key, b, grid$lin[[b]], f)
        else if (tt == "bc") bc_test(key, grid$lin[[b]], f)
        else NULL
      list(bench = b, spec = tt, head = LABELS[[b]],
           es = if (tt == "bc") est_se_bc(f) else est_se(f),
           n = n_obs(key, b, f),
           decline = if (tt == "lin") cost_decline(f) else NULL,
           test = tst)
    })
    # The one summary column, linear tables only: the simple average of the
    # primaries' decline rates (avg_col above). Its content is the cost drop
    # row, which only the linear specification carries here, so the quad and
    # BC tables get no summary column. The pooled-runs and ECI-weighted
    # coefficient-pooling columns it replaced live on in the git history and
    # in fit_pooled_acc (still drawn by the figures).
    if (tt == "lin") cols <- c(cols, list(avg_col(cols)))
    out[[tt]] <- cols
  }
  out
}

# Which term rows actually appear anywhere in this model
active_terms <- function(cols) {
  present <- unique(unlist(lapply(cols, function(cl) cl$es$term)))
  Filter(function(x) x$t %in% present, TERMS)
}

cell <- function(cl, term, which = "est") {
  r <- cl$es[cl$es$term == term, ]
  if (!nrow(r)) return("")
  if (which == "est") fmt(r$est[1])
  else if (is.na(r$se[1])) "" else sprintf("(%s)", fmt(r$se[1]))
}

# Consecutive columns sharing a benchmark, in order: the spans for the header.
# Built by run-length encoding the column heads rather than assuming two columns
# per benchmark, so adding a third specification later needs no change here.
col_groups <- function(cols) {
  r <- rle(vapply(cols, function(cl) cl$head, character(1)))
  data.frame(head = r$values, span = r$lengths, stringsAsFactors = FALSE)
}

# Which data columns begin a new benchmark group (used to draw the separators
# that carry the pairing now that the sub-labels are gone). The first group needs
# no separator, so it is excluded.
group_starts <- function(cols) {
  g <- col_groups(cols)
  cumsum(c(1, head(g$span, -1)))[-1]
}

# One rendered test row per test family (quad vs linear, bc vs linear), each
# with its chi-square cells and p cells. df is folded into the row LABEL when
# every column in the row shares it -- quad always does, adding the same three
# terms everywhere; bc does except where fm13's fixed lambda_time drops its df
# to 2, in which case the df moves into the cells. Returns NULL when no column
# carries a test of this family, and otherwise list(label, cells, ps).
test_row_data <- function(cols, tr) {
  tests <- lapply(cols, function(cl)
    if (!is.null(cl$test) && identical(cl$test$row, tr)) cl$test else NULL)
  live <- !vapply(tests, is.null, logical(1))
  if (!any(live)) return(NULL)
  kind <- unique(vapply(tests[live], `[[`, character(1), "kind"))[1]
  dfs  <- unique(vapply(tests[live], `[[`, numeric(1), "df"))
  nm   <- c(quad = "quadratic", bc = "Box-Cox")[[tr]]
  label <- if (length(dfs) == 1)
    sprintf("%s vs linear, %s (df %d)", nm, kind, dfs) else
      sprintf("%s vs linear, %s", nm, kind)
  cells <- vapply(tests, function(t) {
    if (is.null(t)) "" else if (length(dfs) == 1) fmt(t$stat, 2) else
      sprintf("%s (df %d)", fmt(t$stat, 2), t$df)
  }, character(1))
  ps <- vapply(tests, function(t) if (is.null(t)) "" else fmt_p(t$p),
               character(1))
  list(label = label, cells = cells, ps = ps)
}

TEST_ROWS <- c("quad", "bc")

## ---- HTML ---------------------------------------------------------------------

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  gsub("<", "&lt;", x, fixed = TRUE)
}

# Real minus signs on negative estimates: sprintf leaves an ASCII hyphen,
# which is narrower than the digits it sits against and reads as a dash
# rather than a sign. Applied at render time, one wrapper per format, so the
# shared cell builders stay plain text; anchored to the leading position,
# where a data cell's sign lives (labels, which can carry real hyphens, are
# never wrapped).
html_minus <- function(x) sub("^-", "&minus;", x)
rtf_minus  <- function(x) sub("^-", "\\\\u8722-", x)

html_table <- function(key, label, cols, tt) {
  trm <- active_terms(cols)
  kind <- unique(unlist(lapply(cols, function(cl) cl$test$kind)))
  # the envelope has no standard errors, so it gets no (empty) SE rows either
  has_se <- any(!is.na(unlist(lapply(cols, function(cl) cl$es$se))))

  o <- c('<!DOCTYPE html>', '<html lang="en"><head><meta charset="UTF-8" />',
         sprintf('<title>%s</title>', html_escape(label)),
         '<style>',
         # #fcfcfb throughout: the same surface colour the figures are painted
         # on, and the single background the viewers now use -- a table page
         # inside the viewer's iframe should show no seam against it
         'body{font-family:"Segoe UI",Arial,sans-serif;margin:12px;color:#1d1d1d;background:#fcfcfb}',
         'h1{font-size:1.25em;margin:0 0 4px 0}',
         'p.sub{color:#5e5e5e;margin:0 0 16px 0;font-size:.9em}',
         'table{border-collapse:collapse;background:#fcfcfb;font-size:.86em}',
         'th,td{padding:4px 10px;text-align:right;white-space:nowrap}',
         'th:first-child,td:first-child{text-align:left}',
         'thead th{border-bottom:1px solid #1d1d1d;font-weight:600;vertical-align:bottom}',
         'thead th[colspan]{text-align:center;padding-bottom:3px}',
         # headers wrap inside a width-capped block, so a long benchmark name
         # no longer sets its column's width -- the data cells do; max-width
         # on the th itself is unreliable under auto table layout, a block
         # inside the cell is not
         'thead th .hd{white-space:normal;overflow-wrap:break-word;max-width:6.5em;margin:0 auto}',
         # the pairing is carried by these separators now that the columns are
         # not individually labelled linear / quadratic
         'th.grp,td.grp{border-left:1px solid #e6e6e2}',
         'tr.se td{color:#5e5e5e;padding-top:0}',
         'tr.gap td{border-top:1px solid #d9d9d5}',
         'tfoot td{border-top:1px solid #1d1d1d;font-size:.92em;color:#5e5e5e;',
         '  text-align:left;white-space:normal;padding-top:8px;max-width:900px}',
         '</style></head><body>',
         sprintf('<h1>%s</h1>', html_escape(label)),
         sprintf('<p class="sub">%s</p>',
                 if (has_se) "Standard errors in parentheses."
                 else "Point estimates only; see the notes for why no standard errors."),
         '<table><thead><tr><th></th>')
  # One spanning header per benchmark. The two columns beneath are linear and
  # quadratic in that order, unlabelled: the quadratic is identifiable as the one
  # carrying the second-order rows, and dropping the two words per column is what
  # takes the table down to a readable width.
  grp <- col_groups(cols)
  for (i in seq_len(nrow(grp)))
    o <- c(o, sprintf('<th colspan="%d"%s><div class="hd">%s</div></th>',
                      grp$span[i],
                      if (i > 1) ' class="grp"' else '',
                      html_escape(grp$head[i])))
  o <- c(o, '</tr></thead><tbody>')

  starts <- group_starts(cols)
  cls <- function(j) if (j %in% starts) ' class="grp"' else ''

  for (x in trm) {
    o <- c(o, sprintf('<tr><td>%s</td>', x$h))
    for (j in seq_along(cols))
      o <- c(o, sprintf('<td%s>%s</td>', cls(j),
                        html_minus(cell(cols[[j]], x$t, "est"))))
    o <- c(o, '</tr>')
    if (has_se) {
      o <- c(o, '<tr class="se"><td></td>')
      for (j in seq_along(cols))
        o <- c(o, sprintf('<td%s>%s</td>', cls(j),
                          html_minus(cell(cols[[j]], x$t, "se"))))
      o <- c(o, '</tr>')
    }
  }

  # rate of cost decline, above N; the rule that separated the summary block from
  # the coefficients moves up onto this row
  has_decl <- any(vapply(cols, function(cl) !is.null(cl$decline), logical(1)))
  has_decl_se <- any(vapply(cols, function(cl)
    !is.null(cl$decline) && !is.na(cl$decline$se), logical(1)))
  if (has_decl) {
    o <- c(o, '<tr class="gap"><td>cost drop, %/qtr</td>')
    for (j in seq_along(cols))
      o <- c(o, sprintf('<td%s>%s</td>', cls(j),
                        html_minus(decl_cell(cols[[j]], "est"))))
    o <- c(o, '</tr>')
    if (has_decl_se) {
      o <- c(o, '<tr class="se"><td></td>')
      for (j in seq_along(cols))
        o <- c(o, sprintf('<td%s>%s</td>', cls(j),
                          html_minus(decl_cell(cols[[j]], "se"))))
      o <- c(o, '</tr>')
    }
  }

  # the at/past-SOTA decline rates (cost tables only): point estimates, no
  # SE rows -- the auxiliary SOTA regression sits outside every fit's own
  # error model
  has_sota <- any(vapply(cols, function(cl) !is.null(cl$sota), logical(1)))
  if (has_sota) {
    for (i in seq_along(SOTA_ROW_LABELS)) {
      o <- c(o, sprintf('<tr%s><td>%s</td>',
                        if (i == 1 && !has_decl) ' class="gap"' else '',
                        SOTA_ROW_LABELS[i]))
      for (j in seq_along(cols))
        o <- c(o, sprintf('<td%s>%s</td>', cls(j),
                          html_minus(sota_cell(cols[[j]], i))))
      o <- c(o, '</tr>')
    }
  }

  o <- c(o, sprintf('<tr%s><td>N</td>',
                    if (has_decl || has_sota) '' else ' class="gap"'))
  for (j in seq_along(cols))
    o <- c(o, sprintf('<td%s>%s</td>', cls(j),
                      if (is.na(cols[[j]]$n)) "" else
                        sprintf("%d", cols[[j]]$n)))
  o <- c(o, '</tr>')

  for (tr in TEST_ROWS) {
    trd <- test_row_data(cols, tr)
    if (is.null(trd)) next
    o <- c(o, sprintf('<tr><td>%s, &chi;&sup2;</td>', html_escape(trd$label)))
    for (j in seq_along(cols))
      o <- c(o, sprintf('<td%s>%s</td>', cls(j), trd$cells[j]))
    o <- c(o, '</tr><tr class="se"><td>p</td>')
    # html_escape, not raw: fmt_p can return "<0.0001", and an unescaped "<0"
    # is swallowed by the parser as the start of a tag -- the cell renders empty,
    # which reads as "no test" exactly where the test is most significant.
    for (j in seq_along(cols))
      o <- c(o, sprintf('<td%s>%s</td>', cls(j), html_escape(trd$ps[j])))
    o <- c(o, '</tr>')
  }

  o <- c(o, sprintf('</tbody><tfoot><tr><td colspan="%d">%s</td></tr></tfoot></table>',
                    length(cols) + 1, notes_html(key, kind, tt)),
         '</body></html>')
  o
}

## ---- RTF ----------------------------------------------------------------------

rtf_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\\{", "\\\\{", x)
  gsub("\\}", "\\\\}", x)
}

# One row: `widths` are cumulative twips, `cells` the already-escaped strings.
#
# Every cell must be `\pard\intbl ... \cell`, and the table must be closed with a
# `\pard` after the last `\row`. Omitting `\pard` is not a cosmetic slip: Word
# cannot then tell where the table ends, reports "a table in this document has
# become corrupted", and its recovery path renders the rows a second time -- two
# stacked copies of the same table from a file that contains it once.
#
# `\b` is closed explicitly with `\b0` and `\fs18` restated per cell rather than
# relied on from the document default, because `\pard` resets paragraph but not
# character formatting, so bold would otherwise leak from the header row into
# every row beneath it.
# `merge` marks horizontal spans: "first" on the leading cell of a span, "cont"
# on each cell absorbed into it. A merged span still needs one \cell per column,
# so the \cell / \cellx counts the validator checks stay equal.
# `sep` marks cells that open a new benchmark group and get a left rule.
# `align` is "l", "r" or "c" per cell.
rtf_row <- function(cells, widths, bold = FALSE, top = FALSE, bottom = FALSE,
                    merge = NULL, sep = NULL, align = NULL) {
  n <- length(cells)
  if (is.null(merge)) merge <- rep("", n)
  if (is.null(sep))   sep   <- rep(FALSE, n)
  if (is.null(align)) align <- c("l", rep("r", n - 1))
  defn <- paste0("\\trowd\\trgaph60\\trleft0",
                 paste0(vapply(seq_len(n), function(i) paste0(
                   if (merge[i] == "first") "\\clmgf" else
                     if (merge[i] == "cont") "\\clmrg" else "",
                   if (sep[i]) "\\clbrdrl\\brdrs\\brdrw10" else "",
                   if (top) "\\clbrdrt\\brdrs\\brdrw10" else "",
                   if (bottom) "\\clbrdrb\\brdrs\\brdrw10" else "",
                   sprintf("\\cellx%d", widths[i])), character(1)),
                   collapse = ""))
  body <- vapply(seq_len(n), function(i) {
    sprintf("\\pard\\intbl\\itap1%s\\fs18%s %s\\cell",
            switch(align[i], l = "\\ql", r = "\\qr", c = "\\qc"),
            if (bold) "\\b" else "\\b0", cells[i])
  }, character(1))
  c(defn, paste0(body, collapse = ""), "\\row")
}

rtf_table <- function(key, label, cols, tt) {
  trm <- active_terms(cols)
  kind <- unique(unlist(lapply(cols, function(cl) cl$test$kind)))
  has_se <- any(!is.na(unlist(lapply(cols, function(cl) cl$es$se))))

  # With the per-column "linear"/"quadratic" words gone, the widest data string
  # is an SE like "(11.102)" and the widest label is the test row. The Box-Cox
  # columns take the count to fourteen (four benchmark trios plus the pooled
  # pair), which no longer fits portrait letter at a readable width, so the page
  # is LANDSCAPE: 15840 wide less 2*720 margins = 14400 usable, and
  # 2400 + 14*850 = 14300 fits.
  w1 <- 2400; wc <- 850
  widths <- c(w1, w1 + seq_len(length(cols)) * wc)

  starts <- group_starts(cols)
  # +1 throughout because the row-label column occupies position 1
  sep_data <- c(FALSE, seq_along(cols) %in% starts)

  # \uc1: one fallback character after each \uN unicode escape -- the
  # convention rtf_minus's unicode-minus-with-hyphen-fallback relies on
  o <- c("{\\rtf1\\ansi\\ansicpg1252\\uc1\\deff0",
         "{\\fonttbl{\\f0\\fswiss Calibri;}}",
         "\\paperw15840\\paperh12240\\landscape\\margl720\\margr720\\margt720\\margb720",
         "\\f0\\fs18",
         sprintf("\\pard\\ql{\\b\\fs24 %s}\\par", rtf_escape(label)),
         sprintf("\\pard\\ql{\\i %s}\\par\\par",
                 if (has_se) "Standard errors in parentheses."
                 else "Point estimates only; see the notes for why no standard errors."))

  # Spanning benchmark header: one label per pair, centred over it. The absorbed
  # columns carry empty text; \clmrg does the joining.
  grp <- col_groups(cols)
  head_cells <- c(""); head_merge <- c(""); head_align <- c("l")
  for (i in seq_len(nrow(grp))) {
    head_cells <- c(head_cells, rtf_escape(grp$head[i]),
                    rep("", grp$span[i] - 1))
    head_merge <- c(head_merge, "first", rep("cont", grp$span[i] - 1))
    head_align <- c(head_align, rep("c", grp$span[i]))
  }
  o <- c(o, rtf_row(head_cells, widths, bold = TRUE, bottom = TRUE,
                    merge = head_merge, sep = sep_data, align = head_align))

  for (x in trm) {
    o <- c(o, rtf_row(c(rtf_escape(x$p),
                        vapply(cols, function(cl)
                          rtf_minus(cell(cl, x$t, "est")), character(1))),
                      widths, sep = sep_data))
    if (has_se)
      o <- c(o, rtf_row(c("", vapply(cols, function(cl)
                                       rtf_minus(cell(cl, x$t, "se")),
                                     character(1))), widths, sep = sep_data))
  }
  has_decl <- any(vapply(cols, function(cl) !is.null(cl$decline), logical(1)))
  has_decl_se <- any(vapply(cols, function(cl)
    !is.null(cl$decline) && !is.na(cl$decline$se), logical(1)))
  if (has_decl) {
    o <- c(o, rtf_row(c("cost drop, %/qtr",
                        vapply(cols, function(cl)
                          rtf_minus(decl_cell(cl, "est")), character(1))),
                      widths, top = TRUE, sep = sep_data))
    if (has_decl_se)
      o <- c(o, rtf_row(c("", vapply(cols, function(cl)
                                       rtf_minus(decl_cell(cl, "se")),
                                     character(1))), widths, sep = sep_data))
  }
  has_sota <- any(vapply(cols, function(cl) !is.null(cl$sota), logical(1)))
  if (has_sota) for (i in seq_along(SOTA_ROW_LABELS))
    o <- c(o, rtf_row(c(SOTA_ROW_LABELS[i],
                        vapply(cols, function(cl)
                          rtf_minus(sota_cell(cl, i)), character(1))),
                      widths, top = i == 1 && !has_decl, sep = sep_data))
  o <- c(o, rtf_row(c("N", vapply(cols, function(cl)
                                    if (is.na(cl$n)) "" else as.character(cl$n),
                                  character(1))),
                    widths, top = !has_decl && !has_sota, sep = sep_data))
  for (tr in TEST_ROWS) {
    trd <- test_row_data(cols, tr)
    if (is.null(trd)) next
    o <- c(o, rtf_row(c(rtf_escape(sprintf("%s, chi-sq", trd$label)), trd$cells),
                      widths, sep = sep_data))
    o <- c(o, rtf_row(c("p", trd$ps), widths, sep = sep_data))
  }
  # \pard closes the table: without it the notes paragraph is still inside table
  # context and Word treats the whole run as a malformed table.
  o <- c(o, "\\pard\\ql\\par",
         sprintf("\\pard\\ql{\\fs16 %s}\\par", rtf_escape(notes_plain(key, kind, tt))),
         "}")
  o
}

## ---- the notes that keep each table honest -------------------------------------

# Notes for the cost-direction tables. Shorter than the accuracy notes: no
# quadratic or Box-Cox columns, and the pooling needs no ECI machinery.
notes_cost <- function(key, tt) {
  dir <- paste(
    "This model fits LN COST as the response, linear in logit accuracy and time (years,",
    "centered within benchmark), so the intercept is the fitted log cost of 50% accuracy at the",
    "benchmark's reference date. Runs scoring zero -- that is, no better than guessing, since",
    "accuracy is rescaled from each benchmark's guessing floor -- and runs scoring n/n are both",
    "excluded, logit 0 and logit 1 being unusable as coordinates. Accuracy is a",
    "REGRESSOR here, so its sampling noise is measurement error, which attenuates the",
    "logit-accuracy slope.")
  se <- switch(key,
    costols = paste(
      "Robust (HC1) standard errors. This is the reverse regression of the plain logistic: the",
      "TYPICAL cost of a run scoring a given accuracy, not a frontier."),
    costsfa = , costsfab = paste(
      "Robust (sandwich) standard errors, clustered at the model x effort group level, where the",
      "one-sided inefficiency term lives (u >= 0 is excess cost). Read the time coefficient with",
      "care: it tracks the dense cheap edge of the model-effort cells, which grows dearer as",
      "expensive reasoning configurations arrive at every accuracy level -- the cost-direction",
      "dual of the run-distribution story, not the record's decline."),
    costgridols = paste(
      "No standard errors: the response is the record cost ln C_a(t) sampled on a uniform",
      "(logit accuracy, date) grid -- the mirror of the accuracy-direction fits' (cost, date)",
      "grid -- and grid nodes are not observations. N is the nodes at which a level is defined;",
      "they resample a few dozen staircase corners, so N measures resolution, not information."),
    costgridolsenv = paste(
      "No standard errors: this is the grid OLS's objective -- the record cost ln C_a(t) sampled",
      "on the uniform (logit accuracy, date) grid -- minimised subject to the cost envelope's",
      "constraints: the surface must lie at or below every positive-accuracy run's log cost, rise",
      "with accuracy and never rise with date. Grid nodes are not observations and the constrained",
      "fit has no likelihood, so errors would be fictions twice over. N is the defined grid nodes."))
  extra <- if (key == "costsfab") paste(
    "log sigma_u x time lets the inefficiency spread move with date; the constant-scale variant",
    "is the table before this one.") else ""
  decl <- paste(
    "cost drop, %/qtr is 100*(1 - exp(b_time / 4)): the time coefficient transformed directly --",
    "the estimand is a single coefficient here, not the ratio of two as in the accuracy-direction",
    "tables. Standard error by the delta method on the transform, where a covariance exists.",
    "Linear specification only: with the quadratic's or Box-Cox's extra terms the time slope",
    "moves with accuracy and date, so no single rate describes those columns.")
  sota <- paste(
    "The at/past-SOTA rows ask the fitted surface how much faster cost falls near the",
    "frontier: an auxiliary OLS of the accuracy record's logit on date (one observation per",
    "release date) predicts each level's achievement date, and the surface's instantaneous",
    "d ln cost/dt is read at that date and at each of the following eight quarters -- the decline",
    "rate of a level the moment it becomes achievable, and as it recedes behind the frontier --",
    "averaged over the grid's accuracy levels and expressed quarterly. Each horizon keeps the",
    "levels whose own window lies inside the benchmark's observed dates, so the surface is never",
    "extrapolated; longer horizons therefore average over fewer, earlier-achieved levels, and a",
    "row is blank where none qualify. Read the deep horizons with that composition shift in",
    "mind. On the linear specification every row repeats the cost drop row by construction; the",
    "curved specifications let them differ, which is where faster-near-SOTA shows up. No",
    "standard errors: the auxiliary SOTA regression sits outside every fit's own error model.")
  tst <- switch(key,
    costsfa = , costsfab = paste(
      "Quadratic adds logit accuracy^2, time^2 and logit accuracy x time; the test is a",
      "likelihood ratio test of all three jointly -- a real one, since this model has a genuine",
      "likelihood."),
    costols = paste(
      "Quadratic adds logit accuracy^2, time^2 and logit accuracy x time, tested jointly by a",
      "Wald test on the robust covariance."),
    paste("No quadratic or Box-Cox tests: there is no likelihood behind this fit and no sampling",
          "distribution to appeal to."))
  bc <- paste(
    paste("BC columns are BOX-TIDWELL: only the REGRESSORS are transformed, and LN COST is always the",
          "response. ln cost is linear in phi(odds; lambda_odds), phi(time; lambda_time) and their",
          "product, with phi(x; lambda) = (x^lambda - 1)/lambda (log at lambda = 0) applied to the ODDS",
          "a/(1-a) -- at lambda_odds = 0 phi(odds) IS logit accuracy -- and to years since November 18,",
          "2021, when the GPT-3 API became generally available. The family nests the linear model at (lambda_odds, lambda_time) = (0, 1) with",
          "no product term. A response-side lambda was tried and removed: phi(cost) has to be inverted to",
          "be read back, and that inverse has a pole at a finite index, which the profile walked into --",
          "fitted costs near e^20 per task against a dearest observed run under $1, with 1% of grid nodes",
          "carrying up to 44% of the mean decline rate.",
          if (key %in% c("costsfa", "costsfab"))
            paste("For the SFA duals it was doubly unwelcome: a transform of the dependent variable soaks",
                  "up residual skewness, and the composite residual's skewness is exactly what identifies",
                  "the one-sided inefficiency term. Even the regressor transforms alone absorb the",
                  "frontier term on swe_bench_verified -- read that column's sigma_u accordingly.")
          else ""),
    "lambdas are profiled",
    switch(key,
           costols = paste("against the Gaussian profile likelihood -- the residual sum of squares on the",
                           "transformed scale with the classical Box-Cox Jacobian term, which is what makes",
                           "fits on different response scales comparable,"),
           costgridols = paste("against the grid's Gaussian profile likelihood (transformed-scale residual",
                               "sum of squares with the Box-Cox Jacobian term),"),
           costgridolsenv = paste("against the grid's Gaussian profile likelihood (transformed-scale residual",
                                  "sum of squares with the Box-Cox Jacobian term), subject to the cost",
                                  "envelope's constraints,"),
           "against the likelihood (no Jacobian: the response is untransformed),"),
    "and are reported without standard errors. A lambda_time whose profile runs to the edge of",
    "its search box ([-8, 8]) is treated as unidentified and locked at 1, and fm13's is fixed at",
    "1 a priori (its short span identifies no time curvature); the pooled-runs column's lambda_odds",
    "transforms exp(C/10) in place of the odds, a pure relabeling of lambda.",
    if (key %in% c("costsfa", "costsfab"))
      "The BC LR row tests the nesting restrictions, one df per free lambda plus the product term."
    else "")
  avgn <- paste(
    "Average, primary is the SIMPLE AVERAGE over the primary benchmarks",
    sprintf("(%s)", paste(intersect(PRIMARY_BENCHES, benches),
                          collapse = ", ")),
    "of the decline-rate rows, taken in log space: the mean of the implied annual log-cost",
    "changes, converted back to a quarterly rate. It carries no coefficient estimates -- slopes",
    "sit on benchmark-specific scales -- and no standard error: per-benchmark sampling errors,",
    "gauged by a delete-one-model jackknife, dwarf the differences among candidate weightings,",
    "which is also why this transparent, calibration-free summary replaced the ECI-weighted and",
    "pooled-data columns of earlier versions (see the rate comparison table for that weighted",
    "summary).")
  switch(tt,
         lin  = paste(dir, se, extra, decl, sota, avgn),
         quad = paste(dir, se, extra, tst, sota, avgn),
         bc   = paste(dir, se, extra, bc, sota, avgn))
}

notes_plain <- function(key, kind, tt = "lin") {
  if (key %in% COST_KEYS) return(notes_cost(key, tt))
  se <- switch(key,
    S = "Robust (HC1) standard errors: the quasibinomial objective is a quasi-likelihood, so inverse-Hessian errors would be wrong.",
    paretologit = paste("No standard errors: the response is the empirical Pareto frontier P_t(c) sampled at the nodes of the",
                        "envelope's fixed cost-date grid, and grid nodes are not observations. N is the nodes at which the",
                        "frontier is defined; they resample a few dozen staircase corners, so N measures resolution, not",
                        "information."),
    A = , B = "Robust (sandwich) standard errors, as the SFA objective is a quasi-likelihood in the Papke-Wooldridge sense.",
    paretologitenv = paste("No standard errors: this is the frontier logit's objective -- the empirical Pareto frontier P_t(c)",
                           "sampled on the fixed cost-date grid -- minimised subject to the envelope's constraints: the surface",
                           "must lie at or above every run (perfect scores clipped by their own sample size) and be monotone in",
                           "cost and date. Grid nodes are not observations and the constrained fit has no likelihood, so errors",
                           "would be fictions twice over. N is the nodes at which the frontier is defined."))
  tst <- if (is.null(kind) || !length(kind)) {
    "No quadratic-versus-linear test: there is no likelihood to test and no sampling distribution to appeal to."
  } else if (kind[1] == "LR") {
    "Quadratic adds ln cost^2, time^2 and ln cost x time; the test is a likelihood ratio test of all three jointly."
  } else {
    paste("Quadratic adds ln cost^2, time^2 and ln cost x time. Quasibinomial glms report no logLik, so the joint test",
          "of the three is a Wald test on the robust covariance, not a likelihood ratio test.")
  }
  decl <- paste("cost drop, %/qtr is 100*(1 - exp(-b_time / b_ln cost / 4)): the percentage fall in the cost of a fixed",
                "accuracy level per quarter (time is in years, hence the 4). Linear specification only -- with the",
                "quadratic's ln cost x time term the cost slope moves with date, so no single rate describes the column.",
                if (key %in% c("paretologit", "paretologitenv")) "" else
                  paste("Standard error by the delta method on the same robust covariance, taken on the transformed",
                        "quantity. Being symmetric it can reach past 100% where the estimated drop is near total."))
  pb <- intersect(PRIMARY_BENCHES, benches)
  avgn <- paste(
    "Average, primary is the SIMPLE AVERAGE over the primary benchmarks",
    sprintf("(%s)", paste(pb, collapse = ", ")),
    "of the cost drop row, taken in log space: the mean of the implied annual log-cost changes,",
    "converted back to a quarterly rate. It carries no coefficient estimates -- the slopes are in",
    "each benchmark's own logit units, with no common scale to average on -- and no standard",
    "error: per-benchmark sampling errors, gauged by a delete-one-model jackknife, dwarf the",
    "differences among candidate weightings, which is also why this transparent, calibration-free",
    "summary replaced the ECI-weighted and pooled-data columns of earlier versions (see the rate",
    "comparison table for that weighted summary).")
  bc <- paste(
    "BC columns are a Box-Cox alternative to the quadratic: the index is linear in phi(cost), phi(time) and",
    "their product, with phi(x; lambda) = (x^lambda - 1)/lambda (log at lambda = 0) applied to LEVEL cost per",
    "task and to years since November 18, 2021, when the GPT-3 API became generally available -- so the BC",
    "intercept is the fit at $1 per task on November 18, 2022, where both transforms vanish. phi is increasing whatever lambda is, so the surface is monotone in",
    "cost at every date -- the quadratic's bending back toward the data cannot happen -- while the product term",
    "still allows the cost slope one sign change over time.",
    paste("These are BOX-TIDWELL fits: phi acts on the REGRESSORS only and the index is the plain logit,",
          "so no Jacobian arises. A response-side lambda on the odds was tried and removed -- a transformed",
          "response has to be inverted to be read back, and that inverse has a pole at a finite index",
          "which the profile walked into.",
          if (key %in% c("A", "B"))
            paste("For the SFA families it would have been doubly unwelcome: their one-sided inefficiency",
                  "term is identified off the response's asymmetry, exactly the signal a free link-shape",
                  "parameter would compete for.")
          else ""),
    "The lambdas are estimated per benchmark by profiling",
    switch(key, S = "the quasibinomial deviance,", A = , B = "the likelihood,",
           paretologit = "the probability-scale quasi-likelihood,",
           paretologitenv = "the probability-scale quasi-likelihood, subject to the envelope's constraints,"),
    "and are reported without standard errors: the profile provides none, and",
    if (key %in% c("paretologit", "paretologitenv"))
      "this model reports none anywhere."
    else
      paste("the coefficient standard errors are conditional on the profiled lambdas --",
            "they carry no lambda uncertainty."),
    if (key %in% c("A", "B"))
      paste("The BC specification nests the linear one (lambda_cost = 0, lambda_time = 1, no product term) --",
            "its LR row tests exactly those restrictions -- but not the quadratic.")
    else if (key %in% c("S", "paretologit", "paretologitenv"))
      paste("The BC specification nests the linear one (lambda_cost = 0, lambda_time = 1, no product term)",
            "but not the quadratic.")
    else
      "The BC specification nests the linear one (lambda_cost = 0, lambda_time = 1, no product term) but not the quadratic.",
    "fm13's five months of data sit 5.7-6.1 years from the origin, over which every lambda_time fits alike, so",
    "lambda_time is fixed at 1 there rather than estimated; elsewhere it is weakly identified and best read as",
    "a shape the data tolerates rather than demands -- a lambda_time whose profile runs to the edge of its",
    "search box ([-8, 8]) is treated as unidentified and locked at 1. This table has no summary column:",
    "the average column's content is the cost drop row, which only the linear specification carries.")
  intro <- paste("Time is measured in years and centered within benchmark, so the intercept is the",
                 "frontier at each benchmark's own reference date.")
  # one table per specification now, so each carries only its own paragraphs
  switch(tt,
         lin  = paste(intro, decl, se, avgn),
         quad = paste(intro, se, tst),
         bc   = paste(se, bc))
}

notes_html <- function(key, kind, tt) {
  x <- notes_plain(key, kind, tt)
  x <- html_escape(x)
  x <- gsub("ln cost^2", "ln cost<sup>2</sup>", x, fixed = TRUE)
  x <- gsub("logit accuracy^2", "logit accuracy<sup>2</sup>", x, fixed = TRUE)
  gsub("time^2", "time<sup>2</sup>", x, fixed = TRUE)
}

## ---- emit ----------------------------------------------------------------------

dir.create(out_path("tables"), showWarnings = FALSE, recursive = TRUE)
SPEC_TITLE <- c(lin = "linear", quad = "full quadratic", bc = "Box-Cox")
for (m in MODELS) {
  colsets <- build_model(m$key)
  # the pre-split combined tables, if still present, are stale now
  unlink(out_path("tables", sprintf("regression_%s.%s", m$key, c("html", "rtf"))))
  for (tt in names(colsets)) {
    lab <- sprintf("%s -- %s specification", m$label, SPEC_TITLE[[tt]])
    writeLines(html_table(m$key, lab, colsets[[tt]], tt),
               out_path("tables", sprintf("regression_%s_%s.html", m$key, tt)))
    writeLines(rtf_table(m$key, lab, colsets[[tt]], tt),
               out_path("tables", sprintf("regression_%s_%s.rtf", m$key, tt)))
    cat("wrote regression_", m$key, "_", tt, ".html / .rtf\n", sep = "")
  }
}
