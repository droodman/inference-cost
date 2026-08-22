# Regression tables: one per statistical model, in RTF and HTML.
#
# Columns are benchmark x specification (AIME linear, AIME quadratic, AIME
# Box-Cox, Chess linear, ...), estimates with standard errors in parentheses
# beneath, plus a final Pooled pair combining the benchmarks on the ECI
# capability scale (see the pooling section below; the Box-Cox specification is
# excluded from pooling -- its per-benchmark transforms leave no common scale).
# Model order and naming follow the plot viewer, which is the version of these
# names that has actually been read by a human:
#
#   S            Logistic, all tests
#   A            Stochastic frontier
#   B            Stochastic frontier (time-dependent inefficiency spread)
#   paretologit  Logistic, Pareto points
#   envelope     Strict logistic envelope
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
src_source("fit_specs.R")
src_source("envelope_frontier.R")
src_source("boxcox_frontier.R")

d <- load_runs(drop_gpt4o_chess = FALSE)
benches <- sort(unique(d$benchmark))

## ---- benchmark discriminations for pooling --------------------------------------
#
# alpha_b = estimated_slope_scaled: each benchmark's 2PL discrimination, in
# logits per point of Epoch's ECI (Epoch Capabilities Index) scale -- the scale
# anchored at Claude 3.5 Sonnet = 130 and GPT-5 = 150, on which the `edi`
# column gives the benchmark's difficulty. TASK_LABEL (prepare_data.R) already
# holds each benchmark's exact name in that table, so it doubles as the join
# key.
#
# data/edi_scores.csv was downloaded from https://epoch.ai/data/edi_scores.csv
# on 2026-08-20.
ECI <- read.csv(data_path("edi_scores.csv"), stringsAsFactors = FALSE)
ALPHA <- setNames(ECI$estimated_slope_scaled[match(TASK_LABEL, ECI$benchmark_name)],
                  names(TASK_LABEL))
stopifnot(!anyNA(ALPHA[benches]))

## ---- row layout ---------------------------------------------------------------

# term -> (plain label for RTF, HTML label). Order is the table's row order; a
# term absent from a model is simply skipped, so one layout serves all five.
TERMS <- list(
  list(t = "(Intercept)",      p = "Intercept",          h = "Intercept"),
  list(t = "lncost",           p = "ln cost",            h = "ln cost"),
  list(t = "tc",               p = "time",               h = "time"),
  list(t = "I(lncost^2)",      p = "ln cost^2",          h = "ln cost<sup>2</sup>"),
  list(t = "I(tc^2)",          p = "time^2",             h = "time<sup>2</sup>"),
  list(t = "lncost:tc",        p = "ln cost x time",     h = "ln cost &times; time"),
  # the Box-Cox specification's rows (boxcox_frontier.R): transformed cost and
  # time, their product, and the profiled transform parameters themselves
  list(t = "phic",             p = "BC cost",            h = "BC cost"),
  list(t = "phit",             p = "BC time",            h = "BC time"),
  list(t = "phixt",            p = "BC cost x BC time",  h = "BC cost &times; BC time"),
  list(t = "lambda_cost",      p = "lambda_cost",        h = "&lambda;<sub>cost</sub>"),
  list(t = "lambda_time",      p = "lambda_time",        h = "&lambda;<sub>time</sub>"),
  list(t = "logsig_(Intercept)", p = "log sigma_u",      h = "log &sigma;<sub>u</sub>"),
  list(t = "logsig_tc",        p = "log sigma_u x time", h = "log &sigma;<sub>u</sub> &times; time")
)

MODELS <- list(
  list(key = "S",           label = "Logistic, all tests"),
  list(key = "A",           label = "Stochastic frontier"),
  list(key = "B",           label = "Stochastic frontier (time-dependent inefficiency spread)"),
  list(key = "paretologit", label = "Logistic, Pareto points"),
  list(key = "envelope",    label = "Strict logistic envelope")
)

## ---- fitting and extraction ---------------------------------------------------

# One fit per specification x benchmark, fitted UP FRONT: fit_family() fits
# every benchmark in a single call, so for the likelihood families the loop is
# over specifications only. The predecessor (fit_one) called fit_family() once
# per benchmark-spec cell and kept [[b]], recomputing each SFA fit four times
# for one use of it. The envelope and the frontier logit fit one benchmark at a
# time, so they gain nothing but lose nothing either.
fit_grid <- function(key) {
  lapply(TIME_FORMS, function(form) {
    if (key == "envelope")
      setNames(lapply(benches, function(b)
        fit_envelope(d[d$benchmark == b, ], formula = form)), benches)
    else if (key == "paretologit")
      setNames(lapply(benches, function(b)
        fit_pareto_logit(d[d$benchmark == b, ], formula = form)), benches)
    else fit_family(key, form, d)
  })
}

# estimate/SE table with one row per term. maxLik keeps a beta_/logsig_ prefix;
# strip beta_ only, so the sigma terms stay distinguishable from the frontier's.
est_se <- function(fit) {
  # This branch must precede the glm one: the frontier logit IS a glm, but its
  # rows are grid nodes, so any SE the glm machinery reports is a fiction.
  if (inherits(fit, c("envelope_frontier", "pareto_grid_logit"))) {
    cf <- coef(fit)
    return(data.frame(term = names(cf), est = unname(cf),
                      se = NA_real_, stringsAsFactors = FALSE))
  }
  if (inherits(fit, "glm")) {
    m <- lmtest::coeftest(fit, vcov = sandwich::vcovHC(fit, type = "HC1"))
    return(data.frame(term = rownames(m), est = m[, 1], se = m[, 2],
                      stringsAsFactors = FALSE))
  }
  m <- summary_robust(fit)
  data.frame(term = sub("^beta_", "", rownames(m)), est = m[, 1], se = m[, 2],
             stringsAsFactors = FALSE)
}

# For the frontier logit the "sample" is the grid nodes carrying a defined
# frontier value, read off the fit itself; for everything else it is the runs.
n_obs <- function(key, b, fit = NULL) {
  if (key == "paretologit") attr(fit, "n_grid") else sum(d$benchmark == b)
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
  if (inherits(fit, c("envelope_frontier", "pareto_grid_logit")))
    return(list(b = b, V = NULL))
  V <- if (inherits(fit, "glm")) sandwich::vcovHC(fit, type = "HC1") else
    vcov_robust(fit)
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

# percentage points, so one decimal rather than the coefficients' three
decl_cell <- function(cl, which = "est") {
  if (is.null(cl$decline)) return("")
  if (which == "est") fmt(cl$decline$est, 1)
  else if (is.na(cl$decline$se)) "" else sprintf("(%s)", fmt(cl$decline$se, 1))
}

# Quadratic block test. LR where a likelihood exists, Wald where it does not,
# nothing for the envelope. Returns label + statistic + df + p, plus `row`
# naming which test row of the table the result belongs to.
quad_test <- function(key, b, fit_lin, fit_quad) {
  extra <- c("I(lncost^2)", "I(tc^2)", "lncost:tc")
  # no test for the envelope or the frontier logit: no likelihood, and no
  # sampling distribution for a Wald statistic built on grid nodes
  if (key %in% c("envelope", "paretologit")) return(NULL)
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

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}
fmt_p <- function(p) {
  if (is.na(p)) "" else if (p < 0.0001) "<0.0001" else formatC(p, format = "f", digits = 4)
}

## ---- pooling across benchmarks ---------------------------------------------------
#
# In the 2PL behind Epoch's ECI, logit accuracy on benchmark b is
# alpha_b * (C - D_b) with C a common capability index, so every frontier
# coefficient on the logit scale is alpha_b times the corresponding slope of C.
# Dividing by alpha_b maps each benchmark's estimate onto the SAME capability
# scale, where averaging is legitimate, and the 2PL's own information weights
# are w_b = alpha_b^2, giving per term
#
#     theta_k = sum_b alpha_b beta_kb / sum_b alpha_b^2
#
# in ECI points per unit regressor. The pooled cost-decline ratio
# -theta_t / theta_x is then -sum(alpha beta_t) / sum(alpha beta_x): the pooled
# summary is what the pooled slopes imply, not a separate aggregate. The sigma
# parameters are logit-scale too, but log sigma_u is a LOG of one, so it converts
# by subtracting log alpha_b, and its time slope is already scale-free; both pool
# with the plain w_b weights.
#
# Benchmark fits are independent, so pooled variances are the correspondingly
# weighted sums of per-benchmark pieces, and the decline's delta method runs on
# the pooled 2x2 block. For the envelope and the frontier logit the alpha^2
# weights borrow an information interpretation those fits cannot support -- no
# likelihood -- so their pooled figures are mechanical averages, point estimates
# like the rest of those tables.
pooled_col <- function(key, tt, fitlist) {
  cv <- lapply(fitlist, coef_vcov)
  bs <- names(fitlist)
  disp <- vapply(TERMS, function(x) x$t, character(1))
  present <- disp[disp %in% unique(unlist(lapply(cv, function(x) names(x$b))))]

  rows <- lapply(present, function(k) {
    # a benchmark contributes where the term was estimated; an aliased NA drops
    # out and the weights renormalise over the rest
    use <- bs[vapply(bs, function(b)
      k %in% names(cv[[b]]$b) && is.finite(cv[[b]]$b[[k]]), logical(1))]
    if (!length(use)) return(NULL)
    a  <- ALPHA[use]
    cc <- if (grepl("^logsig_", k)) a^2 / sum(a^2) else a / sum(a^2)
    off <- if (k == "logsig_(Intercept)") -sum(a^2 / sum(a^2) * log(a)) else 0
    est <- sum(cc * vapply(use, function(b) cv[[b]]$b[[k]], numeric(1))) + off
    haveV <- all(vapply(use, function(b)
      !is.null(cv[[b]]$V) && k %in% rownames(cv[[b]]$V), logical(1)))
    se <- if (haveV)
      sqrt(sum(cc^2 * vapply(use, function(b) cv[[b]]$V[k, k], numeric(1))))
    else NA_real_
    data.frame(term = k, est = est, se = se, stringsAsFactors = FALSE)
  })
  es <- do.call(rbind, rows)

  decline <- NULL
  if (tt == "lin" && all(c("lncost", "tc") %in% es$term)) {
    need <- c("lncost", "tc")
    b2 <- setNames(es$est[match(need, es$term)], need)
    V2 <- NULL
    if (all(vapply(bs, function(b) !is.null(cv[[b]]$V), logical(1)))) {
      W  <- sum(ALPHA[bs]^2)
      V2 <- matrix(0, 2, 2, dimnames = list(need, need))
      for (b in bs) V2 <- V2 + (ALPHA[[b]] / W)^2 * cv[[b]]$V[need, need]
    }
    decline <- decline_from(b2, V2)
  }

  list(bench = "pooled", spec = tt, head = "Pooled (ECI pts)", es = es,
       n = as.integer(sum(vapply(bs, function(b)
         n_obs(key, b, fitlist[[b]]), numeric(1)))),
       decline = decline, test = NULL)
}

# The S family's profiled lambdas, cached per benchmark as the starting point
# for the slower families' profiles (S runs first in MODELS, and its glm inner
# loop makes its profile the cheap one). A start, not a constraint: each family
# still profiles its own lambdas.
BC_LSTART <- new.env(parent = emptyenv())

# Build the whole grid for one model: a list of columns, each carrying its
# estimates, N and (for quad and bc) the test; the pooled pair comes last.
build_model <- function(key) {
  cols <- list()
  grid <- fit_grid(key)
  for (b in benches) {
    fits <- lapply(grid, `[[`, b)
    fits$bc <- fit_bc(key, d[d$benchmark == b, ],
                      lambda_start = BC_LSTART[[b]] %||% c(0, 1))
    if (key == "S") BC_LSTART[[b]] <- unname(attr(fits$bc, "bc_lambda"))
    tsts <- list(quad = quad_test(key, b, fits$lin, fits$quad),
                 bc   = bc_test(key, fits$lin, fits$bc))
    for (tt in c(names(TIME_FORMS), "bc")) {
      cols[[length(cols) + 1]] <- list(
        bench = b, spec = tt,
        # the benchmark name only; the specifications are not labelled per column
        # but read off the pattern -- within each trio the quadratic is the one
        # carrying the second-order rows, the Box-Cox the one carrying the BC
        # rows, and the order is always linear, quadratic, Box-Cox
        head = LABELS[[b]],
        es = if (tt == "bc") est_se_bc(fits$bc) else est_se(fits[[tt]]),
        n = n_obs(key, b, fits[[tt]]),
        decline = if (tt == "lin") cost_decline(fits[[tt]]) else NULL,
        test = tsts[[tt]])
    }
  }
  for (tt in names(TIME_FORMS))
    cols[[length(cols) + 1]] <- pooled_col(key, tt, grid[[tt]])
  cols
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

html_table <- function(key, label, cols) {
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
    o <- c(o, sprintf('<th colspan="%d"%s>%s</th>', grp$span[i],
                      if (i > 1) ' class="grp"' else '',
                      html_escape(grp$head[i])))
  o <- c(o, '</tr></thead><tbody>')

  starts <- group_starts(cols)
  cls <- function(j) if (j %in% starts) ' class="grp"' else ''

  for (x in trm) {
    o <- c(o, sprintf('<tr><td>%s</td>', x$h))
    for (j in seq_along(cols))
      o <- c(o, sprintf('<td%s>%s</td>', cls(j), cell(cols[[j]], x$t, "est")))
    o <- c(o, '</tr>')
    if (has_se) {
      o <- c(o, '<tr class="se"><td></td>')
      for (j in seq_along(cols))
        o <- c(o, sprintf('<td%s>%s</td>', cls(j), cell(cols[[j]], x$t, "se")))
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
      o <- c(o, sprintf('<td%s>%s</td>', cls(j), decl_cell(cols[[j]], "est")))
    o <- c(o, '</tr>')
    if (has_decl_se) {
      o <- c(o, '<tr class="se"><td></td>')
      for (j in seq_along(cols))
        o <- c(o, sprintf('<td%s>%s</td>', cls(j), decl_cell(cols[[j]], "se")))
      o <- c(o, '</tr>')
    }
  }

  o <- c(o, sprintf('<tr%s><td>N</td>', if (has_decl) '' else ' class="gap"'))
  for (j in seq_along(cols))
    o <- c(o, sprintf('<td%s>%d</td>', cls(j), cols[[j]]$n))
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
                    length(cols) + 1, notes_html(key, kind)),
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

rtf_table <- function(key, label, cols) {
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

  o <- c("{\\rtf1\\ansi\\ansicpg1252\\deff0",
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
                        vapply(cols, function(cl) cell(cl, x$t, "est"), character(1))),
                      widths, sep = sep_data))
    if (has_se)
      o <- c(o, rtf_row(c("", vapply(cols, function(cl) cell(cl, x$t, "se"),
                                     character(1))), widths, sep = sep_data))
  }
  has_decl <- any(vapply(cols, function(cl) !is.null(cl$decline), logical(1)))
  has_decl_se <- any(vapply(cols, function(cl)
    !is.null(cl$decline) && !is.na(cl$decline$se), logical(1)))
  if (has_decl) {
    o <- c(o, rtf_row(c("cost drop, %/qtr",
                        vapply(cols, function(cl) decl_cell(cl, "est"), character(1))),
                      widths, top = TRUE, sep = sep_data))
    if (has_decl_se)
      o <- c(o, rtf_row(c("", vapply(cols, function(cl) decl_cell(cl, "se"),
                                     character(1))), widths, sep = sep_data))
  }
  o <- c(o, rtf_row(c("N", vapply(cols, function(cl) as.character(cl$n), character(1))),
                    widths, top = !has_decl, sep = sep_data))
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
         sprintf("\\pard\\ql{\\fs16 %s}\\par", rtf_escape(notes_plain(key, kind))),
         "}")
  o
}

## ---- the notes that keep each table honest -------------------------------------

notes_plain <- function(key, kind) {
  se <- switch(key,
    S = "Robust (HC1) standard errors: the quasibinomial objective is a quasi-likelihood, so inverse-Hessian errors would be wrong.",
    paretologit = paste("No standard errors: the response is the empirical Pareto frontier P_t(c) sampled at the nodes of the",
                        "envelope's fixed cost-date grid, and grid nodes are not observations. N is the nodes at which the",
                        "frontier is defined; they resample a few dozen staircase corners, so N measures resolution, not",
                        "information."),
    A = , B = "Robust (sandwich) standard errors, as the SFA objective is a quasi-likelihood in the Papke-Wooldridge sense.",
    envelope = paste("No standard errors: the envelope is the solution of a constrained optimisation with no likelihood behind it,",
                     "so none are reported rather than invented. N is the runs it must clear."))
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
                if (key %in% c("envelope", "paretologit")) "" else
                  paste("Standard error by the delta method on the same robust covariance, taken on the transformed",
                        "quantity. Being symmetric it can reach past 100% where the estimated drop is near total."))
  pool <- paste(
    "Pooled maps the benchmarks onto the common scale of Epoch's ECI (Epoch Capabilities Index) 2PL, in which",
    "logit accuracy on benchmark b is alpha_b (C - D_b): each slope over alpha_b estimates the same",
    "capability-scale slope, and the pooled coefficient is their information-weighted (w = alpha^2) average,",
    "sum(alpha b) / sum(alpha^2), in ECI points per unit regressor. The discriminations (estimated_slope_scaled",
    sprintf("in data/edi_scores.csv, downloaded from https://epoch.ai/data/edi_scores.csv on 2026-08-20) are %s,",
            paste(sprintf("%s %.3f", benches, ALPHA[benches]), collapse = ", ")),
    sprintf("so one logit is worth %.1f-%.1f ECI points and pooled slopes read several times larger than the",
            min(1 / ALPHA[benches]), max(1 / ALPHA[benches])),
    "logit-scale columns beside them.",
    "The pooled cost drop is the same transform of the pooled slopes, -sum(alpha b_time) / sum(alpha b_ln cost)",
    "inside it. The pooled intercept averages capability net of difficulty at each benchmark's own reference",
    "date, so unlike the slopes it carries no clean interpretation, and pooled N sums the benchmark columns.",
    if (key %in% c("A", "B"))
      paste("log sigma_u, the log of a logit-scale spread, converts as log sigma_u - log alpha_b before averaging;",
            "its time slope is scale-free and pools directly. Fits are independent across benchmarks, so pooled",
            "standard errors sum the per-benchmark covariance pieces.")
    else if (key == "S")
      "Fits are independent across benchmarks, so pooled standard errors sum the per-benchmark covariance pieces."
    else
      paste("For this model the alpha^2 weights borrow an information interpretation the fit cannot support --",
            "there is no likelihood behind it -- so the pooled figures are mechanical averages."),
    "One scale caveat: gpqa accuracy is rescaled for its 0.25 guessing floor before the logit (prepare_data.R),",
    "a scale on which Epoch's alpha_gpqa was not necessarily estimated.")
  bc <- paste(
    "BC columns are a Box-Cox alternative to the quadratic: the logit is linear in phi(cost), phi(time) and",
    "their product, with phi(x; lambda) = (x^lambda - 1)/lambda (log at lambda = 0) applied to LEVEL cost per",
    "task and to years since mid-2020, GPT-3's release -- so the BC intercept is the fit at $1 per task in",
    "mid-2021, where both transforms vanish. phi is increasing whatever lambda is, so the surface is monotone in",
    "cost at every date -- the quadratic's bending back toward the data cannot happen -- while the product term",
    "still allows the cost slope one sign change over time.",
    "lambda_cost and lambda_time are estimated per benchmark by profiling",
    switch(key, S = "the quasibinomial deviance,", A = , B = "the likelihood,",
           paretologit = "the grid deviance,",
           envelope = "the envelope's mean fitted height,"),
    "and are reported without standard errors: the profile provides none, and",
    if (key %in% c("envelope", "paretologit"))
      "this model reports none anywhere."
    else
      paste("the coefficient standard errors are conditional on the profiled lambdas --",
            "they carry no lambda uncertainty."),
    if (key %in% c("A", "B"))
      paste("The BC specification nests the linear one (lambda_cost = 0, lambda_time = 1, no product term) --",
            "its LR row tests exactly those restrictions -- but not the quadratic.")
    else
      "The BC specification nests the linear one (lambda_cost = 0, lambda_time = 1, no product term) but not the quadratic.",
    "fm13's five months of data sit 5.7-6.1 years from the origin, over which every lambda_time fits alike, so",
    "lambda_time is fixed at 1 there rather than estimated; elsewhere it is weakly identified and best read as",
    "a shape the data tolerates rather than demands -- a lambda_time of exactly 3 or -2 sits on the edge of the",
    "search box, the flat profile having run to the wall. The pooled pair covers the linear and quadratic",
    "specifications only: per-benchmark lambdas put the BC slopes on different transforms of cost and time,",
    "leaving no common scale for the alpha-weighted average to land on.")
  paste("Time is measured in years and centered within benchmark, so the intercept is the frontier at each benchmark's own",
        "reference date (the BC columns instead use uncentered years since mid-2020).", decl, se, tst, bc, pool)
}

notes_html <- function(key, kind) {
  x <- notes_plain(key, kind)
  x <- html_escape(x)
  gsub("ln cost^2", "ln cost<sup>2</sup>", gsub("time^2", "time<sup>2</sup>", x,
       fixed = TRUE), fixed = TRUE)
}

## ---- emit ----------------------------------------------------------------------

dir.create(out_path("tables"), showWarnings = FALSE, recursive = TRUE)
for (m in MODELS) {
  cols <- build_model(m$key)
  writeLines(html_table(m$key, m$label, cols),
             out_path("tables", sprintf("regression_%s.html", m$key)))
  writeLines(rtf_table(m$key, m$label, cols),
             out_path("tables", sprintf("regression_%s.rtf", m$key)))
  cat("wrote regression_", m$key, ".html / .rtf\n", sep = "")
}
