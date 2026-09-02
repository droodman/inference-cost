# Console report: quarterly cost-decline rates from the cost-direction dual
# fits (cost_frontier.R), beside the accuracy-direction fits they mirror and
# the nonparametric staircase check. All linear specifications; every column
# targets the same estimand, the fall per quarter in the cost of fixed
# frontier performance.
#
# Columns, in the order printed:
#
#   data start the benchmark's earliest run, as YYYY-MM: how much history a
#              row's rates rest on, and the first thing to check when one
#              benchmark disagrees with the rest
#   staircase  pareto_decline_qtr, the model-free average over a QUARTER
#              horizon (blank where the benchmark spans less than that)
#   grid OLS   fit_lncost_grid: ln C_a(t) sampled on the (logit acc, date)
#              grid, OLS through the samples; rate = its tc coefficient
#   BC surf    the Box-Cox grid OLS surface's own INSTANTANEOUS d lnC/dt,
#              averaged over the staircase check's lattice and expressed
#              quarterly -- what the curved model claims the rate is, node by
#              node. No horizon is differenced: the record needs one because
#              it is a step function, a fitted surface does not. On a linear
#              fit this reproduces the grid OLS column exactly
#   grid+env   fit_lncost_grid_env: the grid OLS objective under the cost
#              envelope's constraints; rate = its tc coefficient
#   OLS runs   lm of ln cost on logit acc and date over all positive-accuracy
#              runs -- model S's reverse regression, the typical run's cost
#   SFA cost   fit_cost_sfa: stochastic cost frontier, half-normal per
#              model x effort; rate = its tc coefficient
#   par.logit  fit_pareto_logit (accuracy direction); rate = -b_t/b_x
#   pl env     fit_pareto_logit_env (accuracy direction); rate = -b_t/b_x
#
# Reading the spread: the frontier-per-se columns (staircase, grid OLS, BC
# surf, grid+env, par.logit, pl env) all target the record's decline, and differ by
# loss direction and weighting; if the surface were truly logit-linear they
# would agree. SFA cost and OLS runs instead follow the model-effort CLOUD at
# fixed accuracy, whose dense cheap edge grows dearer as expensive reasoning
# configurations arrive -- so negative "declines" there are the cloud
# drifting up while the record collapses, the sharpest statement in this file
# of how much of the story is frontier movement rather than typical-run
# movement. See cost_frontier.R's fit_cost_sfa block.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("fit_store.R")   # every fit comes from the shared store

d <- load_runs()

# -b_t/b_x of an accuracy-direction fit, as a quarterly percentage
rate_acc <- function(fit) {
  cf <- coef(fit)
  100 * (1 - exp(-cf[["tc"]] / cf[["lncost"]] / 4))
}

# One row per benchmark, collected once so the console print and the saved
# table cannot disagree. `bc_surface` is the smoothed twin of the staircase
# column: pareto_decline_qtr's aggregate on the SAME nodes and horizon, read
# off the Box-Cox Pareto-grid cost fit's surface (surface_decline_qtr,
# cost_frontier.R) -- the smoothed-versus-actual comparison, everything else
# held fixed. The other columns are the linear fits' single coefficients.
rate_rows <- do.call(rbind, lapply(bench_levels(d$benchmark), function(b) {
  s <- d[d$benchmark == b, ]
  np <- pareto_decline_qtr(s)
  ms <- surface_decline_qtr(store_cost_bc("costgridols")[[b]], s)
  data.frame(
    bench      = b,
    benchmark  = unname(LABELS[[b]]),
    # when the benchmark's history starts: how much data each rate rests on,
    # which is the first thing to check when one row disagrees with the rest
    start      = format(min(s$releasedate), "%Y-%m"),
    # alpha_b^2, the 2PL information weight the pooled row uses. Reported for
    # EVERY benchmark, pooled or not: it is a property of the benchmark rather
    # than of the pool, and seeing it for the secondaries is how one judges
    # what adding them would do -- several carry MORE information per
    # observation than any primary does.
    alpha2     = ALPHA[[b]]^2,
    staircase  = if (is.null(np)) NA_real_ else np$pct_qtr,
    bc_surface = if (is.null(ms)) NA_real_ else ms$pct_qtr,
    grid_ols   = cost_decline_qtr(store_cost("costgridols")$lin[[b]]),
    grid_env   = cost_decline_qtr(store_cost("costgridolsenv")$lin[[b]]),
    sfa_cost   = cost_decline_qtr(store_cost("costsfa")$lin[[b]]),
    ols_runs   = cost_decline_qtr(store_cost("costols")$lin[[b]]),
    par_logit  = rate_acc(store_grid("paretologit")$lin[[b]]),
    pl_env     = rate_acc(store_grid("paretologitenv")$lin[[b]]))
}))

## ---- pooling the primary benchmarks onto the ECI scale ---------------------------
#
# Every column here is a rate of decline of the log cost of FIXED performance,
# and all three of the constructions above express it the same way,
# r = 100(1 - exp(g/4)) with g the annual log-cost change. So g = 4 log(1 -
# r/100) recovers the underlying rate from any column, pooling happens in g,
# and the same formula converts back -- no column needs special handling to
# get onto a common footing.
#
# WEIGHTS. Fixing accuracy on benchmark b fixes a capability level, so each g_b
# already answers "how fast does the cost of a fixed capability fall?" and needs
# no ECI rescaling -- what ECI supplies is how much each benchmark should count.
# In the 2PL the Fisher information carried by benchmark b is proportional to
# alpha_b^2 (regression_tables.R), so the pool is the alpha^2-weighted mean.
# Unlike the coefficient pooling there, this cannot be inverse-variance
# weighted: four of these columns are grid or nonparametric quantities with no
# standard error at all.
#
# THE TWO ACCURACY COLUMNS ARE DIFFERENT. Their rate is a RATIO of logit-scale
# slopes, -b_t/b_x, in which alpha_b cancels within a benchmark but not across
# them. regression_tables.R's pooled_col already settles this case: pool the
# slopes onto the ECI scale first, theta_k = sum(alpha b_k)/sum(alpha^2), and
# take the ratio of the pooled slopes, so the pooled summary is what the pooled
# slopes imply. That is used verbatim here, which is also why these two entries
# agree with the pooled decline in the regression tables.
#
# Primary benchmarks only, matching pooled_col; the newer benchmarks keep their
# own rows but stay out of the pool.
PB <- intersect(PRIMARY_BENCHES, bench_levels(d$benchmark))
stopifnot(!anyNA(ALPHA[PB]))

qtr_to_g <- function(r) 4 * log(1 - r / 100)     # %/qtr -> annual log change
g_to_qtr <- function(g) 100 * (1 - exp(g / 4))

pool_qtr <- function(r) {                        # alpha^2-weighted, r over PB
  ok <- !is.na(r)
  if (!any(ok)) return(NA_real_)
  w <- ALPHA[PB][ok]^2
  g_to_qtr(sum(w * qtr_to_g(r[ok])) / sum(w))
}

# -sum(alpha b_t) / sum(alpha b_x): the ratio of the ECI-pooled slopes
pool_ratio <- function(fits) {
  a  <- ALPHA[PB]
  bt <- vapply(PB, function(b) coef(fits[[b]])[["tc"]], 0)
  bx <- vapply(PB, function(b) coef(fits[[b]])[["lncost"]], 0)
  g_to_qtr(-sum(a * bt) / sum(a * bx))
}

pri <- rate_rows[match(PB, rate_rows$bench), ]
pooled <- data.frame(
  bench = "pooled", benchmark = "Pooled, primary (ECI-weighted)", start = "",
  # the pool's denominator: summed over the PRIMARY benchmarks only, so each
  # primary's share of the pooled estimate is its own cell over this one
  alpha2 = sum(ALPHA[PB]^2),
  staircase  = pool_qtr(pri$staircase),
  bc_surface = pool_qtr(pri$bc_surface),
  grid_ols   = pool_qtr(pri$grid_ols),
  grid_env   = pool_qtr(pri$grid_env),
  sfa_cost   = pool_qtr(pri$sfa_cost),
  ols_runs   = pool_qtr(pri$ols_runs),
  par_logit  = pool_ratio(store_grid("paretologit")$lin),
  pl_env     = pool_ratio(store_grid("paretologitenv")$lin))

pc <- function(x) ifelse(is.na(x), "", sprintf("%.1f%%", x))

# The HTML gets a real minus sign. sprintf leaves an ASCII hyphen, which is
# narrower than the digits it sits against and reads as a dash rather than a
# sign. The CSV keeps the hyphen -- it has to parse back as a number -- and
# so does the console, where the entity would print literally.
pc_html <- function(x) sub("^-", "&minus;", pc(x))

cat("cost decline per quarter at fixed accuracy\n")
cat("(positive = cheaper; frontier columns first, then the run-cloud pair)\n")
cat(sprintf("%-6s %8s %8s %10s %9s %9s %9s | %9s %9s | %10s %8s\n",
            "bench", "from", "alpha^2", "staircase", "grid OLS", "BC surf",
            "grid+env", "OLS runs", "SFA cost", "par.logit", "pl env"))
for (i in seq_len(nrow(rate_rows) + 1)) {
  r <- if (i <= nrow(rate_rows)) rate_rows[i, ] else pooled
  if (i > nrow(rate_rows)) cat(strrep("-", 96), "\n")
  cat(sprintf(
    "%-6s %8s %8.4f %10s %9s %9s %9s | %9s %9s | %10s %8s\n",
    r$bench, r$start, r$alpha2, pc(r$staircase), pc(r$grid_ols),
    pc(r$bc_surface),
    pc(r$grid_env), pc(r$ols_runs), pc(r$sfa_cost), pc(r$par_logit),
    pc(r$pl_env)))
}

## ---- the same comparison, saved as a table (HTML + CSV) ------------------------------

dir.create(out_path("tables"), showWarnings = FALSE, recursive = TRUE)

csv <- rbind(rate_rows, pooled)[, c("benchmark", "start", "alpha2",
                     "staircase", "grid_ols", "bc_surface", "grid_env",
                     "ols_runs", "sfa_cost", "par_logit", "pl_env")]
names(csv) <- c("benchmark", "data_start", "alpha_sq", "staircase_pct_qtr",
                "grid_ols_pct_qtr", "bc_surface_pct_qtr",
                "grid_ols_env_pct_qtr", "ols_runs_pct_qtr",
                "sfa_cost_pct_qtr", "pareto_logit_pct_qtr",
                "pareto_logit_env_pct_qtr")
write.csv(csv, out_path("tables", "rate_comparison.csv"), row.names = FALSE)
cat("\nwrote rate_comparison.csv\n")

esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  gsub("<", "&lt;", x, fixed = TRUE)
}
o <- c('<!DOCTYPE html>', '<html lang="en"><head><meta charset="UTF-8" />',
       '<title>Rate comparison</title>',
       '<style>',
       'body{font-family:"Segoe UI",Arial,sans-serif;margin:12px;color:#1d1d1d;background:#fcfcfb}',
       'h1{font-size:1.25em;margin:0 0 4px 0}',
       'p.sub{color:#5e5e5e;margin:0 0 16px 0;font-size:.9em}',
       'table{border-collapse:collapse;background:#fcfcfb;font-size:.9em}',
       # Centred by default; the rules below left-align the label column and
       # the footnote, and the headers restate centring explicitly.
       'th,td{padding:3px 12px;text-align:center;white-space:nowrap}',
       # Only the FIRST header row's leading cell is left-aligned. The second
       # row's leading cell is a spanning group label -- the rowspan cells
       # above belong to row 1 -- so a bare th:first-child rule caught it and
       # shoved it left instead of centring it over its pair of columns.
       'thead tr:first-child th:first-child{text-align:left}',
       'thead th{border-bottom:1px solid #1d1d1d;font-weight:600;vertical-align:bottom;text-align:center}',
       # Header labels wrap inside a width-capped block rather than in the
       # cell itself: th inherits the nowrap above, and max-width on a th is
       # unreliable under auto table layout, while a block inside it is not.
       # Capping the block keeps a long label from setting its column's
       # width -- the data cells do -- and margin auto keeps the block
       # centred under the spanning group headers.
       'thead th .hd{white-space:normal;overflow-wrap:break-word;max-width:7em;margin:0 auto}',
       # No width cap on a SPANNING header: the cap exists to stop a long
       # leaf label from setting its one column's width, while a group label
       # has its whole span to sit in -- "Accuracy models" was wrapping for
       # no reason. Uncapped, the block also fills the cell, which is what
       # the group rule below hangs on.
       'thead th[colspan] .hd{max-width:none}',
       # The first row's group rules ride on that block rather than on the
       # cell, so each stops at the cell's padding instead of running into
       # its neighbour: "Cost models" and "Accuracy models" read as two
       # rules with a gap between them, not one continuous line. The row
       # BELOW keeps its rule on the cells, where it should stay unbroken --
       # it is the header/body separator.
       'thead tr:first-child th[colspan]{border-bottom:0}',
       'thead tr:first-child th[colspan] .hd{border-bottom:1px solid #1d1d1d;padding-bottom:3px}',
       'tbody th:first-child,tbody td:first-child{text-align:left}',
       'tfoot td{border-top:1px solid #1d1d1d;font-size:.92em;color:#5e5e5e;',
       '  text-align:left;white-space:normal;padding-top:8px;max-width:900px}',
       # The order row leads the body; its extra bottom padding is the
       # half-row gap that keeps the specification labels from reading as
       # one more benchmark.
       '.order-row th,.order-row td{color:#4d4d4d;padding-bottom:1em}',
       # The pooled row closes the body: a rule and a half-row gap above it,
       # and an italic label, so it reads as a summary of the rows rather
       # than as another benchmark among them.
       '.pooled-row td{border-top:1px solid #1d1d1d;padding-top:1em}',
       '.pooled-row td:first-child{font-style:italic}',
       '</style></head><body>',
       '<h1>Average quarterly rate of decline in cost of given accuracy</h1>',
       '<table><thead><tr>',
       # Benchmark alone keeps its bare label: it is the left-aligned first
       # column, where a centred block would drift off the column's edge.
       '<th rowspan="2">Benchmark</th>',
       '<th rowspan="2"><div class="hd">Data start</div></th>',
       '<th rowspan="2"><div class="hd">&alpha;<sup>2</sup> weight</div></th>',
       '<th rowspan="2"><div class="hd">Non-parametric</div></th>',
       '<th colspan="5"><div class="hd">Cost models</div></th>',
       '<th colspan="2"><div class="hd">Accuracy models</div></th>',
       '</tr><tr>',
       '<th colspan="2"><div class="hd">Model frontier</div></th>',
       '<th><div class="hd">Model frontier, require envelopment</div></th>',
       '<th><div class="hd">Model all data</div></th>',
       '<th><div class="hd">Stochastic frontier analysis</div></th>',
       '<th><div class="hd">Model frontier</div></th>',
       '<th><div class="hd">Model frontier, require envelopment</div></th>',
       '</tr></thead><tbody>')
order_cells <- c("<i>Model order</i>", "", "", "", "Linear", "Box-Cox",
                 "Linear", "Linear", "Linear", "Linear", "Linear")
o <- c(o, sprintf('<tr class="order-row">%s</tr>',
                  paste0(sprintf('<td>%s</td>', order_cells), collapse = '')))
for (i in seq_len(nrow(rate_rows) + 1)) {
  r <- if (i <= nrow(rate_rows)) rate_rows[i, ] else pooled
  cells <- c(esc(r$start), sprintf("%.4f", r$alpha2),
             pc_html(unlist(r[, c("staircase", "grid_ols", "bc_surface",
                                  "grid_env", "ols_runs", "sfa_cost",
                                  "par_logit", "pl_env")])))
  # the pooled row is a summary, not another benchmark: a rule above it and
  # the label in italics keep it from reading as a twelfth row of data
  o <- c(o, sprintf('<tr%s><td>%s</td>%s</tr>',
                    if (i > nrow(rate_rows)) ' class="pooled-row"' else '',
                    esc(r$benchmark),
                    paste0(sprintf('<td>%s</td>', cells), collapse = '')))
}
o <- c(o, '</tbody><tfoot><tr><td colspan="11">',
       paste("All numbers are estimates of the average quarterly drop in the cost",
             "of a given level of accuracy on a given benchmark, over the years of available data.",
             "The \"non-parametric\" values are averages over even 100&times;100 grids that span the benchmark's release date and accuracy ranges,",
             "with the latter scaled within the state of the art (SOTA) for that benchmark at each given time. At each point, the lowest cost",
             "of performance at least as good one QUARTER later is found and divided by the initial cost. The geometric mean of the ratios is the quarterly decline rate.",
             "Over so short a horizon most levels do not move at all, and each such node enters the geometric mean as a ratio of 1, so the rate is a minority of real drops",
             "averaged in with a majority of unchanged records.",
             "Results in all cost model columns are coefficients on release year as an explanator for log cost, again reexpressed quarterly; except that the nonlinear Box-Cox variant,",
             "having no single such coefficient, is summarized by the average over the same grid of the instantaneous rate of decline the fitted surface implies at each point.",
             "That average differences nothing over time &mdash; the record needs a horizon because it is a step function, a fitted surface does not &mdash; and on a linear fit it reduces",
             "exactly to the coefficient in the column beside it. In the accuracy model columns accuracy",
             "is the dependent variable and log cost as an explanatory variable, so the rates are extracted as -b_t/b_x.",
             "\"Model frontier\" means modeling the empirical frontier as realized at a grid of points, either with ordinary least squares (OLS; for log cost) or with a logit link (for accuracy).",
             "\"Model frontier, require envelopment\" means the same, but with the constraint that the fitted surface is never above any data point (for cost) or below (for accuracy).",
             "\"Model all data\" means modeling all runs, not just the frontier, with OLS. \"Stochastic frontier analysis\" models the frontier and the distribution of runs around it",
             "with a half-normal distribution of inefficiency.",
             "The final row pools the five primary benchmarks &mdash; AIME, Chess Puzzles, FrontierMath tiers 1&ndash;3, GPQA Diamond and Mystery Game Puzzles &mdash; onto the common",
             "capability scale of Epoch's ECI (Epoch Capabilities Index), whose 2PL writes logit accuracy on benchmark b as &alpha;<sub>b</sub>(C &minus; D<sub>b</sub>) for a shared",
             "capability C. Holding accuracy fixed on a benchmark holds C fixed, so each column's rate already answers the same question and needs no rescaling; what ECI supplies is",
             "how much each benchmark counts, and the 2PL's information weight is &alpha;<sub>b</sub><sup>2</sup> &mdash; the &alpha;<sup>2</sup> column, shown for every",
             "benchmark because it is a property of the benchmark and not of the pool. Each primary's share of the pooled estimate is its own &alpha;<sup>2</sup> over the final",
             "row's, which sums &alpha;<sup>2</sup> across the five primaries alone: AIME carries 36.7% of the pool, FrontierMath tiers 1&ndash;3 21.2%, Mystery 20.9%, and Chess and",
             "GPQA 10.6% each. Note that several SECONDARY benchmarks carry a larger &alpha;<sup>2</sup> than any primary, so they would dominate the pool if added &mdash; their",
             "exclusion rests on their short histories, not on low information per observation. Every column but the last two is pooled as the",
             "&alpha;<sup>2</sup>-weighted mean of the annual log-cost change, converted back to a quarterly rate. Inverse-variance weighting is unavailable here because the grid and",
             "non-parametric columns carry no standard errors. The two accuracy-model rates are ratios &minus;b_t/b_x in which &alpha;<sub>b</sub> cancels within a benchmark but not",
             "across them, so they instead take the ratio of the ECI-pooled slopes, &minus;&Sigma;&alpha;b_t / &Sigma;&alpha;b_x, matching the pooled decline in the regression tables.",
             "For the envelope and grid fits the &alpha;<sup>2</sup> weights borrow an information interpretation those fits cannot support &mdash; there is no likelihood behind them",
             "&mdash; so their pooled entries are mechanical averages."),
       '</td></tr></tfoot></table></body></html>')
writeLines(o, out_path("tables", "rate_comparison.html"))
cat("wrote rate_comparison.html\n")

cat("\ncost-direction fit details\n")
cat(sprintf("%-6s %28s %28s\n", "", "grid OLS", "grid OLS + envelope"))
cat(sprintf("%-6s %9s %8s %9s %9s %8s %9s\n", "bench",
            "$/logit", "tc", "nodes", "$/logit", "tc", "touch"))
for (b in bench_levels(d$benchmark)) {
  fo <- store_cost("costgridols")$lin[[b]]
  fe <- store_cost("costgridolsenv")$lin[[b]]
  cat(sprintf("%-6s %9.3f %8.3f %9d %9.3f %8.3f %9.4f\n", b,
              coef(fo)[["la"]], coef(fo)[["tc"]], attr(fo, "n_grid"),
              coef(fe)[["la"]], coef(fe)[["tc"]], fe$slack_envelope))
}
cat("\n$/logit is the la coefficient: log dollars per logit of accuracy,\n")
cat("the surface's steepness in the cost direction. touch is the constrained\n")
cat("fit's minimum run slack in log dollars; 0 means it touches a run.\n")