# Console report: quarterly cost-decline rates from the cost-direction dual
# fits (cost_frontier.R), beside the accuracy-direction fits they mirror and
# the nonparametric staircase check. Every column targets the same estimand,
# the fall per quarter in the cost of fixed frontier performance.
#
# TWO tables (console blocks and HTML alike), both led by the data-start and
# model-free columns:
#
#   data start the benchmark's earliest run, as YYYY-MM: how much history a
#              row's rates rest on, and the first thing to check when one
#              benchmark disagrees with the rest
#   model-free pareto_decline_qtr, the model-free average over a QUARTER
#              horizon (blank where the benchmark spans less than that);
#              printed as "staircase" in the console
#
# Table 1, the frontier and all-data cost models, each as a linear /
# Box-Tidwell pair:
#
#   grid OLS   fit_lncost_grid: ln C_a(t) sampled on the (logit acc, date)
#              grid, OLS through the samples; rate = its tc coefficient
#   BC surf    the Box-Cox grid OLS surface's own INSTANTANEOUS d lnC/dt,
#              averaged over the staircase check's lattice and expressed
#              quarterly -- what the curved model claims the rate is, node by
#              node. No horizon is differenced: the record needs one because
#              it is a step function, a fitted surface does not. On a linear
#              fit this reproduces the grid OLS column exactly
#   grid+env   fit_lncost_grid_env: the grid OLS objective under the cost
#              envelope's constraints; rate = its tc coefficient. Its BC twin
#              is summarized like BC surf
#   OLS runs   lm of ln cost on logit acc and date over all positive-accuracy
#              runs -- model S's reverse regression, the typical run's cost.
#              Its BC twin averages the surface's rate over the RUNS
#              (cloud_decline_qtr), the model's own population
#
# Table 2, the SFA pair and the accuracy models, again as linear / BC pairs:
#
#   SFA        fit_cost_sfa: stochastic cost frontier, half-normal per
#              model x effort; rate = its tc coefficient
#   SFAb       the same with log sigma_u linear in date (costsfab)
#   (BC twins of both average the frontier's rate over the runs, as OLS BC)
#   par.logit  fit_pareto_logit (accuracy direction); rate = -b_t/b_x, and
#              its BC twin the ratio of node-averaged slopes
#              (pareto_bc_decline_qtr)
#   pl env     fit_pareto_logit_env, same pair of summaries
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
    staircase  = if (is.null(np)) NA_real_ else np$pct_qtr,
    bc_surface = if (is.null(ms)) NA_real_ else ms$pct_qtr,
    grid_ols   = cost_decline_qtr(store_cost("costgridols")$lin[[b]]),
    grid_env   = cost_decline_qtr(store_cost("costgridolsenv")$lin[[b]]),
    # the envelope's Box-Cox twin, summarized over the same lattice as the
    # BC surf column
    grid_env_bc = {
      me <- surface_decline_qtr(store_cost_bc("costgridolsenv")[[b]], s)
      if (is.null(me)) NA_real_ else me$pct_qtr
    },
    sfa_cost   = cost_decline_qtr(store_cost("costsfa")$lin[[b]]),
    # the SFA pair's Box-Tidwell twins, and the time-varying-sigma_u SFA:
    # cloud models, so their BC rates average the fitted surface's
    # instantaneous d lnC/dt over the runs, exactly as ols_runs_bc does
    sfa_bc     = cloud_decline_qtr(store_cost_bc("costsfa")[[b]], s),
    sfab       = cost_decline_qtr(store_cost("costsfab")$lin[[b]]),
    sfab_bc    = cloud_decline_qtr(store_cost_bc("costsfab")[[b]], s),
    ols_runs   = cost_decline_qtr(store_cost("costols")$lin[[b]]),
    # the all-data Box-Cox twin: the cloud model's instantaneous rate
    # averaged over the runs themselves (cloud_decline_qtr, cost_frontier.R)
    ols_runs_bc = cloud_decline_qtr(store_cost_bc("costols")[[b]], s),
    par_logit  = rate_acc(store_grid("paretologit")$lin[[b]]),
    # the Box-Tidwell variant of the Pareto-logit fit: no single -b_t/b_x
    # ratio exists, so the surface's instantaneous d lnC/dt is averaged over
    # the accuracy grid's defined nodes (pareto_bc_decline_qtr,
    # boxcox_frontier.R) -- the accuracy mirror of the BC surf column
    par_lg_bc  = pareto_bc_decline_qtr(store_bc("paretologit")[[b]], s),
    pl_env     = rate_acc(store_grid("paretologitenv")$lin[[b]]),
    # the envelope-constrained accuracy fit's Box-Tidwell twin, the same
    # node-averaged-slopes ratio as par_lg_bc
    pl_env_bc  = pareto_bc_decline_qtr(store_bc("paretologitenv")[[b]], s))
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
# the pooling weight: alpha^2 * observed history * the record path's average
# p(1-p) -- the benchmark's integrated Fisher information about the shared
# capability path (pool_weights, envelope_frontier.R)
WPOOL <- pool_weights(d)[PB]

qtr_to_g <- function(r) 4 * log(1 - r / 100)     # %/qtr -> annual log change
g_to_qtr <- function(g) 100 * (1 - exp(g / 4))

pool_qtr <- function(r) {                        # alpha^2*T-weighted, over PB
  ok <- !is.na(r)
  if (!any(ok)) return(NA_real_)
  g_to_qtr(sum(WPOOL[ok] * qtr_to_g(r[ok])) / sum(WPOOL[ok]))
}

# the ratio of the ECI-pooled slopes: each slope converted to the ECI scale
# (b/alpha) and pooled with weight WPOOL, so WPOOL/alpha multiplies each raw
# slope and the common denominator cancels out of the ratio
pool_ratio <- function(fits) {
  wa <- WPOOL / ALPHA[PB]
  bt <- vapply(PB, function(b) coef(fits[[b]])[["tc"]], 0)
  bx <- vapply(PB, function(b) coef(fits[[b]])[["lncost"]], 0)
  g_to_qtr(-sum(wa * bt) / sum(wa * bx))
}

pri <- rate_rows[match(PB, rate_rows$bench), ]
# the primaries' pooling shares, spelled out in the footnote
share_txt <- paste(sprintf("%s %.1f%%", LABELS[PB], 100 * WPOOL / sum(WPOOL)),
                   collapse = ", ")
pooled <- data.frame(
  bench = "pooled", benchmark = "Pooled, primary (ECI-weighted)", start = "",
  staircase  = pool_qtr(pri$staircase),
  bc_surface = pool_qtr(pri$bc_surface),
  grid_ols   = pool_qtr(pri$grid_ols),
  grid_env   = pool_qtr(pri$grid_env),
  grid_env_bc = pool_qtr(pri$grid_env_bc),
  sfa_cost   = pool_qtr(pri$sfa_cost),
  sfa_bc     = pool_qtr(pri$sfa_bc),
  sfab       = pool_qtr(pri$sfab),
  sfab_bc    = pool_qtr(pri$sfab_bc),
  ols_runs   = pool_qtr(pri$ols_runs),
  ols_runs_bc = pool_qtr(pri$ols_runs_bc),
  par_logit  = pool_ratio(store_grid("paretologit")$lin),
  # no likelihood-scale slopes to pool for the BC variants: mechanical
  # weighted means, like the grid and nonparametric columns
  par_lg_bc  = pool_qtr(pri$par_lg_bc),
  pl_env     = pool_ratio(store_grid("paretologitenv")$lin),
  pl_env_bc  = pool_qtr(pri$pl_env_bc))

# The other summary: the UNWEIGHTED mean of the five primaries' annual log
# changes, converted back to a quarterly rate -- the equal-weights limit of
# random-effects pooling, transparent and free of the ECI calibration. Placed
# just above the ECI-weighted row so the two summaries read as a bracket.
# (A pooled-DATA regression row -- one fit over the stacked primaries with
# benchmark fixed effects -- lived here for a while and was dropped: its
# leverage weighting is a defensible precision story, but a delete-one-model
# jackknife put the per-benchmark sampling errors at 0.3-0.8 annual log
# units, within which every summary here agrees, and the row's shared-lambda
# Box-Cox cell answered a different question than the per-benchmark fits.
# The machinery survives in fit_pooled_cost / fit_pooled_acc_bs for the
# figures and any future return.)
savg <- function(r) {
  ok <- !is.na(r)
  if (!any(ok)) return(NA_real_)
  g_to_qtr(mean(qtr_to_g(r[ok])))
}
# the accuracy ratios are averaged as rates like every other column here --
# equal weights need no slope pooling, unlike the ECI row's pool_ratio
simple <- data.frame(
  bench = "s.avg", benchmark = "Simple average, primary", start = "",
  staircase  = savg(pri$staircase),
  bc_surface = savg(pri$bc_surface),
  grid_ols   = savg(pri$grid_ols),
  grid_env   = savg(pri$grid_env),
  grid_env_bc = savg(pri$grid_env_bc),
  sfa_cost   = savg(pri$sfa_cost),
  sfa_bc     = savg(pri$sfa_bc),
  sfab       = savg(pri$sfab),
  sfab_bc    = savg(pri$sfab_bc),
  ols_runs   = savg(pri$ols_runs),
  ols_runs_bc = savg(pri$ols_runs_bc),
  par_logit  = savg(pri$par_logit),
  par_lg_bc  = savg(pri$par_lg_bc),
  pl_env     = savg(pri$pl_env),
  pl_env_bc  = savg(pri$pl_env_bc))
summary_rows <- rbind(simple, pooled)

pc <- function(x) ifelse(is.na(x), "", sprintf("%.1f%%", x))

# The HTML gets a real minus sign. sprintf leaves an ASCII hyphen, which is
# narrower than the digits it sits against and reads as a dash rather than a
# sign. The console keeps the hyphen, where the entity would print literally.
pc_html <- function(x) sub("^-", "&minus;", pc(x))

# Two blocks, mirroring the two saved tables: the frontier/all-data cost
# models, then the SFA pair and the accuracy models.
console_block <- function(hdr, keys) {
  f <- paste0("%-6s %8s ", paste(rep("%10s", length(keys)), collapse = " "),
              "\n")
  cat(do.call(sprintf, c(list(f, "bench", "from"), as.list(hdr))))
  for (i in seq_len(nrow(rate_rows) + nrow(summary_rows))) {
    r <- if (i <= nrow(rate_rows)) rate_rows[i, ] else
      summary_rows[i - nrow(rate_rows), ]
    if (i == nrow(rate_rows) + 1)
      cat(strrep("-", 16 + 11 * length(keys)), "\n")
    cat(do.call(sprintf, c(list(f, r$bench, r$start),
                           lapply(keys, function(k) pc(r[[k]])))))
  }
}
cat("cost decline per quarter at fixed accuracy (positive = cheaper)\n\n")
cat("frontier and all-data cost models\n")
console_block(c("model-free", "grid OLS", "BC surf", "grid+env", "env BC",
                "OLS runs", "OLS BC"),
              c("staircase", "grid_ols", "bc_surface", "grid_env",
                "grid_env_bc", "ols_runs", "ols_runs_bc"))
cat("\nstochastic frontier and accuracy models\n")
console_block(c("model-free", "SFA", "SFA BC", "SFAb", "SFAb BC",
                "par.logit", "par.lg BC", "pl env", "plenv BC"),
              c("staircase", "sfa_cost", "sfa_bc", "sfab", "sfab_bc",
                "par_logit", "par_lg_bc", "pl_env", "pl_env_bc"))

## ---- the same comparison, saved as a table (HTML) ------------------------------

dir.create(out_path("tables"), showWarnings = FALSE, recursive = TRUE)

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
       # The two summary rows close the body: a rule and a half-row gap above
       # the FIRST of them, and italic labels on both, so they read as
       # summaries of the rows rather than as more benchmarks among them.
       '.pooled-row td{border-top:1px solid #1d1d1d;padding-top:1em}',
       '.pooled-row td:first-child,.summary2-row td:first-child{font-style:italic}',
       # subtitles for the two tables the page now holds
       'h2{font-size:1.05em;margin:18px 0 4px 0}',
       '</style></head><body>',
       '<h1>Average quarterly rate of decline in cost of given accuracy</h1>')

# One rendered table: subtitle, its own two header rows over the shared
# Benchmark / Data start / Model-free stub, the model-order row, every
# benchmark and summary row restricted to `keys`, and its own footnote.
html_tbl <- function(subtitle, group_html, order_cells, keys, notes) {
  out <- c(sprintf('<h2>%s</h2>', subtitle),
           '<table><thead><tr>',
           # Benchmark alone keeps its bare label: it is the left-aligned
           # first column, where a centred block would drift off the edge.
           '<th rowspan="2">Benchmark</th>',
           '<th rowspan="2"><div class="hd">Data start</div></th>',
           '<th rowspan="2"><div class="hd">Model-free</div></th>',
           group_html,
           '</tr></thead><tbody>',
           sprintf('<tr class="order-row">%s</tr>',
                   paste0(sprintf('<td>%s</td>',
                                  c("<i>Model order</i>", "", "",
                                    order_cells)), collapse = '')))
  for (i in seq_len(nrow(rate_rows) + nrow(summary_rows))) {
    r <- if (i <= nrow(rate_rows)) rate_rows[i, ] else
      summary_rows[i - nrow(rate_rows), ]
    cells <- c(esc(r$start), pc_html(unlist(r[, c("staircase", keys)])))
    # the summary rows are summaries, not more benchmarks: a rule above the
    # first and italic labels on both keep them apart from the data rows
    out <- c(out, sprintf('<tr%s><td>%s</td>%s</tr>',
                          if (i == nrow(rate_rows) + 1)
                            ' class="pooled-row"' else
                            if (i > nrow(rate_rows))
                              ' class="summary2-row"' else '',
                          esc(r$benchmark),
                          paste0(sprintf('<td>%s</td>', cells),
                                 collapse = '')))
  }
  c(out, sprintf('</tbody><tfoot><tr><td colspan="%d">', length(keys) + 3),
    notes, '</td></tr></tfoot></table>')
}

# Footnote text shared by both tables' openings.
NOTE_ESTIMAND <- paste(
  "All numbers are estimates of the average quarterly drop in the cost",
  "of a given level of accuracy on a given benchmark, over the years of available data.")
NOTE_MODELFREE <- paste(
  "The model-free values are averages over grids with one date node every ~13 days (the same time resolution for every benchmark, so longer",
  "histories carry proportionally more nodes) and 100 accuracy-uniform levels spanning the benchmark's achieved range, truncated at each date to the",
  "state of the art (SOTA) by then &mdash; the same defined-where-achieved rule the fitted models' grids apply. At each point, the lowest cost",
  "of performance at least as good one QUARTER later is found and divided by the initial cost. The geometric mean of the ratios is the quarterly decline rate.",
  "Over so short a horizon most levels do not move at all, and each such node enters the geometric mean as a ratio of 1, so the rate is a minority of real drops",
  "averaged in with a majority of unchanged records.")
NOTE_SUMMARIES <- paste(
  "The two closing rows summarize the five primary benchmarks &mdash; AIME, Chess Puzzles, FrontierMath tiers 1&ndash;3, GPQA Diamond and Mystery Game Puzzles &mdash;",
  "two ways. The first is the SIMPLE AVERAGE: the unweighted mean of the five annual log-cost changes, converted back to a quarterly rate &mdash; transparent, free of any",
  "external calibration, and the equal-weights limit of random-effects pooling, appropriate insofar as the benchmarks are exchangeable draws whose true rates differ.",
  "The final row instead pools the same five estimates with weights &alpha;<sub>b</sub><sup>2</sup> &times; the benchmark's observed history in years &times; the record",
  "path's average of p(1&minus;p), where &alpha;<sub>b</sub> is the benchmark's discrimination in the 2PL behind Epoch's ECI (Epoch Capabilities Index), which writes logit",
  "accuracy as &alpha;<sub>b</sub>(C &minus; D<sub>b</sub>) for a shared capability C. The three factors are the benchmark's integrated Fisher information about the shared",
  "capability path: per-response information is &alpha;<sup>2</sup>p(1&minus;p), and a rate is pinned down over its time base, so short histories and saturated benchmarks",
  "are both discounted automatically.",
  "The two summaries typically differ by a point or two; a delete-one-model jackknife puts each benchmark's own sampling error at several times that,",
  "so the gap between them is a difference of weighting philosophies, not one the data can resolve.")

t1 <- html_tbl(
  "Frontier and all-data cost models",
  c('<th colspan="6"><div class="hd">Cost models</div></th>',
    '</tr><tr>',
    '<th colspan="2"><div class="hd">Model frontier</div></th>',
    '<th colspan="2"><div class="hd">Model frontier, require envelopment</div></th>',
    '<th colspan="2"><div class="hd">Model all data</div></th>'),
  rep(c("Linear", "Box-Cox"), 3),
  c("grid_ols", "bc_surface", "grid_env", "grid_env_bc",
    "ols_runs", "ols_runs_bc"),
  paste(NOTE_ESTIMAND, NOTE_MODELFREE,
        "Results in the Linear columns are coefficients on release year as an explanator for log cost, reexpressed quarterly. The nonlinear Box-Cox variants,",
        "having no single such coefficient, are summarized by the average of the instantaneous rate of decline the fitted surface implies at each point &mdash; over the same",
        "lattice as the model-free column for the two frontier fits, and over the runs themselves for the all-data fit, whose population is the cloud rather than the record.",
        "That average differences nothing over time &mdash; the record needs a horizon because it is a step function, a fitted surface does not &mdash; and on a linear fit it reduces",
        "exactly to the coefficient in the column beside it.",
        "\"Model frontier\" means modeling the empirical frontier as realized at a grid of points with ordinary least squares (OLS).",
        "\"Model frontier, require envelopment\" means the same, but with the constraint that the fitted surface is never above any data point.",
        "\"Model all data\" means modeling all runs, not just the frontier, with OLS.",
        NOTE_SUMMARIES,
        "Every column here is pooled as the weighted (or unweighted) mean of the annual log-cost change, converted back to a quarterly rate; inverse-variance weighting is",
        "unavailable because the grid and model-free columns carry no standard errors, and for those fits the &alpha;<sup>2</sup> weights borrow an information interpretation",
        "they cannot support &mdash; there is no likelihood behind them &mdash; so their pooled entries are mechanical averages. The five primaries'",
        "shares under the weight are",
        sprintf("%s.", share_txt),
        "Several SECONDARY benchmarks carry a larger &alpha;<sup>2</sup> than any primary; under this weighting their short histories discount them automatically,",
        "so their continued exclusion is a choice about data maturity rather than a necessity.",
        "(A third summary &mdash; one pooled regression over the five primaries' stacked runs with benchmark fixed effects &mdash; appeared in earlier versions of this",
        "table; it weights benchmarks by data volume and leverage and landed within the same sampling error, and its Box-Cox cell imposed a shared curvature the",
        "benchmarks reject, so it was removed.)"))

t2 <- html_tbl(
  "Stochastic frontier and accuracy models",
  c('<th colspan="4"><div class="hd">Cost models</div></th>',
    '<th colspan="4"><div class="hd">Accuracy models</div></th>',
    '</tr><tr>',
    '<th colspan="2"><div class="hd">Stochastic frontier analysis</div></th>',
    '<th colspan="2"><div class="hd">Stochastic frontier analysis, time-varying inefficiency spread</div></th>',
    '<th colspan="2"><div class="hd">Model frontier</div></th>',
    '<th colspan="2"><div class="hd">Model frontier, require envelopment</div></th>'),
  rep(c("Linear", "Box-Cox"), 4),
  c("sfa_cost", "sfa_bc", "sfab", "sfab_bc",
    "par_logit", "par_lg_bc", "pl_env", "pl_env_bc"),
  paste(NOTE_ESTIMAND,
        "The Data start and Model-free columns repeat the first table's; see its notes for their construction.",
        "Stochastic frontier analysis models ALL runs: log cost is a frontier in accuracy and date plus a half-normal inefficiency term (one draw per model &times;",
        "reasoning-effort group) plus noise; the second variant lets the log inefficiency spread vary linearly with date. Read these columns with care: they track",
        "the dense cheap edge of the run cloud, which grows dearer as expensive reasoning configurations arrive at every accuracy level, so a negative decline here",
        "is the cloud drifting up while the record falls. Linear SFA rates are the time coefficient reexpressed quarterly; the Box-Cox variants average the fitted",
        "frontier's instantaneous rate over the runs, as the first table's all-data fit does. A blank Box-Cox cell means that average is not finite: an extreme",
        "profiled lambda_time makes the surface's rate overflow at the earliest runs.",
        "In the accuracy model columns accuracy is the dependent variable and log cost an explanatory variable, fit at a grid of points with a logit link",
        "(\"Model frontier\"), optionally constrained so the fitted surface is never below any data point (\"require envelopment\"). The linear rates are extracted",
        "as -b_t/b_x. The Box-Cox (Box-Tidwell) columns have no single such ratio, so they take the grid analogue: with z the fitted logit index, the ratio of the",
        "node-averaged slopes, &minus;mean(&part;z/&part;t) / mean(&part;z/&part;ln c) over the accuracy grid's defined nodes where the surface rises in cost &mdash;",
        "equivalently a &part;z/&part;ln c-weighted mean of each node's implied d ln cost/dt, so near-flat nodes whose pointwise ratio diverges carry no weight &mdash;",
        "which reduces exactly to -b_t/b_x on a linear fit.",
        NOTE_SUMMARIES,
        "The two LINEAR accuracy-model rates are ratios &minus;b_t/b_x in which &alpha;<sub>b</sub> cancels within a benchmark but not across them, so the weighted row",
        "instead takes the ratio of the ECI-pooled slopes, &minus;&Sigma;&alpha;Tb_t / &Sigma;&alpha;Tb_x, matching the pooled decline in the regression tables, while",
        "the simple average averages the five ratios directly. Every other column in this table is pooled as the weighted (or unweighted) mean of the annual log-cost",
        "change, converted back to a quarterly rate."))

o <- c(o, t1, t2, '</body></html>')
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