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
#   staircase  pareto_decline_qtr, the model-free average (blank where the
#              benchmark spans less than the one-year horizon)
#   grid OLS   fit_lncost_grid: ln C_a(t) sampled on the (logit acc, date)
#              grid, OLS through the samples; rate = its tc coefficient
#   BC surf    the same statistic as the staircase column -- the one-year
#              decline averaged over the check's nodes -- but read off the
#              Box-Cox grid OLS surface instead of the records themselves,
#              so it asks whether smoothing changes the aggregate
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
    staircase  = if (is.null(np)) NA_real_ else np$pct_qtr,
    bc_surface = if (is.null(ms)) NA_real_ else ms$pct_qtr,
    grid_ols   = cost_decline_qtr(store_cost("costgridols")$lin[[b]]),
    grid_env   = cost_decline_qtr(store_cost("costgridolsenv")$lin[[b]]),
    sfa_cost   = cost_decline_qtr(store_cost("costsfa")$lin[[b]]),
    ols_runs   = cost_decline_qtr(store_cost("costols")$lin[[b]]),
    par_logit  = rate_acc(store_grid("paretologit")$lin[[b]]),
    pl_env     = rate_acc(store_grid("paretologitenv")$lin[[b]]))
}))

pc <- function(x) ifelse(is.na(x), "", sprintf("%.1f%%", x))

# The HTML gets a real minus sign. sprintf leaves an ASCII hyphen, which is
# narrower than the digits it sits against and reads as a dash rather than a
# sign. The CSV keeps the hyphen -- it has to parse back as a number -- and
# so does the console, where the entity would print literally.
pc_html <- function(x) sub("^-", "&minus;", pc(x))

cat("cost decline per quarter at fixed accuracy\n")
cat("(positive = cheaper; frontier columns first, then the run-cloud pair)\n")
cat(sprintf("%-6s %8s %10s %9s %9s %9s | %9s %9s | %10s %8s\n",
            "bench", "from", "staircase", "grid OLS", "BC surf", "grid+env",
            "OLS runs", "SFA cost", "par.logit", "pl env"))
for (i in seq_len(nrow(rate_rows))) {
  r <- rate_rows[i, ]
  cat(sprintf(
    "%-6s %8s %10s %9s %9s %9s | %9s %9s | %10s %8s\n",
    r$bench, r$start, pc(r$staircase), pc(r$grid_ols), pc(r$bc_surface),
    pc(r$grid_env), pc(r$ols_runs), pc(r$sfa_cost), pc(r$par_logit),
    pc(r$pl_env)))
}

## ---- the same comparison, saved as a table (HTML + CSV) ------------------------------

dir.create(out_path("tables"), showWarnings = FALSE, recursive = TRUE)

csv <- rate_rows[, c("benchmark", "start", "staircase", "grid_ols",
                     "bc_surface", "grid_env", "ols_runs", "sfa_cost",
                     "par_logit", "pl_env")]
names(csv) <- c("benchmark", "data_start", "staircase_pct_qtr",
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
       '</style></head><body>',
       '<h1>Average quarterly rate of decline in cost of given accuracy</h1>',
       '<table><thead><tr>',
       # Benchmark alone keeps its bare label: it is the left-aligned first
       # column, where a centred block would drift off the column's edge.
       '<th rowspan="2">Benchmark</th>',
       '<th rowspan="2"><div class="hd">Data start</div></th>',
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
order_cells <- c("<i>Model order</i>", "", "", "Linear", "Box-Cox", "Linear",
                 "Linear", "Linear", "Linear", "Linear")
o <- c(o, sprintf('<tr class="order-row">%s</tr>',
                  paste0(sprintf('<td>%s</td>', order_cells), collapse = '')))
for (i in seq_len(nrow(rate_rows))) {
  r <- rate_rows[i, ]
  cells <- c(esc(r$start),
             pc_html(unlist(r[, c("staircase", "grid_ols", "bc_surface",
                                  "grid_env", "ols_runs", "sfa_cost",
                                  "par_logit", "pl_env")])))
  o <- c(o, sprintf('<tr><td>%s</td>%s</tr>', esc(r$benchmark),
                    paste0(sprintf('<td>%s</td>', cells), collapse = '')))
}
o <- c(o, '</tbody><tfoot><tr><td colspan="10">',
       paste("All numbers are estimates of the average quarterly drop in the cost",
             "of a given level of accuracy on a given benchmark, over the years of available data.",
             "The \"non-parametric\" values are averages over even 100&times;100 grids that span the benchmark's release date and accuracy ranges,",
             "with the latter scaled within the state of the art (SOTA) for that benchmark at each given time. At each point, the lowest cost",
             "of performance at least as good one year later is found and divided by the initial cost. The geometric mean of the ratios is reexpressed as a quarterly decline rate.",          
             "Results in all cost model columns are coefficients on release year as an explanator for log cost, again reexpressed quarterly; except that in the nonlinear Box-Cox variant,",
             "the statistic is computed as in the non-parametric column, using the model's best-fit surface. In the accuracy model columns accuracy",
             "is the dependent variable and log cost as an explanatory variable, so the rates are extracted as -b_t/b_x.",
             "\"Model frontier\" means modeling the empirical frontier as realized at a grid of points, either with ordinary least squares (OLS; for log cost) or with a logit link (for accuracy).",
             "\"Model frontier, require envelopment\" means the same, but with the constraint that the fitted surface is never above any data point (for cost) or below (for accuracy).",
             "\"Model all data\" means modeling all runs, not just the frontier, with OLS. \"Stochastic frontier analysis\" models the frontier and the distribution of runs around it", 
             "with a half-normal distribution of inefficiency."),
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