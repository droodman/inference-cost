# Console report: quarterly cost-decline rates from the cost-direction dual
# fits (cost_frontier.R), beside the accuracy-direction fits they mirror and
# the nonparametric staircase check. All linear specifications; every column
# targets the same estimand, the fall per quarter in the cost of fixed
# frontier performance.
#
#   staircase  pareto_decline_qtr, the model-free average (blank where the
#              benchmark spans less than the one-year horizon)
#   grid OLS   fit_lncost_grid: ln C_a(t) sampled on the (logit acc, date)
#              grid, OLS through the samples; rate = its tc coefficient
#   grid+env   fit_lncost_grid_env: the grid OLS objective under the cost
#              envelope's constraints; rate = its tc coefficient
#   SFA cost   fit_cost_sfa: stochastic cost frontier, half-normal per
#              model x effort; rate = its tc coefficient
#   OLS runs   lm of ln cost on logit acc and date over all positive-accuracy
#              runs -- model S's reverse regression, the typical run's cost
#   par.logit  fit_pareto_logit (accuracy direction); rate = -b_t/b_x
#   pl env     fit_pareto_logit_env (accuracy direction); rate = -b_t/b_x
#
# Reading the spread: the frontier-per-se columns (staircase, grid OLS,
# grid+env, par.logit, pl env) all target the record's decline, and differ by
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

cat("cost decline per quarter at fixed accuracy, all estimates\n")
cat("(positive = cheaper; frontier columns first, then the run-cloud pair)\n")
cat(sprintf("%-6s %10s %9s %9s %9s | %9s %9s | %10s %8s\n",
            "bench", "staircase", "BC surf", "grid OLS", "grid+env",
            "SFA cost", "OLS runs", "par.logit", "pl env"))
for (i in seq_len(nrow(rate_rows))) {
  r <- rate_rows[i, ]
  cat(sprintf(
    "%-6s %10s %9s %9s %9s | %9s %9s | %10s %8s\n",
    r$bench, pc(r$staircase), pc(r$bc_surface), pc(r$grid_ols),
    pc(r$grid_env), pc(r$sfa_cost), pc(r$ols_runs), pc(r$par_logit),
    pc(r$pl_env)))
}

## ---- the same comparison, saved as a table (HTML + CSV) ------------------------------

dir.create(out_path("tables"), showWarnings = FALSE, recursive = TRUE)

csv <- rate_rows[, -1]
names(csv) <- c("benchmark", "staircase_pct_qtr", "bc_surface_pct_qtr",
                "grid_ols_pct_qtr", "grid_ols_env_pct_qtr",
                "sfa_cost_pct_qtr", "ols_runs_pct_qtr",
                "pareto_logit_pct_qtr", "pareto_logit_env_pct_qtr")
write.csv(csv, out_path("tables", "rate_comparison.csv"), row.names = FALSE)
cat("\nwrote rate_comparison.csv\n")

esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  gsub("<", "&lt;", x, fixed = TRUE)
}
hd <- c("Benchmark", "Staircase", "BC surface", "Grid OLS",
        "Grid OLS + envelope", "SFA cost", "OLS runs", "Pareto logit",
        "Pareto logit + envelope")
o <- c('<!DOCTYPE html>', '<html lang="en"><head><meta charset="UTF-8" />',
       '<title>Rate comparison</title>',
       '<style>',
       'body{font-family:"Segoe UI",Arial,sans-serif;margin:12px;color:#1d1d1d;background:#fcfcfb}',
       'h1{font-size:1.25em;margin:0 0 4px 0}',
       'p.sub{color:#5e5e5e;margin:0 0 16px 0;font-size:.9em}',
       'table{border-collapse:collapse;background:#fcfcfb;font-size:.9em}',
       'th,td{padding:3px 12px;text-align:right;white-space:nowrap}',
       'th:first-child,td:first-child{text-align:left}',
       'thead th{border-bottom:1px solid #1d1d1d;font-weight:600;vertical-align:bottom}',
       'thead th .hd{white-space:normal;overflow-wrap:break-word;max-width:7em}',
       'th.grp,td.grp{border-left:1px solid #e6e6e2}',
       'tfoot td{border-top:1px solid #1d1d1d;font-size:.92em;color:#5e5e5e;',
       '  text-align:left;white-space:normal;padding-top:8px;max-width:900px}',
       '</style></head><body>',
       '<h1>Cost decline per quarter at fixed accuracy, all estimates</h1>',
       '<p class="sub">Positive = cheaper. Nonparametric and smoothed-surface',
       'aggregates first, then the frontier fits, then the run-cloud pair.</p>',
       '<table><thead><tr>',
       paste0(sprintf('<th%s><div class="hd">%s</div></th>',
                      ifelse(seq_along(hd) %in% c(4, 6, 8), ' class="grp"', ''),
                      hd), collapse = ''),
       '</tr></thead><tbody>')
for (i in seq_len(nrow(rate_rows))) {
  r <- rate_rows[i, ]
  cells <- pc(unlist(r[, 3:10]))
  o <- c(o, sprintf('<tr><td>%s</td>%s</tr>', esc(r$benchmark),
                    paste0(sprintf('<td%s>%s</td>',
                                   ifelse(seq_along(cells) %in% c(2, 4, 6),
                                          ' class="grp"', ''),
                                   cells), collapse = '')))
}
o <- c(o, '</tbody><tfoot><tr><td colspan="9">',
       paste("Every column targets the fall per quarter in the cost of fixed",
             "frontier performance. Staircase is the nonparametric check",
             "(pareto_frontiers.R, saved as staircase_check.html): record",
             "costs read off the empirical staircase at every grid node, over",
             "a one-year horizon, geometric-mean change compounded to a",
             "quarterly rate. BC surface is the SAME aggregate -- same nodes,",
             "horizon and compounding -- with the empirical record replaced",
             "by the Box-Cox Pareto-grid cost fit's surface, so the gap",
             "between the two columns is purely what smoothing does (levels",
             "of exactly 1 are excluded there; the node sets differ",
             "slightly). The remaining columns are the LINEAR fits' single",
             "coefficients: the cost-direction fits' time coefficient",
             "transformed directly, the accuracy-direction pair's -b_t/b_x",
             "ratio. SFA cost and OLS runs follow the model-effort run cloud",
             "rather than the record; negative values there are the cloud",
             "drifting dearer while the record collapses."),
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
