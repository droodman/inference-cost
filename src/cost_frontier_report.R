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
#   cost env   fit_cost_envelope: highest plane under every run's log cost;
#              rate = its tc coefficient
#   SFA cost   fit_cost_sfa: stochastic cost frontier, half-normal per
#              model x effort; rate = its tc coefficient
#   OLS runs   lm of ln cost on logit acc and date over all positive-accuracy
#              runs -- model S's reverse regression, the typical run's cost
#   par.logit  fit_pareto_logit (accuracy direction); rate = -b_t/b_x
#   acc env    fit_envelope (accuracy direction); rate = -b_t/b_x
#
# Reading the spread: the frontier-per-se columns (staircase, grid OLS, cost
# env, par.logit, acc env) all target the record's decline, and differ by
# loss direction and weighting; if the surface were truly logit-linear they
# would agree. SFA cost and OLS runs instead follow the model-effort CLOUD at
# fixed accuracy, whose dense cheap edge grows dearer as expensive reasoning
# configurations arrive -- so negative "declines" there are the cloud
# drifting up while the record collapses, the sharpest statement in this file
# of how much of the story is frontier movement rather than typical-run
# movement. See cost_frontier.R's fit_cost_sfa block.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("cost_frontier.R")

d <- load_runs()

# -b_t/b_x of an accuracy-direction fit, as a quarterly percentage
rate_acc <- function(fit) {
  cf <- coef(fit)
  100 * (1 - exp(-cf[["tc"]] / cf[["lncost"]] / 4))
}

cat("cost decline per quarter at fixed accuracy, all estimates\n")
cat("(positive = cheaper; frontier columns first, then the run-cloud pair)\n")
cat(sprintf("%-6s %10s %9s %9s | %9s %9s | %10s %8s\n",
            "bench", "staircase", "grid OLS", "cost env",
            "SFA cost", "OLS runs", "par.logit", "acc env"))
for (b in sort(unique(d$benchmark))) {
  s <- d[d$benchmark == b, ]
  np <- pareto_decline_qtr(s)
  cat(sprintf("%-6s %10s %8.1f%% %8.1f%% | %8.1f%% %8.1f%% | %9.1f%% %7.1f%%\n",
              b,
              if (is.null(np)) "" else sprintf("%.1f%%", np$pct_qtr),
              cost_decline_qtr(fit_lncost_grid(s)),
              cost_decline_qtr(fit_cost_envelope(s)),
              cost_decline_qtr(fit_cost_sfa(s)),
              cost_decline_qtr(lm(lncost ~ la + tc, data = iso_runs(s))),
              rate_acc(fit_pareto_logit(s)),
              rate_acc(fit_envelope(s))))
}

cat("\ncost-direction fit details\n")
cat(sprintf("%-6s %28s %28s\n", "", "grid OLS", "cost envelope"))
cat(sprintf("%-6s %9s %8s %9s %9s %8s %9s\n", "bench",
            "$/logit", "tc", "nodes", "$/logit", "tc", "touch"))
for (b in sort(unique(d$benchmark))) {
  s <- d[d$benchmark == b, ]
  fo <- fit_lncost_grid(s)
  fe <- fit_cost_envelope(s)
  cat(sprintf("%-6s %9.3f %8.3f %9d %9.3f %8.3f %9.4f\n", b,
              coef(fo)[["la"]], coef(fo)[["tc"]], attr(fo, "n_grid"),
              coef(fe)[["la"]], coef(fe)[["tc"]], fe$slack_envelope))
}
cat("\n$/logit is the la coefficient: log dollars per logit of accuracy,\n")
cat("the surface's steepness in the cost direction. touch is the cost\n")
cat("envelope's minimum run slack in log dollars; 0 means it touches a run.\n")
