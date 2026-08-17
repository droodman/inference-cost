# Empirical Pareto frontier, the nonparametric counterpart to plot_frontiers.R.
#
#   f_t(c) = max { a_i : t_i <= t, c_i <= c }
#
# A running maximum over everything released by date t and costing no more than
# c. Non-decreasing in BOTH arguments by construction, so its curves can never
# cross and it can never run backwards -- the failure modes the parametric
# frontiers had. What it cannot do is separate signal from luck: a single
# fortunate run sets the frontier permanently, and with 20-46% of runs scoring
# exactly 0 and n_samples as low as 30, one lucky draw is a real possibility.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("frontier_viz.R")

## ---- data ---------------------------------------------------------------------------

# gpt-4o is retained here (unlike plot_frontiers.R). The frontier is a maximum,
# so dropping a model can only lower it, and there is no fit for an outlier to
# distort -- flip this to TRUE to match the parametric figures exactly.
DROP_GPT4O <- FALSE

d <- load_runs(drop_gpt4o_chess = DROP_GPT4O)
dates <- bench_dates(d)

## ---- the frontier ---------------------------------------------------------------------

# Jump points of the running maximum, extended flat to the benchmark's largest
# cost -- correct, not cosmetic: for c above every observed cost the max is
# unchanged, so the frontier really is flat out there.
pareto_steps <- function(sub, q) {
  s <- sub[sub$releasedate <= q, ]
  if (!nrow(s)) return(NULL)
  s <- s[order(s$cost), ]
  m <- cummax(s$acc)
  keep <- c(TRUE, diff(m) > 0)
  data.frame(cost = c(s$cost[keep], max(sub$cost)),
             value = c(m[keep], m[length(m)]))
}

pareto_curves <- function(data, dates_by_bench) {
  do.call(rbind, lapply(names(dates_by_bench), function(b) {
    sub <- data[data$benchmark == b, ]
    do.call(rbind, lapply(dates_by_bench[[b]], function(q) {
      st <- pareto_steps(sub, q)
      if (is.null(st)) return(NULL)
      st$qdate <- q
      st$year <- 2023 + as_t(q)
      st$benchmark <- b
      st
    }))
  }))
}

curves <- pareto_curves(d, dates)

## ---- figure ----------------------------------------------------------------------------

p <- frontier_plot(
  curves, d, step = TRUE,
  title = "Empirical Pareto frontier of accuracy by cost per task",
  subtitle = "f_t(c) = best accuracy among all runs released by t costing no more than c",
  ylab = "Best accuracy achieved",
  notes = c(
    paste("Non-decreasing in cost and in time by construction, so curves cannot",
          "cross and a step persists once set -- a released model stays available."),
    paste("A running maximum, so one lucky run fixes the frontier permanently;",
          "read closely-spaced late steps with that in mind."),
    if (DROP_GPT4O) "Chess excludes model gpt-4o; gpt-4o-mini is retained."
    else "All models retained, including gpt-4o on chess."))

ggsave(out_path("pareto_frontier.png"), p, width = 10, height = 7.5, dpi = 200,
       device = ragg::agg_png)
cat("wrote pareto_frontier.png\n")

## ---- how much of the frontier rests on a single run? -------------------------------------

cat("\nsteps in the final Pareto frontier, and the best accuracy reached\n")
for (b in sort(unique(d$benchmark))) {
  sub <- d[d$benchmark == b, ]
  st <- pareto_steps(sub, max(sub$releasedate))
  cat(sprintf("%-6s %3d steps over %5d runs   best acc %.3f at $%.4f\n",
              b, nrow(st) - 1, nrow(sub), max(st$value),
              st$cost[which.max(st$value)]))
}
