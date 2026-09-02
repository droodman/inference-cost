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
src_source("envelope_frontier.R")   # pareto_decline_qtr(), for the rate check below

## ---- data ---------------------------------------------------------------------------

d <- load_runs()
dates <- bench_dates(d)
benches <- bench_levels(d$benchmark)

## ---- the frontier ---------------------------------------------------------------------

# pareto_steps() and pareto_curves() are in frontier_viz.R: this figure draws the
# staircase on its own, while plot_paretologit.R and plot_paretologitenv.R draw the
# same object underneath their fitted curves, and all three must agree exactly.
curves <- pareto_curves(d, dates)

## ---- figure ----------------------------------------------------------------------------

p <- frontier_plot(
  curves, d, step = TRUE,
  ylab = "Best accuracy achieved",
  notes = c(
    paste("Non-decreasing in cost and in time by construction, so curves cannot",
          "cross and a step persists once set -- a released model stays available."),
    paste("A running maximum, so one lucky run fixes the frontier permanently;",
          "read closely-spaced late steps with that in mind.")))

ggsave(out_path("pareto_frontier.png"), p, width = 10, height = fig_height(length(benches)), dpi = 200,
       device = ragg::agg_png)
cat("wrote pareto_frontier.png\n")

## ---- how much of the frontier rests on a single run? -------------------------------------

cat("\nsteps in the final Pareto frontier, and the best accuracy reached\n")
for (b in bench_levels(d$benchmark)) {
  sub <- d[d$benchmark == b, ]
  st <- pareto_steps(sub, max(sub$releasedate))
  cat(sprintf("%-6s %3d steps over %5d runs   best acc %.3f at $%.4f\n",
              b, nrow(st) - 1, nrow(sub), max(st$value),
              st$cost[which.max(st$value)]))
}

## ---- each all-time record's own cost history ----------------------------------------------

# The iso-accuracy view with the levels chosen by the DATA: a contour is born
# the moment a model's best score beats every earlier model's, at that
# model's (release date, record-run cost), and traces C_a(t) -- the lowest
# cost at which that accuracy had since been matched -- from then on.
# iso_pareto_steps() computes exactly that curve; the succession of ceilings
# replaces the fixed ladder of round levels.
#
# The point: records are born DEAR -- set by a frontier model near the top of
# its own budget curve, where accuracy is saturating in compute -- and
# cheapen when later models pass the level. Each contour is one traversal of
# the cliff-then-shelf profile that defeats the additive fitted cost
# surfaces: steep early drops as the frontier moves past the level, then
# commodity pricing. Levels are nested by construction (each record beats the
# last), so the staircases cannot cross.

# One row per record event: the model, its release date, the record level,
# and the cost of the record run (= C_a at birth, since no earlier run
# reaches the level).
ceiling_events <- function(sub) {
  mb <- aggregate(acc ~ model, data = sub, FUN = max)
  mb$releasedate <- sub$releasedate[match(mb$model, sub$model)]
  mb <- mb[order(mb$releasedate, mb$acc), ]
  cm <- cummax(mb$acc)
  keep <- c(TRUE, mb$acc[-1] > cm[-nrow(mb)]) & mb$acc > 0
  mb[keep, , drop = FALSE]
}

rec_curves <- list(); rec_births <- list()
for (b in bench_levels(d$benchmark)) {
  sub <- d[d$benchmark == b, ]
  ev <- ceiling_events(sub)
  for (i in seq_len(nrow(ev))) {
    st <- iso_pareto_steps(sub, ev$acc[i])
    st$acc <- ev$acc[i]
    st$benchmark <- b
    rec_curves[[length(rec_curves) + 1]] <- st
    rec_births[[length(rec_births) + 1]] <- data.frame(
      benchmark = b, model = ev$model[i], date = st$date[1], acc = ev$acc[i],
      cost = st$cost[1], last = st$cost[nrow(st)])
  }
}
rec_curves <- do.call(rbind, rec_curves)
rec_births <- do.call(rbind, rec_births)

iso_ranges <- do.call(rbind, lapply(bench_levels(d$benchmark), function(b) {
  s <- d[d$benchmark == b, ]
  data.frame(benchmark = b, date = range(s$releasedate), cost = range(s$cost))
}))

# iso_acc_plot supplies the chrome (same panels, scales and cost floor as
# every other iso figure); the curves themselves go in as step layers, since
# a running minimum is a staircase, not a path through samples.
rc <- rec_curves; rc$benchmark <- factor(LABELS[rc$benchmark], levels = LABELS)
rb <- rec_births; rb$benchmark <- factor(LABELS[rb$benchmark], levels = LABELS)
empty <- data.frame(date = as.Date(character()), cost = numeric(),
                    acc = numeric(), seg = integer(), benchmark = character())
p_rec <- iso_acc_plot(
  empty, d, ranges = iso_ranges,
  notes = c(
    paste("Solid: the cost history of each all-time accuracy record. A",
          "contour is born when a model first beats every earlier score",
          "(dot) and traces the lowest cost at which that accuracy had since",
          "been matched."),
    paste("Colour is the record's accuracy, on the same scale as the",
          "observed runs behind it."),
    paste("Records are born dear -- set by a frontier model near the top of",
          "its own budget curve -- and cheapen in a few large drops as later",
          "models pass the level, rather than at a steady rate."),
    "Nonparametric: nothing on this figure is fitted.")) +
  geom_step(data = rc, aes(date, cost, group = interaction(benchmark, acc),
                           colour = acc),
            direction = "hv", linewidth = 0.6, inherit.aes = FALSE) +
  geom_point(data = rb, aes(date, cost, colour = acc), size = 1.6,
             inherit.aes = FALSE)

ggsave(out_path("isoaccuracy_records.png"), p_rec, width = 10, height = fig_height(length(benches)),
       dpi = 200, device = ragg::agg_png)
cat("wrote isoaccuracy_records.png\n")

cat("\nall-time accuracy records: born dear, then commoditized\n")
cat(sprintf("%-6s %-32s %-11s %6s %10s %10s %7s\n", "bench", "model",
            "date", "acc", "born at", "now", "drop"))
for (i in seq_len(nrow(rec_births))) {
  e <- rec_births[i, ]
  cat(sprintf("%-6s %-32s %-11s %5.1f%% %10.4f %10.4f %6.0fx\n",
              e$benchmark, substr(e$model, 1, 32),
              format(e$date, "%Y-%m-%d"), 100 * e$acc, e$cost, e$last,
              e$cost / e$last))
}

## ---- cost decline at fixed frontier performance, read off the staircase alone ------------

# The tables' "cost drop, %/qtr" (regression_tables.R) with no model at all:
# at every (accuracy, date) node -- accuracy levels uniform WITHIN each
# date's state of the art, dates uniform (decline_nodes,
# envelope_frontier.R) -- the record cost of the node's level one QUARTER
# later against the record now, averaged geometrically. A benchmark observed
# for less than a quarter has no measurable horizon and drops out of the
# table. `moved` is the share of nodes whose record changed at all over the
# quarter; over so short a horizon most levels do not move, and an unmoved
# record enters the geometric mean as a ratio of 1, so the rate is that
# minority's real drops averaged in with a majority of unchanged records,
# and the two columns have to be read together.
cat("\ncost decline at fixed frontier performance, from the staircase alone",
    "\n(measured over a quarter, expressed as a quarterly rate)\n")
cat(sprintf("%-6s %8s %8s %7s\n", "bench", "%/qtr", "moved", "nodes"))
check <- list()
for (b in bench_levels(d$benchmark)) {
  r <- pareto_decline_qtr(d[d$benchmark == b, ])
  if (is.null(r)) next
  cat(sprintf("%-6s %7.1f%% %7.1f%% %7d\n",
              b, r$pct_qtr, 100 * r$share_moved, r$n_nodes))
  check[[b]] <- r
}

## ---- the same check, saved as a table (HTML + CSV) -----------------------------------

dir.create(out_path("tables"), showWarnings = FALSE, recursive = TRUE)

write.csv(data.frame(benchmark = unname(LABELS[names(check)]),
                     decline_pct_qtr = vapply(check, `[[`, 0, "pct_qtr"),
                     share_records_moved = vapply(check, `[[`, 0,
                                                  "share_moved"),
                     n_nodes = vapply(check, `[[`, 0L, "n_nodes"),
                     row.names = NULL),
          out_path("tables", "staircase_check.csv"), row.names = FALSE)
cat("\nwrote staircase_check.csv\n")

esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  gsub("<", "&lt;", x, fixed = TRUE)
}
o <- c('<!DOCTYPE html>', '<html lang="en"><head><meta charset="UTF-8" />',
       '<title>Staircase check</title>',
       '<style>',
       # the regression tables' palette, so this page sits beside them
       'body{font-family:"Segoe UI",Arial,sans-serif;margin:12px;color:#1d1d1d;background:#fcfcfb}',
       'h1{font-size:1.25em;margin:0 0 4px 0}',
       'p.sub{color:#5e5e5e;margin:0 0 16px 0;font-size:.9em}',
       'table{border-collapse:collapse;background:#fcfcfb;font-size:.9em}',
       'th,td{padding:3px 14px;text-align:right;white-space:nowrap}',
       'th:first-child,td:first-child{text-align:left}',
       'thead th{border-bottom:1px solid #1d1d1d;font-weight:600}',
       'tfoot td{border-top:1px solid #1d1d1d;font-size:.92em;color:#5e5e5e;',
       '  text-align:left;white-space:normal;padding-top:8px;max-width:760px}',
       '</style></head><body>',
       '<h1>Cost decline at fixed frontier performance, staircase alone</h1>',
       '<p class="sub">The nonparametric sense check: no model anywhere.</p>',
       '<table><thead><tr><th>Benchmark</th><th>Decline, %/qtr</th>',
       '<th>Records moved</th><th>Grid nodes</th></tr></thead><tbody>')
for (b in names(check)) {
  r <- check[[b]]
  o <- c(o, sprintf('<tr><td>%s</td><td>%.1f%%</td><td>%.1f%%</td><td>%d</td></tr>',
                    esc(LABELS[[b]]), r$pct_qtr, 100 * r$share_moved,
                    r$n_nodes))
}
o <- c(o, '</tbody><tfoot><tr><td colspan="4">',
       paste("At every node of an (accuracy, date) lattice -- 100 accuracy",
             "levels placed uniformly WITHIN each date's state of the art",
             "(a midpoint lattice on (0, SOTA], shrunk per date rather than",
             "clipped, so every date contributes equally and every node is",
             "well-defined), 100 dates, only dates with a full quarter of",
             "data ahead of them -- the record cost of the node's level one",
             "QUARTER later is compared with the record now, and the",
             "geometric-mean log change is the quarterly rate the regression",
             "tables print, so the columns are directly comparable. Records",
             "moved is the share of nodes whose record changed at all over",
             "the quarter: over so short a horizon most levels do not move,",
             "and an unmoved record enters the geometric mean as a ratio of",
             "1 -- an unchanged cost, not a zero one -- so the rate is that",
             "minority's real drops averaged in with a majority of unchanged",
             "records, and a low rate on a low share is sparse jumps rather",
             "than a slow frontier. A benchmark observed for less than a",
             "quarter has no measurable horizon and no row. Because each date's lattice",
             "spans its own achieved range, a fixed node index is a rising",
             "absolute level as the SOTA climbs; each node's change is still",
             "a fixed-level quantity."),
       '</td></tr></tfoot></table></body></html>')
writeLines(o, out_path("tables", "staircase_check.html"))
cat("wrote staircase_check.html\n")
