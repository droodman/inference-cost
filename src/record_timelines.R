# Record timelines: for each PRIMARY benchmark and each new absolute
# performance record, the sequence of models that delivered at least that
# accuracy at lower cost than any predecessor -- the raw material of the
# staircases and the record-cost fits, listed model by model.
#
# A timeline opens when a run sets a new best-ever accuracy on its benchmark
# (per release date: the day's best run, when it beats the running record;
# records of exactly zero are excluded -- every run "achieves" at least
# nothing, so its cost record would track cheap-model arrival, not the
# frontier, the same reasoning that drops a = 0 from pareto_decline_qtr).
# The record-setting model is the timeline's first row; each later row is
# the model that next held the level's cost record -- accuracy at or above
# the record, cost strictly below every earlier qualifying model's. At most
# one model per release date can enter (the day's cheapest qualifier).
#
# One table, timelines stacked: benchmark and accuracy record repeat within
# a timeline, exactly the two columns that identify it. Cost per task is
# appended so the size of each cost drop is visible, and effort is folded
# into the model name where it is informative. Written as HTML and CSV to
# output/tables/, and printed to the console.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("frontier_viz.R")   # load_runs, LABELS, PRIMARY_BENCHES

d <- load_runs()

# "default" and "none" are the no-choice effort values; anything else is a
# configuration worth naming (reasoning levels, token budgets).
model_label <- function(model, effort) {
  ifelse(effort %in% c("default", "none"), model,
         paste0(model, " (", effort, ")"))
}

timelines <- do.call(rbind, lapply(
  intersect(PRIMARY_BENCHES, bench_levels(d$benchmark)), function(b) {
    s <- d[d$benchmark == b, ]
    dates <- sort(unique(s$releasedate))

    # the accuracy records: each date's best run, where it beats the running
    # record and scores above zero
    day_max <- vapply(seq_along(dates), function(i)
      max(s$acc[s$releasedate == dates[i]]), numeric(1))
    is_rec <- day_max > cummax(c(0, head(day_max, -1))) & day_max > 0

    do.call(rbind, lapply(which(is_rec), function(i0) {
      a <- day_max[i0]
      rec <- Inf
      do.call(rbind, lapply(seq(i0, length(dates)), function(i) {
        q <- s[s$releasedate == dates[i] & s$acc >= a, , drop = FALSE]
        if (!nrow(q)) return(NULL)
        q <- q[which.min(q$cost), ]
        if (q$cost >= rec) return(NULL)
        rec <<- q$cost
        data.frame(bench = b,
                   benchmark = LABELS[[b]],
                   level = a,
                   model = model_label(q$model, q$effort),
                   acc = q$acc,
                   date = q$releasedate,
                   cost = q$cost,
                   stringsAsFactors = FALSE)
      }))
    }))
  }))

fmt_cost <- function(x) paste0("$", formatC(x, format = "fg", digits = 3))
fmt_lev  <- function(a) sprintf("%.1f%%", 100 * a)

## ---- console -------------------------------------------------------------------

cat("cost records at each accuracy record, primary benchmarks\n")
prev <- ""
for (i in seq_len(nrow(timelines))) {
  r <- timelines[i, ]
  id <- paste(r$benchmark, r$level)
  if (id != prev) cat("\n")
  prev <- id
  cat(sprintf("%-24s %7s  %-42s %7s %s %12s\n", r$benchmark, fmt_lev(r$level),
              r$model, fmt_lev(r$acc), format(r$date, "%Y-%m-%d"),
              fmt_cost(r$cost)))
}

## ---- files ---------------------------------------------------------------------

dir.create(out_path("tables"), showWarnings = FALSE, recursive = TRUE)

csv <- data.frame(benchmark = timelines$benchmark,
                  accuracy_record = timelines$level,
                  cheapest_model = timelines$model,
                  accuracy = timelines$acc,
                  release_date = format(timelines$date, "%Y-%m-%d"),
                  cost_per_task_usd = timelines$cost)
write.csv(csv, out_path("tables", "record_timelines.csv"), row.names = FALSE)
cat("\nwrote record_timelines.csv\n")

esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  gsub("<", "&lt;", x, fixed = TRUE)
}
o <- c('<!DOCTYPE html>', '<html lang="en"><head><meta charset="UTF-8" />',
       '<title>Record timelines</title>',
       '<style>',
       # the regression tables' palette, so this page sits beside them
       'body{font-family:"Segoe UI",Arial,sans-serif;margin:12px;color:#1d1d1d;background:#fcfcfb}',
       'h1{font-size:1.25em;margin:0 0 4px 0}',
       'table{border-collapse:collapse;background:#fcfcfb;font-size:.9em}',
       'th,td{padding:3px 12px;text-align:left;white-space:nowrap}',
       'td.num{text-align:right}',
       'thead th{border-bottom:1px solid #1d1d1d;font-weight:600}',
       'tr.gap td{border-top:1px solid #d9d9d5}',
       'tfoot td{border-top:1px solid #1d1d1d;font-size:.92em;color:#5e5e5e;',
       '  white-space:normal;padding-top:8px;max-width:760px}',
       '</style></head><body>',
       '<h1>Record timelines</h1>',
       '<table><thead><tr><th>Benchmark</th><th>Accuracy record</th>',
       '<th>Cheapest model</th><th>Accuracy</th><th>Release date</th>',
       '<th>Cost per task</th></tr></thead><tbody>')
prev <- ""
for (i in seq_len(nrow(timelines))) {
  r <- timelines[i, ]
  id <- paste(r$benchmark, r$level)
  o <- c(o, sprintf('<tr%s><td>%s</td><td class="num">%s</td><td>%s</td><td class="num">%s</td><td>%s</td><td class="num">%s</td></tr>',
                    if (id != prev && prev != "") ' class="gap"' else '',
                    esc(r$benchmark), fmt_lev(r$level), esc(r$model),
                    fmt_lev(r$acc), format(r$date, "%Y-%m-%d"),
                    fmt_cost(r$cost)))
  prev <- id
}
o <- c(o, '</tbody><tfoot><tr><td colspan="6">',
       paste("Each timeline opens when a model sets a new best-ever accuracy on the",
             "benchmark (the day's best run, records of exactly zero excluded); its",
             "later rows are the models that next held that level's cost record --",
             "accuracy at or above the record at strictly lower cost than every",
             "earlier qualifying model, at most one model per release date. Accuracy",
             "is the qualifying run's own score, which can exceed the record it is",
             "matching. Primary benchmarks only. Costs are per task, as in the rest",
             "of the analysis."),
       '</td></tr></tfoot></table></body></html>')
writeLines(o, out_path("tables", "record_timelines.html"))
cat("wrote record_timelines.html\n")

## ---- figures: one per benchmark --------------------------------------------------
#
# The same timelines drawn on the (release date, cost) plane, one figure per
# benchmark: each accuracy record's cost trace as a connected colored line
# marching right and down. Labels name the models, placed by ggrepel; a run
# that holds the cost record for several levels at once (common at the cheap
# end) is labeled once.

for (b in unique(timelines$bench)) {
  tl <- timelines[timelines$bench == b, ]
  tl$lev <- factor(fmt_lev(tl$level), levels = fmt_lev(sort(unique(tl$level))))

  # one label per model per date: near-identical runs of the same model often
  # hold several levels' records at once, and labeling each dot double-prints
  # the name
  labs <- tl[!duplicated(tl[c("model", "date")]), ]

  n_lev <- nlevels(tl$lev)
  p <- ggplot(tl, aes(date, cost)) +
    geom_line(aes(group = lev, colour = lev), linewidth = 0.55) +
    geom_point(aes(colour = lev), size = 1.7) +
    ggrepel::geom_text_repel(
      data = labs, aes(label = model), colour = INK_SECOND, size = 2.5,
      segment.colour = INK_MUTED, segment.size = 0.25, min.segment.length = 0.3,
      box.padding = 0.3, point.padding = 0.35, max.overlaps = Inf, seed = 1) +
    scale_y_log10(breaks = 10^(-5:2), labels = dollar_log) +
    scale_x_date(expand = expansion(mult = c(0.05, 0.05))) +
    scale_colour_manual(name = "accuracy record",
                        values = colorRampPalette(PALETTE)(n_lev)) +
    guides(colour = guide_legend(nrow = 1)) +
    labs(title = paste0(LABELS[[b]], ": cost records at each accuracy record"),
         x = NULL, y = "Cost per task (log scale)",
         caption = paste(
           "Each line follows one accuracy record's cost record over time:",
           "models scoring at least the record at lower cost than every",
           "predecessor.\nA model holding several levels' records at once is",
           "labeled once.")) +
    frontier_theme()

  f <- sprintf("record_timeline_%s.png", b)
  ggsave(out_path(f), p, width = 11, height = 7, dpi = 200,
         device = ragg::agg_png)
  cat("wrote", f, "\n")
}
