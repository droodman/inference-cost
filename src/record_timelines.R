# Record timelines at fixed accuracy levels: for each PRIMARY benchmark and
# each of the levels 10%, 50% and 90%, the sequence of models that delivered
# at least that accuracy at lower cost than any predecessor -- the raw
# material of the staircases and the record-cost fits at three canonical
# levels, listed model by model.
#
# A timeline opens on the first release date any run attains its level (on
# the rescaled scale, guessing floor at 0 -- so 50% is midway between
# guessing and perfect, the level at which the 2PL puts capability exactly
# at the benchmark's ECI difficulty D_b). Its first row is that day's
# cheapest qualifying run; each later row is the model that next held the
# level's cost record -- accuracy at or above the level, cost strictly below
# every earlier qualifying model's. At most one model per release date can
# enter (the day's cheapest qualifier). A benchmark contributes only the
# timelines whose levels it has reached: Mystery Game Puzzles appears at 10%
# but not 50%.
#
# One table, timelines stacked, and ONE figure: every trace on a single
# (release date, cost) plane -- log cost is comparable across benchmarks, so
# the plate that gave each benchmark its own panel is retired -- with each
# record holder COLOURED by the capability its own accuracy implies,
# C = logit(a)/alpha_b + D_b, and each trace named -- "GPQA Diamond, 10%" --
# by a standalone label at its opening dot. Models go by the registry's
# display names (model_display, prepare_data.R), with effort folded in where
# it is informative. Written as HTML to output/tables/, and printed to the
# console.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("frontier_viz.R")   # load_runs, LABELS, ALPHA, EDI, PRIMARY_BENCHES

d <- load_runs()

LEVELS <- c(0.1, 0.5, 0.9)   # the accuracy levels the timelines track

# "default" and "none" are the no-choice effort values; anything else is a
# configuration worth naming (reasoning levels, token budgets).
model_label <- function(model, effort) {
  ifelse(effort %in% c("default", "none"), model,
         paste0(model, " (", effort, ")"))
}

timelines <- do.call(rbind, lapply(
  intersect(PRIMARY_BENCHES, bench_levels(d$benchmark)), function(b) {
    s <- d[d$benchmark == b, ]
    do.call(rbind, lapply(LEVELS, function(lev) {
      rec <- Inf
      do.call(rbind, lapply(sort(unique(s$releasedate)), function(dt) {
        q <- s[s$releasedate == dt & s$acc >= lev, , drop = FALSE]
        if (!nrow(q)) return(NULL)
        q <- q[which.min(q$cost), ]
        if (q$cost >= rec) return(NULL)
        rec <<- q$cost
        data.frame(bench = b,
                   benchmark = LABELS[[b]],
                   level = lev,
                   model = model_label(q$model_display, q$effort),
                   acc = q$acc,
                   # the run's own capability on the anchored ECI scale
                   eci = qlogis(q$acc) / ALPHA[[b]] + EDI[[b]],
                   date = q$releasedate,
                   cost = q$cost,
                   stringsAsFactors = FALSE)
      }))
    }))
  }))

fmt_cost <- function(x) paste0("$", formatC(x, format = "fg", digits = 3))
fmt_acc  <- function(a) sprintf("%.1f%%", 100 * a)
fmt_lev  <- function(a) sprintf("%.0f%%", 100 * a)
fmt_eci  <- function(e) sprintf("%.1f", e)

## ---- console -------------------------------------------------------------------

cat("cost records at 10%/50%/90% accuracy, primary benchmarks\n")
prev <- ""
for (i in seq_len(nrow(timelines))) {
  r <- timelines[i, ]
  id <- paste(r$benchmark, r$level)
  if (id != prev) cat("\n")
  prev <- id
  cat(sprintf("%-24s %4s %-42s %7s %7s %s %12s\n", r$benchmark,
              fmt_lev(r$level), r$model, fmt_acc(r$acc), fmt_eci(r$eci),
              format(r$date, "%Y-%m-%d"), fmt_cost(r$cost)))
}

## ---- files ---------------------------------------------------------------------

dir.create(out_path("tables"), showWarnings = FALSE, recursive = TRUE)

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
       '<h1>Record timelines at 10%, 50%, and 90% performance</h1>',
       '<table><thead><tr><th>Benchmark</th><th>Level</th>',
       '<th>Cheapest model</th><th>Accuracy</th><th>Equivalent ECI</th>',
       '<th>Release date</th><th>Cost per task</th></tr></thead><tbody>')
prev <- ""
for (i in seq_len(nrow(timelines))) {
  r <- timelines[i, ]
  id <- paste(r$benchmark, r$level)
  o <- c(o, sprintf('<tr%s><td>%s</td><td class="num">%s</td><td>%s</td><td class="num">%s</td><td class="num">%s</td><td>%s</td><td class="num">%s</td></tr>',
                    if (id != prev && prev != "") ' class="gap"' else '',
                    esc(r$benchmark), fmt_lev(r$level), esc(r$model),
                    fmt_acc(r$acc), fmt_eci(r$eci),
                    format(r$date, "%Y-%m-%d"), fmt_cost(r$cost)))
  prev <- id
}
o <- c(o, '</tbody><tfoot><tr><td colspan="7">',
       paste("Each timeline opens on the first release date any run attains its",
             "accuracy level -- 10%, 50%, or 90%, rescaled so the guessing floor",
             "is 0, which makes 50% the level at which the 2PL puts capability",
             "at the benchmark's ECI difficulty; its later rows are the models",
             "that next held that level's cost record -- accuracy at or above",
             "the level at strictly lower cost than every earlier qualifying",
             "model, at most one model per release date. Accuracy is the",
             "qualifying run's own score, which can exceed the level it is",
             "matching, and Equivalent ECI converts it to the anchored ECI",
             "capability scale, logit(a)/&alpha;<sub>b</sub> + D<sub>b</sub>.",
             "Primary benchmarks only; a benchmark contributes only the levels",
             "it has reached, so Mystery Game Puzzles appears at 10% but not",
             "50%. Costs are per task, as in the rest of the analysis."),
       '</td></tr></tfoot></table></body></html>')
writeLines(o, out_path("tables", "record_timelines.html"))
cat("\nwrote record_timelines.html\n")

## ---- figure: one graph, every trace on the cost scale -----------------------------
#
# The timelines drawn on the (release date, cost) plane: one connected trace
# per benchmark x level marching right and down, each cost-record holder
# coloured by the equivalent ECI its own accuracy implies. Every trace is
# named by a standalone label -- "GPQA Diamond, 10%" -- at its opening dot,
# larger and darker than the model labels so the two kinds of text read as
# different layers. The same run often holds several of its benchmark's
# levels at once, and is labeled once. Labels name the models, horizontal,
# placed by ggrepel.
#
# LIGHT MODE, whatever frontier_viz.R's DARK toggle says: this figure is
# destined for a light context, so its chrome uses the light branch's
# constants literally rather than the globals, which would repaint it if the
# toggle ever moves. The blues are frontier_viz.R's BLUE ramp -- the original
# light-theme palette, dark = high capability on a light surface.
LABEL_SIZE <- 3.4   # these names are the point of the figure
LT <- list(ink = "#0b0b0b", second = "#52514e", muted = "#898781",
           gridline = "#e1e0d9", axis = "#c3c2b7", surface = "#fcfcfb")

light_theme <- theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    plot.background  = element_rect(fill = LT$surface, colour = NA),
    panel.background = element_rect(fill = LT$surface, colour = NA),
    panel.grid.major = element_line(colour = LT$gridline, linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.line  = element_line(colour = LT$axis, linewidth = 0.3),
    axis.text  = element_text(colour = LT$ink, size = 8),
    axis.title = element_text(colour = LT$ink, size = 9),
    plot.title    = element_text(colour = LT$ink, face = "bold", size = 13),
    plot.caption  = element_text(colour = LT$muted, size = 7.5, hjust = 0),
    legend.position = "top", legend.justification = "left",
    legend.text  = element_text(colour = LT$second, size = 8),
    legend.title = element_text(colour = LT$second, size = 8),
    plot.margin = margin(12, 16, 10, 12))

TRACE_LABEL_SIZE <- 4.1   # the trace names, a notch above the model labels

tl <- timelines

# One combined label frame -- model names AND the standalone trace names --
# so a single ggrepel layer places every piece of text: two layers would
# repel within themselves but overprint each other. Per-label size and
# colour ride identity scales, which is why the DOTS take the ECI gradient
# on fill rather than colour.
#
# Model names: one per (benchmark, model, date) -- a run that takes several
# of its benchmark's levels at once would double- or triple-print.
labs <- tl[!duplicated(tl[c("bench", "model", "date")]), ]
labs$lab  <- labs$model
labs$sz   <- LABEL_SIZE
labs$col  <- LT$second
# The trace names, one at each trace's opening dot.
tr <- tl[!duplicated(tl[c("bench", "level")]), ]
tr$lab <- paste0(fmt_lev(tr$level), " on ", tr$benchmark)
tr$sz  <- TRACE_LABEL_SIZE
tr$col <- LT$ink
all_labs <- rbind(labs[c("date", "cost", "lab", "sz", "col")],
                  tr[c("date", "cost", "lab", "sz", "col")])

# One fillable marker (shapes 21-24, so the ECI gradient stays on fill) per
# benchmark -- the diamond goes to GPQA Diamond, naturally. Mystery Game
# Puzzles gets a literal question mark instead: its down-triangle beside
# FrontierMath's up-triangle read as a pair that isn't one. A character
# marker draws in the point COLOUR and ignores fill, so its dots ride the
# identity colour scale the labels already use, their hexes computed from
# the same BLUE ramp over the same domain the fill scale is pinned to.
SHAPES <- c(aime = 21, chess = 22, fm13 = 24, gpqa = 23)
names(SHAPES) <- LABELS[names(SHAPES)]
eci_rng <- range(tl$eci)
myst <- tl[tl$bench == "mystery", ]
myst$dotcol <- scales::gradient_n_pal(BLUE)(
  scales::rescale(myst$eci, from = eci_rng))

p <- ggplot(tl, aes(date, cost)) +
  # heavier than the dots' stroke, mid-grey rather than black: the traces
  # outrank the dots but must not upstage the blues or the text (a full-black
  # experiment did exactly that)
  geom_line(aes(group = interaction(bench, level)),
            colour = LT$muted, linewidth = 0.8) +
  geom_point(data = tl[tl$bench != "mystery", ],
             aes(fill = eci, shape = benchmark), size = 2.4,
             colour = LT$surface, stroke = 0.3) +
  # geom_text rather than a "?" point shape: a character pch cannot be
  # bolded, and the plain glyph sat too faint against the trace lines
  geom_text(data = myst, aes(colour = dotcol), label = "?", size = 3.4,
            fontface = "bold") +
  ggrepel::geom_text_repel(
    data = all_labs, aes(label = lab, size = sz, colour = col),
    segment.colour = LT$muted, segment.size = 0.25, min.segment.length = 0.3,
    box.padding = 0.3, point.padding = 0.35, max.overlaps = Inf, seed = 1) +
  scale_size_identity() +
  scale_colour_identity() +
  # no key: the trace names on the plot already say which benchmark is
  # which, and the shapes only need to separate traces where they cross
  scale_shape_manual(values = SHAPES, guide = "none") +
  scale_y_log10(breaks = 10^(-5:2), labels = dollar_log) +
  scale_x_date(expand = expansion(mult = c(0.05, 0.05))) +
  scale_fill_gradientn(
    # limits pinned to the FULL eci range: mystery's dots colour off-scale
    # (the ? marker is drawn in colour, not fill) and must agree with the bar
    colours = BLUE, name = "Equivalent ECI score", limits = eci_rng,
    guide = guide_colourbar(barheight = grid::unit(0.35, "cm"),
                            barwidth = grid::unit(7, "cm"),
                            direction = "horizontal",
                            ticks.colour = LT$surface)) +
  labs(title = "Cost records at 10%, 50%, and 90% performance, primary benchmarks",
       x = NULL, y = "Cost per task (log scale)",
       caption = paste(
         "Each trace follows one benchmark's cost record at a fixed accuracy",
         "level -- 10%, 50%, or 90%, guessing floor rescaled to 0 -- from the",
         "date the level was first attained:\nmodels scoring at least the",
         "level at lower cost than every predecessor, the trace named at its",
         "opening dot. Each dot is coloured by the capability its run's own",
         "accuracy implies\non the ECI scale, logit(a)/alpha + D. A run",
         "holding several of its benchmark's levels at once is labeled once.",
         "A benchmark contributes only the levels it has reached:",
         "Mystery Game Puzzles appears at 10% only.")) +
  light_theme

f <- "record_timelines.png"
ggsave(out_path(f), p, width = 12, height = 8, dpi = 200,
       device = ragg::agg_png)
cat(sprintf("wrote %s (%d traces, one graph)\n", f,
            nrow(unique(tl[c("bench", "level")]))))
