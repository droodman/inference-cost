# Record timelines at fixed accuracy levels: for each PRIMARY benchmark and
# each of the levels 25% and 75%, the sequence of models that delivered at
# least that accuracy at lower cost than any predecessor -- the raw material
# of the staircases and the record-cost fits at two canonical levels, listed
# model by model.
#
# A timeline opens on the first release date any run attains its level (on
# the rescaled scale, guessing floor at 0 -- so 50% would be midway between
# guessing and perfect, the level at which the 2PL puts capability exactly
# at the benchmark's ECI difficulty D_b; 25% and 75% straddle it). Its first
# row is that day's cheapest qualifying run; each later row is the model
# that next held the level's cost record -- accuracy at or above the level,
# cost strictly below every earlier qualifying model's. At most one model
# per release date can enter (the day's cheapest qualifier). A benchmark
# contributes only the timelines whose levels it has reached: Mystery Game
# Puzzles appears at 25% but not 75%.
#
# One table, timelines stacked, and ONE figure: every trace on a single
# (release date, cost) plane -- log cost is comparable across benchmarks, so
# the plate that gave each benchmark its own panel is retired -- with each
# trace named -- "25% on GPQA Diamond"
# -- by a standalone label at its opening dot. Models go by the registry's
# display names (model_display, prepare_data.R), with effort folded in where
# it is informative. Written as HTML to output/tables/, and printed to the
# console.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("frontier_viz.R")   # load_runs, LABELS, ALPHA, EDI, PRIMARY_BENCHES

d <- load_runs()

# The accuracy levels the timelines track.
LEVELS <- c(0.25, 0.75)

# Hand placement for the TRACE NAMES only: the six standalone labels that name
# each staircase at its opening dot. One nudge() call per label, in panel units
# (days, log10 dollars) from the label's own dot, landed exactly (force_pull =
# 0 below).
#
# The MODEL labels are left entirely to ggrepel. They carry no leader line, so
# a hand placement can only help if it lands the label against its own dot --
# which is what repel already does, and it re-solves the whole layout when the
# data move, where a fixed nudge goes stale. The trace names are the exception
# because they are thrown deliberately far, with a leader to follow.
#
# `lab` must match the rendered text exactly: a trace label is
# "<level>% on <benchmark>". `bench` disambiguates a label that appears on two
# traces; trace names are unique, so nothing needs it today -- it is kept
# because a model label, should one ever be pinned again, would.
#
# Scale, at the current 12x8in figure: the panel runs about 103 days to the
# inch across and 0.88 log decades down it, so one trace-label line height is
# y = 0.17. These rows do NOT get the main layer's baseline lift (nudge_y
# below), so y = 0 sits ON the dot, not above it. The placements were tuned by
# eye and will need retuning if the data, the figure size, or the label sizes
# change.
nudge <- function(lab, x, y, bench = NA_character_)
  data.frame(lab = lab, bench = bench, x = x, y = y, stringsAsFactors = FALSE)

TRACE_NUDGES <- rbind(
  # the column beside the o3 dot cluster, where the openers of GPQA-75,
  # Chess-25 and AIME-75 all sit at 2025-04: GPQA-75 at its diamond's height
  # (horizontal leader, so it crosses no trace), then Chess and AIME-75
  nudge("75% on GPQA Diamond",            -100, 0),
  nudge("25% on Chess Puzzles",           -130, 0.1),
  nudge("25% on AIME (OTIS Mock)",         -90, -.01),
  nudge("75% on AIME (OTIS Mock)",        -145, 0.0),
  # GPQA-25 opens at the panel's left edge (GPT-4 Turbo): well above its
  # dot, in the empty band around $0.15 over early 2024, with a leader down
  # to it -- the rows nearer the dot are taken by AIME-75 and Chess-25.
  # Hand-placed since the house face: Inter runs wider than the Cambria this
  # plate used to be set in, and repel -- which cannot see the hand-placed
  # labels -- ran the name into AIME-75.
  nudge("25% on GPQA Diamond",              90, 0.85),
  # both FrontierMath names straight above their triangles, in the open
  # space. The en dash matches the display name (LABELS,
  # frontier_viz.R); the match is on the rendered text, so it must agree.
  # 0.32, not 0.25: Inter's taller line put FM-25 onto the "o3 (medium)"
  # label under it.
  nudge("25% on FrontierMath, tiers 1–3",    0, 0.32),
  nudge("75% on FrontierMath, tiers 1–3",    0, 0.25))

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

# the level set as prose, for titles and notes: "10%, 50%, and 90%"
lev_txt <- {
  v <- fmt_lev(LEVELS)
  if (length(v) == 1) v else if (length(v) == 2) paste(v[1], "and", v[2]) else
    paste0(paste(head(v, -1), collapse = ", "), ", and ", tail(v, 1))
}

## ---- console -------------------------------------------------------------------

cat(sprintf("cost records at %s accuracy, primary benchmarks\n",
            paste(fmt_lev(LEVELS), collapse = "/")))
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
       sprintf('<h1>Record timelines at %s performance</h1>', lev_txt),
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
       paste(sprintf("Each timeline opens on the first release date any run attains its accuracy level -- %s, rescaled so the guessing floor",
                     lev_txt),
             "is 0, on which scale 50% is the level at which the 2PL puts",
             "capability at the benchmark's ECI difficulty; its later rows are",
             "the models",
             "that next held that level's cost record -- accuracy at or above",
             "the level at strictly lower cost than every earlier qualifying",
             "model, at most one model per release date. Accuracy is the",
             "qualifying run's own score, which can exceed the level it is",
             "matching, and Equivalent ECI converts it to the anchored ECI",
             "capability scale, logit(a)/&alpha;<sub>b</sub> + D<sub>b</sub>.",
             "Primary benchmarks only; a benchmark contributes only the levels",
             "it has reached. Costs are per task, as in the rest of the",
             "analysis."),
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
# The house style, from frontier_viz.R: this figure used to carry a private
# light theme against the dark experiment, and now shares the one theme with
# every other plate. Chrome at 0.65 of the guide's sizes puts the ticks at 13
# pt, matching the trace names (4.7 mm = 13.4 pt), so no text on the plate
# outsizes the names -- which are the point of the figure.
LABEL_SIZE <- 3.9   # these names are the point of the figure
TRACE_LABEL_SIZE <- 4.7   # the trace names, a notch above the model labels
TIMELINE_TYPE_SCALE <- 0.65
# The text geoms do NOT inherit the theme's base_family -- geom_text's
# default family is the device's, not the theme's -- so every text layer
# below names FONT (frontier_viz.R) explicitly.
# ONE dot colour. The dots used to be shaded by the capability each run's own
# accuracy implies, which asked the reader to track a third variable on a
# plate already carrying dates, costs, two levels and four benchmarks; the
# ECI column survives in the HTML table for anyone who wants it. The first
# house categorical, teal: one series, one colour.
DOT <- CAT[["teal"]]

# The FIGURE alone drops Mystery Game Puzzles; the table and the console
# listing above still carry it. Everything downstream reads tl, so the
# colour bar's range, the shape key and the label frames all follow.
FIG_OMIT <- "mystery"
tl <- timelines[!timelines$bench %in% FIG_OMIT, ]

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
labs$col  <- INK_SECOND
# NA segment colour draws no leader at all: a model label must sit close
# enough to its own dot to be read off it, never tethered by a line.
labs$seg  <- NA_character_
# The trace names, one at each trace's opening dot.
tr <- tl[!duplicated(tl[c("bench", "level")]), ]
tr$lab <- paste0(fmt_lev(tr$level), " on ", tr$benchmark)
tr$sz  <- TRACE_LABEL_SIZE
tr$col <- INK_PRIMARY
# the trace names DO get leaders -- they are placed well off their dots
tr$seg <- INK_MUTED
all_labs <- rbind(labs[c("bench", "date", "cost", "lab", "sz", "col", "seg")],
                  tr[c("bench", "date", "cost", "lab", "sz", "col", "seg")])
# hand-placed rows leave the main layer BY ROW, not by label: the same
# model name can label two traces' dots, and only the first (the one
# match() finds) is hand-placed -- unless TRACE_NUDGES carries a bench
# column naming which benchmark's dot is meant (NA: first match). Only the
# trace names are hand-placed today, and those are unique.
nudge_idx <- match(TRACE_NUDGES$lab, all_labs$lab)
if ("bench" %in% names(TRACE_NUDGES)) {
  spec <- !is.na(TRACE_NUDGES$bench)
  nudge_idx[spec] <- match(paste(TRACE_NUDGES$bench[spec], TRACE_NUDGES$lab[spec]),
                           paste(all_labs$bench, all_labs$lab))
}
# A nudge whose label matches nothing is a silent no-op -- the label just
# falls back to repel -- so say so: the match is on the rendered text, and a
# display-name change (LABELS, frontier_viz.R) can break it without a trace.
if (any(is.na(nudge_idx)))
  warning("TRACE_NUDGES labels not found on the plate: ",
          paste(TRACE_NUDGES$lab[is.na(nudge_idx)], collapse = "; "))
main_labs <- if (any(!is.na(nudge_idx)))
  all_labs[-nudge_idx[!is.na(nudge_idx)], ] else all_labs

# One fillable marker (shapes 21-24, so the ECI gradient stays on fill) per
# benchmark -- the diamond goes to GPQA Diamond, naturally.
SHAPES <- c(aime = 21, chess = 22, fm13 = 24, gpqa = 23)
names(SHAPES) <- LABELS[names(SHAPES)]

p <- ggplot(tl, aes(date, cost)) +
  # heavier than the dots' stroke, mid-grey rather than black: the traces
  # outrank the dots but must not upstage the teal or the text (a full-black
  # experiment did exactly that)
  geom_line(aes(group = interaction(bench, level)),
            colour = INK_MUTED, linewidth = 0.8) +
  geom_point(aes(shape = benchmark), fill = DOT, size = 2.4,
             colour = SURFACE, stroke = 0.3) +
  # nudge_y (in log10-dollar panel units) starts every label a step ABOVE
  # its dot, so the placements read consistently up-from-the-point; repel
  # still resolves collisions from there. Trace names get a LARGER lift than
  # model labels, so at a trace's opening dot the stack reads dot, model,
  # trace name on top.
  ggrepel::geom_text_repel(
    data = main_labs,
    aes(label = lab, size = sz, colour = col, segment.colour = seg),
    family = FONT,
    nudge_y = ifelse(grepl("^[0-9]+% on ", main_labs$lab), 0.5, 0.12),
    segment.size = 0.25, min.segment.length = 0.3,
    box.padding = 0.3, point.padding = 0.35, max.overlaps = Inf, seed = 1) +
  scale_size_identity() +
  scale_colour_identity() +
  # no key: the trace names on the plot already say which benchmark is
  # which, and the shapes only need to separate traces where they cross
  scale_shape_manual(values = SHAPES, guide = "none") +
  # The breaks already run to $100, but the panel stopped just above the
  # dearest run (~$2), so $10 never rendered. Pinning the ceiling AT $10
  # brings its tick and gridline in and gives the top trace some air; nothing
  # is clipped, the costliest record being well below it.
  # Zero expansion on top, so the panel ENDS at $10 rather than floating a
  # log-scale 5% above it; the bottom keeps its expansion so the cheapest
  # dots are not clipped against the axis.
  scale_y_log10(breaks = 10^(-5:2), labels = dollar_log,
                limits = c(NA, 10),
                expand = expansion(mult = c(0.05, 0))) +
  # Pinned on the RIGHT only: the panel ends on 2027 with no expansion, so
  # that break lands exactly at the edge and is labelled. The left end is
  # left to the data (plus the usual 5%), since forcing it back to 2023
  # opened a margin with nothing in it. Yearly breaks are given explicitly --
  # the default picks its own and skipped the endpoint.
  # as.Date(c(NA, ...)), not c(NA, as.Date(...)): c() dispatches on its first
  # argument, so a leading logical NA strips the Date class and the scale
  # rejects the limits.
  scale_x_date(limits = as.Date(c(NA, "2027-01-01")),
               date_breaks = "1 year", date_labels = "%Y",
               expand = expansion(mult = c(0.05, 0))) +
  # No title and no notes: the document and the deck both caption this plate
  # themselves, and stripping them here means the PNG and Figure 2.svg are
  # the same picture rather than the SVG being a trimmed copy.
  labs(x = "Release date of AI model", y = "Cost per task (log scale)") +
  frontier_theme(TIMELINE_TYPE_SCALE) +
  # the 2027 tick sits exactly on the panel's right edge (no expansion
  # there, above), so half its label hangs outside the panel; the theme's
  # slim outer margin cannot hold it and clipped it to "202"
  theme(plot.margin = margin(PLOT_PAD, 24, PLOT_PAD, PLOT_PAD))

# The hand-placed trace names (TRACE_NUDGES above): ONE repel layer per label,
# with scalar nudges, so no vector-to-row alignment can go wrong (a single
# layer with vector nudges once swapped two labels' placements).
# force_pull = 0: by default repel drags each label back toward its dot,
# honouring a long nudge only partway; a hand placement is meant literally.
for (i in which(!is.na(nudge_idx))) p <- p +
  ggrepel::geom_text_repel(
    data = all_labs[nudge_idx[i], ],
    aes(label = lab, size = sz, colour = col, segment.colour = seg),
    family = FONT,
    nudge_x = TRACE_NUDGES$x[i], nudge_y = TRACE_NUDGES$y[i],
    force_pull = 0,
    segment.size = 0.25, min.segment.length = 0.3,
    box.padding = 0.3, point.padding = 0.35, max.overlaps = Inf, seed = 1)

save_png(out_path("record_timelines.png"), p, width = 12, height = 8)

# Figure 2 of the report: one plate, so nothing is subset -- only the title
# and the notes come off. Kept at the figure's own 12 x 8 rather than the
# quartet's shape.
report_figure(p, 2, height = 8, width = 12)
cat(sprintf("wrote record_timelines.png (%d traces, one graph)\n",
            nrow(unique(tl[c("bench", "level")]))))

## ---- the figure's data, as CSV ------------------------------------------------
#
# One row per PLOTTED dot, in the order the traces are drawn, with the text
# the figure puts beside it. The HTML table above is the same series dressed
# for reading ($ and % formatting, Mystery Game Puzzles included); this is
# the series the FIGURE shows, raw, for reuse -- so it follows tl, not
# timelines, and Mystery is absent here exactly as it is from the plate.
#
# Two columns the table has no reason to carry, both about the labelling:
# label_drawn marks the row whose model name is actually rendered, since a
# run holding two of its benchmark's levels is plotted twice and labelled
# once; trace_label carries the standalone trace name, which appears only at
# the trace's opening dot.
csv <- tl[c("bench", "benchmark", "level", "model", "acc", "date", "cost")]
names(csv)[names(csv) == "model"] <- "model_label"
csv$label_drawn <- !duplicated(tl[c("bench", "model", "date")])
csv$trace_label <- ifelse(!duplicated(tl[c("bench", "level")]),
                          paste0(fmt_lev(tl$level), " on ", tl$benchmark), "")
# beside the PNG, not with the HTML table: this is the figure's data, and
# the two travel together
write.csv(csv, out_path("slides/Figure 2.csv"), row.names = FALSE)
cat("wrote slides/Figure 2.csv", "
", sep = "")
