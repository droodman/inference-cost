# Accuracy against release date, every run, one panel per benchmark.
#
# The plainest possible view of the raw material behind every fit in the
# analysis: no frontier, no model, no smoothing -- just where the runs are and
# which model produced them. It answers the questions the fitted figures
# cannot, like which model is responsible for an odd corner of a surface, and
# how thin a benchmark's early history really is.
#
# WHY THE DOTS FORM VERTICAL STRIPS. A release date is a property of the model,
# not of the run: every one of a model's runs -- each effort setting at each
# token budget, up to about sixty of them -- shares one date. So a model is a
# vertical column of dots whose height is the accuracy range its budget sweep
# covers, and the panel is a series of such columns. Reading a column bottom to
# top reads a model's truncation curve from starved to unconstrained.
#
# WHY ONE LABEL PER MODEL, NOT PER DOT. Given that structure, labeling dots
# would print the same name up to sixty times in a single vertical line. The
# label goes on the model's BEST run instead -- one name per column, at the top
# of what that model achieved -- which is both the useful anchor and the only
# one that leaves the panel readable. It still leaves up to 166 names in a
# panel (GPQA), so:
#
# WHY NO ggrepel. At GPQA's density the labels' own ink exceeds the area of the
# panel, so no placement algorithm can separate them -- repulsion would only
# spend minutes converging to a layout that is still overlapping, while
# dragging names arbitrarily far from the dots they name. Deterministic
# placement is the honest choice: each name is rotated to a vertical sliver
# (dates cluster hard in 2024-25, so thin-in-x is what the crowding demands)
# and anchored just above its own dot, so a name is always found where its
# model is. Names spread over accuracy as well as date, which resolves most
# collisions; what remains is genuine crowding, and the figure is drawn large
# enough to zoom into rather than pretending otherwise.
#
# The canvas is deliberately bigger than the analysis figures' 10 x fig_height:
# this is a reference plate to be opened and zoomed, not a report figure to be
# read at a glance.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("frontier_viz.R")   # load_runs, LABELS, bench_levels, theme, colors

d <- load_runs()
d$panel <- factor(LABELS[d$benchmark], levels = LABELS)

# One label per benchmark x model, on that model's best run: sort accuracy
# descending within the pair and keep the first of each. Named `peaks`, not
# `labs`, which is ggplot's own.
o     <- order(d$benchmark, d$model, -d$acc)
ds    <- d[o, ]
peaks <- ds[!duplicated(ds[c("benchmark", "model")]), ]

# A rotated name's length runs vertically, so it competes with the panel's
# HEIGHT -- and once width and aspect ratio are fixed, so is that height: 10
# inches at 1:2.04 leaves 2.26in of panel, while the longest name in the data
# (38 characters) is 2.33in at this size. Sizing every label so the worst case
# fits would mean dropping to 1.4, below what this figure has ever used, since
# the pinch is long names on NEAR-PERFECT runs: they start high on the axis
# and still need their full length above them.
#
# So the axis is sized for the 90th percentile instead of the maximum, and the
# panel edge clips the rest. Ten percent of names lose their tails; a name
# reads upward from its dot, so what survives is the START of the string --
# "claude-sonnet-4-5-2025..." is still recognisable. Trading a tenth of the
# names' tails keeps the other 800 at full size.
LABEL_SIZE <- 2.6
LABEL_FIT  <- 0.90   # share of names that must fit whole inside the panel
CHAR_EM    <- 0.60   # width of a character as a fraction of the font size

# Canvas width in inches. This, not the pixel count, is what decides whether
# the type survives a word processor: a document scales an image to its text
# column (about 6.5in), shrinking every font by FIG_W/6.5. At 16in that was a
# 2.5x reduction -- panel headings at 4pt, unreadable however many pixels back
# them -- so the figure sits at 10in, where the reduction is 1.5x.
FIG_W  <- 10

# Held at the 16 x 32.7in plate's proportions, so narrowing the canvas shortens
# it in step rather than leaving a 1:3.3 ribbon. Height follows from the width,
# and the row height from the height -- which is what squeezes LABEL_SIZE above.
ASPECT <- 32.7 / 16
n_row  <- ceiling(nlevels(droplevels(d$panel)) / 2)
FIG_H  <- FIG_W * ASPECT
ROW_H  <- (FIG_H - 1.5) / n_row

# How high the axis must reach for a given share of names to fit whole. A
# name's length runs VERTICALLY, L inches; on a panel of height h spanning
# [0, M] it covers (L/h)*M in accuracy units, so a model scoring a needs
#
#     a + (L/h)*M <= M   <=>   M >= a / (1 - L/h)
#
# Taking the LABEL_FIT quantile of that rather than the maximum is what lets
# the labels stay large: the top decile overruns the panel and is clipped.
# Never going below 1 keeps the whole 0-100% band on every benchmark, so
# headroom is added only where long names sit on high scores -- a benchmark
# topping out at 26% (Mystery) needs none and gets none, instead of carrying
# an empty band because some other benchmark's names are long.
y_top <- function(pk, plot_h) {
  L <- nchar(pk$model) * LABEL_SIZE * .pt * CHAR_EM / 72.27
  need <- pk$acc / pmax(0.02, 1 - L / plot_h)
  max(1, unname(quantile(need, LABEL_FIT)))
}

# One y scale across all panels, so the headroom is the largest any
# benchmark needs.
p <- ggplot(d, aes(releasedate, acc)) +
  geom_point(colour = INK_MUTED, size = 0.5, alpha = 0.35) +
  geom_point(data = peaks, colour = PALETTE[length(PALETTE)], size = 0.7) +
  geom_text(data = peaks, aes(label = model), angle = 90, hjust = -0.07,
            size = LABEL_SIZE, colour = INK_SECOND, lineheight = 0.85) +
  facet_wrap(~panel, scales = "free_x", ncol = 2, drop = FALSE) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y",
               expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(limits = c(0, y_top(peaks, ROW_H - 0.9)),
                     breaks = seq(0, 1, 0.25),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Accuracy by model release date: every run, every benchmark",
       x = NULL, y = "Accuracy",
       # pad_caption pads to CAPTION_LINES but never WRAPS, so a paragraph
       # longer than the canvas runs off the right edge. Wrap each to what
       # FIG_W inches of 7.5pt caption holds; three paragraphs at two lines
       # each is exactly the six-line budget.
       caption = pad_caption(unlist(lapply(c(
         paste("Every run in the analysis. A model's release date is the same",
               "for all its runs, so each model is a VERTICAL COLUMN of dots",
               "spanning the accuracy range its effort and token-budget sweep",
               "covers; the darker dot tops each column."),
         paste("Names label one dot per model -- its best run -- not one per",
               "run, which would repeat each name up to sixty times. They read",
               "upward from the dot they name; where a panel is crowded they",
               "overlap, and the longest tenth are cut off at the panel top,",
               "leaving the start of the name. Open the file at full size."),
         paste("Accuracy is rescaled from each benchmark's guessing floor to 1",
               "(prepare_data.R), so 0 means no better than chance. The band",
               "above 100% is label space, not attainable accuracy.")),
         strwrap, width = round(FIG_W * 18))))) +
  frontier_theme() +
  theme(axis.text = element_text(colour = "white"))

f <- "accuracy_scatter.png"
ggsave(out_path(f), p, width = FIG_W, height = FIG_H, dpi = 200,
       limitsize = FALSE, device = ragg::agg_png)
cat("wrote", f, "\n")
cat(sprintf("  %d runs, %d model-benchmark labels over %d panels\n",
            nrow(d), nrow(peaks), nlevels(droplevels(d$panel))))
# say how many names the panel edge cut, so the trade stays visible rather
# than being something a reader has to notice for themselves
{
  ph   <- ROW_H - 0.9
  L    <- nchar(peaks$model) * LABEL_SIZE * .pt * CHAR_EM / 72.27
  need <- peaks$acc / pmax(0.02, 1 - L / ph)
  yt   <- y_top(peaks, ph)
  cat(sprintf("  %.1f x %.1f in (1:%.2f), axis to %.2f, %d of %d names clipped (%.0f%%)\n",
              FIG_W, FIG_H, ASPECT, yt, sum(need > yt), nrow(peaks),
              100 * mean(need > yt)))
}
for (b in bench_levels(d$benchmark))
  cat(sprintf("  %-26s %4d runs  %3d models\n", b, sum(d$benchmark == b),
              sum(peaks$benchmark == b)))

