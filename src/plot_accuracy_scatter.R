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

# Headroom for the rotated names, which run upward from their dots: the axis
# is still labeled only to 100%, so the empty band above it reads as margin
# rather than as accuracy above perfect. It has to grow with LABEL_SIZE, since
# a rotated name's length is vertical: at this size the longest name in the
# data (38 characters) spans about 0.8 of the accuracy range. The binding
# cases are long names on near-perfect runs -- claude-sonnet-4-5 and the
# gemini-2.5-pro previews on MATH level 5 -- which reach about 1.54; below
# that they are truncated at the panel edge.
LABEL_SIZE <- 2.6
Y_MAX      <- 1.58

# Per-row height in inches. Also grown with LABEL_SIZE: the extra headroom
# above would otherwise be taken out of the 0-100% band the dots live in,
# which would undo half of the size increase.
ROW_H <- 5.2

p <- ggplot(d, aes(releasedate, acc)) +
  geom_point(colour = INK_MUTED, size = 0.5, alpha = 0.35) +
  geom_point(data = peaks, colour = PALETTE[length(PALETTE)], size = 0.7) +
  geom_text(data = peaks, aes(label = model), angle = 90, hjust = -0.07,
            size = LABEL_SIZE, colour = INK_SECOND, lineheight = 0.85) +
  facet_wrap(~panel, scales = "free_x", ncol = 2, drop = FALSE) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y",
               expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(limits = c(0, Y_MAX), breaks = seq(0, 1, 0.25),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Accuracy by model release date: every run, every benchmark",
       x = NULL, y = "Accuracy",
       caption = pad_caption(c(
         paste("Every run in the analysis. A model's release date is the same",
               "for all its runs, so each model is a VERTICAL COLUMN of dots",
               "spanning the accuracy range its effort and token-budget sweep",
               "covers; the darker dot tops each column."),
         paste("Names label one dot per model -- its best run -- not one per",
               "run, which would repeat each name up to sixty times. They are",
               "placed deterministically, reading upward from the dot they",
               "name; where a panel is crowded the names overlap, and the",
               "plate is sized to be zoomed."),
         paste("Accuracy is rescaled from each benchmark's guessing floor to",
               "1 (prepare_data.R), so 0 means no better than chance. The band",
               "above 100% is label space, not attainable accuracy.")))) +
  frontier_theme()

n_row <- ceiling(nlevels(droplevels(d$panel)) / 2)
f <- "accuracy_scatter.png"
ggsave(out_path(f), p, width = 16, height = 1.5 + ROW_H * n_row, dpi = 200,
       limitsize = FALSE, device = ragg::agg_png)
cat("wrote", f, "\n")
cat(sprintf("  %d runs, %d model-benchmark labels over %d panels\n",
            nrow(d), nrow(peaks), nlevels(droplevels(d$panel))))
for (b in bench_levels(d$benchmark))
  cat(sprintf("  %-26s %4d runs  %3d models\n", b, sum(d$benchmark == b),
              sum(peaks$benchmark == b)))
