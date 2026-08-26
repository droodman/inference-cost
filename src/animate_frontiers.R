# Animated quartets: one frame per grid date (semiannual), frontiers and runs
# accumulating.
#
# gganimate is deliberately not used. It is built for TWEENING between states,
# interpolating positions frame to frame; here the content is discrete -- a
# frontier appears, dots arrive -- and nothing should slide. Rendering one PNG
# per grid date and encoding them is simpler, exact, and reuses frontier_plot()
# unchanged, so frames are identical in style to the static figures.
#
# av bundles ffmpeg's libraries in its CRAN Windows binary, so no external
# ffmpeg install is needed.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("fit_specs.R")   # the shared spec grid; brings the other files with it

FPS       <- 1      # grid dates per second
HOLD_LAST <- 6      # extra frames on the final date, so it can be read
WIDTH     <- 1600
HEIGHT    <- 1200

d <- load_runs()
tbar  <- bench_tbar(d)     # from t - tc, so it matches whatever the fit used
dates <- bench_dates(d)
benches <- sort(unique(d$benchmark))

# Which specifications to animate. A_lin and A_quad are the pair worth watching
# side by side: they differ only in the quadratic time term and the lncost:tc
# interaction, and that difference is what drives the cheap-end reversal.
ANIMATE <- c("A_lin", "A_quad", "S_lin")

# One global frame timeline: the union of every benchmark's grid dates, so all
# four panels advance on the same clock even though their spans differ. One extra
# date is prepended, six months before the first, so the animation opens on
# empty axes and the viewer sees the grid before anything lands on it. Stepping
# back by one grid interval keeps the cadence uniform.
timeline <- sort(unique(do.call(c, unname(dates))))
timeline <- c(seq(min(timeline), by = "-6 months", length.out = 2)[2], timeline)

# Fixed axis ranges, so every frame -- including the empty opening one -- shares
# the panel geometry of the last. Without this the free_x scales would be read
# off whatever data a frame contains and the axes would crawl as runs accumulate.
axis_ranges <- do.call(rbind, lapply(benches, function(b) {
  s <- d[d$benchmark == b, ]
  data.frame(benchmark = b, cost = range(s$cost), value = c(0, 1))
}))

EMPTY_CURVES <- data.frame(cost = numeric(0), value = numeric(0),
                           qdate = as.Date(character(0)),
                           benchmark = character(0), year = numeric(0))

## ---- fits (once, on all data -- the animation reveals, it does not refit) ----------
#
# Refitting per frame would answer a different question ("what would we have
# believed at the time"), which is worth doing but is not what this shows. Here
# the fitted frontier is fixed and the animation walks along its time axis.

specs <- fit_all_specs(d)

## ---- frame rendering ------------------------------------------------------------------

# `curve_fn(upto)` returns the curves visible at date `upto`; points are always
# those released by `upto`. Colour limits are pinned to the full span so a given
# date keeps its colour as frames advance.
render_movie <- function(curve_fn, title, subtitle, ylab, notes, outfile,
                         step = FALSE) {
  clim <- range(2023 + as_t(timeline))
  dir <- file.path(tempdir(), paste0("frames_", tools::file_path_sans_ext(outfile)))
  dir.create(dir, showWarnings = FALSE)
  files <- character(0)

  for (i in seq_along(timeline)) {
    upto <- timeline[i]
    cur <- curve_fn(upto)
    pts <- d[d$releasedate <= upto, ]
    # Empty frames are rendered, not skipped: the opening frame is deliberately
    # bare, and axis_ranges keeps its panels identical to every later frame.
    if (is.null(cur) || !nrow(cur)) cur <- EMPTY_CURVES
    p <- frontier_plot(cur, pts, title = title,
                       subtitle = sprintf("%s   |   through %s", subtitle,
                                          format(upto, "%b %Y")),
                       ylab = ylab, notes = notes, step = step,
                       colour_limits = clim, ranges = axis_ranges)
    f <- file.path(dir, sprintf("f%04d.png", i))
    ggsave(f, p, width = WIDTH / 160, height = HEIGHT / 160, dpi = 160,
           device = ragg::agg_png)
    files <- c(files, f)
  }
  files <- c(files, rep(files[length(files)], HOLD_LAST))   # hold the final frame
  av::av_encode_video(files, output = out_path(outfile), framerate = FPS)
  cat("wrote", outfile, sprintf("(%d frames)\n", length(files)))
}

## ---- the movies -------------------------------------------------------------------------

# Parametric: curves for every grid date up to `upto`, so they accumulate.
param_curves <- function(fitset) function(upto) {
  dd <- lapply(dates, function(x) x[x <= upto])
  dd <- dd[vapply(dd, length, 1L) > 0]
  if (!length(dd)) return(NULL)
  frontier_curves(fitset[names(dd)], d, dd, tbar)
}

# Pareto: recomputed at each visible grid date, which is the honest thing -- the
# empirical frontier at date q genuinely only knows runs released by q.
# pareto_curves() itself is shared (frontier_viz.R); all this adds is trimming the
# date grid to the frame being drawn, so an animation frame and the corresponding
# static figure cannot disagree about where the staircase sits.
pareto_curves_upto <- function(upto) {
  dd <- lapply(dates, function(qs) qs[qs <= upto])
  dd <- dd[lengths(dd) > 0]
  if (!length(dd)) return(NULL)
  pareto_curves(d, dd)
}

NOTES <- c("Frontiers and runs accumulate semiannually; models are fitted once on the full sample.",
           "All models retained; data prepared once in prepare_data.R, which also drops duplicates.")

for (k in ANIMATE) {
  sp <- specs[[k]]
  render_movie(param_curves(sp$fits),
               sprintf("Fitted accuracy frontier by cost per task -- model %s (%s)",
                       sp$family, sp$time),
               sp$subtitle,
               if (sp$family == "S") "Fitted accuracy" else "Frontier accuracy",
               NOTES, sprintf("frontier_anim_%s.mp4", k))
}

render_movie(pareto_curves_upto,
             "Empirical Pareto frontier of accuracy by cost per task",
             "running maximum over runs released by each date",
             "Best accuracy achieved",
             c("Frontiers and runs accumulate semiannually; the frontier at each date uses only runs released by then.",
               "All models retained; data prepared once in prepare_data.R."),
             "frontier_anim_pareto.mp4", step = TRUE)
