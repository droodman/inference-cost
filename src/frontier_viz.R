# Shared plotting machinery for the frontier figures.
#
# Sourced by plot_frontiers.R (parametric SFA + plain fractional logit),
# pareto_frontiers.R (nonparametric running maximum) and animate_frontiers.R.
# Each caller supplies curves and points; everything about how a frontier figure
# LOOKS lives here, so the figures cannot drift apart.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
library(ggplot2)
src_source("prepare_data.R")   # build_runs(): the single source of analysis data

## ---- palette (sequential blue, ordinal steps 250-700; light surface) -----------

BLUE <- c("#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
          "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b")
INK_PRIMARY <- "#0b0b0b"
INK_SECOND  <- "#52514e"
INK_MUTED   <- "#898781"
GRIDLINE    <- "#e1e0d9"
AXIS        <- "#c3c2b7"
SURFACE     <- "#fcfcfb"

LABELS <- c(aime = "AIME (OTIS Mock)", chess = "Chess Puzzles",
            fm13 = "FrontierMath, tiers 1-3", gpqa = "GPQA Diamond")

EPOCH <- as.Date("2023-01-01")          # t = 0
as_t  <- function(date) as.numeric(as.Date(date) - EPOCH) / 365.25

dollar_log <- function(x) {
  ifelse(is.na(x), "", paste0("$", formatC(x, format = "fg", drop0trailing = TRUE)))
}

## ---- date grid ------------------------------------------------------------------

# Quarterly cadence anchored to the LAST observation and stepped backwards, not
# to calendar quarter starts. The final curve then lands exactly on the last run
# (so no observation sits outside every frontier) while consecutive curves stay
# exactly three months apart -- the spacing between curves is read as elapsed
# time, so it has to be uniform.
quarter_dates <- function(first, last) {
  n <- ceiling(as.numeric(as.Date(last) - as.Date(first)) / 365.25 * 4) + 1
  qs <- seq(as.Date(last), by = "-3 months", length.out = n)
  sort(qs[qs >= as.Date(first)])
}

## ---- model-agnostic frontier index ------------------------------------------------

# maxLik fits name coefficients beta_*, glm fits do not; strip the prefix so one
# accessor serves both. Absent terms return 0, so a fit without the quadratic
# evaluates as if its coefficient were zero rather than erroring.
frontier_coefs <- function(fit) {
  cf <- coef(fit)
  names(cf) <- sub("^beta_", "", names(cf))
  get1 <- function(nm) if (nm %in% names(cf)) unname(cf[[nm]]) else 0
  c(b0 = get1("(Intercept)"), bx = get1("lncost"),
    bt = get1("tc"), btt = get1("I(tc^2)"),
    bxt = get1("lncost:tc"), bxx = get1("I(lncost^2)"))
}

# Every term any specification in the repo can carry. A coefficient vector
# missing one contributes zero for it (see the `get*` accessors), so linear,
# quadratic and interaction fits all evaluate through this one function --
# which is the point: a term present in the FIT but absent HERE would be
# silently dropped, and the surface scored would not be the surface fitted.
frontier_index <- function(co, lncost, tc) {
  gete <- function(nm) if (nm %in% names(co)) co[[nm]] else 0
  gete("b0") + gete("bx") * lncost + gete("bt") * tc +
    gete("btt") * tc^2 + gete("bxt") * lncost * tc +
    gete("bxx") * lncost^2
}

# Slope of the frontier in log cost at a given date. Constant without the
# interaction; with it, this is what must stay positive for the frontier to obey
# free disposal -- more spend cannot buy less accuracy.
cost_slope <- function(co, tc) co[["bx"]] + co[["bxt"]] * tc

# Fitted frontier (u = 0) over a cost grid, for each date in `dates`.
frontier_curves <- function(fitset, data, dates_by_bench, tbar, n_cost = 200) {
  do.call(rbind, lapply(names(fitset), function(b) {
    co <- frontier_coefs(fitset[[b]])
    sub <- data[data$benchmark == b, ]
    cost <- exp(seq(log(min(sub$cost)), log(max(sub$cost)), length.out = n_cost))
    g <- expand.grid(cost = cost, qdate = dates_by_bench[[b]])
    g$value <- plogis(frontier_index(co, log(g$cost), as_t(g$qdate) - tbar[[b]]))
    g$benchmark <- b
    g$year <- 2023 + as_t(g$qdate)
    g
  }))
}

## ---- caption block ----------------------------------------------------------------

# The caption sits below the panels, so its height comes out of them: a figure
# with two notes gets a taller plotting area than one with five, and flipping
# between figures in the viewer makes the panels jump. Every caption is padded to
# the same number of lines instead, so the drawing area is identical everywhere
# and only the notes themselves change.
#
# Five is the current maximum (envelope and Pareto figures); adding a sixth note
# anywhere means raising this, and the warning below says so rather than letting
# one figure quietly grow its caption and shrink its panels.
CAPTION_LINES <- 5

pad_caption <- function(notes) {
  lines <- unlist(strsplit(as.character(notes), "\n", fixed = TRUE))
  if (length(lines) > CAPTION_LINES) {
    warning(sprintf(paste("caption has %d lines but CAPTION_LINES is %d;",
                          "this figure's panels will be shorter than the rest"),
                    length(lines), CAPTION_LINES))
  }
  # Padding with a space, not "", because a trailing empty line contributes no
  # height -- the block would collapse back to the number of real notes.
  paste(c(lines, rep(" ", max(0, CAPTION_LINES - length(lines)))),
        collapse = "\n")
}

## ---- the figure -------------------------------------------------------------------

# `curves` needs cost, value, qdate, benchmark, year; `pts` needs cost, acc, year,
# benchmark. `step = TRUE` draws staircases (Pareto), FALSE draws lines (fitted).
# `ranges` (benchmark, cost, value) pins each panel's axes via an invisible layer.
# Needed whenever a frame may hold little or no data -- an animation's opening
# frame has none at all, and with scales = "free_x" a panel's range would
# otherwise be read off whatever happens to be present, so the axes would drift
# from frame to frame.
frontier_plot <- function(curves, pts, title, subtitle, ylab,
                          notes = character(0), step = FALSE,
                          colour_limits = NULL, ranges = NULL) {
  curves$benchmark <- factor(LABELS[curves$benchmark], levels = LABELS)
  pts$benchmark    <- factor(LABELS[pts$benchmark],    levels = LABELS)
  if (is.null(colour_limits)) colour_limits <- range(curves$year)
  blank_layer <- NULL
  if (!is.null(ranges)) {
    ranges$benchmark <- factor(LABELS[ranges$benchmark], levels = LABELS)
    blank_layer <- geom_blank(data = ranges, aes(cost, value), inherit.aes = FALSE)
  }

  geom_curve_layer <- if (step) {
    geom_step(aes(group = qdate, colour = year), direction = "hv", linewidth = 0.6)
  } else {
    geom_line(aes(group = qdate, colour = year), linewidth = 0.6)
  }

  base_notes <- c(
    paste("Curves are quarterly, stepped back from each benchmark's last run so",
          "spacing is exactly three months and the final curve lands on the last",
          "observation."),
    paste("Dots are observed runs, on the same time scale as the curves."))

  ggplot(curves, aes(cost, value)) +
    blank_layer +
    geom_point(data = pts, aes(cost, acc, colour = year), size = 0.35, alpha = 0.3,
               inherit.aes = FALSE) +
    geom_curve_layer +
    facet_wrap(~benchmark, scales = "free_x", nrow = 2, drop = FALSE) +
    scale_x_log10(breaks = 10^(-5:1), labels = dollar_log) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_colour_gradientn(
      colours = BLUE, name = NULL, breaks = 2023:2026, limits = colour_limits,
      guide = guide_colourbar(barheight = grid::unit(0.35, "cm"),
                              barwidth = grid::unit(7, "cm"),
                              direction = "horizontal", ticks.colour = SURFACE)) +
    labs(title = title, subtitle = subtitle,
         x = "Cost per task (log scale)", y = ylab,
         caption = pad_caption(c(base_notes, notes))) +
    frontier_theme()
}

# Shared chrome, so every figure in the repo reads as one family.
frontier_theme <- function() {
  theme_minimal(base_size = 11, base_family = "sans") +
    theme(
      plot.background  = element_rect(fill = SURFACE, colour = NA),
      panel.background = element_rect(fill = SURFACE, colour = NA),
      panel.grid.major = element_line(colour = GRIDLINE, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line  = element_line(colour = AXIS, linewidth = 0.3),
      axis.text  = element_text(colour = INK_MUTED, size = 8),
      axis.title = element_text(colour = INK_SECOND, size = 9),
      strip.text = element_text(colour = INK_PRIMARY, face = "bold", size = 10,
                                hjust = 0),
      plot.title    = element_text(colour = INK_PRIMARY, face = "bold", size = 13),
      plot.subtitle = element_text(colour = INK_SECOND, size = 9.5),
      plot.caption  = element_text(colour = INK_MUTED, size = 7.5, hjust = 0),
      legend.position = "top", legend.justification = "left",
      legend.text = element_text(colour = INK_SECOND, size = 8),
      plot.margin = margin(12, 16, 10, 12))
}

## ---- iso-accuracy cost contours ------------------------------------------------
#
# The frontier inverted for cost at a fixed accuracy target:
#   ln c(t) = [logit(y) - b0 - b_tc*tc - b_tc2*tc^2] / b_lncost
# so each contour traces what a given performance level costs over time. Colour
# now carries accuracy rather than date -- still a magnitude, so the same
# sequential blue ramp applies, light (low) to dark (high).
#
# Contours are blanked outside each benchmark's observed cost range: the
# inversion is happy to report that 95% accuracy costs $10^4, but that is the
# functional form extrapolating, not a finding. Lines simply stop instead.
# With an lncost:tc interaction the inversion gains a date-dependent denominator:
#   ln c(t) = [logit(y) - b0 - b_t*tc - b_tt*tc^2] / (b_x + b_xt*tc)
# which is singular where the cost slope crosses zero and carries the wrong sign
# beyond it. Those dates are blanked rather than drawn: a frontier sloping the
# wrong way in cost is a functional-form artefact, not a cheaper frontier.
# Solving z(ln c, tc) = logit(y) for ln c. Without a (ln c)^2 term this is a
# division; with one it is a quadratic
#     bxx*u^2 + (bx + bxt*tc)*u + (b0 + bt*tc + btt*tc^2 - logit(y)) = 0
# whose two roots straddle the surface's turning point. Monotonicity is imposed
# over the observed cost range, so at most one root lies inside it -- take that
# one, and return NA when neither does (the target is unreachable there).
iso_cost_curves <- function(fitset, data, tbar,
                            levels = seq(0.05, 0.95, by = 0.10), n_date = 300,
                            min_slope = 0.05) {
  do.call(rbind, lapply(names(fitset), function(b) {
    co  <- frontier_coefs(fitset[[b]])
    sub <- data[data$benchmark == b, ]
    urng <- range(sub$lncost)
    dts <- seq(min(sub$releasedate), max(sub$releasedate), length.out = n_date)
    g   <- expand.grid(date = dts, acc = levels)
    tc  <- as_t(g$date) - tbar[[b]]
    aa <- co[["bxx"]]
    bb <- cost_slope(co, tc)                       # bx + bxt*tc
    cc <- co[["b0"]] + co[["bt"]] * tc + co[["btt"]] * tc^2 - qlogis(g$acc)

    if (abs(aa) < 1e-10) {
      u <- -cc / bb
      u[bb < min_slope] <- NA_real_                # slope too flat to invert
    } else {
      disc <- bb^2 - 4 * aa * cc
      r1 <- (-bb + sqrt(pmax(disc, 0))) / (2 * aa)
      r2 <- (-bb - sqrt(pmax(disc, 0))) / (2 * aa)
      inr <- function(x) is.finite(x) & x >= urng[1] & x <= urng[2]
      u <- ifelse(inr(r1), r1, ifelse(inr(r2), r2, NA_real_))
      u[disc < 0] <- NA_real_
    }
    g$cost <- exp(u)
    g$cost[g$cost < min(sub$cost) | g$cost > max(sub$cost)] <- NA_real_
    g$benchmark <- b
    g
  }))
}

# `pts` are the observed runs, placed at (release date, cost) and coloured by the
# accuracy they achieved -- the same scale as the contours, so a run's colour can
# be read against the contour it sits on.
iso_cost_plot <- function(curves, pts, title, subtitle, notes = character(0),
                          ranges = NULL) {
  curves$benchmark <- factor(LABELS[curves$benchmark], levels = LABELS)
  pts$benchmark    <- factor(LABELS[pts$benchmark],    levels = LABELS)
  blank_layer <- NULL
  if (!is.null(ranges)) {
    ranges$benchmark <- factor(LABELS[ranges$benchmark], levels = LABELS)
    blank_layer <- geom_blank(data = ranges, aes(date, cost), inherit.aes = FALSE)
  }

  ggplot(curves, aes(date, cost)) +
    blank_layer +
    geom_point(data = pts, aes(releasedate, cost, colour = acc),
               size = 0.35, alpha = 0.3, inherit.aes = FALSE) +
    geom_line(aes(group = acc, colour = acc), linewidth = 0.6, na.rm = TRUE) +
    facet_wrap(~benchmark, scales = "free_y", nrow = 2, drop = FALSE) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    scale_y_log10(breaks = 10^(-5:1), labels = dollar_log) +
    scale_colour_gradientn(
      colours = BLUE, name = NULL, limits = c(0, 1),
      breaks = seq(0, 1, 0.25), labels = scales::percent_format(accuracy = 1),
      guide = guide_colourbar(barheight = grid::unit(0.35, "cm"),
                              barwidth = grid::unit(7, "cm"),
                              direction = "horizontal", ticks.colour = SURFACE)) +
    labs(title = title, subtitle = subtitle,
         x = "Model release date", y = "Cost per task (log scale)",
         caption = pad_caption(notes)) +
    frontier_theme()
}

## ---- shared data prep ---------------------------------------------------------------

# Analysis sample: build_runs() supplies the data (already deduplicated there,
# mirroring the .do file's `duplicates drop`), and this adds only the derived
# time variables. `drop_gpt4o_chess` is explicit rather than assumed -- the
# figures have differed on it before, and a silent default is how they drifted.
load_runs <- function(drop_gpt4o_chess = FALSE) {
  d <- build_runs()
  d$t <- as_t(d$releasedate)
  d <- d[stats::complete.cases(d[c("acc", "lncost", "t", "model", "effort")]), ]
  if (drop_gpt4o_chess)
    d <- d[!(d$benchmark == "chess" & d$model == "gpt-4o"), ]
  # demeaned within benchmark: beta_tc is then the improvement rate at each
  # benchmark's own reference date, and tc is decorrelated from its square
  d$tc <- ave(d$t, d$benchmark, FUN = function(x) x - mean(x))
  d$year <- 2023 + d$t

  # prepare_data.R already dropped duplicates; the SFA fits are therefore called
  # with dedup = FALSE. Warn rather than silently re-dropping if any survive, so
  # a change upstream shows up instead of being papered over here.
  dup <- duplicated(d[c("benchmark", "model", "effort", "acc", "cost", "releasedate")])
  if (any(dup))
    warning(sprintf("%d duplicate row(s) survived prepare_data.R", sum(dup)))
  d
}

# Reference date per benchmark, recovered from the data rather than recomputed:
# t - tc is exactly the demeaning constant, whatever sample tc was built on. A
# separately-computed mean(t) would silently disagree if rows were dropped after
# tc was assigned, shifting the plotted curves off the model that produced them.
bench_tbar <- function(d) {
  tapply(d$t - d$tc, d$benchmark, function(x) x[1])
}

# ONE global quarterly grid, anchored to the latest run across all benchmarks;
# each benchmark takes the dates at or after its own first run. Anchoring per
# benchmark instead would give fm13 its own grid (its last run is 3 days before
# the others'), putting its curves on dates no other panel has -- frames where
# three panels sit frozen while one advances. A single anchor makes the panels
# share a clock exactly.
#
# The upper end is deliberately NOT capped per benchmark: the final date covers
# every benchmark's last run, so no data falls outside the last frame. That holds
# because all four end within days of each other. A benchmark ending much earlier
# would have its fitted curve extrapolated forward -- cap it here if that arises.
bench_dates <- function(d) {
  benches <- sort(unique(d$benchmark))
  grid <- quarter_dates(min(d$releasedate), max(d$releasedate))
  setNames(lapply(benches, function(b) {
    grid[grid >= min(d$releasedate[d$benchmark == b])]
  }), benches)
}
