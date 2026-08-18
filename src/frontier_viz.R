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

# The two partial derivatives of the frontier index.
#
#   d z / d ln c = b_x  + 2*b_xx*ln c + b_xt*tc      free disposal
#   d z / d tc   = b_t  + 2*b_tt*tc   + b_xt*ln c    models stay available
#
# Both must stay >= 0 for the fitted surface to be a frontier at all. Note that
# each depends on BOTH coordinates once the full quadratic is in play: the cost
# slope moves with date through b_xt and with cost itself through b_xx. An
# earlier version took tc alone and returned b_x + b_xt*tc, which silently
# stopped being the cost slope the moment a (ln c)^2 term was added.
frontier_dcost <- function(co, lncost, tc) {
  co[["bx"]] + 2 * co[["bxx"]] * lncost + co[["bxt"]] * tc
}

frontier_dtime <- function(co, lncost, tc) {
  co[["bt"]] + 2 * co[["btt"]] * tc + co[["bxt"]] * lncost
}

# Each derivative is linear in (ln c, tc) jointly, so it is non-negative
# everywhere on a rectangle exactly when it is non-negative at all four corners.
# That is what makes the check below cheap and exact rather than a grid search,
# and it is the same fact fit_envelope() uses to impose monotonicity.
corner_grid <- function(lnc_range, tc_range) {
  expand.grid(lncost = range(lnc_range), tc = range(tc_range))
}

# Worst case of both derivatives over the observed rectangle.
frontier_slope_bounds <- function(co, lncost, tc) {
  g <- corner_grid(lncost, tc)
  c(dcost = min(frontier_dcost(co, g$lncost, g$tc)),
    dtime = min(frontier_dtime(co, g$lncost, g$tc)))
}

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
# Six is the current maximum (the envelope figure, which must explain the solid
# curve, the dashed staircase, monotonicity and the fixed grid). Adding a seventh
# note anywhere means raising this, and the warning below says so rather than
# letting one figure quietly grow its caption and shrink its panels.
CAPTION_LINES <- 6

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

## ---- the empirical Pareto staircase, as an overlay ---------------------------------

# P_t(c) = max{ a_i : c_i <= c, t_i <= t }, the nonparametric frontier that
# pareto_frontier.png draws on its own and that the fitted figures draw underneath
# themselves for comparison. It lives here, not in the scripts, because three
# figures now need it and two of them had grown identical private copies.
#
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

# One staircase per benchmark per date, in the column layout frontier_plot()
# expects. Dates are taken by index so `q` keeps its Date class whatever the
# iteration style: lapply over a Date vector preserves it, but `for (q in dts)`
# strips it and hands the body a bare number for as_t() to reinterpret, which is
# what the callers' own loops do.
pareto_curves <- function(data, dates_by_bench) {
  do.call(rbind, lapply(names(dates_by_bench), function(b) {
    sub <- data[data$benchmark == b, ]
    qs <- dates_by_bench[[b]]
    do.call(rbind, lapply(seq_along(qs), function(k) {
      q <- qs[k]
      st <- pareto_steps(sub, q)
      if (is.null(st)) return(NULL)
      st$qdate <- q
      st$year <- 2023 + as_t(q)
      st$benchmark <- b
      st
    }))
  }))
}

# The dashed staircase drawn OVER a fitted frontier figure. Returned as a layer
# rather than described in prose in each script, so every figure that shows a
# fitted curve against the empirical frontier shows it identically -- same dash,
# same weight, same colour scale.
pareto_step_layer <- function(steps) {
  steps$benchmark <- factor(LABELS[steps$benchmark], levels = LABELS)
  geom_step(data = steps, aes(cost, value, group = qdate, colour = year),
            direction = "hv", linewidth = 0.4, linetype = "22",
            inherit.aes = FALSE)
}

# The note that must accompany that layer, so the caption cannot drift from what
# is drawn.
# Says what the staircase IS, and no longer that the fit "approximates" it. That
# phrasing invited the reader to expect the curve to hug the staircase at each
# date, which no fit here promises: the staircase is a running maximum carrying
# old records forward, while the curve is the surface at one instant.
PARETO_STEP_NOTE <- paste(
  "Dashed: the empirical Pareto staircase P_t(c) = max{a_i : c_i <= c, t_i <= t}",
  "at the same date; a running maximum, so it carries older records forward.")

# Belongs with any iso-accuracy figure whose fit can bend in cost. The contour is
# the whole level set, not a function of date: where the surface bends in cost a
# target is met at two costs, and the curve simply turns back on itself.
ISO_BRANCH_NOTE <- paste(
  "Contours are complete iso-accuracy sets. Where the surface bends in cost a",
  "target is met at two costs, so a contour can turn back on itself rather than",
  "run left to right.")

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
#     aa*u^2 + bb*u + cc = 0,  aa = bxx, bb = bx + bxt*tc,
#                              cc = b0 + bt*tc + btt*tc^2 - logit(y)
# whose two roots straddle the surface's turning point in cost.
#
# WHICH ROOT is the whole difficulty, and getting it wrong is visible: the code
# here used to take whichever root happened to fall inside the observed cost
# range, on the reasoning that monotonicity meant only one ever would. That holds
# for the envelope, where monotonicity is imposed -- but the SFA and logit fits
# only have it checked, and under the full quadratic they fail it on every
# benchmark. Both roots then sit inside the range, the preference silently
# switched from one to the other as the date moved, and the contour jumped
# between branches: a near-vertical segment spanning three orders of magnitude in
# cost, which looks like a finding and is an artefact of root bookkeeping.
#
# BOTH roots are drawn, each as its own contour, and they are told apart rather
# than mixed. Differentiating, dz/du = bb + 2*aa*u, so at the two roots
#     u+ = (-bb + sqrt(disc)) / (2*aa)   ->   dz/du = +sqrt(disc)
#     u- = (-bb - sqrt(disc)) / (2*aa)   ->   dz/du = -sqrt(disc)
# whatever the sign of aa -- the aa cancels. So u+ is ALWAYS the branch where
# accuracy rises with spend ("rising") and u- always the branch where it falls
# ("falling"), which is what makes them separable at all.
#
# Drawing only u+ would be defensible for a frontier, but these fits describe the
# whole distribution of runs, where a region in which more spend buys less is a
# pattern worth seeing rather than a violation worth hiding. So both go on the
# figure, distinguished by linetype, and they meet at the turning point where
# sqrt(disc) = 0 -- the contour reads as a fold, which is what it is.
#
# Blanked in two cases, both honest gaps rather than lines:
#   disc < 0        the target is not attained at any cost on that date
#   u outside urng  the answer would be extrapolation beyond observed costs
# Grouping each branch's contiguous stretches separately (iso_segments) is what
# keeps the two from ever being joined into one line, which is how the old
# single-branch code produced its near-vertical jumps.
iso_acc_curves <- function(fitset, data, tbar,
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
    # NOT the cost slope: this is the coefficient on u in the quadratic below.
    # The two coincide only when bxx = 0, where the quadratic degenerates to the
    # division and the u-coefficient IS the derivative.
    bb <- co[["bx"]] + co[["bxt"]] * tc
    cc <- co[["b0"]] + co[["bt"]] * tc + co[["btt"]] * tc^2 - qlogis(g$acc)

    if (abs(aa) < 1e-10) {
      # No curvature in cost: one root, and the min_slope guard is needed here
      # because the inversion is a division that blows up as bb -> 0.
      u <- -cc / bb
      u[bb < min_slope] <- NA_real_
      roots <- list(rising = u)
      disc <- rep(NA_real_, nrow(g))
    } else {
      disc  <- bb^2 - 4 * aa * cc
      slope <- sqrt(pmax(disc, 0))                 # |dz/du| at either root
      roots <- list(rising  = (-bb + slope) / (2 * aa),
                    falling = (-bb - slope) / (2 * aa))
      # disc < 0 means the target is not attained at ANY cost on that date
      roots <- lapply(roots, function(u) { u[disc < 0] <- NA_real_; u })
      # No min_slope guard on this branch: nothing blows up here (the roots stay
      # bounded, and a near-linear surface sends the spurious root outside the
      # cost range, where the clip below removes it), and gapping at sqrt(disc)
      # ~ 0 would punch a hole exactly at the turning point where the two
      # branches meet -- the one place the contour is genuinely continuous.
    }
    both <- do.call(rbind, lapply(names(roots), function(br) {
      u <- roots[[br]]
      # "to the extent they are in the range of actual data": the observed cost
      # range is the only clip either branch gets.
      u[is.finite(u) & (u < urng[1] | u > urng[2])] <- NA_real_
      h <- g
      h$cost <- exp(u)
      h$branch <- br
      h$disc <- disc
      h
    }))
    both$benchmark <- b
    # ordered as ONE curve per accuracy level, not two branches
    iso_segments(both)
  }))
}

# Contours are blanked wherever the inversion has no honest answer, and the gaps
# can fall in the MIDDLE of a contour. geom_line drops NA rows and joins whatever
# is left, so an interior gap gets drawn as a straight chord across it -- a
# segment that looks like a fitted contour but is purely an artefact of the hole
# it spans. Numbering the contiguous stretches and grouping on that number makes
# the line break at the gap instead, which is what a gap should look like.
iso_segments <- function(g) {
  # nominal spacing of the date grid: the step between neighbouring samples of
  # the same branch. Anything larger between consecutive drawn points means
  # something was skipped, and the line must break rather than bridge it.
  dd <- sort(unique(as.numeric(g$date)))
  spacing <- if (length(dd) > 1) min(diff(dd)) else 1

  do.call(rbind, lapply(split(g, g$acc, drop = TRUE), function(s) {
    dnum <- as.numeric(s$date)

    # Do the two roots actually MEET inside the drawn window? They coincide where
    # disc = 0, so the fold is in the window only if disc dips below zero
    # somewhere in it AND is non-negative somewhere else -- a level with disc < 0
    # at every date is unreachable throughout, every cost is already NA, and
    # there is nothing to orient. If disc stays positive throughout, the target
    # is met at two costs at every date and the turning point lies outside the
    # panel: the branches are then separate pieces of one curve, and joining them
    # would draw a chord up the edge of the panel rather than a fold.
    ok <- if (is.null(s$disc)) integer(0) else which(s$disc >= 0)
    fold_inside <- length(ok) > 0 && any(s$disc < 0, na.rm = TRUE)

    if (fold_inside) {
      # Orient so the path JOINS at the fold: in along one branch, round the
      # turning point, back out along the other. The fold sits at whichever end
      # of the feasible range disc is smaller.
      fold_at_start <- s$disc[ok[which.min(dnum[ok])]] <=
        s$disc[ok[which.max(dnum[ok])]]
      sgn_fall <- if (fold_at_start) -1 else 1
      key <- ifelse(s$branch == "falling", 1L, 2L)
      sgn <- ifelse(s$branch == "falling", sgn_fall, -sgn_fall)
      s <- s[order(key, sgn * dnum), ]
    } else {
      s <- s[order(s$branch, dnum), ]
    }

    s <- s[!is.na(s$cost), ]
    if (!nrow(s)) return(s)
    # One break rule covers everything: consecutive drawn points that are not
    # neighbours on the date grid start a new segment. At the fold the two
    # branches share a date, so the step is zero and the curve stays whole; a
    # cost-range gap, or the far ends of two branches that never meet, step much
    # further and break.
    s$seg <- cumsum(c(FALSE,
                      abs(diff(as.numeric(s$date))) > 1.5 * spacing)) + 1L
    s
  }))
}

# `pts` are the observed runs, placed at (release date, cost) and coloured by the
# accuracy they achieved -- the same scale as the contours, so a run's colour can
# be read against the contour it sits on.
iso_acc_plot <- function(curves, pts, title, subtitle, notes = character(0),
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
    # geom_path, NOT geom_line: geom_line reorders by x, which would undo the
    # path ordering and reduce the fold back to two overlapping branches. Grouped
    # on the contiguous stretch so a blanked interior still breaks the curve;
    # branch is deliberately NOT in the grouping, because the two roots are one
    # isoquant and are drawn as one line.
    geom_path(aes(group = interaction(benchmark, acc, seg), colour = acc),
              linewidth = 0.6, na.rm = TRUE) +
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
