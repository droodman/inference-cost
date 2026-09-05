# Shared plotting machinery for the frontier figures.
#
# Sourced by plot_frontiers.R (parametric SFA + plain fractional logit) and
# pareto_frontiers.R (nonparametric running maximum).
# Each caller supplies curves and points; everything about how a frontier figure
# LOOKS lives here, so the figures cannot drift apart.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
library(ggplot2)
src_source("prepare_data.R")   # build_runs(): the single source of analysis data

## ---- theme + palette ------------------------------------------------------------
#
# EXPERIMENT (round 2): dark surface with plasma. One flag switches the whole
# look, because palette and chrome cannot be chosen separately: which END of a
# palette is legible depends on the background. On the light surface viridis had
# to be REVERSED (its yellow end vanished into white, so dark took the high
# values); on black the failure mode is mirrored -- the dark end vanishes -- so
# the palette runs UNREVERSED and the high values glow instead. Plasma over
# inferno/magma because its low end (deep blue-purple) is still separable from
# the background, so early curves recede without disappearing outright.
#
# DARK <- FALSE restores the light viridis experiment exactly; PALETTE <- BLUE
# under DARK <- FALSE restores the original blues.
DARK <- TRUE

BLUE <- c("#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
          "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b")

if (DARK) {
  # unreversed (bright = late/high), with the first 20% clipped off: plasma's
  # deepest blue-purples sat too close to the black surface, so the low end now
  # starts at a violet that reads as data rather than background
  PALETTE     <- viridisLite::plasma(10, begin = 0.2)
  INK_PRIMARY <- "#f2f1ec"
  INK_SECOND  <- "#c6c5bf"
  INK_MUTED   <- "#8f8e88"
  GRIDLINE    <- "#26262c"
  AXIS        <- "#44444c"
  SURFACE     <- "#101014"
} else {
  PALETTE     <- rev(viridisLite::viridis(10))  # reversed: dark = late/high
  INK_PRIMARY <- "#0b0b0b"
  INK_SECOND  <- "#52514e"
  INK_MUTED   <- "#898781"
  GRIDLINE    <- "#e1e0d9"
  AXIS        <- "#c3c2b7"
  SURFACE     <- "#fcfcfb"
}

# Axis and tick ink, deliberately the brightest thing on the panel after the
# data: pure white in the dark theme, where INK_MUTED left the scales reading
# as background furniture rather than as the frame you navigate the plot by.
# Not a literal "white" in the theme itself -- that would vanish on the light
# surface if DARK is turned off -- so it follows the toggle like every other
# colour here.
AXIS_INK <- if (DARK) "#ffffff" else INK_PRIMARY

# Display names, in the panels' order: the PRIMARY benchmarks first
# (PRIMARY_BENCHES, prepare_data.R), the rest alphabetical -- matching
# bench_levels(), which orders every script's loops the same way. Keys follow
# the benchmark ids in the analysis data (fm13 is fm_tiers1_3_v2's legacy
# key, prepare_data.R).
LABELS <- c(
  aime                     = "AIME (OTIS Mock)",
  chess                    = "Chess Puzzles",
  fm13                     = "FrontierMath, tiers 1-3",
  gpqa                     = "GPQA Diamond",
  mystery                  = "Mystery Game Puzzles",
  fm_2025_02_private       = "FrontierMath 2025-02",
  fm_tier4_2025_07_private = "FrontierMath, tier 4 (2025-07)",
  fm_tier4_v2              = "FrontierMath, tier 4 v2",
  math_lvl5                = "MATH, level 5",
  simpleqa                 = "SimpleQA Verified",
  swe_bench_verified       = "SWE-Bench Verified")

# The panel set extended by the pooled pseudo-benchmark, SIXTH -- right
# after the primaries it stacks. Deliberately a SEPARATE constant: adding
# pooled to LABELS itself would draw an empty sixth panel in every figure
# that never fits it (facet_wrap drop = FALSE). Scripts that draw a pooled
# panel pass this as frontier_plot()/iso_acc_plot()'s `labels`.
LABELS_POOLED <- append(LABELS, c(pooled = "Pooled primaries (ECI scale)"),
                        after = match("mystery", names(LABELS)))

# Benchmark facets are laid out two across, however many benchmarks the data
# carries; figure height follows the row count. 3 inches per row plus 1.5 for
# legend and caption reproduces the historical 10 x 7.5 quartet exactly at
# n = 4.
fig_height <- function(n_bench) 1.5 + 3 * ceiling(n_bench / 2)

EPOCH <- as.Date("2023-01-01")          # t = 0
as_t  <- function(date) as.numeric(as.Date(date) - EPOCH) / 365.25

dollar_log <- function(x) {
  ifelse(is.na(x), "", paste0("$", formatC(x, format = "fg", drop0trailing = TRUE)))
}

## ---- date grid ------------------------------------------------------------------

# Semiannual cadence anchored to the LAST observation and stepped backwards, not
# to calendar half-year starts. The final curve then lands exactly on the last run
# (so no observation sits outside every frontier) while consecutive curves stay
# exactly six months apart -- the spacing between curves is read as elapsed
# time, so it has to be uniform. Was quarterly; at that cadence successive
# Pareto staircases often coincided, so the figures showed fewer staircases
# than smooth curves and the counts looked mismatched.
grid_dates <- function(first, last) {
  n <- ceiling(as.numeric(as.Date(last) - as.Date(first)) / 365.25 * 2) + 1
  qs <- seq(as.Date(last), by = "-6 months", length.out = n)
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
# and it is the same fact envelope_constraints() uses to impose monotonicity.
corner_grid <- function(lnc_range, tc_range) {
  expand.grid(lncost = range(lnc_range), tc = range(tc_range))
}

# Worst case of both derivatives over the observed rectangle.
frontier_slope_bounds <- function(co, lncost, tc) {
  g <- corner_grid(lncost, tc)
  c(dcost = min(frontier_dcost(co, g$lncost, g$tc)),
    dtime = min(frontier_dtime(co, g$lncost, g$tc)))
}

## ---- Box-Cox primitives ------------------------------------------------------------

# The transform pair the Box-Cox specification is built on. They live HERE, not
# in boxcox_frontier.R, because the model-agnostic curve builders below must
# evaluate and invert BC fits, and boxcox_frontier.R already sources this file
# (via fit_specs.R) -- defining them there would close a source cycle.
# The origin for tau: October 1, 2020, when OpenAI started charging for
# GPT-3. Expressed through as_t so it converts exactly the way bc_tau()
# converts every other date (2020.7495 in decimal years).
BC_T0 <- 2023 + as_t(as.Date("2020-10-01"))

bc_tf <- function(y, l) if (abs(l) < 1e-6) log(y) else (y^l - 1) / l

# Inverse transform; NA where 1 + l*phi <= 0, i.e. where no positive argument
# attains phi -- the BC analogue of a negative discriminant.
bc_inv <- function(phi, l) {
  if (abs(l) < 1e-6) return(exp(phi))
  base <- 1 + l * phi
  ifelse(base > 0, base^(1 / l), NA_real_)
}

# THE RESPONSE IS NEVER TRANSFORMED. Every Box-Cox fit in this analysis is a
# Box-TIDWELL fit: phi() is applied to the REGRESSORS (cost, time, and in the
# cost direction odds), never to the dependent variable. The accuracy models
# therefore always use the plain logit link and the cost models always model
# ln cost, so nothing here needs a parametric link family.
#
# There used to be one: bc_mu(), bc_mu_eta() and bc_link() implemented
# phi(odds; lambda_odds) as a glm link, with ln_bc_inv() (cost_frontier.R)
# the cost direction's counterpart. They are deleted because a transformed
# response has to be INVERTED to be read back, and phi(y; l) inverts as
# ln y = log(1 + l*index)/l, which has a pole at index = -1/l. Any l != 0
# therefore puts a vertical asymptote at a finite index, and the profile
# walked the surface up to it: gpqa and aime reached fitted costs near e^20
# per task against a dearest observed run of e^0.6, and 1% of grid nodes
# carried up to 44% of the mean decline rate. l = 0 is the unique value whose
# inverse is the identity, hence the unique one with no pole -- so rather
# than keep a family whose interior members are unusable, the response
# lambda is gone and the identity is hard-coded. Regressor lambdas are never
# inverted, carry no pole, and stay profiled.

# The Bernoulli quasi-log-likelihood on the PROBABILITY scale, the objective
# the Pareto-grid fits share (fit_pareto_logit_env in envelope_frontier.R,
# the Box-Cox profiles in boxcox_frontier.R). Not -deviance/2: the two differ
# by the saturated log-likelihood, which is constant in the lambdas. mu is
# clamped away from 0/1 so a fitted probability at the boundary cannot emit
# -Inf.
bc_qll <- function(y, mu) {
  # which() + index assignment, not pmin/pmax: this sits in the Box-Cox
  # profiles' innermost loop, where the temporaries pmin/pmax allocate were a
  # measured fraction of the whole profile's time. Same values.
  lo <- which(mu < 1e-12)
  hi <- which(mu > 1 - 1e-12)
  mu[lo] <- 1e-12
  mu[hi] <- 1 - 1e-12
  sum(y * log(mu) + (1 - y) * log(1 - mu))
}

# A BC fit announces itself by the lambda attribute fit_bc() stamps on it;
# every fit the older specifications produce lacks it.
is_bc_fit <- function(fit) !is.null(attr(fit, "bc_lambda"))

# Named pieces of a BC fit: coefficients (beta_ prefix stripped, as in
# frontier_coefs) plus the profiled lambdas. lambda_odds is the doubly-
# transformed family's response-side parameter (fit_bc for the envelope and
# Pareto-grid fits); fits without it -- the run-level families keep the logit
# link -- read as 0, the logit.
bc_pieces <- function(fit) {
  cf <- coef(fit)
  names(cf) <- sub("^beta_", "", names(cf))
  lam <- attr(fit, "bc_lambda")
  list(b0 = unname(cf[["(Intercept)"]]), bx = unname(cf[["phic"]]),
       bt = unname(cf[["phit"]]), bxt = unname(cf[["phixt"]]),
       lc = unname(lam[["lambda_cost"]]), lt = unname(lam[["lambda_time"]]))
}

# tau, the BC time coordinate, from a Date: years since BC_T0. Uncentered --
# the BC specification does not use tbar (load_runs: year = 2023 + t).
bc_tau <- function(date) as_t(date) + 2023 - BC_T0

# Fitted frontier (u = 0) over a cost grid, for each date in `dates`.
frontier_curves <- function(fitset, data, dates_by_bench, tbar, n_cost = 200) {
  do.call(rbind, lapply(names(fitset), function(b) {
    fit <- fitset[[b]]
    sub <- data[data$benchmark == b, ]
    cost <- exp(seq(log(min(sub$cost)), log(max(sub$cost)), length.out = n_cost))
    g <- expand.grid(cost = cost, qdate = dates_by_bench[[b]])
    if (is_bc_fit(fit)) {
      p <- bc_pieces(fit)
      phic <- bc_tf(g$cost, p$lc)
      phit <- bc_tf(bc_tau(g$qdate), p$lt)
      # plogis: the index is on the logit scale, the response never being
      # transformed. Only the regressors carry a lambda.
      g$value <- plogis(p$b0 + p$bx * phic + p$bt * phit + p$bxt * phic * phit)
    } else {
      co <- frontier_coefs(fit)
      g$value <- plogis(frontier_index(co, log(g$cost), as_t(g$qdate) - tbar[[b]]))
    }
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
# `ranges` (benchmark, cost, value) pins each panel's axes via an invisible
# layer, so a panel with sparse data cannot have its axes read off whatever
# happens to be present.
# `title`/`subtitle` default to NULL and the static quartets pass neither: the
# viewer's controls and the filename already say which figure this is, and the
# vertical space goes to the panels instead.
# `labels` lets a caller extend the panel set beyond the data's benchmarks --
# plot_cost_frontier.R appends the pooled pseudo-benchmark in the sixth slot
# -- without touching the global LABELS, whose unused level would otherwise
# draw an empty panel in every other figure (drop = FALSE).
# `free_value` frees the value axis per panel: the pooled panel's value is an
# ECI capability score (~60-165), not a share, so a fixed 0-1 percent axis
# cannot hold both. Each real panel is still pinned to 0-1 through `ranges`,
# and the axis labeller formats each panel's own breaks -- percents at or
# below 1, plain ECI points above.
frontier_plot <- function(curves, pts, title = NULL, subtitle = NULL, ylab,
                          notes = character(0), step = FALSE,
                          colour_limits = NULL, ranges = NULL,
                          labels = LABELS, free_value = FALSE) {
  curves$benchmark <- factor(labels[curves$benchmark], levels = labels)
  pts$benchmark    <- factor(labels[pts$benchmark],    levels = labels)
  if (is.null(colour_limits)) colour_limits <- range(curves$year)
  blank_layer <- NULL
  if (!is.null(ranges)) {
    ranges$benchmark <- factor(labels[ranges$benchmark], levels = labels)
    blank_layer <- geom_blank(data = ranges, aes(cost, value), inherit.aes = FALSE)
  }

  geom_curve_layer <- if (step) {
    geom_step(aes(group = qdate, colour = year), direction = "hv", linewidth = 0.6)
  } else {
    geom_line(aes(group = qdate, colour = year), linewidth = 0.6)
  }

  base_notes <- c(
    paste("Curves are semiannual, stepped back from each benchmark's last run so",
          "spacing is exactly six months and the final curve lands on the last",
          "observation."),
    paste("Dots are observed runs, on the same time scale as the curves."))

  ggplot(curves, aes(cost, value)) +
    blank_layer +
    geom_point(data = pts, aes(cost, acc, colour = year), size = 0.35, alpha = 0.3,
               inherit.aes = FALSE) +
    geom_curve_layer +
    facet_wrap(~benchmark, scales = if (free_value) "free" else "free_x",
               ncol = 2, drop = FALSE) +
    scale_x_log10(breaks = 10^(-5:1), labels = dollar_log) +
    (if (free_value)
      scale_y_continuous(labels = function(v)
        ifelse(v > 1.5, sprintf("%.0f", v),
               scales::percent(v, accuracy = 1)))
     else scale_y_continuous(limits = c(0, 1),
                             labels = scales::percent_format(accuracy = 1))) +
    scale_colour_gradientn(
      colours = PALETTE, name = NULL, breaks = 2023:2026, limits = colour_limits,
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
      axis.text  = element_text(colour = AXIS_INK, size = 8),
      axis.title = element_text(colour = AXIS_INK, size = 9),
      strip.text = element_text(colour = INK_PRIMARY, face = "bold", size = 10,
                                hjust = 0),
      plot.title    = element_text(colour = INK_PRIMARY, face = "bold", size = 13),
      plot.subtitle = element_text(colour = INK_SECOND, size = 9.5),
      plot.caption  = element_text(colour = INK_MUTED, size = 7.5, hjust = 0),
      legend.position = "top", legend.justification = "left",
      legend.text = element_text(colour = INK_SECOND, size = 8),
      # legend.title matched to legend.text rather than left to
      # theme_minimal's default, which is near-black -- invisible against the
      # dark surface. Only record_timelines.R titles a legend today (every
      # other scale passes name = NULL), but the default was wrong for any
      # figure that ever does.
      legend.title = element_text(colour = INK_SECOND, size = 8),
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
pareto_step_layer <- function(steps, labels = LABELS) {
  steps$benchmark <- factor(labels[steps$benchmark], levels = labels)
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

## ---- the empirical staircase in iso-accuracy space ------------------------------

# The same Pareto staircase the frontier figures overlay, read the
# other way round: the minimum cost at which accuracy AT LEAST `level` had been
# achieved by each date,
#
#     C_a(t) = min { c_i : t_i <= t, acc_i >= a }
#
# which is exactly the level set of P_t(c): P_t(c) >= a iff c >= C_a(t). It is
# non-increasing in t by construction -- a released model stays available -- and
# it is a running MINIMUM, so like the staircase it carries old records forward.
#
#   sub    DATA FRAME of one benchmark's runs; needs cost, acc, releasedate
#   level  SCALAR in (0, 1), the accuracy target
#
# Returns a DATA FRAME with columns date and cost -- the jump points of the
# running minimum, extended flat to the benchmark's last run date (correct, not
# cosmetic: past the last release nothing changes, so the record really does
# hold out there) -- or NULL when no run ever reaches the level, which is a
# missing curve rather than a curve at infinite cost.
iso_pareto_steps <- function(sub, level) {
  s <- sub[sub$acc >= level, ]
  if (!nrow(s)) return(NULL)
  s <- s[order(s$releasedate), ]
  m <- cummin(s$cost)
  keep <- c(TRUE, diff(m) < 0)
  data.frame(date = c(s$releasedate[keep], max(sub$releasedate)),
             cost = c(m[keep], m[length(m)]))
}

# One staircase per benchmark per accuracy level, in the column layout
# iso_acc_plot()'s overlay expects; `levels` should be the same vector the
# fitted contours use, so each dashed curve pairs with a solid one shade for
# shade.
#
#   data     DATA FRAME of runs across benchmarks (needs benchmark, cost, acc,
#            releasedate)
#   levels   numeric VECTOR of accuracy targets in (0, 1)
#
# Returns a DATA FRAME with columns date, cost, acc, benchmark.
iso_pareto_curves <- function(data, levels) {
  do.call(rbind, lapply(sort(unique(data$benchmark)), function(b) {
    sub <- data[data$benchmark == b, ]
    do.call(rbind, lapply(levels, function(a) {
      st <- iso_pareto_steps(sub, a)
      if (is.null(st)) return(NULL)
      st$acc <- a
      st$benchmark <- b
      st
    }))
  }))
}

# Drawn over iso_acc_plot() exactly as pareto_step_layer() is drawn over
# frontier_plot(): same dash, same weight, and coloured by the same accuracy
# scale as the fitted contours, so each staircase reads against the contour at
# its own level.
iso_pareto_layer <- function(steps, labels = LABELS) {
  steps$benchmark <- factor(labels[steps$benchmark], levels = labels)
  geom_step(data = steps, aes(date, cost, group = interaction(benchmark, acc),
                              colour = acc),
            direction = "hv", linewidth = 0.4, linetype = "22",
            inherit.aes = FALSE)
}

ISO_PARETO_NOTE <- paste(
  "Dashed: the minimum cost at which accuracy at least each contour's level had",
  "been achieved by each date -- the Pareto staircase read as cost against date.")

# The pooled panel's iso-view layers: contours at ECI capability levels,
# SHADED like the accuracy contours -- each level mapped to its position
# within the pooled capability range `crng` and coloured on the figures'
# shared 0-1 ramp (the percent legend describes the accuracy panels; the
# pooled mapping is positional and each contour carries its ECI value as an
# in-place label). `steps`, when given, adds the pooled record staircases,
# shaded the same way.
pooled_iso_layers <- function(pool_iso, crng, labels = LABELS_POOLED,
                              steps = NULL) {
  shade <- function(x) {
    x$benchmark <- factor(labels[x$benchmark], levels = labels)
    x$lvl <- pmin(pmax((x$acc - crng[1]) / diff(crng), 0), 1)
    x
  }
  # a fit can blank every pooled contour point (the envelope BC has); the
  # label frame must then stay an empty data.frame, not collapse to NULL
  ok <- pool_iso[is.finite(pool_iso$cost), , drop = FALSE]
  lab <- if (nrow(ok)) do.call(rbind, lapply(split(ok, ok$acc), function(s)
    s[which.max(s$date), ])) else ok
  c(if (!is.null(steps)) list(
      geom_step(data = shade(steps),
                aes(date, cost, group = acc, colour = lvl),
                direction = "hv", linewidth = 0.4, linetype = "22",
                inherit.aes = FALSE)),
    list(
      geom_path(data = shade(pool_iso),
                aes(date, cost, group = interaction(acc, seg), colour = lvl),
                linewidth = 0.6, na.rm = TRUE, inherit.aes = FALSE),
      geom_text(data = shade(lab),
                aes(date, cost, label = sprintf("%g", acc)),
                colour = INK_SECOND, size = 2.6, hjust = 1, vjust = -0.6,
                inherit.aes = FALSE)))
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
#
# `cost_cap` (optional) tightens the second blank per LEVEL rather than per
# benchmark: a DATA FRAME with columns benchmark, acc and cost -- in practice
# the iso_pareto_curves() staircases -- whose per-(benchmark, acc) maximum cost
# becomes that contour's ceiling. A contour is then never drawn dearer than the
# empirical record for its own level, and a level with no rows at all (never
# achieved by any run) gets no contour rather than a purely extrapolated one.
# Matching is on the acc values themselves, so pass the SAME `levels` vector to
# both this function and iso_pareto_curves(). NULL (the default) keeps the
# benchmark-wide clip, which is what the S/A/B figures use.
iso_acc_curves <- function(fitset, data, tbar,
                            levels = seq(0.10, 0.90, by = 0.20), n_date = 300,
                            min_slope = 0.05, cost_cap = NULL) {
  do.call(rbind, lapply(names(fitset), function(b) {
    fit <- fitset[[b]]
    sub <- data[data$benchmark == b, ]
    urng <- range(sub$lncost)
    dts <- seq(min(sub$releasedate), max(sub$releasedate), length.out = n_date)
    g   <- expand.grid(date = dts, acc = levels)

    if (is_bc_fit(fit)) {
      # Still closed form: holding z fixed, phi_c = (z - b0 - bt*phit) /
      # (bx + bxt*phit), and the transform inverts analytically. phi is
      # monotone, so there is ONE root -- no falling branch, no fold. The
      # inversion blanks in two honest ways: a non-positive slope in phi_c
      # (free disposal failing at that date), and 1 + lambda*phi_c <= 0, where
      # NO positive cost attains the target (bc_inv returns NA) -- the BC
      # analogue of a negative discriminant. Near-zero positive slopes need no
      # guard of their own: they send phi_c to +/-Inf, which lands outside the
      # observed cost range and is removed by the clip below.
      p <- bc_pieces(fit)
      phit <- bc_tf(bc_tau(g$date), p$lt)
      bb <- p$bx + p$bxt * phit
      # the index a target accuracy must reach: its logit, the response being
      # untransformed
      phic <- (qlogis(g$acc) - p$b0 - p$bt * phit) / bb
      phic[bb <= 0] <- NA_real_
      roots <- list(rising = log(bc_inv(phic, p$lc)))
      disc <- rep(NA_real_, nrow(g))
    } else {
      co <- frontier_coefs(fit)
      tc <- as_t(g$date) - tbar[[b]]
      aa <- co[["bxx"]]
      # NOT the cost slope: this is the coefficient on u in the quadratic below.
      # The two coincide only when bxx = 0, where the quadratic degenerates to
      # the division and the u-coefficient IS the derivative.
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
        slope <- sqrt(pmax(disc, 0))               # |dz/du| at either root
        roots <- list(rising  = (-bb + slope) / (2 * aa),
                      falling = (-bb - slope) / (2 * aa))
        # disc < 0 means the target is not attained at ANY cost on that date
        roots <- lapply(roots, function(u) { u[disc < 0] <- NA_real_; u })
        # No min_slope guard on this branch: nothing blows up here (the roots
        # stay bounded, and a near-linear surface sends the spurious root
        # outside the cost range, where the clip below removes it), and gapping
        # at sqrt(disc) ~ 0 would punch a hole exactly at the turning point
        # where the two branches meet -- the one place the contour is genuinely
        # continuous.
      }
    }
    # Per-row admissibility, beyond the benchmark's observed cost range.
    # Where cost_cap supplies a level's record staircase, a contour point is
    # blanked only where it is BOTH earlier than that level's first record
    # AND dearer than its dearest record -- so each contour extends back to
    # its record's start or up to its cost ceiling, whichever is the more
    # generous, rather than being cut by the ceiling alone. A level cost_cap
    # covers the benchmark for but omits was never achieved by any run:
    # no contour at all.
    cap_u <- rep(Inf, nrow(g))
    birth <- rep(-Inf, nrow(g))
    if (!is.null(cost_cap)) {
      cb <- cost_cap[cost_cap$benchmark == b, ]
      lv <- unique(cb$acc)
      # match on the shared `levels` values, exact because both data frames
      # descend from the same numeric vector -- no character round trip, which
      # would corrupt values like seq()'s 0.45000000000000007
      i <- match(g$acc, lv)
      cap_u <- log(vapply(lv, function(a) max(cb$cost[cb$acc == a]),
                          numeric(1)))[i]
      birth <- vapply(lv, function(a) as.numeric(min(cb$date[cb$acc == a])),
                      numeric(1))[i]
      cap_u[is.na(cap_u)] <- -Inf   # never achieved
      birth[is.na(birth)] <- Inf
    }

    both <- do.call(rbind, lapply(names(roots), function(br) {
      u <- roots[[br]]
      # the observed cost range, then the union rule above, are the only
      # clips either branch gets
      u[is.finite(u) &
          (u < urng[1] | u > urng[2] |
             !(u <= cap_u | as.numeric(g$date) >= birth))] <- NA_real_
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
# `labels` as in frontier_plot(): the caller may extend the panel set (the
# pooled pseudo-benchmark) without touching the global LABELS.
iso_acc_plot <- function(curves, pts, title = NULL, subtitle = NULL,
                          notes = character(0), ranges = NULL,
                          labels = LABELS) {
  curves$benchmark <- factor(labels[curves$benchmark], levels = labels)
  pts$benchmark    <- factor(labels[pts$benchmark],    levels = labels)
  blank_layer <- NULL
  if (!is.null(ranges)) {
    ranges$benchmark <- factor(labels[ranges$benchmark], levels = labels)
    blank_layer <- geom_blank(data = ranges, aes(date, cost), inherit.aes = FALSE)
  }

  # The cost floor: the decade at or below the cheapest run anywhere in `pts`
  # ($0.00001 with the current data). `pts` is the full sample in every caller,
  # so every iso figure computes the SAME floor -- which, with the shared facet
  # scales, is what puts all iso panels on one identical cost axis.
  ymin <- 10^floor(log10(min(pts$cost, na.rm = TRUE)))

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
    # Shared scales, not free_y: cost is in DOLLARS on every panel, so one axis
    # serves all four, the right-hand column drops its y labels exactly as the
    # top row already drops its x labels, and the horizontal space goes to the
    # panels. The price is that each panel spans the union of the cost ranges
    # rather than its own -- on a log scale spanning five decades anyway, cheap.
    facet_wrap(~benchmark, ncol = 2, drop = FALSE) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    # The panel stops exactly at the floor: zero expansion below, so the axis
    # corner sits on a labelled round break rather than leaving dead margin
    # under lines that hug the cheapest costs. A whisker of expansion stays on
    # top so the dearest dots are not clipped in half.
    # mult expansion applies to the WHOLE 5.5-decade range, so even 3% on top
    # put the ceiling 46% above the dearest run; additive 0.04 decades (~10%)
    # is enough that the top dots are not clipped in half.
    scale_y_log10(breaks = 10^(-5:1), labels = dollar_log,
                  limits = c(ymin, NA),
                  expand = expansion(add = c(0, 0.04))) +
    scale_colour_gradientn(
      colours = PALETTE, name = NULL, limits = c(0, 1),
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
# time variables.
#
# READ ONCE PER PROCESS. Every script called this independently and
# fit_store.R kept a snapshot of its own -- eleven reads minutes apart in a
# full run -- so a data file that changed mid-run (OneDrive finishing a sync,
# in the case that prompted this) yielded outputs mixing two datasets: the
# nonparametric columns computed from one, the fitted columns from the other,
# with nothing to reveal the split. One read makes that impossible, and lets
# fit_store.R drop its own copy. Guarded like the other process-level state
# (.fit_store, .fit_cluster_env): library files are re-sourced many times per
# process and an unguarded assignment would empty the cache on each pass.
if (!exists(".runs_cache", inherits = FALSE))
  .runs_cache <- new.env(parent = emptyenv())

load_runs <- function() {
  if (!is.null(.runs_cache$d)) return(.runs_cache$d)
  d <- build_runs()
  d$t <- as_t(d$releasedate)
  d <- d[stats::complete.cases(d[c("acc", "lncost", "t", "model", "effort")]), ]
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
  .runs_cache$d <- d
  d
}

# Reference date per benchmark, recovered from the data rather than recomputed:
# t - tc is exactly the demeaning constant, whatever sample tc was built on. A
# separately-computed mean(t) would silently disagree if rows were dropped after
# tc was assigned, shifting the plotted curves off the model that produced them.
bench_tbar <- function(d) {
  tapply(d$t - d$tc, d$benchmark, function(x) x[1])
}

# Observed history per benchmark, in years: first run to last -- the spans,
# not "years since start", so a benchmark that stopped being tested does not
# keep gaining weight. One input to the pooling weights (pool_weights,
# envelope_frontier.R).
bench_spans <- function(d) {
  tapply(d$t, d$benchmark, function(x) diff(range(x)))
}

# ONE global semiannual grid, anchored to the latest run across all benchmarks;
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
  benches <- bench_levels(d$benchmark)
  grid <- grid_dates(min(d$releasedate), max(d$releasedate))
  setNames(lapply(benches, function(b) {
    grid[grid >= min(d$releasedate[d$benchmark == b])]
  }), benches)
}
