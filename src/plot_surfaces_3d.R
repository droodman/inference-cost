# Interactive 3-D renderings of the Pareto-grid fits, for the plot viewer's
# 2-D / 3-D toggle: each page shows, per benchmark, the EMPIRICAL surface
# (the staircase P_t(c) or the record cost ln C_a(t), NA holes where
# undefined) with the FITTED surface overlaid at half opacity -- rotate the
# scene and the shelf-and-cliff geometry, and how each specification does or
# does not follow it, is directly visible.
#
# The frontier-per-se models only -- the Pareto-grid fits, unconstrained and
# envelope-constrained, in both directions -- because those are the models
# whose empirical reference IS the staircase/record surface drawn beneath
# them (the run-level models' honest reference would be the run cloud, a
# different figure). The constrained fits gain the most from rotation: where
# the wireframe touches the data surface is the visual form of their slack
# diagnostics. Each model is rendered in three specifications and BOTH
# orientations:
#
#   frontier     x = ln cost, y = year, z = accuracy      (performance surface)
#   isoaccuracy  x = accuracy, y = year, z = ln cost      (cost surface)
#
# Accuracy axes are PLAIN accuracy, not logit: the logit stretch made 0.5%
# vs 1% look as different as 50% vs 73%, an emphasis the accuracy-direction
# fits (probability-scale quasi-likelihoods) do not share. The level lattice
# is therefore uniform in accuracy; internally the fits are still evaluated
# at its logit coordinates, so nothing about any FIT changes -- the display
# does. (The cost-direction objectives' own la-uniform node weighting, a
# documented tail-heavy hazard, is unaffected; the display just no longer
# mirrors it.)
#
# so 4 model keys x 3 specifications x 2 views = 24 HTML files,
# surface3d_<key>_<spec>_<view>.html. Each direction is NATIVE in one view and
# INVERTED into the other; the inversions take the rising branch and leave NA
# where the surface bends back or the transform has no preimage -- holes, not
# folds, exactly as the 2-D figures blank.
#
# A third view, surface3d_<key>_<spec>_decline.html (12 more files), plots the
# INSTANTANEOUS RATE OF COST DECLINE at fixed accuracy -- d ln C / dt read off
# each fitted surface, as a quarterly percentage -- over (accuracy, date); see
# the decline section at the end of this file.
#
# Every page also gets a heatmap twin, heatmap_<key>_<spec>_<view>.png: the
# fitted surface seen from directly above, as a static faceted figure (see
# the heatmap section below).
#
# The widgets are static HTML + client-side WebGL (htmlwidgets/plotly): no
# server, so they serve from the repo or GitHub Pages as-is. Written with
# selfcontained = FALSE and a SHARED output/lib directory, so plotly.js is
# committed once (~4 MB) and each page is small.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("fit_store.R")   # pulls the whole fitting stack
suppressMessages(library(plotly))

N_X <- 100   # lattice density along the non-time axis, matching the fits' grids
N_T <- 100

d <- load_runs()
benches <- bench_levels(d$benchmark)
tbar <- bench_tbar(d)

## ---- fits: both frontier-per-se pairs, three specifications each ------------------

# which keys are accuracy-direction; the cost keys below are their duals
ACC_KEYS <- c("paretologit", "paretologitenv")

# All from the shared store (fit_store.R): under run_all.R these are the same
# objects the 2-D figure scripts and the tables already fitted.
fits <- list()
for (key in c(ACC_KEYS, "costgridols", "costgridolsenv")) {
  grid <- if (key %in% ACC_KEYS) store_grid(key) else store_cost(key)
  fits[[key]] <- c(grid, list(
    bc = if (key %in% ACC_KEYS) store_bc(key) else store_cost_bc(key)))
}

## ---- per-benchmark lattices and empirical surfaces --------------------------------
#
# One bundle per benchmark: the two lattices (matching the fits' own grids)
# and the two empirical matrices with NA where undefined. z matrices are
# [y, x] as plotly expects (rows indexed by y).

bundle <- function(b) {
  s  <- d[d$benchmark == b, ]
  su <- iso_runs(s)
  off <- (s$year - s$tc)[1] - BC_T0

  lnc <- seq(min(s$lncost), max(s$lncost), length.out = N_X)
  # uniform in ACCURACY (so heatmap rasters stay regular and the axis reads
  # plainly), carried as logit values because that is the scale the fits and
  # the record are evaluated on
  la  <- qlogis(seq(plogis(min(su$la)), plogis(max(su$la)),
                    length.out = N_X))
  tc  <- seq(min(s$tc), max(s$tc), length.out = N_T)

  # P_t(c) over (lnc, tc), raw accuracy, NA where no run is cheaper & earlier
  so <- s[order(s$lncost), ]
  P <- matrix(NA_real_, N_T, N_X)
  for (k in seq_len(N_T)) {
    el <- so[so$tc <= tc[k], ]
    if (!nrow(el)) next
    m <- cummax(el$acc)
    idx <- findInterval(lnc, el$lncost)
    P[k, ] <- ifelse(idx >= 1, m[pmax(idx, 1)], NA_real_)
  }

  # ln C_a(t) over (la, tc), positive-accuracy runs, NA where unachieved
  suo <- su[order(su$lncost), ]
  R <- matrix(NA_real_, N_T, N_X)
  for (k in seq_len(N_T)) {
    el <- suo[suo$tc <= tc[k], ]
    if (!nrow(el)) next
    m <- cummax(el$la)
    R[k, ] <- vapply(la, function(a) {
      j <- which(m >= a)[1]
      if (is.na(j)) NA_real_ else el$lncost[j]
    }, numeric(1))
  }

  list(s = s, off = off, lnc = lnc, la = la, tc = tc,
       year = tc + tbar[[b]] + 2023, P = P, R = R,
       zrng_cost = range(s$lncost) + c(-0.5, 0.5))
}
bundles <- setNames(lapply(benches, bundle), benches)

## ---- fitted surfaces on each lattice ----------------------------------------------
#
# Native evaluations are direct; the cross evaluations invert the fitted
# surface, taking the RISING branch of a quadratic (NA past the fold) and the
# single root of the monotone BC transforms (NA where the transform has no
# preimage). Everything returns a [N_T, N_X] matrix.

grid_lt <- function(xv, tcv) list(x = rep(xv, times = length(tcv)),
                                  t = rep(tcv, each = length(xv)))
as_z <- function(v) t(matrix(v, nrow = N_X))

# accuracy-direction fit as accuracy over (lnc, tc) -- native
zacc_native <- function(fit, bl) {
  g <- grid_lt(bl$lnc, bl$tc)
  if (is_bc_fit(fit)) {
    p <- bc_pieces(fit)
    phic <- bc_tf(exp(g$x), p$lc)
    phit <- bc_tf(g$t + bl$off, p$lt)
    as_z(plogis(p$b0 + p$bx * phic + p$bt * phit + p$bxt * phic * phit))
  } else {
    as_z(plogis(frontier_index(frontier_coefs(fit), g$x, g$t)))
  }
}

# accuracy-direction fit inverted: ln cost over (la, tc)
zacc_inverted <- function(fit, bl) {
  g <- grid_lt(bl$la, bl$tc)
  if (is_bc_fit(fit)) {
    p <- bc_pieces(fit)
    phit <- bc_tf(g$t + bl$off, p$lt)
    bb <- p$bx + p$bxt * phit
    # g$x IS the logit, and the response is untransformed, so the index a
    # target accuracy must reach is g$x itself
    phic <- (g$x - p$b0 - p$bt * phit) / bb
    phic[bb <= 0] <- NA_real_
    as_z(log(bc_inv(phic, p$lc)))
  } else {
    co <- frontier_coefs(fit)
    aa <- co[["bxx"]]
    bb <- co[["bx"]] + co[["bxt"]] * g$t
    cc <- co[["b0"]] + co[["bt"]] * g$t + co[["btt"]] * g$t^2 - g$x
    u <- if (abs(aa) < 1e-10) {
      out <- -cc / bb; out[bb < 0.05] <- NA_real_; out
    } else {
      disc <- bb^2 - 4 * aa * cc
      out <- (-bb + sqrt(pmax(disc, 0))) / (2 * aa)   # rising branch
      out[disc < 0] <- NA_real_; out
    }
    as_z(u)
  }
}

# cost-direction fit as ln cost over (la, tc) -- native
zcost_native <- function(fit, bl) {
  srf <- cost_surface(fit, bl$s)
  g <- grid_lt(bl$la, bl$tc)
  as_z(srf$f(g$x, g$t))
}

# cost-direction fit inverted: accuracy over (lnc, tc)
zcost_inverted <- function(fit, bl) {
  g <- grid_lt(bl$lnc, bl$tc)
  if (is_cost_bc(fit)) {
    cf <- coef(fit); names(cf) <- sub("^beta_", "", names(cf))
    lam <- attr(fit, "bc_lambda")
    phit <- bc_tf(g$t + bl$off, lam[["lambda_time"]])
    bb <- cf[["phia"]] + cf[["phiat"]] * phit
    # g$x IS ln cost and the response is untransformed, so the index a
    # target cost corresponds to is g$x itself
    phia <- (g$x - cf[["(Intercept)"]] - cf[["phit"]] * phit) / bb
    phia[bb <= 0] <- NA_real_
    odds <- bc_inv(phia, lam[["lambda_odds"]])
    as_z(odds / (1 + odds))
  } else {
    co <- cost_coefs(fit)
    aa <- co[["gaa"]]
    bb <- co[["ga"]] + co[["gat"]] * g$t
    cc <- co[["g0"]] + co[["gt"]] * g$t + co[["gtt"]] * g$t^2 - g$x
    u <- if (abs(aa) < 1e-10) {
      out <- -cc / bb; out[bb < 0.05] <- NA_real_; out
    } else {
      disc <- bb^2 - 4 * aa * cc
      out <- (-bb + sqrt(pmax(disc, 0))) / (2 * aa)   # rising branch
      out[disc < 0] <- NA_real_; out
    }
    as_z(plogis(u))
  }
}

## ---- assembly ----------------------------------------------------------------------

SURFACE  <- "#101014"; INK <- "#f2f1ec"; MUTED <- "#8f8e88"; GRID <- "#26262c"

# The 2-D figures' own ramp (PALETTE, frontier_viz.R -- plasma with the
# darkest 20% clipped off under the dark theme), rebuilt as a plotly
# colorscale so the 3-D surfaces cannot drift from the 2-D figures, clip
# included. viridisLite returns 8-digit hex (trailing alpha), which plotly's
# WebGL parser does not accept -- strip to 6.
COLORSCALE <- Map(function(p, col) list(p, substr(col, 1, 7)),
                  seq(0, 1, length.out = length(PALETTE)), PALETTE)

ax <- function(title) list(title = title, backgroundcolor = SURFACE,
                           gridcolor = GRID, color = MUTED, showbackground = TRUE)

# A cost axis: internally ln cost, labelled in DOLLARS at decadal intervals,
# matching the 2-D figures' scale_*_log10(breaks = 10^(-5:1), dollar_log).
ax_cost <- function(rng) {
  k <- -5:1
  keep <- k * log(10) >= rng[1] & k * log(10) <= rng[2]
  c(ax("cost per task"), list(tickvals = k[keep] * log(10),
                              ticktext = dollar_log(10^k[keep])))
}

# A plain accuracy axis, labelled in percent.
ax_prob <- function(rng) {
  a <- seq(0, 1, 0.25)
  keep <- a >= rng[1] - 0.02 & a <= rng[2] + 0.02
  c(ax("accuracy"), list(tickvals = a[keep],
                         ticktext = sprintf("%d%%", 100 * a[keep])))
}

# Initial viewpoint: the box's near lower corner is the one with the LOWEST
# date and LOWEST cost -- the eye sits in the (min x, min y) quadrant, so the
# frontier view leads with cheap-and-early and the iso view with the low end
# of the level axis at the earliest date (cost is vertical there, so its
# minimum is the bottom face by construction). The eye's DISTANCE sets how
# large the cube renders inside its quadrant (closer = larger); this norm
# (~1.8, vs plotly's default ~2.2) is the zoom that fills the cell without
# clipping axis labels.
CAMERA <- list(eye = list(x = -1.2, y = -1.2, z = 0.7))

# One page: four scenes (2 x 2 benchmark quadrants), each with the empirical
# surface (filled, plasma ramp) and the fitted one as a SEE-THROUGH WIREFRAME:
# hidesurface drops the fill entirely and the x/y contour lines draw the mesh,
# so the fitted shape is legible without occluding the data surface behind it.
# Ink-colored, not black -- black segments would vanish on the dark backdrop.
# The "decline" view has no empirical counterpart (pass zs_emp = NULL): the
# fitted rate surface is drawn filled, over (accuracy, year).
page <- function(zs_emp, zs_fit, xs, view) {
  # Full-bleed cells, two across and however many rows the benchmark count
  # needs (mirroring the 2-D figures' N x 2 facets): scenes carry generous
  # internal padding of their own, so explicit gutters between the domains
  # only compound the dead space; the benchmark titles sit on the seams. The
  # page's pixel height scales with the row count so each scene keeps a
  # usable size, and the viewer's iframe scrolls.
  n_rows <- ceiling(length(benches) / 2)
  doms <- lapply(seq_along(benches), function(i) {
    r <- ceiling(i / 2)
    c0 <- (i - 1) %% 2
    list(x = c(0.5 * c0, 0.5 * c0 + 0.5),
         y = c((n_rows - r) / n_rows, (n_rows - r + 1) / n_rows))
  })
  p <- plot_ly(height = 420 * n_rows + 30)
  lay <- list(paper_bgcolor = SURFACE, font = list(color = INK),
              showlegend = FALSE, margin = list(l = 0, r = 0, t = 22, b = 0))
  ann <- list()
  for (i in seq_along(benches)) {
    b <- benches[i]
    bl <- bundles[[b]]
    sc <- if (i == 1) "scene" else paste0("scene", i)
    if (is.null(zs_emp)) {
      p <- add_surface(p, x = xs[[b]], y = bl$year, z = zs_fit[[b]],
                       scene = sc, colorscale = COLORSCALE, showscale = FALSE,
                       name = "fitted")
    } else {
      p <- add_surface(p, x = xs[[b]], y = bl$year, z = zs_emp[[b]],
                       scene = sc, colorscale = COLORSCALE, showscale = FALSE,
                       name = "empirical")
      xr <- range(xs[[b]]); yr <- range(bl$year)
      mesh <- function(rng, n = 14) list(show = TRUE, color = INK,
                                         width = 2, start = rng[1],
                                         end = rng[2], size = diff(rng) / n)
      p <- add_surface(p, x = xs[[b]], y = bl$year, z = zs_fit[[b]],
                       scene = sc, hidesurface = TRUE, showscale = FALSE,
                       contours = list(x = mesh(xr), y = mesh(yr)),
                       name = "fitted")
    }
    zx <- if (view == "frontier")
      list(xaxis = ax_cost(range(xs[[b]])),
           zaxis = c(ax("accuracy"), list(range = c(0, 1))))
    else if (view == "isoaccuracy")
      list(xaxis = ax_prob(range(xs[[b]])),
           zaxis = c(ax_cost(bl$zrng_cost), list(range = bl$zrng_cost)))
    else
      list(xaxis = ax_prob(range(xs[[b]])),
           zaxis = ax("cost drop, %/qtr"))
    lay[[sc]] <- c(list(domain = doms[[i]], yaxis = ax("year"),
                        aspectmode = "cube", camera = CAMERA), zx)
    # In each quadrant's UPPER-LEFT corner, hanging below the top edge --
    # centred titles straddled the 0.5 seam onto the upper scenes, and even
    # anchored ones read as captions of the plot above; the corner is where
    # the neighbouring scene's centred cube leaves the most clearance, and it
    # matches the 2-D facet labels (left-aligned, frontier_theme).
    ann[[i]] <- list(text = LABELS[[b]], x = doms[[i]]$x[1] + 0.01,
                     y = doms[[i]]$y[2] - 0.01, xanchor = "left",
                     yanchor = "top", xref = "paper", yref = "paper",
                     showarrow = FALSE, font = list(color = INK, size = 14))
  }
  lay$annotations <- ann
  do.call(layout, c(list(p), lay))
}

# plotly keys its internal data references by tempfile-derived RANDOM ids
# (visdat / cur_data / attrs), and they land verbatim in the saved JSON -- the
# second source of meaningless rebuild-to-rebuild churn after saveWidget's
# random element id. They are only cross-references among themselves, so
# renaming them consistently to a stable stem makes identical pages serialize
# identically.
stable_plotly_ids <- function(w, stem) {
  ids <- unique(c(names(w$x$visdat), names(w$x$attrs), w$x$cur_data))
  map <- setNames(paste0(stem, "-", seq_along(ids)), ids)
  names(w$x$visdat) <- unname(map[names(w$x$visdat)])
  names(w$x$attrs)  <- unname(map[names(w$x$attrs)])
  w$x$cur_data      <- unname(map[w$x$cur_data])
  w
}

## ---- heatmap companions ---------------------------------------------------------------
#
# Every 3-D page gets a HEATMAP twin, heatmap_<key>_<spec>_<view>.png: the
# FITTED surface seen from directly above, faceted per benchmark like the
# other 2-D figures. The fill is the 3-D page's z on the same plasma ramp;
# for the overlay views that means the fitted surface alone (the 3-D page's
# wireframe, filled in), the empirical staircase being already available in
# the 2-D figure sets.
#
# The heatmaps are NEVER masked by the data: the fitted surface fills its
# whole domain, and the empirical frontier is DRAWN over it instead, as a
# PAIR of 50%-black step lines: the running extremes of the observed values
# of the plotted variable -- best and worst (positive) accuracy so far in the
# (accuracy, date) views, cheapest and dearest run so far in the (cost, date)
# view. Between the lines, models with those values actually existed by that
# date; outside them, the surface is the functional form extrapolating,
# visible but labeled as such by the band, rather than hidden by a mask.
# (The decline 3-D pages keep their mask: they draw no lines to carry the
# distinction.) The only remaining holes are where the fit itself has no
# value -- a quadratic's fold, a Box-Cox index off its transform's range.

# Each panel's fill is normalized to ITS OWN anchor range -- the range of the
# same benchmark's colored surface on the 3-D page (the empirical staircase P
# for the frontier view, the record R for the iso view, the masked decline
# surface for the decline view, via zs_ref). plotly normalizes every scene's
# colors to that scene's own values, so this is what makes a heatmap panel
# and its 3-D scene agree color-for-color; a single shared scale let one
# extreme benchmark compress everyone else into a corner of the ramp. The
# bracket in each strip label carries the panel's mapping (dark -> bright),
# and values outside it saturate at the endpoints.
heat_anchor <- function(view, b, zs_ref = NULL) {
  bl <- bundles[[b]]
  r <- switch(view,
              frontier    = range(bl$P, na.rm = TRUE),
              isoaccuracy = range(bl$R, na.rm = TRUE),
              range(zs_ref[[b]], na.rm = TRUE))
  if (r[1] == r[2]) r + c(-0.5, 0.5) else r
}

fmt_anchor <- function(view, r) {
  switch(view,
         frontier    = sprintf("[%.0f%%, %.0f%%]", 100 * r[1], 100 * r[2]),
         isoaccuracy = sprintf("[%s, %s]", dollar_log(exp(r[1])),
                               dollar_log(exp(r[2]))),
         sprintf("[%.0f, %.0f %%/qtr]", r[1], r[2]))
}

heat_labels <- function(view, zs_ref = NULL) {
  setNames(vapply(benches, function(b)
    sprintf("%s  %s", LABELS[[b]],
            fmt_anchor(view, heat_anchor(view, b, zs_ref))), ""), benches)
}

# The long data frame geom_raster wants, from a benchmark-keyed list of
# [N_T, N_X] matrices and their x lattices, fills normalized per panel.
# Column-major as.vector makes the date index vary fastest, matching
# rep(x, each = N_T).
heat_df <- function(zs, xs, view, labs_b, zs_ref = NULL) {
  do.call(rbind, lapply(benches, function(b) {
    bl <- bundles[[b]]
    a <- heat_anchor(view, b, zs_ref)
    zn <- (zs[[b]] - a[1]) / (a[2] - a[1])
    zn[zn < 0] <- 0
    zn[zn > 1] <- 1
    data.frame(benchmark = factor(labs_b[[b]], levels = unname(labs_b)),
               x = rep(xs[[b]], each = N_T),
               year = rep(bl$year, times = length(xs[[b]])),
               z = as.vector(zn))
  }))
}

# The empirical frontiers as step paths over (year, x), built from the runs
# themselves rather than the lattice: a running extreme enters at the run
# that sets it and holds until the next record, ending at the panel's last
# date. Both extremes are returned -- the observed RANGE of the plotted
# variable at each date. For the frontier view that is the cheapest and
# dearest run so far; for the other two views the best and worst CLIPPED
# logit accuracy so far (zeros excluded by iso_runs, as everywhere on that
# scale). `side` keeps geom_path from joining the two staircases; `labs_b`
# is the per-panel strip-label vector, so the paths land in the right facets.
frontier_steps <- function(view, labs_b) {
  one_run <- function(v, yr, ymax, extreme) {
    run <- extreme(v)
    new <- !duplicated(run)
    xk <- run[new]
    yk <- yr[new]
    data.frame(x = rep(xk, each = 2),
               year = as.vector(rbind(yk, c(yk[-1], ymax))))
  }
  do.call(rbind, lapply(benches, function(b) {
    bl <- bundles[[b]]
    if (view == "frontier") {
      s <- bl$s[order(bl$s$tc), ]
      v <- s$lncost
    } else {
      s <- iso_runs(bl$s)
      s <- s[order(s$tc), ]
      v <- plogis(s$la)   # the accuracy axes plot plain accuracy
    }
    ymax <- max(bl$year)
    out <- rbind(cbind(one_run(v, s$year, ymax, cummin), side = "lo"),
                 cbind(one_run(v, s$year, ymax, cummax), side = "hi"))
    out$benchmark <- factor(labs_b[[b]], levels = unname(labs_b))
    out
  }))
}

# Axis pieces shared with the 3-D pages' labeling tricks: log cost labeled
# in dollars.
COST_BRK <- log(10^(-5:1))
COST_LAB <- dollar_log(10^(-5:1))

# `fill_limits`, when given, pins the color scale (values beyond it squish
# into the endpoint colors) -- the decline maps calibrate it to the region
# inside the frontier, so extrapolated rates cannot stretch the ramp.
# Time on the HORIZONTAL axis, the modeled variable on the vertical -- the
# reading convention of every other 2-D figure in the repo. No colorbar: the
# fill is per-panel normalized (see heat_anchor), so the mapping lives in
# each strip label's bracket and a caption states the rule once.
# The caption, spelling out what the two black staircases ARE. They are the
# one element of these figures that is not model output, and without saying so
# a reader has no way to tell them from a fitted contour -- nor to know that
# the surface outside them is extrapolating rather than describing. Named for
# the vertical variable, which differs by view: cost in the frontier view,
# accuracy in the other two. Line breaks are manual, sized for the 10-inch
# canvas at the theme's 7.5pt caption.
heat_caption <- function(view) {
  frontier <- view == "frontier"
  extremes <- if (frontier) "cheapest and dearest run"
              else "lowest- and highest-scoring run"
  vert     <- if (frontier) "cost" else "accuracy"
  paste0(
    "Fill runs dark to bright over each panel's own bracketed range -- the ",
    "same per-panel normalization the 3-D pages use, so equal color means ",
    "equal value between a panel and its 3-D scene.\n",
    "The two 50%-black staircases are NOT fitted: they trace the ", extremes,
    " observed up to each date, each a running extreme that steps out at a ",
    "new record and holds until the next.\n",
    "The band between them is the range the data actually cover in ", vert,
    "; outside it the surface is extrapolating, and values beyond the ",
    "bracket saturate at its endpoints.")
}

heat_plot <- function(zs, xs, view, zs_ref = NULL) {
  labs_b <- heat_labels(view, zs_ref)
  hd <- heat_df(zs, xs, view, labs_b, zs_ref)
  p <- ggplot(hd, aes(year, x, fill = z)) +
    geom_raster(na.rm = TRUE) +
    geom_path(data = frontier_steps(view, labs_b), aes(year, x, group = side),
              inherit.aes = FALSE, colour = "black", alpha = 0.5,
              linewidth = 0.7) +
    facet_wrap(~ benchmark, ncol = 2, scales = "free") +
    scale_fill_gradientn(colours = PALETTE, limits = c(0, 1),
                         na.value = SURFACE, guide = "none") +
    labs(x = NULL, caption = heat_caption(view)) +
    frontier_theme() +
    theme(panel.grid.major = element_blank())
  if (view == "frontier") {
    p + scale_y_continuous(name = "Cost per task (log scale)",
                           breaks = COST_BRK, labels = COST_LAB)
  } else {
    p + scale_y_continuous(name = "Accuracy",
                           breaks = seq(0, 1, 0.25),
                           labels = sprintf("%d%%", seq(0, 100, 25)))
  }
}

save_heatmap <- function(zs, xs, view, key, spec, zs_ref = NULL) {
  fh <- sprintf("heatmap_%s_%s_%s.png", key, spec, view)
  ggsave(out_path(fh), heat_plot(zs, xs, view, zs_ref), width = 10,
         height = fig_height(length(benches)), dpi = 200,
         device = ragg::agg_png)
  cat("wrote", fh, "\n")
}

VIEWS <- c("frontier", "isoaccuracy")
for (spec in c("lin", "quad", "bc")) {
  for (key in names(fits)) {
    fset <- fits[[key]][[spec]]
    for (view in VIEWS) {
      xs <- lapply(bundles, function(bl)
        if (view == "frontier") bl$lnc else plogis(bl$la))
      zs_emp <- lapply(bundles, function(bl) if (view == "frontier") bl$P else bl$R)
      zs_fit <- setNames(lapply(benches, function(b) {
        bl <- bundles[[b]]
        if (key %in% ACC_KEYS) {
          if (view == "frontier") zacc_native(fset[[b]], bl)
          else zacc_inverted(fset[[b]], bl)
        } else {
          if (view == "frontier") zcost_inverted(fset[[b]], bl)
          else zcost_native(fset[[b]], bl)
        }
      }), benches)
      w <- page(zs_emp, zs_fit, xs, view)
      f <- sprintf("surface3d_%s_%s_%s.html", key, spec, view)
      # a DETERMINISTIC element id: saveWidget otherwise stamps a random
      # htmlwidget-xxxx id into every page, so identical rebuilds diff in git
      w$elementId <- sub("\\.html$", "", f)
      w <- stable_plotly_ids(w, w$elementId)
      # selfcontained = FALSE with one shared lib/: plotly.js lands in
      # output/lib once, each page stays small, and GitHub Pages serves both
      htmlwidgets::saveWidget(w, out_path(f), selfcontained = FALSE,
                              libdir = "lib", title = f)
      cat("wrote", f, "\n")
      save_heatmap(zs_fit, xs, view, key, spec)
    }
  }
}

## ---- rate-of-decline surfaces --------------------------------------------------------
#
# The estimand every "cost drop, %/qtr" cell summarises, un-collapsed: the
# instantaneous rate of cost decline at fixed accuracy, over (accuracy, date).
# For the linear specification this surface is flat at the tables' single
# number; for the quadratic and Box-Cox specifications it is not, and its
# shape is the answer to "declining faster at the frontier's top or bottom,
# early or late".
#
# Accuracy-direction fits: holding the index z fixed,
#   d ln c / dt = -(dz/dt) / (dz/d ln c),
# evaluated ALONG the surface -- ln c pinned by inverting the fit at each
# (accuracy, date) node (zacc_inverted's rising branch). Cost-direction fits
# model ln C(la, t) directly, so the rate is just d lnC/dt at the node.
# Everything is reported as the tables' transform 100*(1 - exp(rate/4)):
# percent cheaper per quarter, positive = falling.
#
# Blanked, not extrapolated: where the accuracy level exceeded the STATE OF
# THE ART at that date (the bundle's record matrix R is NA exactly there --
# no run had achieved the level yet, so a "decline in its cost" is a claim
# about something that did not exist, and unlike the fitted-surface views
# these pages draw no empirical layer to mark where the fit stops being
# disciplined by data); where the inversion has no admissible root; and where
# the fitted cost leaves the OBSERVED cost range -- the same clip the 2-D
# iso-accuracy contours apply. For the accuracy direction the near-fold
# region (cost slope under 0.05 logits per log dollar) blanks too, mirroring
# zacc_inverted's own guard: the ratio explodes exactly where the surface
# stops being invertible.

# Both return a [N_T, N_X] matrix of quarterly percentage declines on the
# bundle's (la, tc) lattice, NA where blanked. mask = FALSE (the heatmaps)
# skips the DATA masks -- the SOTA record and the observed cost range -- and
# keeps only the fit-domain blanks: there the empirical frontier is drawn
# over the map instead, and the color scale is pinned to the in-frontier
# range, so the unmasked region can be shown without being mistaken for, or
# distorting, the anchored one.
decline_acc <- function(fit, bl, mask = TRUE) {
  u <- zacc_inverted(fit, bl)               # ln cost over (la, tc)
  if (mask) {
    urng <- range(bl$s$lncost)
    u[u < urng[1] | u > urng[2]] <- NA_real_
  }
  tcm <- matrix(bl$tc, N_T, N_X)
  if (is_bc_fit(fit)) {
    # dz/d ln c = (bx + bxt*phit) * lc-transform slope * c = (...) * c^lc;
    # dz/dt     = (bt + bxt*phic) * tau^(lt - 1)
    p <- bc_pieces(fit)
    cost <- exp(u)
    tau <- tcm + bl$off
    den <- (p$bx + p$bxt * bc_tf(tau, p$lt)) * cost^p$lc
    num <- (p$bt + p$bxt * bc_tf(cost, p$lc)) * tau^(p$lt - 1)
    rate <- -num / ifelse(den > 0, den, NA_real_)
  } else {
    co <- frontier_coefs(fit)
    den <- frontier_dcost(co, u, tcm)
    rate <- -frontier_dtime(co, u, tcm) / ifelse(den > 0.05, den, NA_real_)
  }
  if (mask) rate[is.na(bl$R)] <- NA_real_   # beyond the state of the art then
  z <- 100 * (1 - exp(rate / 4))
  z[!is.finite(z)] <- NA_real_
  z
}

decline_cost <- function(fit, bl, mask = TRUE) {
  g <- grid_lt(bl$la, bl$tc)
  urng <- range(bl$s$lncost)
  if (is_cost_bc(fit)) {
    # d lnC/dt = (gt + gat*phia) * tau^(lt - 1): the index derivative, and
    # the index IS ln cost, so no response-transform factor divides it
    cf <- coef(fit)
    names(cf) <- sub("^beta_", "", names(cf))
    lam <- attr(fit, "bc_lambda")
    tau <- g$t + bl$off
    phia <- bc_tf(exp(g$x), lam[["lambda_odds"]])
    phit <- bc_tf(tau, lam[["lambda_time"]])
    eta <- cf[["(Intercept)"]] + cf[["phia"]] * phia + cf[["phit"]] * phit +
      cf[["phiat"]] * phia * phit
    lnC <- eta
    rate <- (cf[["phit"]] + cf[["phiat"]] * phia) *
      tau^(lam[["lambda_time"]] - 1)
  } else {
    co <- cost_coefs(fit)
    lnC <- cost_index(co, g$x, g$t)
    rate <- cost_dtime(co, g$x, g$t)
  }
  rate[is.na(lnC)] <- NA_real_              # fit domain: no cost attains it
  if (mask) rate[lnC < urng[1] | lnC > urng[2]] <- NA_real_
  z <- as_z(100 * (1 - exp(rate / 4)))
  if (mask) z[is.na(bl$R)] <- NA_real_      # beyond the state of the art then
  z[!is.finite(z)] <- NA_real_
  z
}

for (spec in c("lin", "quad", "bc")) {
  for (key in names(fits)) {
    fset <- fits[[key]][[spec]]
    zs <- setNames(lapply(benches, function(b) {
      bl <- bundles[[b]]
      if (key %in% ACC_KEYS) decline_acc(fset[[b]], bl)
      else decline_cost(fset[[b]], bl)
    }), benches)
    xs <- lapply(bundles, function(bl) plogis(bl$la))
    w <- page(NULL, zs, xs, "decline")
    f <- sprintf("surface3d_%s_%s_decline.html", key, spec)
    w$elementId <- sub("\\.html$", "", f)   # deterministic ids, as above
    w <- stable_plotly_ids(w, w$elementId)
    htmlwidgets::saveWidget(w, out_path(f), selfcontained = FALSE,
                            libdir = "lib", title = f)
    cat("wrote", f, "\n")
    # the heatmap twin: UNMASKED (the frontier band is drawn over it
    # instead), with each panel's colors anchored to its masked variant's
    # range -- exactly the range its 3-D scene colors span
    zs_full <- setNames(lapply(benches, function(b) {
      bl <- bundles[[b]]
      if (key %in% ACC_KEYS) decline_acc(fset[[b]], bl, mask = FALSE)
      else decline_cost(fset[[b]], bl, mask = FALSE)
    }), benches)
    save_heatmap(zs_full, xs, "decline", key, spec, zs_ref = zs)
  }
}
