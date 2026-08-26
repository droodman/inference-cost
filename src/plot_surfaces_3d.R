# Interactive 3-D renderings of the Pareto-grid fits, for the plot viewer's
# 2-D / 3-D toggle: each page shows, per benchmark, the EMPIRICAL surface
# (the staircase P_t(c) or the record cost ln C_a(t), NA holes where
# undefined) with the FITTED surface overlaid at half opacity -- rotate the
# scene and the shelf-and-cliff geometry, and how each specification does or
# does not follow it, is directly visible.
#
# The frontier-per-se models only -- the Pareto-grid pair and the envelope
# pair -- because those are the models whose empirical reference IS the
# staircase/record surface drawn beneath them (the run-level models' honest
# reference would be the run cloud, a different figure). The envelope gains
# the most from rotation: its whole claim is one-sided, so where the
# wireframe touches the data surface and how far it floats elsewhere is the
# visual form of its slack diagnostics. Each model is rendered in three
# specifications and BOTH orientations:
#
#   frontier     x = ln cost, y = year, z = accuracy      (performance surface)
#   isoaccuracy  x = logit accuracy, y = year, z = ln cost (cost surface)
#
# so 4 model keys x 3 specifications x 2 views = 24 HTML files,
# surface3d_<key>_<spec>_<view>.html. Each direction is NATIVE in one view and
# INVERTED into the other; the inversions take the rising branch and leave NA
# where the surface bends back or the transform has no preimage -- holes, not
# folds, exactly as the 2-D figures blank.
#
# The widgets are static HTML + client-side WebGL (htmlwidgets/plotly): no
# server, so they serve from the repo or GitHub Pages as-is. Written with
# selfcontained = FALSE and a SHARED output/lib directory, so plotly.js is
# committed once (~4 MB) and each page is small.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("cost_frontier.R")   # pulls the whole fitting stack
suppressMessages(library(plotly))

N_X <- 100   # lattice density along the non-time axis, matching the fits' grids
N_T <- 40

d <- load_runs()
benches <- sort(unique(d$benchmark))
tbar <- bench_tbar(d)

## ---- fits: both frontier-per-se pairs, three specifications each ------------------

ACC_KEYS <- c("paretologit", "envelope")   # which keys are accuracy-direction

fit_one <- function(key, b, spec) {
  s <- d[d$benchmark == b, ]
  switch(key,
    paretologit  = fit_pareto_logit(s, TIME_FORMS[[spec]]),
    envelope     = fit_envelope(s, TIME_FORMS[[spec]]),
    costgridols  = fit_lncost_grid(s, COST_FORMS[[spec]]),
    costenvelope = fit_cost_envelope(s, COST_FORMS[[spec]]))
}

fits <- list()
for (key in c(ACC_KEYS, "costgridols", "costenvelope")) {
  for (spec in c("lin", "quad"))
    fits[[key]][[spec]] <- setNames(lapply(benches, function(b)
      fit_one(key, b, spec)), benches)
  fits[[key]]$bc <- if (key %in% ACC_KEYS) fit_bc_by(key, d) else
    setNames(lapply(benches, function(b)
      fit_cost_bc(key, d[d$benchmark == b, ])), benches)
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
  la  <- seq(min(su$la), max(su$la), length.out = N_X)
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
    as_z(bc_mu(p$b0 + p$bx * phic + p$bt * phit + p$bxt * phic * phit, p$lo))
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
    phic <- (bc_tf(exp(g$x), p$lo) - p$b0 - p$bt * phit) / bb
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
    phia <- (bc_tf(exp(g$x), lam[["lambda_cost"]]) - cf[["(Intercept)"]] -
               cf[["phit"]] * phit) / bb
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

# Initial viewpoint: the box's near lower corner is the one with the LOWEST
# date and LOWEST cost -- the eye sits in the (min x, min y) quadrant, so the
# frontier view leads with cheap-and-early and the iso view with the low end
# of the level axis at the earliest date (cost is vertical there, so its
# minimum is the bottom face by construction).
CAMERA <- list(eye = list(x = -1.5, y = -1.5, z = 0.9))

# One page: four scenes (2 x 2 benchmark quadrants), each with the empirical
# surface (filled, plasma ramp) and the fitted one as a SEE-THROUGH WIREFRAME:
# hidesurface drops the fill entirely and the x/y contour lines draw the mesh,
# so the fitted shape is legible without occluding the data surface behind it.
# Ink-colored, not black -- black segments would vanish on the dark backdrop.
page <- function(zs_emp, zs_fit, xs, view) {
  doms <- list(list(x = c(0, .48), y = c(.55, 1)), list(x = c(.52, 1), y = c(.55, 1)),
               list(x = c(0, .48), y = c(0, .45)), list(x = c(.52, 1), y = c(0, .45)))
  p <- plot_ly()
  lay <- list(paper_bgcolor = SURFACE, font = list(color = INK),
              showlegend = FALSE, margin = list(l = 0, r = 0, t = 30, b = 0))
  ann <- list()
  for (i in seq_along(benches)) {
    b <- benches[i]
    bl <- bundles[[b]]
    sc <- if (i == 1) "scene" else paste0("scene", i)
    p <- add_surface(p, x = xs[[b]], y = bl$year, z = zs_emp[[b]],
                     scene = sc, colorscale = "Plasma", showscale = FALSE,
                     name = "empirical")
    xr <- range(xs[[b]]); yr <- range(bl$year)
    mesh <- function(rng, n = 14) list(show = TRUE, color = INK,
                                       width = 2, start = rng[1],
                                       end = rng[2], size = diff(rng) / n)
    p <- add_surface(p, x = xs[[b]], y = bl$year, z = zs_fit[[b]],
                     scene = sc, hidesurface = TRUE, showscale = FALSE,
                     contours = list(x = mesh(xr), y = mesh(yr)),
                     name = "fitted")
    zx <- if (view == "frontier")
      list(xaxis = ax_cost(range(xs[[b]])),
           zaxis = c(ax("accuracy"), list(range = c(0, 1))))
    else
      list(xaxis = ax("logit accuracy"),
           zaxis = c(ax_cost(bl$zrng_cost), list(range = bl$zrng_cost)))
    lay[[sc]] <- c(list(domain = doms[[i]], yaxis = ax("year"),
                        aspectmode = "cube", camera = CAMERA), zx)
    ann[[i]] <- list(text = LABELS[[b]], x = mean(doms[[i]]$x),
                     y = doms[[i]]$y[2], xref = "paper", yref = "paper",
                     showarrow = FALSE, font = list(color = INK, size = 14))
  }
  lay$annotations <- ann
  do.call(layout, c(list(p), lay))
}

VIEWS <- c("frontier", "isoaccuracy")
for (spec in c("lin", "quad", "bc")) {
  for (key in names(fits)) {
    fset <- fits[[key]][[spec]]
    for (view in VIEWS) {
      xs <- lapply(bundles, function(bl) if (view == "frontier") bl$lnc else bl$la)
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
      # selfcontained = FALSE with one shared lib/: plotly.js lands in
      # output/lib once, each page stays small, and GitHub Pages serves both
      htmlwidgets::saveWidget(w, out_path(f), selfcontained = FALSE,
                              libdir = "lib", title = f)
      cat("wrote", f, "\n")
    }
  }
}
