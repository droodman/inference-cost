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
#   isoaccuracy  x = logit accuracy, y = year, z = ln cost (cost surface)
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
# The widgets are static HTML + client-side WebGL (htmlwidgets/plotly): no
# server, so they serve from the repo or GitHub Pages as-is. Written with
# selfcontained = FALSE and a SHARED output/lib directory, so plotly.js is
# committed once (~4 MB) and each page is small.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("fit_store.R")   # pulls the whole fitting stack
suppressMessages(library(plotly))

N_X <- 100   # lattice density along the non-time axis, matching the fits' grids
N_T <- 40

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

# An accuracy axis: internally logit accuracy (the lattice the fits live on),
# labelled in PERCENT, the same trick ax_cost plays with dollars.
ax_acc <- function(rng) {
  a <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
  keep <- qlogis(a) >= rng[1] & qlogis(a) <= rng[2]
  c(ax("accuracy"), list(tickvals = qlogis(a[keep]),
                         ticktext = sprintf("%g%%", 100 * a[keep])))
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
      list(xaxis = ax("logit accuracy"),
           zaxis = c(ax_cost(bl$zrng_cost), list(range = bl$zrng_cost)))
    else
      list(xaxis = ax_acc(range(xs[[b]])),
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
# Blanked, not extrapolated: where the inversion has no admissible root, and
# where the fitted cost leaves the OBSERVED cost range -- the same clip the
# 2-D iso-accuracy contours apply, because out there the "decline" is the
# functional form talking to itself. For the accuracy direction the near-fold
# region (cost slope under 0.05 logits per log dollar) blanks too, mirroring
# zacc_inverted's own guard: the ratio explodes exactly where the surface
# stops being invertible.

# Both return a [N_T, N_X] matrix of quarterly percentage declines on the
# bundle's (la, tc) lattice, NA where blanked.
decline_acc <- function(fit, bl) {
  u <- zacc_inverted(fit, bl)               # ln cost over (la, tc)
  urng <- range(bl$s$lncost)
  u[u < urng[1] | u > urng[2]] <- NA_real_
  tcm <- matrix(bl$tc, N_T, N_X)
  if (is_bc_fit(fit)) {
    # dz/d ln c = (bx + bxt*phit) * lc-transform slope * c = (...) * c^lc;
    # dz/dt     = (bt + bxt*phic) * tau^(lt - 1)
    p <- bc_pieces(fit)
    cost <- exp(u)
    tau <- tcm + bl$off
    den <- (p$bx + p$bxt * bc_tf(tau, p$lt)) * cost^p$lc
    num <- (p$bt + p$bxt * bc_tf(cost, p$lc)) * tau^(p$lt - 1)
    g <- -num / ifelse(den > 0, den, NA_real_)
  } else {
    co <- frontier_coefs(fit)
    den <- frontier_dcost(co, u, tcm)
    g <- -frontier_dtime(co, u, tcm) / ifelse(den > 0.05, den, NA_real_)
  }
  100 * (1 - exp(g / 4))
}

decline_cost <- function(fit, bl) {
  g <- grid_lt(bl$la, bl$tc)
  urng <- range(bl$s$lncost)
  if (is_cost_bc(fit)) {
    # d lnC/dt = (gt + gat*phia) * tau^(lt - 1) / C^lC  (the index derivative
    # over d phi(C)/d lnC)
    cf <- coef(fit)
    names(cf) <- sub("^beta_", "", names(cf))
    lam <- attr(fit, "bc_lambda")
    tau <- g$t + bl$off
    phia <- bc_tf(exp(g$x), lam[["lambda_odds"]])
    phit <- bc_tf(tau, lam[["lambda_time"]])
    eta <- cf[["(Intercept)"]] + cf[["phia"]] * phia + cf[["phit"]] * phit +
      cf[["phiat"]] * phia * phit
    lnC <- ln_bc_inv(eta, lam[["lambda_cost"]])
    rate <- (cf[["phit"]] + cf[["phiat"]] * phia) * tau^(lam[["lambda_time"]] - 1) /
      exp(lam[["lambda_cost"]] * lnC)
  } else {
    co <- cost_coefs(fit)
    lnC <- cost_index(co, g$x, g$t)
    rate <- cost_dtime(co, g$x, g$t)
  }
  rate[is.na(lnC) | lnC < urng[1] | lnC > urng[2]] <- NA_real_
  as_z(100 * (1 - exp(rate / 4)))
}

for (spec in c("lin", "quad", "bc")) {
  for (key in names(fits)) {
    fset <- fits[[key]][[spec]]
    zs <- setNames(lapply(benches, function(b) {
      bl <- bundles[[b]]
      if (key %in% ACC_KEYS) decline_acc(fset[[b]], bl)
      else decline_cost(fset[[b]], bl)
    }), benches)
    xs <- lapply(bundles, function(bl) bl$la)
    w <- page(NULL, zs, xs, "decline")
    f <- sprintf("surface3d_%s_%s_decline.html", key, spec)
    htmlwidgets::saveWidget(w, out_path(f), selfcontained = FALSE,
                            libdir = "lib", title = f)
    cat("wrote", f, "\n")
  }
}
