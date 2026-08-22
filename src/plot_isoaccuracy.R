# Iso-accuracy cost contours: what a given performance level costs over time.
#
# The companion to plot_frontiers.R with the axes swapped -- release date across,
# cost up the log scale, one contour per accuracy target. Same six parametric
# specifications (three inefficiency families x linear/quadratic time), same
# palette and chrome.
#
# Each contour holds ACCURACY fixed and lets cost vary, so these are isoquants,
# not isocosts. They were called isocost figures until Aug 2026; that name says
# the opposite of what is drawn, hence the rename.
#
# Colour carries accuracy here rather than date. It is still a magnitude, so the
# same sequential blue ramp is right, and the observed runs are plotted at their
# (date, cost) coloured by the accuracy they actually reached -- so a run's shade
# can be compared directly against the contour it lies on.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("fit_specs.R")
src_source("boxcox_frontier.R")   # fit_bc(), for the Box-Cox figures

LEVELS <- seq(0.05, 0.95, by = 0.10)

d <- load_runs(drop_gpt4o_chess = FALSE)
tbar <- bench_tbar(d)
benches <- sort(unique(d$benchmark))

specs <- fit_all_specs(d)

# Pin the panels to the observed date and cost ranges, so blanked contour
# segments cannot pull the axes around.
axis_ranges <- do.call(rbind, lapply(benches, function(b) {
  s <- d[d$benchmark == b, ]
  data.frame(benchmark = b, date = range(s$releasedate), cost = range(s$cost))
}))

NOTES_BASE <- c(
  paste("Contours are accuracy targets from 5% to 95%; a falling contour means",
        "the same performance costs less over time. Dots are observed runs,",
        "coloured by the accuracy they reached."),
  paste("Contours are cut where they leave the benchmark's observed cost range,",
        "so a missing contour means that target lies outside the data rather",
        "than that it costs nothing."))
NOTES_QUAD <- paste("Quadratic adds all three second-order terms (cost^2, time^2,",
                    "cost x time) jointly, so curvature here is the whole block's",
                    "doing, not any one term's.")

for (k in names(specs)) {
  sp <- specs[[k]]
  curves <- iso_acc_curves(sp$fits, d, tbar, levels = LEVELS)
  p <- iso_acc_plot(
    curves, d,
    notes = c(NOTES_BASE,
              if (sp$time == "quad") c(NOTES_QUAD, ISO_BRANCH_NOTE)),
    ranges = axis_ranges)
  f <- sprintf("isoaccuracy_%s.png", k)
  ggsave(out_path(f), p, width = 10, height = 7.5, dpi = 200,
         device = ragg::agg_png)
  cat("wrote", f, "\n")
}

## ---- the Box-Cox specification ---------------------------------------------------
# Same fits as plot_frontiers.R's Box-Cox section, recomputed here as every
# other specification's are. The contours stay closed form: phi(cost) =
# (logit(target) - b0 - b_t phi(t)) / (b_x + b_xt phi(t)), inverted through the
# transform; phi is monotone, so each contour has ONE branch and no fold --
# ISO_BRANCH_NOTE does not apply.
NOTES_BC <- paste("Box-Cox: the logit is linear in phi(cost), phi(years since",
                  "mid-2020) and their product, with the transform parameters",
                  "profiled per benchmark; monotone in cost at every date.")
bc_lam <- list()
for (fam in c("S", "A", "B")) {
  fits <- setNames(lapply(benches, function(b) {
    fit_bc(fam, d[d$benchmark == b, ],
           lambda_start = if (is.null(bc_lam[[b]])) c(0, 1) else bc_lam[[b]])
  }), benches)
  if (fam == "S")
    for (b in benches) bc_lam[[b]] <- unname(attr(fits[[b]], "bc_lambda"))
  curves <- iso_acc_curves(fits, d, tbar, levels = LEVELS)
  p <- iso_acc_plot(curves, d, notes = c(NOTES_BASE, NOTES_BC),
                    ranges = axis_ranges)
  f <- sprintf("isoaccuracy_%s_bc.png", fam)
  ggsave(out_path(f), p, width = 10, height = 7.5, dpi = 200,
         device = ragg::agg_png)
  cat("wrote", f, "\n")
}

## ---- how much of each contour is inside the data? ------------------------------------
# A contour that is mostly blank is the model extrapolating, not a cost decline.

# Rising branch only. iso_acc_curves() now returns both branches, so averaging
# over all its rows would halve every figure here for reasons that have nothing
# to do with coverage.
cat("\nshare of each contour inside the observed cost range (%), model A linear\n")
cv <- iso_acc_curves(specs$A_lin$fits, d, tbar, levels = LEVELS)
cvr <- cv[cv$branch == "rising", ]
print(round(100 * tapply(!is.na(cvr$cost), list(cvr$benchmark, cvr$acc), mean)))

## ---- implied cost decline at 50% accuracy --------------------------------------------
# Cost at the first vs last observed date. Under the linear specification the
# ratio is the same at every accuracy level (the contours are parallel in logs);
# the quadratic breaks that, so the two columns can differ.

cat("\ncost of 50% accuracy, first vs last observed date\n")
cat(sprintf("%-6s %-6s %12s %12s %8s %10s\n",
            "spec", "bench", "at first", "at last", "ratio", "%/yr"))
for (k in names(specs)) {
  for (b in benches) {
    co <- frontier_coefs(specs[[k]]$fits[[b]])
    s  <- d[d$benchmark == b, ]
    ends <- range(s$releasedate)
    tc <- as_t(ends) - tbar[[b]]
    cst <- exp((qlogis(0.5) - co[["b0"]] - co[["bt"]] * tc - co[["btt"]] * tc^2) /
                 co[["bx"]])
    yrs <- diff(as_t(ends))
    cat(sprintf("%-6s %-6s %12.4f %12.4f %8.3f %9.1f%%\n", k, b, cst[1], cst[2],
                cst[2] / cst[1], 100 * (1 - (cst[2] / cst[1])^(1 / yrs))))
  }
}
