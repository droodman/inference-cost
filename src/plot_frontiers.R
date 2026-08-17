# Fitted frontier by cost, one curve per quarter, colour showing elapsed time.
#
# Six figures: three inefficiency families crossed with linear and quadratic
# time. See fit_specs.R for the grid. S (no inefficiency term) is the reference
# case: with no u its "frontier" IS the conditional mean, so roughly half the
# runs sit above it by construction. A and B are worth having only insofar as
# they beat that.
#
# marginaleffects is not used -- it dispatches on classes with predict() methods
# and cannot handle a raw maxLik fit -- but the prediction is one line.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("fit_specs.R")   # brings panel_frontier.R and frontier_viz.R with it

d <- load_runs(drop_gpt4o_chess = FALSE)
tbar  <- bench_tbar(d)     # from t - tc, so it matches whatever the fit used
dates <- bench_dates(d)
benches <- sort(unique(d$benchmark))

specs <- fit_all_specs(d)

NOTES_BASE <- paste("All models retained. Data prepared once in prepare_data.R,",
                    "which also drops the duplicate rows, so every specification",
                    "sees the same sample.")
NOTES_QUAD <- paste("The quadratic time term is not significant on any benchmark",
                    "-- see the LR table printed by this script -- so curvature",
                    "here should not be read as established.")

for (k in names(specs)) {
  sp <- specs[[k]]
  curves <- frontier_curves(sp$fits, d, dates, tbar)
  p <- frontier_plot(
    curves, d,
    title = sprintf("Fitted accuracy frontier by cost per task -- model %s (%s)",
                    sp$family, sp$time),
    subtitle = sp$subtitle,
    ylab = if (sp$family == "S") "Fitted accuracy" else "Frontier accuracy",
    notes = c(NOTES_BASE, if (sp$time == "quad") NOTES_QUAD))
  f <- sprintf("frontier_progression_%s.png", k)
  ggsave(out_path(f), p, width = 10, height = 7.5, dpi = 200,
         device = ragg::agg_png)
  cat("wrote", f, "\n")
}

## ---- does each fitted curve actually envelope the data? -------------------------------
# For S this is a sanity check, not a criticism: a conditional mean SHOULD have
# about half the runs above it. For A and B a high share means the inefficiency
# term is not buying a frontier.

cat("\nshare of observed runs lying ABOVE the fitted curve\n")
cat(sprintf("%-6s %s\n", "bench",
            paste(sprintf("%10s", names(specs)), collapse = "")))
for (b in benches) {
  s <- d[d$benchmark == b, ]
  sh <- vapply(specs, function(sp) {
    co <- frontier_coefs(sp$fits[[b]])
    mean(s$acc > plogis(frontier_index(co, s$lncost, s$tc)))
  }, numeric(1))
  cat(sprintf("%-6s %s\n", b, paste(sprintf("%9.1f%%", 100 * sh), collapse = "")))
}

report_slope_monotonicity(specs, d)
report_scale_lr(specs, d)
report_quad_lr(specs)
