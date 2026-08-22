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
src_source("boxcox_frontier.R")   # fit_bc(), for the Box-Cox figures

d <- load_runs(drop_gpt4o_chess = FALSE)
tbar  <- bench_tbar(d)     # from t - tc, so it matches whatever the fit used
dates <- bench_dates(d)
benches <- sort(unique(d$benchmark))

specs <- fit_all_specs(d)

NOTES_BASE <- paste("All models retained. Data prepared once in prepare_data.R,",
                    "which also drops the duplicate rows, so every specification",
                    "sees the same sample.")
# Describes the test rather than asserting its outcome. The previous wording
# ("the quadratic time term is not significant on any benchmark") both went stale
# -- the LR table it cited reports p < 1e-7 on aime and gpqa -- and named a
# single term the test never isolated: quad adds all three second-order terms at
# once, so the LR is a joint test of the block.
NOTES_QUAD <- paste("Quadratic adds all three second-order terms (cost^2, time^2,",
                    "cost x time) jointly; see this script's LR table and",
                    "monotonicity check for whether the block earns its keep.")

for (k in names(specs)) {
  sp <- specs[[k]]
  curves <- frontier_curves(sp$fits, d, dates, tbar)
  p <- frontier_plot(
    curves, d,
    ylab = if (sp$family == "S") "Fitted accuracy" else "Frontier accuracy",
    notes = c(NOTES_BASE, if (sp$time == "quad") NOTES_QUAD))
  f <- sprintf("frontier_progression_%s.png", k)
  ggsave(out_path(f), p, width = 10, height = 7.5, dpi = 200,
         device = ragg::agg_png)
  cat("wrote", f, "\n")
}

## ---- the Box-Cox specification ---------------------------------------------------
# The third specification (boxcox_frontier.R): profiled per benchmark and per
# family, so it does not flow through fit_all_specs()' formula grid. Refitted
# here rather than cached from regression_tables.R, matching how every other
# fit in this repo is recomputed from the data each run. S is profiled first --
# its glm inner loop is cheap -- and its lambdas seed the slower SFA profiles.
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
  curves <- frontier_curves(fits, d, dates, tbar)
  p <- frontier_plot(
    curves, d,
    ylab = if (fam == "S") "Fitted accuracy" else "Frontier accuracy",
    notes = c(NOTES_BASE, NOTES_BC))
  f <- sprintf("frontier_progression_%s_bc.png", fam)
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

report_quadratic_pathologies(specs, d)
report_scale_lr(specs, d)
report_quad_lr(specs)
