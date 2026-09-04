# Build the analysis dataset from source CSVs. 
#
#   cost_truncated_curves.csv        one row per model x benchmark x effort x budget run
#   Model versions-Grid view.csv  release dates, matched by id prefix
#
# Steps, in the .do file's order (order matters: the duplicate drop keeps the
# first row in collapse's sort order, so sorting must match Stata's):
#   1  strip -maas and trailing -YYYY-MM-DD suffixes from model names
#   2  collapse: MEAN of acc/tokens/cost/full_run_acc, SUM of n_samples,
#      by model-benchmark-effort-budget. After step 1 some names merge, so this
#      really does aggregate: 8840 rows -> 8660.
#   3  release date = earliest Version.release.date among registry ids that start
#      with the model name (case-insensitive)
#   4  cost, lncost, benchmarkid
#   5  duplicates drop on benchmark-model-effort-acc-cost-releasedate. This is
#      the deduplication the analysis relies on: repeated budget steps that
#      elicit the same effort and score collapse to one row.
#
# DELIBERATE DIVERGENCE FROM THE .do FILE: line 39 there keys the duplicate drop
# on model-effort-acc-cost-releasedate WITHOUT benchmark, so a row could be
# dropped because an identical (model, effort, acc, cost) combination appeared
# under a DIFFERENT benchmark -- most easily acc = 0 at the same cost. Adding
# benchmark keeps such rows. If the .do file is re-run, update its line 39 to
# match or the two pipelines will disagree.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
suppressMessages(library(dplyr))

STATA_EPOCH_OFFSET <- 3653   # days from 1960-01-01 (Stata) to 1970-01-01 (R)

# The PRIMARY benchmarks: first in every figure, table and console report,
# and the only ones entering the ECI pooling. The rest are reported alongside
# but pooled nowhere. fm13 is fm_tiers1_3_v2's legacy key (see build_runs).
PRIMARY_BENCHES <- c("aime", "chess", "fm13", "gpqa", "mystery")

# Canonical benchmark order: primaries first (in PRIMARY_BENCHES order), the
# rest alphabetical. Every script's `benches` vector and loop runs through
# this, so results ordering cannot drift between outputs.
bench_levels <- function(x) {
  u <- unique(x)
  c(intersect(PRIMARY_BENCHES, u), sort(setdiff(u, PRIMARY_BENCHES)))
}

TASK_LABEL <- c(
  aime                     = "OTIS Mock AIME 2024-2025",
  chess                    = "Chess Puzzles",
  fm13                     = "FrontierMath-Tiers-1-3-v2-Private",
  fm_2025_02_private       = "FrontierMath-2025-02-28-Private",
  fm_tier4_2025_07_private = "FrontierMath-Tier-4-2025-07-01-Private",
  fm_tier4_v2              = "FrontierMath-Tier-4-v2-Private",
  gpqa                     = "GPQA diamond",
  math_lvl5                = "MATH level 5",
  mystery                  = "Mystery Game Puzzles",
  simpleqa                 = "SimpleQA Verified",
  swe_bench_verified       = "SWE-Bench verified")

# alpha_b = estimated_slope_scaled: each benchmark's 2PL discrimination, in
# logits per point of Epoch's ECI (Epoch Capabilities Index) scale -- the scale
# anchored at Claude 3.5 Sonnet = 130 and GPT-5 = 150, on which the `edi`
# column gives the benchmark's difficulty. TASK_LABEL above holds each
# benchmark's exact name in that table, so it doubles as the join key.
#
# It lives here, beside the join key, because BOTH table scripts pool with it:
# regression_tables.R across the per-benchmark coefficient columns, and
# cost_frontier_report.R for the rate table's pooled row. Two copies of this
# join could silently disagree about which discriminations are in force.
#
# data/edi_scores.csv was downloaded from https://epoch.ai/data/edi_scores.csv
# on 2026-08-20.
ECI <- read.csv(data_path("edi_scores.csv"), stringsAsFactors = FALSE)
ALPHA <- setNames(ECI$estimated_slope_scaled[match(TASK_LABEL, ECI$benchmark_name)],
                  names(TASK_LABEL))
# D_b, the same table's `edi` difficulty column, in ECI points. The 2PL writes
# logit(a) = alpha_b*(C - D_b), so logit(a)/alpha_b + D_b recovers the model's
# capability C on the anchored ECI scale (Claude 3.5 Sonnet = 130, GPT-5 =
# 150). pooled_cost_runs (cost_frontier.R) uses it to place pooled logits on
# that scale; its benchmark fixed effects absorb the shift, so the fits are
# invariant to including D_b -- it buys interpretable units, nothing else.
EDI <- setNames(ECI$edi[match(TASK_LABEL, ECI$benchmark_name)],
                names(TASK_LABEL))

# Greedy .* so the LAST occurrence wins, matching Stata's regexcapturenamed on
# "(?<stub>.*)-maas": the captured stub is everything before the final marker.
clean_model_name <- function(m) {
  m <- sub("^(.*)-maas.*$", "\\1", m)
  sub("^(.*)-[0-9]{4}-[0-9]{2}-[0-9]{2}.*$", "\\1", m)
}

build_runs <- function(curves_csv = data_path("cost_truncated_curves.csv"),
                       models_csv = data_path("Model versions-Grid view.csv")) {

  x <- read.csv(curves_csv, stringsAsFactors = FALSE)
  x$model <- clean_model_name(x$model)

  agg <- x %>%
    group_by(model, benchmark, effort, budget) %>%
    summarise(across(c(acc, mean_tokens_used, cost_per_task_usd, full_run_acc),
                     ~ mean(.x, na.rm = TRUE)),
              n_samples = sum(n_samples, na.rm = TRUE),
              .groups = "drop") %>%
    as.data.frame()

  # The *_public FrontierMath variants are excluded: tiny samples (their
  # staircases have as few as 5 corners and their ECI file carries no
  # discriminations for them), and their private counterparts are in the
  # sample. fm_tiers1_3_v2 keeps its old fm13 key so history, labels and
  # joins carry over, and mystery_agent shortens to mystery -- the key
  # PRIMARY_BENCHES, TASK_LABEL and LABELS all use.
  agg$benchmark[agg$benchmark == "fm_tiers1_3_v2"] <- "fm13"
  agg$benchmark[agg$benchmark == "mystery_agent"]  <- "mystery"
  agg <- agg[!grepl("_public$", agg$benchmark), , drop = FALSE]

  # release date: earliest registry entry whose id starts with the model name
  mv <- read.csv(models_csv, stringsAsFactors = FALSE)
  mv$releasedate <- as.Date(mv$Version.release.date)
  id_lower <- tolower(mv$id)
  first_release <- function(m) {
    hit <- startsWith(id_lower, tolower(m)) & !is.na(mv$releasedate)
    if (!any(hit)) return(as.Date(NA))
    min(mv$releasedate[hit])
  }
  lut <- vapply(unique(agg$model), first_release, as.Date(NA))
  agg$releasedate <- as.Date(unname(lut[agg$model]), origin = "1970-01-01")

  # display-friendly model name, from the same registry: the Model column of
  # the matching row, preferring an exact id match and otherwise the SHORTEST
  # matching id -- the closest thing to exact, where the prefix rule above
  # would let "gpt-5" match gpt-5-mini's row too. A model absent from the
  # registry (or with a blank Model cell) keeps its raw name.
  display_name <- function(m) {
    hit <- which(startsWith(id_lower, tolower(m)))
    if (!length(hit)) return(m)
    ex <- hit[id_lower[hit] == tolower(m)]
    i <- if (length(ex)) ex[1] else hit[which.min(nchar(mv$id[hit]))]
    if (is.na(mv$Model[i]) || mv$Model[i] == "") m else mv$Model[i]
  }
  dlut <- vapply(unique(agg$model), display_name, "")
  agg$model_display <- unname(dlut[agg$model])

  agg$releaseyear  <- (as.numeric(agg$releasedate) + STATA_EPOCH_OFFSET) / 365.25
  agg$cost         <- agg$cost_per_task_usd
  agg$lncost       <- log(agg$cost)
  agg$lncost[!is.finite(agg$lncost)] <- NA_real_   # Stata's ln(0) is missing
  agg$benchmarkid  <- as.integer(factor(agg$benchmark))

  # Rescale each benchmark from its GUESSING FLOOR to 1, onto [0, 1].
  #
  # Unanswered questions score at the probability of guessing correctly, so
  # every benchmark's accuracy has a floor: 25% for GPQA's four-way choices,
  # 0.1% for AIME's integer answers, ~9.2% for the mystery puzzles, and 0
  # wherever a wrong answer and no answer score alike (the FrontierMath set,
  # chess, MATH, SimpleQA, SWE-Bench). Rescaled, accuracy means the same thing
  # everywhere -- the share of the way from guessing to perfect -- which is
  # what makes the logit, the pooling and the cross-benchmark figures
  # comparable. A floor of 0 leaves the column untouched, so this generalises
  # the GPQA-only censoring it replaces rather than adding a step.
  #
  # The floor is READ OFF THE DATA rather than kept as a hand-written list
  # that would silently rot as benchmarks are added: the truncation curves
  # give every benchmark many fully-truncated rows that answer nothing and so
  # score exactly the floor, a point mass of 169-1294 rows here. Taking the
  # most common value in the bottom half of the distribution finds it and
  # cannot be fooled by a mass of perfect scores at the top. A run BELOW its
  # floor (it answered, and answered worse than chance) clips to 0, as GPQA's
  # scores always did.
  guess_floor <- function(a) {
    lo <- a[a <= median(a)]
    as.numeric(names(sort(table(lo), decreasing = TRUE))[1])
  }
  floors <- vapply(split(agg$acc, agg$benchmark), guess_floor, numeric(1))
  nz <- floors[floors > 0]
  if (length(nz))
    message("guessing floors rescaled to 0: ",
            paste(sprintf("%s %.4g", names(nz), nz), collapse = ", "))
  for (b in names(floors)) {
    g <- floors[[b]]
    i <- agg$benchmark == b
    agg$acc[i] <- (pmax(g, agg$acc[i]) - g) / (1 - g)
  }

  # Stata's collapse leaves data sorted by its by() variables, and duplicates
  # drop keeps the FIRST row in that order -- so sort the same way, with radix
  # (byte) collation to match Stata rather than the R locale.
  agg <- agg[order(agg$model, agg$benchmark, agg$effort, agg$budget,
                   method = "radix"), ]
  dup <- duplicated(agg[c("benchmark", "model", "effort", "acc", "cost",
                          "releasedate")])
  agg <- agg[!dup, , drop = FALSE]
  rownames(agg) <- NULL
  agg
}
