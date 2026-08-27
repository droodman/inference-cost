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

  agg$releaseyear  <- (as.numeric(agg$releasedate) + STATA_EPOCH_OFFSET) / 365.25
  agg$cost         <- agg$cost_per_task_usd
  agg$lncost       <- log(agg$cost)
  agg$lncost[!is.finite(agg$lncost)] <- NA_real_   # Stata's ln(0) is missing
  agg$benchmarkid  <- as.integer(factor(agg$benchmark))

  # censor GPQA scores, for which 25% is the expected minimum, to [25%,100%]
  # the full_run_acc values are almost all in this range, but the method of synthesizing
  # rows by budget level put them in [0%,100%]
  agg$acc[agg$benchmark=="gpqa"] <- (pmax(.25, agg$acc[agg$benchmark=="gpqa"]) - .25) / .75

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
