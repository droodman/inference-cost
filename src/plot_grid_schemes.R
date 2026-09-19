# Four ways one COULD lay the node lattice over the (accuracy, date) plane,
# drawn on one benchmark as an illustration -- BENCH below picks which, and
# the filename follows it. The MODEL-FREE decline rate -- the
# "Model-free" column of rate_comparison -- reads the Pareto staircase
# directly, with no functional form anywhere, but it still has to decide
# WHERE to read it. That choice is this figure's subject: the lattice is what
# the average is taken over, so it fixes how the (level, date) plane is
# weighted.
#
# Two binary choices cross:
#
#   SPACING   even in ACCURACY, or even in LOGIT accuracy. Logit is the scale
#             the models' formulas use; accuracy is the scale the reader
#             thinks in. Evenly spaced logits pile levels up near 0 and 1 and
#             thin them through the middle, because the logistic saturates.
#
#   CEILING   CLIPPED to SOTA, or SCALED within it. Clipped keeps one lattice
#             for all dates and drops the levels above the day's best score,
#             so early dates carry few nodes and late dates many. Scaled
#             restretches the same node COUNT under each day's ceiling, so
#             every date carries equal weight and the levels track the
#             frontier upward.
#
# decline_nodes() (envelope_frontier.R) implements the top-left cell:
# accuracy-even midpoints, clipped. Its own comment records that the scaled
# variant -- a_j(t) = SOTA(t)*(j - 1/2)/n -- was tried there and replaced, so
# the bottom-left cell is not hypothetical but the predecessor.
#
# The DRAWN lattice is deliberately coarse -- 10 x 10 -- because the point is
# where nodes fall, not how many. The RATES are not: each is
# pareto_decline_qtr()'s exact arithmetic (record cost at both ends of a
# quarter, geometric mean, 100*(1 - exp(mean(dln)/(4dt)))) at the production
# resolution, so the top-left cell reproduces the published model-free figure
# for this benchmark rather than approximating it.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("cost_frontier.R")

BENCH     <- "gpqa"
N_ACC     <- 10     # drawn levels per date
N_DATE    <- 10     # drawn dates
N_ACC_FIT <- 100    # decline_nodes()' n_level, the resolution the rates use
DT        <- 0.25   # pareto_decline_qtr()'s horizon, in years

d  <- load_runs()
dd <- d[d$benchmark == BENCH, ]     # the full frame, as decline_nodes uses

# SOTA(t): the best accuracy anyone had reached by t -- a running maximum over
# date, ignoring cost. The ceiling every scheme below has to respect.
ss <- dd[order(dd$releasedate), ]
m  <- cummax(ss$acc)
keep <- c(TRUE, diff(m) > 0)
sota_steps <- rbind(
  data.frame(date = ss$releasedate[keep], acc = m[keep]),
  data.frame(date = max(ss$releasedate), acc = max(m)))

A_LO <- min(dd$acc[dd$acc > 0])     # decline_nodes' floor: smallest positive
A_HI <- max(dd$acc)
# The logit lattice needs a FINITE ceiling. This benchmark has runs at exactly
# 1.0, and qlogis(1) is Inf -- which silently collapses every logit rung onto
# 1.0 and yields a lattice of identical levels whose record never moves, i.e.
# a spurious 0%/qtr. The accuracy-even schemes are untouched by this (they
# never take a logit) and must keep A_HI = max(acc) to reproduce
# decline_nodes exactly, so the clamp applies to the logit branch alone.
A_HI_L <- max(dd$acc[dd$acc < 1])
TBAR <- bench_tbar(d)[[BENCH]]

# Midpoints of the range, even in accuracy or in logit -- no node at zero or
# exactly at the maximum, as decline_nodes requires.
rungs <- function(lo, hi, logit, n) {
  j <- (seq_len(n) - 0.5) / n
  if (logit) {
    hi <- min(hi, A_HI_L)           # no node at 1: logit cannot place one
    plogis(qlogis(lo) + (qlogis(hi) - qlogis(lo)) * j)
  } else lo + (hi - lo) * j
}

sota_tc <- function(tk) max(dd$acc[dd$tc <= tk])

# One scheme's nodes over a given set of dates and resolution: the coarse pair
# for drawing, the production pair for the rate.
nodes <- function(logit, scaled, tcs, n_acc) {
  do.call(rbind, lapply(tcs, function(tk) {
    ceil <- sota_tc(tk)
    a <- if (scaled) rungs(A_LO, ceil, logit, n_acc)
         else rungs(A_LO, A_HI, logit, n_acc)
    # clipped: one lattice for every date, the levels above the day's ceiling
    # dropped -- decline_nodes' `lev <= max(acc[tc <= tk])`
    if (!scaled) a <- a[a <= ceil]
    if (!length(a)) return(NULL)
    data.frame(tc = tk, a = a)
  }))
}

## ---- the model-free rate, pareto_decline_qtr's arithmetic ------------------------

srt <- dd[order(dd$lncost), ]
# ln C_a(t) = cheapest run no later than t scoring at least a; over cost-sorted
# runs cummax(acc) is the running best, so the first index reaching the level
# is the record.
rec <- function(t, levels) {
  el <- srt[srt$tc <= t, ]
  mm <- cummax(el$acc)
  vapply(levels, function(a) {
    j <- which(mm >= a)[1]
    if (is.na(j)) NA_real_ else el$lncost[j]
  }, numeric(1))
}

# Record at BOTH ends of the horizon, never the cost of a run that merely
# achieves the level: the staircase is flat between jumps, so a run's own cost
# typically overpays and measuring from it would count that slack as decline.
free_rate <- function(nd) {
  dln <- unlist(lapply(unique(nd$tc), function(tk) {
    a <- nd$a[nd$tc == tk]
    rec(tk + DT, a) - rec(tk, a)
  }))
  dln <- dln[is.finite(dln)]
  if (!length(dln)) return(NA_real_)
  100 * (1 - exp(mean(dln) / (4 * DT)))
}

# production dates: GRID_TIME_STEP cadence, keeping only those with a full
# horizon still ahead -- decline_nodes' own filter
tc_fit <- grid_tc_seq(range(dd$tc))
tc_fit <- tc_fit[tc_fit + DT <= max(dd$tc)]

## ---- the drawn lattice -----------------------------------------------------------

dts    <- seq(min(dd$releasedate), max(dd$releasedate), length.out = N_DATE)
tc_draw <- as_t(dts) - TBAR

SPACING <- c("Even in accuracy", "Even in logit accuracy")
CEILING <- c("Clipped to SOTA", "Scaled within SOTA")

grid <- do.call(rbind, lapply(c(FALSE, TRUE), function(lg)
  do.call(rbind, lapply(c(FALSE, TRUE), function(sc) {
    g <- nodes(lg, sc, tc_draw, N_ACC)
    g$date    <- as.Date((g$tc + TBAR) * 365.25, origin = EPOCH)
    g$acc     <- g$a
    g$spacing <- factor(SPACING[lg + 1], levels = SPACING)
    g$ceiling <- factor(CEILING[sc + 1], levels = CEILING)
    g$rate    <- free_rate(nodes(lg, sc, tc_fit, N_ACC_FIT))
    g
  }))))

rates <- unique(grid[c("spacing", "ceiling", "rate")])
rates$date <- min(dts)
rates$acc  <- 0.99
# the same decline stated two ways: the compound quarterly percentage the
# tables carry, and the annual factor it implies -- (1 - r)^4 is the year's
# cost multiplier, so its reciprocal is the fold reduction
rates$fold <- 1 / (1 - rates$rate / 100)^4
# "×", the multiplication sign, not the letter x: house style for factors
rates$lab  <- sprintf("Cost decline %.1f%%/qtr,
%.1f×/year",
                      rates$rate, rates$fold)

p <- ggplot(grid, aes(date, acc)) +
  geom_step(data = sota_steps, aes(date, acc), direction = "hv",
            colour = PALETTE[length(PALETTE)], linewidth = 0.7,
            inherit.aes = FALSE) +
  geom_point(colour = INK_SECOND, size = 1.1, alpha = 0.9) +
  geom_text(data = rates, aes(date, acc, label = lab), hjust = 0, vjust = 1,
            colour = INK_PRIMARY, size = 4.4, lineheight = 0.95,
            family = FONT, inherit.aes = FALSE) +
  facet_grid(ceiling ~ spacing) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  # no title: the slide and the document both caption this plate themselves
  labs(x = "Model release date", y = "Accuracy") +
  # slide sizing: 0.6 puts the ticks at 12 pt on the 10-inch slide, where the
  # facet plates' 0.5 (10 pt) does not project. strip.text.x is CENTRED over
  # each column, overriding the theme's left-aligned default, which was set
  # for the benchmark-name strips of the faceted plates and puts a
  # two-condition header off to one side here.
  frontier_theme(0.6) +
  theme(strip.text.x = element_text(hjust = 0.5))

dir.create(out_path("slides"), showWarnings = FALSE, recursive = TRUE)
FIG <- sprintf("grid_schemes_%s.png", BENCH)
save_png(out_path("slides", FIG), p, width = 10, height = 5.625)
cat("wrote slides/", FIG, "\n", sep = "")

# Figure 7 of the report: the same plate as vector, its title dropped (the
# document captions it) but its in-panel rate labels kept. The y-axis title
# goes too: set horizontal above the axis, as the house style has it, the
# lone word "Accuracy" over a titleless plate reads as a heading, and the
# column strips already say what the axis is.
report_figure(p + labs(y = NULL), 7, height = 5.625)
print(rates[c("ceiling", "spacing", "rate", "fold")], row.names = FALSE)
cat(sprintf("\npublished model-free for %s: %+.2f%%\n", BENCH,
            -pareto_decline_qtr(dd)$pct_qtr))
