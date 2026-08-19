# Regression tables: one per statistical model, in RTF and HTML.
#
# Columns are benchmark x specification (AIME linear, AIME quadratic, Chess
# linear, ...), estimates with standard errors in parentheses beneath. Model
# order and naming follow the plot viewer, which is the version of these names
# that has actually been read by a human:
#
#   S            Logistic, all tests
#   A            Stochastic frontier
#   B            Stochastic frontier (time-dependent inefficiency spread)
#   paretologit  Logistic, Pareto points
#   envelope     Strict logistic envelope
#
# Three things differ across models and the tables must not paper over them:
#
#   * Standard errors. A and B are maxLik fits with robust (sandwich) errors --
#     the Papke-Wooldridge point that the objective is a QUASI-likelihood, so
#     inverse-Hessian errors are wrong. S is a quasibinomial glm, given HC1
#     errors for the same reason. The envelope and the Pareto-frontier logit
#     have NO standard errors at all: the envelope is a constrained optimisation
#     with no likelihood behind it, and the frontier logit is fitted to P_t(c)
#     sampled on a fixed grid, whose nodes are not observations. Inventing
#     errors for either would be worse than leaving the cells empty.
#
#   * The quadratic-vs-linear test. A and B have real likelihoods, so this is a
#     likelihood ratio test on 3 df (quad adds lncost^2, tc^2 and lncost:tc).
#     S's quasibinomial glm reports no logLik, so it gets the Wald analogue on
#     the same three terms, LABELLED as Wald rather than passed off as LR. The
#     envelope and the frontier logit get neither, for the same reason they get
#     no standard errors.
#
#   * Sample. The Pareto-frontier logit's N is the number of grid nodes at
#     which P_t(c) is defined -- nodes cheaper and earlier than every run have
#     no frontier value. Those nodes resample a few dozen staircase corners, so
#     N there measures the sampling resolution, not independent information.

source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
src_source("fit_specs.R")
src_source("envelope_frontier.R")

d <- load_runs(drop_gpt4o_chess = FALSE)
benches <- sort(unique(d$benchmark))

## ---- row layout ---------------------------------------------------------------

# term -> (plain label for RTF, HTML label). Order is the table's row order; a
# term absent from a model is simply skipped, so one layout serves all five.
TERMS <- list(
  list(t = "(Intercept)",      p = "Intercept",          h = "Intercept"),
  list(t = "lncost",           p = "ln cost",            h = "ln cost"),
  list(t = "tc",               p = "time",               h = "time"),
  list(t = "I(lncost^2)",      p = "ln cost^2",          h = "ln cost<sup>2</sup>"),
  list(t = "I(tc^2)",          p = "time^2",             h = "time<sup>2</sup>"),
  list(t = "lncost:tc",        p = "ln cost x time",     h = "ln cost &times; time"),
  list(t = "logsig_(Intercept)", p = "log sigma_u",      h = "log &sigma;<sub>u</sub>"),
  list(t = "logsig_tc",        p = "log sigma_u x time", h = "log &sigma;<sub>u</sub> &times; time")
)

MODELS <- list(
  list(key = "S",           label = "Logistic, all tests"),
  list(key = "A",           label = "Stochastic frontier"),
  list(key = "B",           label = "Stochastic frontier (time-dependent inefficiency spread)"),
  list(key = "paretologit", label = "Logistic, Pareto points"),
  list(key = "envelope",    label = "Strict logistic envelope")
)

## ---- fitting and extraction ---------------------------------------------------

fit_one <- function(key, form, b) {
  if (key == "envelope")
    return(fit_envelope(d[d$benchmark == b, ], formula = form))
  if (key == "paretologit")
    return(fit_pareto_logit(d[d$benchmark == b, ], formula = form))
  fit_family(key, form, d)[[b]]
}

# estimate/SE table with one row per term. maxLik keeps a beta_/logsig_ prefix;
# strip beta_ only, so the sigma terms stay distinguishable from the frontier's.
est_se <- function(fit) {
  # This branch must precede the glm one: the frontier logit IS a glm, but its
  # rows are grid nodes, so any SE the glm machinery reports is a fiction.
  if (inherits(fit, c("envelope_frontier", "pareto_grid_logit"))) {
    cf <- coef(fit)
    return(data.frame(term = names(cf), est = unname(cf),
                      se = NA_real_, stringsAsFactors = FALSE))
  }
  if (inherits(fit, "glm")) {
    m <- lmtest::coeftest(fit, vcov = sandwich::vcovHC(fit, type = "HC1"))
    return(data.frame(term = rownames(m), est = m[, 1], se = m[, 2],
                      stringsAsFactors = FALSE))
  }
  m <- summary_robust(fit)
  data.frame(term = sub("^beta_", "", rownames(m)), est = m[, 1], se = m[, 2],
             stringsAsFactors = FALSE)
}

# For the frontier logit the "sample" is the grid nodes carrying a defined
# frontier value, read off the fit itself; for everything else it is the runs.
n_obs <- function(key, b, fit = NULL) {
  if (key == "paretologit") attr(fit, "n_grid") else sum(d$benchmark == b)
}

## ---- rate of cost decline ------------------------------------------------------
#
# Holding accuracy fixed in z = b0 + b_x*ln c + b_t*tc:
#
#     0 = b_x d(ln c) + b_t d(tc)   =>   d ln c / d tc = -b_t / b_x
#
# so the summary is minus the TIME coefficient over the LN COST one, in log
# dollars per year. Units settle the direction: b_t/b_x is (logits/yr) over
# (logits/log-$) = log-$/yr, while b_x/b_t would be years per log-dollar, which
# is a cost-time tradeoff rather than a rate of decline.
#
# Reported PER QUARTER: tc is in years, so the annual log decline -b_t/b_x is
# divided by 4 before the transform x -> 100*(1 - exp(x)), which turns log
# points into a percentage drop and flips the sign so a genuine decline reads
# positive. (Compounding, not division, relates it to the annual figure
# plot_isoaccuracy.R prints: 1 - drop_yr = (1 - drop_qtr)^4.)
#
# The transformation is applied INSIDE the function differentiated, not to the
# ratio afterwards: the delta method is not invariant to reparametrisation, and
# an SE computed on -b_t/b_x is not the SE of 1 - exp(-b_t/b_x). (Composing by
# hand would give 100*exp(x)*SE_x, which the accompanying test confirms.)
#
# Linear specifications only. With the quadratic's lncost:tc interaction the cost
# slope is b_x + b_xt*tc, so the ratio is no longer a single number for the
# column and reporting one would be a fiction.
#
# The standard error is the delta method. car::deltaMethod() is the usual package
# for this and msm::deltamethod() the other; neither is installed here, so the
# gradient comes from numDeriv (which is) and is checked against the closed form
#     d/db_x (-b_t/b_x) =  b_t/b_x^2
#     d/db_t (-b_t/b_x) = -1/b_x
# in the accompanying test rather than trusted blind.
DECLINE <- function(p) 100 * (1 - exp(-p[["tc"]] / p[["lncost"]] / 4))

cost_decline <- function(fit) {
  b <- coef(fit)
  names(b) <- sub("^beta_", "", names(b))
  need <- c("lncost", "tc")
  if (!all(need %in% names(b))) return(NULL)
  est <- DECLINE(b)

  # the envelope has no covariance matrix, and the frontier logit's would be a
  # fiction built on grid nodes, so both point estimates stand alone
  if (inherits(fit, c("envelope_frontier", "pareto_grid_logit")))
    return(list(est = est, se = NA_real_))

  V <- if (inherits(fit, "glm")) sandwich::vcovHC(fit, type = "HC1") else
    vcov_robust(fit)
  if (is.null(rownames(V))) dimnames(V) <- list(names(coef(fit)), names(coef(fit)))
  dimnames(V) <- list(sub("^beta_", "", rownames(V)),
                      sub("^beta_", "", colnames(V)))

  g <- function(p) { names(p) <- need; DECLINE(p) }
  J <- numDeriv::grad(g, b[need])
  list(est = est, se = sqrt(drop(t(J) %*% V[need, need] %*% J)))
}

# percentage points, so one decimal rather than the coefficients' three
decl_cell <- function(cl, which = "est") {
  if (is.null(cl$decline)) return("")
  if (which == "est") fmt(cl$decline$est, 1)
  else if (is.na(cl$decline$se)) "" else sprintf("(%s)", fmt(cl$decline$se, 1))
}

# Quadratic block test. LR where a likelihood exists, Wald where it does not,
# nothing for the envelope. Returns label + statistic + df + p.
quad_test <- function(key, b, fit_lin, fit_quad) {
  extra <- c("I(lncost^2)", "I(tc^2)", "lncost:tc")
  # no test for the envelope or the frontier logit: no likelihood, and no
  # sampling distribution for a Wald statistic built on grid nodes
  if (key %in% c("envelope", "paretologit")) return(NULL)
  if (key %in% c("A", "B")) {
    df <- sum(activePar(fit_quad)) - sum(activePar(fit_lin))
    stat <- 2 * (as.numeric(logLik(fit_quad)) - as.numeric(logLik(fit_lin)))
    return(list(kind = "LR", stat = max(stat, 0), df = df,
                p = pchisq(max(stat, 0), df, lower.tail = FALSE)))
  }
  V <- sandwich::vcovHC(fit_quad, type = "HC1")
  bq <- coef(fit_quad)
  i <- intersect(extra, names(bq))
  i <- i[!is.na(bq[i])]
  if (!length(i)) return(NULL)
  stat <- drop(t(bq[i]) %*% solve(V[i, i, drop = FALSE]) %*% bq[i])
  list(kind = "Wald", stat = stat, df = length(i),
       p = pchisq(stat, length(i), lower.tail = FALSE))
}

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}
fmt_p <- function(p) {
  if (is.na(p)) "" else if (p < 0.0001) "<0.0001" else formatC(p, format = "f", digits = 4)
}

# Build the whole grid for one model: a list of columns, each carrying its
# estimates, N and (for quad) the test.
build_model <- function(key) {
  cols <- list()
  for (b in benches) {
    fits <- lapply(TIME_FORMS, function(f) fit_one(key, f, b))
    tst  <- quad_test(key, b, fits$lin, fits$quad)
    for (tt in names(TIME_FORMS)) {
      cols[[length(cols) + 1]] <- list(
        bench = b, spec = tt,
        # the benchmark name only; linear vs quadratic is not labelled per column
        # but read off the pattern -- the quadratic column is the one carrying the
        # second-order rows, and it is always the right-hand member of the pair
        head = LABELS[[b]],
        es = est_se(fits[[tt]]), n = n_obs(key, b, fits[[tt]]),
        decline = if (tt == "lin") cost_decline(fits[[tt]]) else NULL,
        test = if (tt == "quad") tst else NULL)
    }
  }
  cols
}

# Which term rows actually appear anywhere in this model
active_terms <- function(cols) {
  present <- unique(unlist(lapply(cols, function(cl) cl$es$term)))
  Filter(function(x) x$t %in% present, TERMS)
}

cell <- function(cl, term, which = "est") {
  r <- cl$es[cl$es$term == term, ]
  if (!nrow(r)) return("")
  if (which == "est") fmt(r$est[1])
  else if (is.na(r$se[1])) "" else sprintf("(%s)", fmt(r$se[1]))
}

# Consecutive columns sharing a benchmark, in order: the spans for the header.
# Built by run-length encoding the column heads rather than assuming two columns
# per benchmark, so adding a third specification later needs no change here.
col_groups <- function(cols) {
  r <- rle(vapply(cols, function(cl) cl$head, character(1)))
  data.frame(head = r$values, span = r$lengths, stringsAsFactors = FALSE)
}

# Which data columns begin a new benchmark group (used to draw the separators
# that carry the pairing now that the sub-labels are gone). The first group needs
# no separator, so it is excluded.
group_starts <- function(cols) {
  g <- col_groups(cols)
  cumsum(c(1, head(g$span, -1)))[-1]
}

# df is folded into the row LABEL when every column shares it -- which they all
# do here, quad adding the same three terms everywhere. That removes "(df 3)"
# from eight cells, and the width of that string was setting the column width for
# the whole table.
test_label <- function(cols, kind) {
  dfs <- unique(unlist(lapply(cols, function(cl) cl$test$df)))
  if (length(dfs) == 1) sprintf("%s vs linear (df %d)", kind[1], dfs) else
    sprintf("%s vs linear", kind[1])
}

test_cell <- function(cl, cols) {
  if (is.null(cl$test)) return("")
  dfs <- unique(unlist(lapply(cols, function(c2) c2$test$df)))
  if (length(dfs) == 1) fmt(cl$test$stat, 2) else
    sprintf("%s (df %d)", fmt(cl$test$stat, 2), cl$test$df)
}

## ---- HTML ---------------------------------------------------------------------

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  gsub("<", "&lt;", x, fixed = TRUE)
}

html_table <- function(key, label, cols) {
  trm <- active_terms(cols)
  has_test <- any(vapply(cols, function(cl) !is.null(cl$test), logical(1)))
  kind <- unique(unlist(lapply(cols, function(cl) cl$test$kind)))
  # the envelope has no standard errors, so it gets no (empty) SE rows either
  has_se <- any(!is.na(unlist(lapply(cols, function(cl) cl$es$se))))

  o <- c('<!DOCTYPE html>', '<html lang="en"><head><meta charset="UTF-8" />',
         sprintf('<title>%s</title>', html_escape(label)),
         '<style>',
         # #fcfcfb throughout: the same surface colour the figures are painted
         # on, and the single background the viewers now use -- a table page
         # inside the viewer's iframe should show no seam against it
         'body{font-family:"Segoe UI",Arial,sans-serif;margin:12px;color:#1d1d1d;background:#fcfcfb}',
         'h1{font-size:1.25em;margin:0 0 4px 0}',
         'p.sub{color:#5e5e5e;margin:0 0 16px 0;font-size:.9em}',
         'table{border-collapse:collapse;background:#fcfcfb;font-size:.86em}',
         'th,td{padding:4px 10px;text-align:right;white-space:nowrap}',
         'th:first-child,td:first-child{text-align:left}',
         'thead th{border-bottom:1px solid #1d1d1d;font-weight:600;vertical-align:bottom}',
         'thead th[colspan]{text-align:center;padding-bottom:3px}',
         # the pairing is carried by these separators now that the columns are
         # not individually labelled linear / quadratic
         'th.grp,td.grp{border-left:1px solid #e6e6e2}',
         'tr.se td{color:#5e5e5e;padding-top:0}',
         'tr.gap td{border-top:1px solid #d9d9d5}',
         'tfoot td{border-top:1px solid #1d1d1d;font-size:.92em;color:#5e5e5e;',
         '  text-align:left;white-space:normal;padding-top:8px;max-width:900px}',
         '</style></head><body>',
         sprintf('<h1>%s</h1>', html_escape(label)),
         sprintf('<p class="sub">%s</p>',
                 if (has_se) "Standard errors in parentheses."
                 else "Point estimates only; see the notes for why no standard errors."),
         '<table><thead><tr><th></th>')
  # One spanning header per benchmark. The two columns beneath are linear and
  # quadratic in that order, unlabelled: the quadratic is identifiable as the one
  # carrying the second-order rows, and dropping the two words per column is what
  # takes the table down to a readable width.
  grp <- col_groups(cols)
  for (i in seq_len(nrow(grp)))
    o <- c(o, sprintf('<th colspan="%d"%s>%s</th>', grp$span[i],
                      if (i > 1) ' class="grp"' else '',
                      html_escape(grp$head[i])))
  o <- c(o, '</tr></thead><tbody>')

  starts <- group_starts(cols)
  cls <- function(j) if (j %in% starts) ' class="grp"' else ''

  for (x in trm) {
    o <- c(o, sprintf('<tr><td>%s</td>', x$h))
    for (j in seq_along(cols))
      o <- c(o, sprintf('<td%s>%s</td>', cls(j), cell(cols[[j]], x$t, "est")))
    o <- c(o, '</tr>')
    if (has_se) {
      o <- c(o, '<tr class="se"><td></td>')
      for (j in seq_along(cols))
        o <- c(o, sprintf('<td%s>%s</td>', cls(j), cell(cols[[j]], x$t, "se")))
      o <- c(o, '</tr>')
    }
  }

  # rate of cost decline, above N; the rule that separated the summary block from
  # the coefficients moves up onto this row
  has_decl <- any(vapply(cols, function(cl) !is.null(cl$decline), logical(1)))
  has_decl_se <- any(vapply(cols, function(cl)
    !is.null(cl$decline) && !is.na(cl$decline$se), logical(1)))
  if (has_decl) {
    o <- c(o, '<tr class="gap"><td>cost drop, %/qtr</td>')
    for (j in seq_along(cols))
      o <- c(o, sprintf('<td%s>%s</td>', cls(j), decl_cell(cols[[j]], "est")))
    o <- c(o, '</tr>')
    if (has_decl_se) {
      o <- c(o, '<tr class="se"><td></td>')
      for (j in seq_along(cols))
        o <- c(o, sprintf('<td%s>%s</td>', cls(j), decl_cell(cols[[j]], "se")))
      o <- c(o, '</tr>')
    }
  }

  o <- c(o, sprintf('<tr%s><td>N</td>', if (has_decl) '' else ' class="gap"'))
  for (j in seq_along(cols))
    o <- c(o, sprintf('<td%s>%d</td>', cls(j), cols[[j]]$n))
  o <- c(o, '</tr>')

  if (has_test) {
    o <- c(o, sprintf('<tr><td>%s, &chi;&sup2;</td>',
                      html_escape(test_label(cols, kind))))
    for (j in seq_along(cols))
      o <- c(o, sprintf('<td%s>%s</td>', cls(j), test_cell(cols[[j]], cols)))
    o <- c(o, '</tr><tr class="se"><td>p</td>')
    # html_escape, not raw: fmt_p can return "<0.0001", and an unescaped "<0"
    # is swallowed by the parser as the start of a tag -- the cell renders empty,
    # which reads as "no test" exactly where the test is most significant.
    for (j in seq_along(cols))
      o <- c(o, sprintf('<td%s>%s</td>', cls(j),
                        if (is.null(cols[[j]]$test)) "" else
                          html_escape(fmt_p(cols[[j]]$test$p))))
    o <- c(o, '</tr>')
  }

  o <- c(o, sprintf('</tbody><tfoot><tr><td colspan="%d">%s</td></tr></tfoot></table>',
                    length(cols) + 1, notes_html(key, kind)),
         '</body></html>')
  o
}

## ---- RTF ----------------------------------------------------------------------

rtf_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\\{", "\\\\{", x)
  gsub("\\}", "\\\\}", x)
}

# One row: `widths` are cumulative twips, `cells` the already-escaped strings.
#
# Every cell must be `\pard\intbl ... \cell`, and the table must be closed with a
# `\pard` after the last `\row`. Omitting `\pard` is not a cosmetic slip: Word
# cannot then tell where the table ends, reports "a table in this document has
# become corrupted", and its recovery path renders the rows a second time -- two
# stacked copies of the same table from a file that contains it once.
#
# `\b` is closed explicitly with `\b0` and `\fs18` restated per cell rather than
# relied on from the document default, because `\pard` resets paragraph but not
# character formatting, so bold would otherwise leak from the header row into
# every row beneath it.
# `merge` marks horizontal spans: "first" on the leading cell of a span, "cont"
# on each cell absorbed into it. A merged span still needs one \cell per column,
# so the \cell / \cellx counts the validator checks stay equal.
# `sep` marks cells that open a new benchmark group and get a left rule.
# `align` is "l", "r" or "c" per cell.
rtf_row <- function(cells, widths, bold = FALSE, top = FALSE, bottom = FALSE,
                    merge = NULL, sep = NULL, align = NULL) {
  n <- length(cells)
  if (is.null(merge)) merge <- rep("", n)
  if (is.null(sep))   sep   <- rep(FALSE, n)
  if (is.null(align)) align <- c("l", rep("r", n - 1))
  defn <- paste0("\\trowd\\trgaph60\\trleft0",
                 paste0(vapply(seq_len(n), function(i) paste0(
                   if (merge[i] == "first") "\\clmgf" else
                     if (merge[i] == "cont") "\\clmrg" else "",
                   if (sep[i]) "\\clbrdrl\\brdrs\\brdrw10" else "",
                   if (top) "\\clbrdrt\\brdrs\\brdrw10" else "",
                   if (bottom) "\\clbrdrb\\brdrs\\brdrw10" else "",
                   sprintf("\\cellx%d", widths[i])), character(1)),
                   collapse = ""))
  body <- vapply(seq_len(n), function(i) {
    sprintf("\\pard\\intbl\\itap1%s\\fs18%s %s\\cell",
            switch(align[i], l = "\\ql", r = "\\qr", c = "\\qc"),
            if (bold) "\\b" else "\\b0", cells[i])
  }, character(1))
  c(defn, paste0(body, collapse = ""), "\\row")
}

rtf_table <- function(key, label, cols) {
  trm <- active_terms(cols)
  has_test <- any(vapply(cols, function(cl) !is.null(cl$test), logical(1)))
  kind <- unique(unlist(lapply(cols, function(cl) cl$test$kind)))
  has_se <- any(!is.na(unlist(lapply(cols, function(cl) cl$es$se))))

  # With the per-column "linear"/"quadratic" words gone, the widest data string
  # is an SE like "(11.102)" and the widest label is the test row, so the columns
  # can be much narrower: 2400 + 8*950 = 10000 twips. That fits PORTRAIT letter
  # (12240 wide less 2*720 margins = 10800 usable), so the tables no longer need a
  # landscape page to be pasted into.
  w1 <- 2400; wc <- 950
  widths <- c(w1, w1 + seq_len(length(cols)) * wc)

  starts <- group_starts(cols)
  # +1 throughout because the row-label column occupies position 1
  sep_data <- c(FALSE, seq_along(cols) %in% starts)

  o <- c("{\\rtf1\\ansi\\ansicpg1252\\deff0",
         "{\\fonttbl{\\f0\\fswiss Calibri;}}",
         "\\paperw12240\\paperh15840\\margl720\\margr720\\margt720\\margb720",
         "\\f0\\fs18",
         sprintf("\\pard\\ql{\\b\\fs24 %s}\\par", rtf_escape(label)),
         sprintf("\\pard\\ql{\\i %s}\\par\\par",
                 if (has_se) "Standard errors in parentheses."
                 else "Point estimates only; see the notes for why no standard errors."))

  # Spanning benchmark header: one label per pair, centred over it. The absorbed
  # columns carry empty text; \clmrg does the joining.
  grp <- col_groups(cols)
  head_cells <- c(""); head_merge <- c(""); head_align <- c("l")
  for (i in seq_len(nrow(grp))) {
    head_cells <- c(head_cells, rtf_escape(grp$head[i]),
                    rep("", grp$span[i] - 1))
    head_merge <- c(head_merge, "first", rep("cont", grp$span[i] - 1))
    head_align <- c(head_align, rep("c", grp$span[i]))
  }
  o <- c(o, rtf_row(head_cells, widths, bold = TRUE, bottom = TRUE,
                    merge = head_merge, sep = sep_data, align = head_align))

  for (x in trm) {
    o <- c(o, rtf_row(c(rtf_escape(x$p),
                        vapply(cols, function(cl) cell(cl, x$t, "est"), character(1))),
                      widths, sep = sep_data))
    if (has_se)
      o <- c(o, rtf_row(c("", vapply(cols, function(cl) cell(cl, x$t, "se"),
                                     character(1))), widths, sep = sep_data))
  }
  has_decl <- any(vapply(cols, function(cl) !is.null(cl$decline), logical(1)))
  has_decl_se <- any(vapply(cols, function(cl)
    !is.null(cl$decline) && !is.na(cl$decline$se), logical(1)))
  if (has_decl) {
    o <- c(o, rtf_row(c("cost drop, %/qtr",
                        vapply(cols, function(cl) decl_cell(cl, "est"), character(1))),
                      widths, top = TRUE, sep = sep_data))
    if (has_decl_se)
      o <- c(o, rtf_row(c("", vapply(cols, function(cl) decl_cell(cl, "se"),
                                     character(1))), widths, sep = sep_data))
  }
  o <- c(o, rtf_row(c("N", vapply(cols, function(cl) as.character(cl$n), character(1))),
                    widths, top = !has_decl, sep = sep_data))
  if (has_test) {
    o <- c(o, rtf_row(c(rtf_escape(sprintf("%s, chi-sq", test_label(cols, kind))),
                        vapply(cols, function(cl) test_cell(cl, cols), character(1))),
                      widths, sep = sep_data))
    o <- c(o, rtf_row(c("p vs linear", vapply(cols, function(cl) if (is.null(cl$test)) "" else
      fmt_p(cl$test$p), character(1))), widths, sep = sep_data))
  }
  # \pard closes the table: without it the notes paragraph is still inside table
  # context and Word treats the whole run as a malformed table.
  o <- c(o, "\\pard\\ql\\par",
         sprintf("\\pard\\ql{\\fs16 %s}\\par", rtf_escape(notes_plain(key, kind))),
         "}")
  o
}

## ---- the notes that keep each table honest -------------------------------------

notes_plain <- function(key, kind) {
  se <- switch(key,
    S = "Robust (HC1) standard errors: the quasibinomial objective is a quasi-likelihood, so inverse-Hessian errors would be wrong.",
    paretologit = paste("No standard errors: the response is the empirical Pareto frontier P_t(c) sampled at the nodes of the",
                        "envelope's fixed cost-date grid, and grid nodes are not observations. N is the nodes at which the",
                        "frontier is defined; they resample a few dozen staircase corners, so N measures resolution, not",
                        "information."),
    A = , B = "Robust (sandwich) standard errors, as the SFA objective is a quasi-likelihood in the Papke-Wooldridge sense.",
    envelope = paste("No standard errors: the envelope is the solution of a constrained optimisation with no likelihood behind it,",
                     "so none are reported rather than invented. N is the runs it must clear."))
  tst <- if (is.null(kind) || !length(kind)) {
    "No quadratic-versus-linear test: there is no likelihood to test and no sampling distribution to appeal to."
  } else if (kind[1] == "LR") {
    "Quadratic adds ln cost^2, time^2 and ln cost x time; the test is a likelihood ratio test of all three jointly."
  } else {
    paste("Quadratic adds ln cost^2, time^2 and ln cost x time. Quasibinomial glms report no logLik, so the joint test",
          "of the three is a Wald test on the robust covariance, not a likelihood ratio test.")
  }
  decl <- paste("cost drop, %/qtr is 100*(1 - exp(-b_time / b_ln cost / 4)): the percentage fall in the cost of a fixed",
                "accuracy level per quarter (time is in years, hence the 4). Linear specification only -- with the",
                "quadratic's ln cost x time term the cost slope moves with date, so no single rate describes the column.",
                if (key %in% c("envelope", "paretologit")) "" else
                  paste("Standard error by the delta method on the same robust covariance, taken on the transformed",
                        "quantity. Being symmetric it can reach past 100% where the estimated drop is near total."))
  paste("Time is measured in years and centered within benchmark, so the intercept is the frontier at each benchmark's own",
        "reference date.", decl, se, tst)
}

notes_html <- function(key, kind) {
  x <- notes_plain(key, kind)
  x <- html_escape(x)
  gsub("ln cost^2", "ln cost<sup>2</sup>", gsub("time^2", "time<sup>2</sup>", x,
       fixed = TRUE), fixed = TRUE)
}

## ---- emit ----------------------------------------------------------------------

dir.create(out_path("tables"), showWarnings = FALSE, recursive = TRUE)
for (m in MODELS) {
  cols <- build_model(m$key)
  writeLines(html_table(m$key, m$label, cols),
             out_path("tables", sprintf("regression_%s.html", m$key)))
  writeLines(rtf_table(m$key, m$label, cols),
             out_path("tables", sprintf("regression_%s.rtf", m$key)))
  cat("wrote regression_", m$key, ".html / .rtf\n", sep = "")
}
