## Viewers
* [Plots of model fits](https://droodman.github.io/inference-cost/viewers/frontier_plot_viewer.html)
* [Tables of model fits](https://droodman.github.io/inference-cost/viewers/regression_table_viewer.html)
* [Cost record timelines](https://droodman.github.io/inference-cost/viewers/record_timeline_viewer.html)
(Claude has put long notes under the tables, which I need to clean up.)

## Models
Given a set of time, cost, accuracy triples $\left(t_i, c_i, a_i\right)$, define the Pareto accuracy frontier by

$$P_t(c) = \max\limits_{i | c_i \le c, t_i \le t} a_i$$
<br>
and the Pareto cost frontier by

$$C_t(a) = \min\limits_{i | a_i \ge a, t_i \le t} c_i$$

<br>

* The "Logistic, Pareto points" model fits a logistic to these points, with controls based on time and cost.
* The "string logistic model" finds the lowest logistic surface that stays above the Pareto frontier throughout.
* The "logistic, all points" model simply fits the logistic to all triples.
* The "stochastic frontier" model is also fit to all points. The asymmetric error term is half-normal.
* The "Stochastic frontier (time-dependent inefficiency spread)" differs only in allowing the log of the variance of the half-normal term to depend on time. The idea is to capture the widening spread of results trailing behind the frontier as it progresses.

## Controls sets
* Linear: Time (in years) and log cost.
* Quadratic: Adds time², (log cost)² and time × log cost.
* Box-Tidwell: Drops linear and quadratic terms in favor of [Box-Cox transforms](https://en.wikipedia.org/wiki/Power_transform#Box%E2%80%93Cox_transformation) of time and cost, with time expressed as years since the release of GPT3 in mid-2020. Each is allowed its own exponent (where 0=log and 1=linear). The product of these two is also included, with the same exponents. 

## Code
Everything is in `src/`. `run_all.R` sources the output scripts into one process so the heavy fits are computed once and shared.

Driver:
* `run_all.R` — runs every output script below in order, except `plot_grid_schemes.R`.

Libraries (sourced, never run directly):
* `paths.R` — path resolution and `src_source()`, so scripts work from the repo root or `src/`.
* `prepare_data.R` — builds the analysis dataset from the source CSVs.
* `frontier_viz.R` — shared figure machinery: the Epoch AI house style (palette, type, theme, export helpers; `DARK <- TRUE` restores the dark experiment) and the figure builders every plate uses (`frontier_plot()`, `iso_acc_plot()`, `isocost_plot()`, `report_figure()`); also loads the runs and defines the time coordinate.
* `fit_specs.R` — the parametric specification grid: inefficiency (A/B/S) × controls (lin/quad/bc).
* `fractional_frontier.R` — stochastic frontier with a fractional-logit response and half-normal or truncated-normal inefficiency.
* `panel_frontier.R` — the same, with inefficiency as a group-level (model × effort) effect.
* `envelope_frontier.R` — accuracy-direction frontier fits: the Pareto staircase sampled on a grid, logit fitted through it, with and without envelope constraints.
* `cost_frontier.R` — cost-direction duals of those fits, pricing misfit in cost rather than accuracy.
* `boxcox_frontier.R` — the Box-Cox specification and the λ profile searches.
* `fit_store.R` — computes each heavy fit once and hands the same object to every consumer.

Output scripts:
* `plot_accuracy_scatter.R` — raw accuracy-vs-date scatter plate, no fits.
* `pareto_frontiers.R` — the nonparametric staircase figure, the cost history of every all-time accuracy record, and a staircase check table.
* `record_timelines.R` — cost records at 25% and 75% performance, per primary benchmark: an HTML table, one figure with every trace on a single (date, cost) plane, and Figure 2 of the report with its data as CSV.
* `plot_frontiers.R` — stochastic-frontier and plain-logit figures in all three views, plus fit diagnostics on the console.
* `plot_paretologit.R` — Pareto-frontier logit figures, all three views; the construction slides and Figure 4.
* `plot_paretologitenv.R` — the envelope-constrained variant, all three views.
* `plot_cost_frontier.R` — cost-direction dual figures, all three views; construction slides for the Pareto-grid cost fit.
* `plot_surfaces_3d.R` — interactive 3-D surfaces and their heatmap twins for the Pareto-grid pair, including the rate-of-decline surfaces; Figure 6.
* `plot_grid_schemes.R` — four ways of laying the node lattice under SOTA, on one benchmark; a slide and Figure 7. Standalone, not in `run_all.R`.
* `cost_frontier_report.R` — cost-decline rate comparison table.
* `regression_tables.R` — HTML and RTF regression tables.

Each fit is drawn in **three** 2-D views, the three non-redundant ways to read one surface z(cost, date). One is an isochrone plate, holding the date fixed and sweeping cost: accuracy against cost (`<model>_<spec>.png`, `frontier_progression_*`). Two are contour plates, each holding a different variable fixed and running date across: cost against date at fixed accuracy (`isoaccuracy_*`) and accuracy against date at fixed budget (`isocost_*`). Transposing the isochrone plate — cost up, accuracy across — is the same curve objects remapped and shows nothing new, so it is not drawn. The staircases under the two contour plates are the same empirical Pareto frontier sliced the two ways — C_a(t) = min{c : a_i ≥ a, t_i ≤ t} and A_c(t) = max{a_i : c_i ≤ c, t_i ≤ t}.

## Outputs
Everything the scripts write lands in `output/`. Figures are PNGs at the top level, 3-D pages are HTML beside them, and two subfolders hold the report and slide material (`slides/`) and the tables (`tables/`). Every faceted figure lays its benchmarks out two across, primaries first, then the rest alphabetically; figures of pooled fits add a sixth panel, "Pooled primaries (ECI scale)", after the primaries.

### Fitted-model figures: `<view>_<model>_<spec>.png`
The three parts of the name:

* `<view>` — which of the three views (above) the plate shows. No prefix (or `frontier_progression_` for the A/B/S families) is the isochrone plate, accuracy against cost; `isoaccuracy_` is cost against date at fixed accuracy; `isocost_` is accuracy against date at fixed budget; `heatmap_` is the fitted surface seen from directly above, with the empirical frontier drawn over it (heatmaps also take a fourth view, `decline`, the instantaneous rate of cost decline over (accuracy, date)).
* `<model>` — the fit. Accuracy direction: `S` (logistic, all tests), `A` (stochastic frontier, half-normal inefficiency), `B` (stochastic frontier, time-dependent inefficiency spread), `paretologit` (logistic through the Pareto points), `paretologitenv` (the same, envelope-constrained). Cost direction, in the same order: `costols`, `costsfa`, `costsfab`, `costgridols`, `costgridolsenv`. The 3-D pages and heatmaps exist for the Pareto-grid four only (`paretologit`, `paretologitenv`, `costgridols`, `costgridolsenv`).
* `<spec>` — the controls set: `lin`, `quad` or `bc` (Box-Tidwell).

So `isoaccuracy_costgridols_bc.png` is the cost-against-date view of the least-squares Pareto-grid cost fit with Box-Tidwell controls.

### Nonparametric and reference figures
* `pareto_frontier.png` — the empirical staircase P_t(c), one curve per half-year.
* `isoaccuracy_records.png` — the cost history of every all-time accuracy record.
* `record_timelines.png` — the cost-record holders at 25% and 75% accuracy, every trace on one (date, cost) plane, models named.
* `accuracy_scatter.png` — every run, accuracy against release date, one name per model. A reference plate meant to be opened at full size.

### 3-D pages: `surface3d_<model>_<spec>_<view>.html`
Interactive plotly renderings of the empirical surface with the fitted one as a wireframe over it, one scene per benchmark; `<view>` is `frontier`, `isoaccuracy` or `decline`. Static HTML plus client-side WebGL: `lib/` holds the one shared copy of plotly.js, so each page is small and both serve from GitHub Pages as they are. The plot viewer's Rendering control switches between a model's 2-D plate, its 3-D page and its heatmap twin.

### `slides/`
* `Figure <N>.svg` — the plates as the report carries them, numbered by the document, not by any script: vector, with the embedded title and notes stripped because the document supplies both. Figures 2, 4, 6 and 7 are written by `record_timelines.R`, `plot_paretologit.R`, `plot_surfaces_3d.R` and `plot_grid_schemes.R`. Figure 1 is not produced by any script.
* `Figure 2.csv` — the data behind Figure 2, one row per plotted dot. `Figure 1.csv` likewise, hand-made.
* `<plate>_<model>_<spec>_{dots,steps,fit}.png` — construction slides at 16:9 for a deck: the observations alone, the record staircases laid over them, then the fitted surface over both; the third stage is the finished plate restricted to the first four benchmarks. Built for `paretologit_bc` (both the frontier and isocost plates) and `isoaccuracy_costgridols_bc`.
* `grid_schemes_<benchmark>.png` — the node-lattice illustration, one benchmark.

### `tables/`
* `regression_<model>_<spec>.html` / `.rtf` — one regression table per fit, with notes; the model and spec keys are the figures'.
* `rate_comparison.html` — cost-decline rates across the cost-direction fits.
* `record_timelines.html` — the 25% / 75% cost-record timelines, model by model, with equivalent ECI.
* `staircase_check.html` / `.csv` — the model-free decline rate read off the Pareto staircase, per benchmark.

[GitHub repo](https://github.com/droodman/inference-cost)
