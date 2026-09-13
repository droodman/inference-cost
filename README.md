## Viewers
* [Plots of model fits](https://droodman.github.io/inference-cost/viewers/frontier_plot_viewer.html)
* [Tables of model fits](https://droodman.github.io/inference-cost/viewers/regression_table_viewer.html)
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
* `run_all.R` — runs every output script in order.

Libraries (sourced, never run directly):
* `paths.R` — path resolution and `src_source()`, so scripts work from the repo root or `src/`.
* `prepare_data.R` — builds the analysis dataset from the source CSVs.
* `frontier_viz.R` — shared figure machinery; also loads the runs and defines the time coordinate.
* `fit_specs.R` — the parametric specification grid: inefficiency (A/B/S) × controls (lin/quad/bc).
* `fractional_frontier.R` — stochastic frontier with a fractional-logit response and half-normal or truncated-normal inefficiency.
* `panel_frontier.R` — the same, with inefficiency as a group-level (model × effort) effect.
* `envelope_frontier.R` — accuracy-direction frontier fits: the Pareto staircase sampled on a grid, logit fitted through it, with and without envelope constraints.
* `cost_frontier.R` — cost-direction duals of those fits, pricing misfit in cost rather than accuracy.
* `boxcox_frontier.R` — the Box-Cox specification and the λ profile searches.
* `fit_store.R` — computes each heavy fit once and hands the same object to every consumer.

Output scripts:
* `plot_accuracy_scatter.R` — raw accuracy-vs-date scatter plate, no fits.
* `pareto_frontiers.R` — nonparametric staircase figure.
* `record_timelines.R` — cost records at 25% and 75% performance, per primary benchmark (HTML table + single figure, dots coloured by equivalent ECI).
* `plot_frontiers.R` — stochastic-frontier and plain-logit figures, both views, plus fit diagnostics.

Each fit is drawn in **three** 2-D views, the three non-redundant ways to read one surface z(cost, date). One is an isochrone plate, holding the date fixed and sweeping cost: accuracy against cost (`<model>_<order>.png`, `frontier_progression_*`). Two are contour plates, each holding a different variable fixed and running date across: cost against date at fixed accuracy (`isoaccuracy_*`) and accuracy against date at fixed budget (`isocost_*`). Transposing the isochrone plate — cost up, accuracy across — is the same curve objects remapped and shows nothing new, so it is not drawn. The dashed staircases under the two contour plates are the same empirical Pareto frontier sliced the two ways — C_a(t) = min{c : a_i ≥ a, t_i ≤ t} and A_c(t) = max{a_i : c_i ≤ c, t_i ≤ t}.
* `plot_paretologit.R` — Pareto-frontier logit figures, both views.
* `plot_paretologitenv.R` — the envelope-constrained variant, both views.
* `plot_cost_frontier.R` — cost-direction dual figures, both views.
* `plot_surfaces_3d.R` — interactive 3-D surfaces for the Pareto-grid pair.
* `cost_frontier_report.R` — cost-decline rate comparison table.
* `regression_tables.R` — HTML and RTF regression tables.

[GitHub repo](https://github.com/droodman/inference-cost)
