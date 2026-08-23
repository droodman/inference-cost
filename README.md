## Viewers (best accessed via [GitHub Pages rendering of this README](https://droodman.github.io/inference-cost))
* [Plots](https://droodman.github.io/inference-cost/viewers/frontier_plot_viewer.html)
* [Tables](https://droodman.github.io/inference-cost/viewers/regression_table_viewer.html)

## Models
Given a set of time, cost, accuracy triples $\left(t_i, c_i, a_i\right)$, define the Pareto frontier by

$$P_t(c) = \max\limits_{i | c_i \le c, t_i \le t} a_i$$
<br>

* The "Logistic, Pareto points" model fits a logistic to these points, with controls based on time and cost.
* The "string logistic model" finds the lowest logistic surface that stays above the Pareto frontier throughout.
* The "logistic, all points" model simply fits the logistic to all triples.
* The "stochastic frontier" model is also fit to all points. The asymmetric error term is half-normal.
* The "Stochastic frontier (time-dependent inefficiency spread)" differs only in allowing the log of the variance of the half-normal term to depend on time. The idea is to capture the widening spread of results trailing behind the frontier as it progresses.

## Controls sets
* Linear: Time (in years) and log cost.
* Quadratic: Adds time², (log cost)² and time × log cost.
* Box-Cox: Drops linear and quadratic terms in favor of [Box-Cox transforms](https://en.wikipedia.org/wiki/Power_transform#Box%E2%80%93Cox_transformation) of time and cost, with time expressed as years since the release of GPT3 in mid-2020. Each is allowed its own exponent (where 0=log and 1=linear). The product of these two is also included, with the same exponents. 

