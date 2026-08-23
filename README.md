## Viewers
* [Plots](viewers/frontier_plot_viewer.html)
* [Tables](viewers/regression_table_viewer.html)

## Models
Given a set of time, cost, accuracy triples (t_i, c_i, a_i), define the Pareto frontier by $P_t(c) = \max\limits_{i such that c_i \le c, t_i \le t} a_i$.

* The "Logistic, Pareto points" model fits a logistic to these points with time and cost as controls, along with higher-order terms in some variants.
* The "string logistic model" finds the lowest logistic surface that stays above the Pareto frontier throughout.
* The "logistic, all points" model simply fits the logistic to all triples.
* The "stochastic frontier" model is also fit to all points. The asymmetric error term is half-normal.
* The "Stochastic frontier (time-dependent inefficiency spread)" differs only in allowing the log of the variance of the half-normal term to depend on time. The idea is to capture the widening spread of results trailing behind the frontier as it progresses.

## Controls sets
* Linear: Just time (in years) and log cost.
* Quadratic: Adds time^2, (log cost)^2 and time * log cost.
* Box-Cox: Drops linear and quadratic terms in favor of [Box-Cox transforms](https://en.wikipedia.org/wiki/Power_transform#Box%E2%80%93Cox_transformation) of time and cost, with time expressed as years since the release of GPT3 in mid-2020. Each is allowed its own exponent (where 0=log and 1=linear). The product of these two is also included, with the same exponents. 

