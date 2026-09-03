using Pkg; Pkg.activate(".")

using CSV, DataFrames, Statistics, StatsFuns, Dates, GLM, GLMakie, CairoMakie, FixedEffectModels, TensorCore

date2year(t) = (Dates.value(t) - Dates.value(Date(2020,1,1))) / 365.25 + 2020

df = CSV.read("data/cost_truncated_curves.csv", DataFrame)
rename!(df, :cost_per_task_usd=>:cost)

map!(x -> match(r"(.*?)(-maas)?$"                                      , x)[1], df.model, df.model)  # remove any -maas suffixes
map!(x -> match(r"(.*?)(-[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9])?$", x)[1], df.model, df.model)  # remove any release date suffixes

# distill instances of multiple runs for same model-benchmark-effort-budget 
df = combine(groupby(df, [:model, :benchmark, :effort, :budget]), 
                          :acc=>mean=>:acc, :cost=>mean=>:cost, :n_samples=>sum=>:n_samples)

df.lncost = log.(df.cost); replace!(df.lncost, -Inf=>missing)

models = CSV.read("data/Model versions-Grid view.csv", DataFrame)  # https://airtable.com/appDFXXgaG1xLtXGL/tblNoPbI37OaCgVKo/viw42jmOv5ayC2n3M
dropmissing!(models, ["id", "Version release date"])
models.id = lowercase.(models.id)

df.releaseyear = missings(Float64, nrow(df))
for r ∈ eachrow(df)
  I = findall(lowercase(r.model) .== first.(models.id, length(r.model)))
  iszero(length(I)) || (r.releaseyear = date2year.(minimum(models."Version release date"[I])))
end
dropmissing!(df, [:releaseyear, :acc, :lncost])
unique!(df, [:model, :effort, :acc, :cost, :releaseyear, :benchmark])  # duplicate rows arise when the budget exceeds what's fruitful

df.logitacc = logit.(min.(df.acc, 1 .− .5 ./ df.n_samples))  # censor perfect scores to 1-1/2n
deleteat!(df, isinf.(df.logitacc))

# Pareto cost frontier over an accuracy-time grid
cost_frontier(df, accgrid, tgrid) =
  [minimum(r.lncost for r ∈ eachrow(df) if r.logitacc≥a && r.releaseyear≤t; init=Inf) |> (x->isinf(x) ? missing : x) for a∈accgrid, t∈tgrid]

SOTA_frontier(df, tgrid) =
  [maximum(r.acc for r ∈ eachrow(df) if r.releaseyear≤t; init=-Inf) for t∈tgrid]

# return dataframe with cost frontier that ranges between extrema of logitacc and release year in a provided dataframe
function cost_frontier_df(df)
  accgrid = range(extrema(skipmissing(df.logitacc))...; length=100)
  tgrid = range(extrema(df.releaseyear)...; length=40)
  DataFrame(logitacc = vec(accgrid ⊗ ones(size(tgrid,1))), 
            releaseyear = vec(ones(size(accgrid,1)) ⊗ tgrid),  
            lncost = vec(cost_frontier(df, accgrid, tgrid)))
end

gdf = groupby(df, :benchmark)
frontierplots = [cost_frontier_df(subdf) for subdf ∈ gdf]

_b(m,var) = coef(m)[findfirst(coefnames(m).==var)]

GLMakie.activate!()
f = Figure(size=(1500,1000))
for (g,(c,r)) ∈ enumerate(Iterators.product(1:2,1:2))
  _df = frontierplots[g]

  Axis3(f[r,c]; limits=(nothing, (2024, 2026.5), nothing), title=keys(gdf)[g].benchmark, xlabel="logit accuracy", zlabel="ln cost")
  surface!(eachcol(_df)...)

  cost_ols = reg(_df, @formula(lncost ~ logitacc + releaseyear))
  surface!(_df.logitacc, _df.releaseyear, 
            _b(cost_ols,"releaseyear") .* _df.releaseyear .+ _b(cost_ols,"logitacc") .* _df.logitacc .+ _b(cost_ols,"(Intercept)"); colormap=(:viridis,.5))
end
f |> display

cost_ols = reg(frontierplots[1], @formula(lncost ~ logitacc + releaseyear))
# cost_ols = reg(frontierplots[1], @formula(lncost ~ logitacc*logitacc + logitacc & releaseyear + releaseyear*releaseyear))

set_theme!(theme_dark())
f = Figure(size=(750,500))
Axis(f[1,1])
tbounds = collect(extrema(df.releaseyear))
tgrid = range(tbounds...; length=1000)
for a ∈ .01:.05:.99
  lines!(tbounds, _b(cost_ols,"releaseyear") .* tbounds .+ _b(cost_ols,"logitacc") * logit(a) .+ _b(cost_ols,"(Intercept)"); color=a, colorrange=(0,1), colormap=:viridis, linewidth=2)
  lines!(collect(tgrid), vec(cost_frontier(df, [logit(a)], tgrid)); color=a, colorrange=(0,1), colormap=:viridis, linewidth=2)
end
f |> display
