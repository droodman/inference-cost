using Pkg; Pkg.activate(".")

using CSV, DataFrames, Statistics, StatsFuns, Dates, GLM, GLMakie, CairoMakie, FixedEffectModels, TensorCore

date2year(t) = (Dates.value(t) - Dates.value(Date(2020,1,1))) / 365.25 + 2020

df = CSV.read("data/caisi_curves_all.csv", DataFrame)
rename!(df, :cost_per_task_usd=>:cost)

map!(x -> match(r"(.*?)(-maas)?$"                                      , x)[1], df.model, df.model)  # remove any -maas suffixes
map!(x -> match(r"(.*?)(-[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9])?$", x)[1], df.model, df.model)  # remove any release date suffixes

# distill instances of multiple runs for same model-benchmark-effort-budget 
df = combine(groupby(df, [:model, :benchmark, :effort, :budget]), 
                          :acc=>mean=>:acc, :mean_tokens_used=>mean=>:mean_tokens_used, :cost=>mean=>:cost, :full_run_acc=>mean=>:full_run_acc, :n_samples=>sum=>:n_samples)
df.lncost = log.(df.cost)
replace!(df.lncost, -Inf=>missing)
models = CSV.read("data/Model versions-Grid view.csv", DataFrame)  # https://airtable.com/appDFXXgaG1xLtXGL/tblNoPbI37OaCgVKo/viw42jmOv5ayC2n3M
dropmissing!(models, :id)
models.id = lowercase.(models.id)

df.releaseyear = missings(Float64, nrow(df))
for r ∈ eachrow(df)
  I = findfirst(lowercase(r.model) .== first.(models.id, length(r.model)))
  isnothing(I) || (r.releaseyear = date2year.(models."Version release date"[I]))
end
dropmissing!(df, [:releaseyear, :acc, :lncost])
unique!(df, [:model, :effort, :acc, :cost, :releaseyear, :benchmark])  # duplicate rows arise when increasing the budget elicits the same effort and score

reg(df, @formula(acc ~ releaseyear + lncost))

cost_frontier(df, accgrid, tgrid) =
  [minimum(r.lncost for r ∈ eachrow(df) if r.acc≥a && r.releaseyear≤t; init=Inf) for a∈accgrid, t∈tgrid]

function cost_frontier_df(df)
  accgrid =  range(.01, .99; length=100)
  tgrid = range(extrema(df.releaseyear)...; length=100)
  DataFrame(acc = vec(logit.(accgrid) ⊗ ones(size(tgrid,1))), 
            releaseyear = vec(ones(size(accgrid,1)) ⊗ tgrid),  
            lncost = vec(cost_frontier(df, accgrid, tgrid)))
end

gdf = groupby(df, :benchmark)
frontierplots = [cost_frontier_df(subdf) for subdf ∈ gdf]

f = Figure(size=(1500,1000))
for (g,(c,r)) ∈ enumerate(Iterators.product(1:2,1:2))
  Axis3(f[r,c]; limits=(nothing, (2024, 2026.5), nothing), title=keys(gdf)[g].benchmark, ylabel="logit accuracy", zlabel="ln cost")
  surface!(eachcol(frontierplots[g])...)
end
f |> display

df.logitacc = logit.(df.acc)
allowmissing!(df, :logitacc)
replace!(df.logitacc, -Inf=>missing, Inf=>missing)
cost_ols = reg(groupby(df, :benchmark)[1], @formula(lncost ~ logitacc + releaseyear))
