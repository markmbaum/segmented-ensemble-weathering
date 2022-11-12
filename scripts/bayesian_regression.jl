using Turing
using Distributions
using CSV
using Arrow
using DataFrames
using TableTransforms
using StatsBase

using Plots, StatsPlots
using PyPlot
using Seaborn

true |> pygui
plt.style.use("dark_background")
plt.rc("axes", grid=false)

##

function ylabel(label::String, ax)::Nothing
    ax.set_ylabel(label, ha="right", rotation=0)
    ax.yaxis.set_label_coords(0, 1.02)
    return nothing
end

##

ens = CSV.File(
    joinpath(
        "..",
        "data",
        "segmented_weathering_categorized.csv"
    )
) |> DataFrame

ens = filter(r -> (r.lat <= 20) & (r["Area Range"] == 0), ens)

sort!(ens, :lat)

## regression coefficients from the non-bayesian fit

const β_mle = [-3.97796052e-01, -4.88751669e+02, 9.11326384e+01, 2.72619880e+00]
const b_mle = -3.4161702721927476

## check it looks identical

fig, ax = plt.subplots(1,1)
Seaborn.scatterplot(
    x=ens.lat, 
    y=ens.whak .|> log,
    size=ens.A,
    hue=ens.lat,
    palette="flare",
    edgecolor="k",
    legend=false,
    ax=ax
)
X = ens[:,[:lat,:A]] |> Matrix
X = hcat(X, sqrt.(ens.A), sqrt.(ens.A) .* ens.lat)
for ra ∈ LinRange(sqrt.(extrema(ens.A))..., 5)
    X_ = copy(X)
    X_[:,2] .= ra^2
    X_[:,3] .= ra
    X_[:,4] .= ra .* X_[:,1]
    ax.plot(X_[:,1], X_*β_mle .+ b_mle, color="white", alpha=0.5)
end

##

@model function regression(X, y)
    
    b ~ Normal(b_mle, 1)
    β ~ filldist(Normal(0, 10), 4)
    μ = X*β .+ b

    g ~ Normal(0, 5)
    γ ~ Normal(0, 5)
    ν = view(X,:,1)*γ .+ g
    σ = @. log(1 + exp(ν))
    
    y ~ MvNormal(μ, σ)
end

##

X = ens[:,[:lat,:A]] |> Matrix
mₛ = maximum(X, dims=1)
X ./= mₛ
X = hcat(X, sqrt.(X[:,2]), X[:,1] .* sqrt.(X[:,2]))

y = ens.whak .|> log;

##

model = regression(X, y)
chain = sample(model, NUTS(), 10_000)

##

df = DataFrame(chain)[:,chain.name_map.parameters] .|> Float32
Arrow.write(joinpath("..", "data", "chain.feather"), df)

##

fig, axs = plt.subplots(1, 3, figsize=(12,4), sharey=true)

for (ax,a) ∈ zip(axs, LinRange(extrema(ens.A)..., length(axs)))

    Seaborn.scatterplot(
        x=X[:,1], 
        y=y,
        size=ens.A,
        hue=ens.lat,
        palette="flare",
        edgecolor="k",
        legend=false,
        ax=ax
    )

    X_ = copy(X)
    a /= mₛ[2]
    X_[:,2] .= a
    X_[:,3] .= sqrt(a)
    X_[:,4] .= sqrt(a) .* X_[:,1]
    for row ∈ eachrow(df[1:100,:])
        β = values(row[2:5]) |> collect
        b = row[:b]
        μ = X_*β .+ b
        ax.plot(X_[:,1], μ, color="white", alpha=0.1)
        γ = row[:γ]
        g = row[:g]
        σ = log.(1 .+ exp.(view(X,:,1)*γ .+ g))
        q = quantile.(Normal.(μ,σ), 0.975)
        ax.plot(X_[:,1], q, color="gray", alpha=0.1)
        q = quantile.(Normal.(μ,σ), 0.025)
        ax.plot(X_[:,1], q, color="gray", alpha=0.1)
    end

    ax.set_xticks([])
    ylabel("log WHAK\$_T\$\n[Tmole/yr]", ax)
    ax.set_xlabel("\$\\theta\$")
end

fig.tight_layout()
fig.savefig(joinpath("..", "plots", "bayesian_regressions"), dpi=500)