using Images
using GEOCLIM
using UnPack
using StaticArrays
using DataFrames
using CSV
using Random: shuffle!

using PyPlot
true |> pygui
plt.style.use("dark_background")

include("EnsembleTools.jl") 
using .EnsembleTools

##

idx2cart(i::Int, n::Int) = CartesianIndex((i-1) % n + 1, ((i-1) ÷ n) + 1)

idx2cart(I::Vector, n::Int) = map(i -> idx2cart(i,n), I)

const connectivity_matrix = SMatrix{3, 3, Bool, 9}(false, true, false, true, true, true, false, true, false)

function landmasses(mask::BitMatrix)
    #assign unique labels to isolated landmasses
    labels = label_components(mask, connectivity_matrix)
    #handle the longitudinal wrapping
    n, m = labels |> size
    for i ∈ 1:n
        if (labels[i,1] > 0) & (labels[i,m]  > 0)
            #landmasses are in contact across the wrap
            if labels[i,1] != labels[i,m]
                labels[labels .== labels[i,m]] .= labels[i,1]
            end
        end
    end
    #take out label gaps just in case
    for (i,label) ∈ enumerate(labels |> unique |> sort)
        labels[labels .== label] .= i-1
    end
    #indices of each landmass, excluding the first one which is ocean
    map(I -> idx2cart(I, n), component_indices(labels))[2:end]
end

landmasses(c::Climatology) = landmasses(c.f .== 1)

𝓉(x::Int) = x
𝓉(x::AbstractFloat) = round(x, sigdigits=4)

## load GCM results

df = joinpath(
    "..",
    "data",
    "large-ens-convex"
) |> EnsembleTools.ensembleframe

## compute weathering rates on grid cels

CO₂ = 1000
df[!,:whak] = map(c -> whak.(c.r, c.T, CO₂*1e-6, 0.18, 11.1, 288.15, 285e-6) .|> Float32, df.clim)
df[!,:mac] = map(c -> mac.(c.r, c.T, CO₂*1e-6, 11.1, 288.15, 285e-6, Λ=0.0084, L=1.5) .|> Float32, df.clim);

## plot isolated landmasses for a few example maps

fig_s = plt.figure(figsize=(12,7))
fig_w = plt.figure(figsize=(12,7))
nrows = 4
ncols = 4
idx = 1:size(df,1) |> collect
shuffle!(idx)
idx = idx[1:nrows*ncols]
Δϕ = 2π/96
for n ∈ 1:nrows*ncols
    ax_s = fig_s.add_subplot(nrows, ncols, n, projection="mollweide")
    ax_w = fig_w.add_subplot(nrows, ncols, n, projection="mollweide")

    c = df[idx[n],:clim]
    z = fill(NaN, size(c))
    for land ∈ landmasses(c)
        z[land] .= rand()
    end
    ax_s.pcolormesh(
        LinRange(-π, π, 96),
        c.lat * (π/180),
        z,
        cmap="Pastel2",
        vmin=0,
        vmax=2
    )
    ax_s.set_xticks([])
    ax_s.set_yticks([])

    ax_w.pcolormesh(
        LinRange(-π, π, 96),
        df[idx[n],:clim].lat * (π/180),
        df[idx[n],:mac],
        cmap="viridis",
        shading="gouraud"
    )
    ax_w.set_xticks([])
    ax_w.set_yticks([])
end
fig_s.tight_layout()
fig_s.savefig(joinpath("..", "plots", "segmentation.png"), dpi=500)
fig_w.tight_layout()
fig_w.savefig(joinpath("..", "plots", "weathering_maps.png"), dpi=500)

## segment the land masses

df[!,:lands] = map(landmasses, df[!,:clim]);

## compute weathering for individual landmasses

N = map(sum ∘ length, df[!,:lands]) |> sum

res = Dict()
for k ∈ [:n, :p]
    res[k] = zeros(Int, N)
end
for k ∈ [:A, :whak, :meanwhak, :mac, :meanmac, :q, :r, :T, :O, :meanlat, :minlat, :perimeter]
    res[k] = zeros(N)
end

n = 1
for row ∈ eachrow(df)
    for land ∈ landmasses(row[:clim])
        #trial number and consolidation parameter are direct
        res[:n][n] = row[:n]
        res[:p][n] = row[:p]    
        #modify a copy of the climatology to contain the land-sea mask for a single land mass
        c = deepcopy(row[:clim])
        c.mask .= 0
        c.mask[land] .= 1
        c.f .*= c.mask
        #use GEOCLIM functions to compute everything
        res[:A][n] = sum(c.A .* c.f)
        res[:whak][n] = GEOCLIM.landsum(row[:whak], c)
        res[:mac][n] = GEOCLIM.landsum(row[:mac], c)
        res[:meanwhak][n] = GEOCLIM.landmean(row[:whak], c)
        res[:meanmac][n] = GEOCLIM.landmean(row[:mac], c)
        res[:q][n] = GEOCLIM.landmean(c.r, c)
        res[:r][n] = GEOCLIM.landsum(c.r, c)
        res[:T][n] = GEOCLIM.landmean(c.T, c)
        res[:O][n] = meanoceandistance(c)
        res[:meanlat][n] = meanlandlatitude(c)
        res[:minlat][n] = c.lat[sum(c.mask, dims=2) .|> !iszero |> vec] .|> abs |> minimum
        res[:perimeter][n] = perimeter(c)
        #increment counter
        n += 1
    end
end

res = res |> DataFrame

## create a final data frame

CSV.write(
    joinpath(
        "..",
        "data",
        "segmented_weathering.csv"
    ),
    res,
    transform=(col,val)->𝓉(val)
)
