##

#===
Almost all of this is taken from:
    https://github.com/markmbaum/random-volcanic-climate/blob/main/src/RandomVolcanicClimate.jl
with some modifications/simplifications, mostly concerning function arguments.
A full discussion of the equations and parameters is in the associated research paper:
    https://arxiv.org/abs/2208.02793
===#

using GEOCLIM
using Roots
using CSV
using DataFrames
using StatsBase

using Seaborn
using PyPlot

true |> pygui
plt.rc("axes", grid=false)
plt.style.use("dark_background")

##

function ylabel(label::String, ax)::Nothing
    ax.set_ylabel(label, ha="right", rotation=0)
    ax.yaxis.set_label_coords(0, 1.02)
    return nothing
end

##

#seconds in a day
const day = 60*60*24
#seconds in a year
const yr = 365.25*day
#surface gravity [m/s^2]
const 𝐠 = 9.8
#molar mass of CO2 [kg/mole]
const 𝛍 = 0.044
#the Earth's mean radius [m]
const 𝐑ₑ = 6.371e6
#the Earth's surface area [m^2]
const 𝐒ₑ = 4π*𝐑ₑ^2

##

#modern insolation [W/m^2]
const F₀ = 1366.0
#temporal scaling of insolation [Gyr]
const t₀ = 4.5
#preindustrial temperature [K]
const T₀ = 288.0
#planetary albedo [dimensionless]
const α₀ = 0.3
#OLR temperature scaling [W/m^2/K]
const a₀ = 2.0
#OLR CO2 scaling [W/m^2]
const b₀ = 5.0
#equilibrium outgoing longwave radiation [W/m^2]
const OLR₀ = 239.05
#preindustrial atmospheric CO2 concentration
const f₀ = 285.0
#reference atmospheric pressure excluding CO2
const P₀ = 1e5 - 28.5
#default parameter for CO2 ocean-atmosphere partioning [teramole]
const h₀ = 2.3269250670587494e20/1e12
#default weathering rate [teramole/yr]
const W₀ = 7.0
#reference total carbon in ocean-atmosphere [teramole]
const C₀ = 3.193e18/1e12
#reference precipitation rate [m/s]
const p₀ = 1.0/yr
#fractional change in precipitation per K change in temperature
const ϵ = 0.03
#proportionality between precipitation and runoff
const Γ = 0.2

##

#instellation over time [W/m^2]
𝒻F(t) = F₀/(1 + (2/5)*(1 - t/t₀))

#surface area averaged radiation [W/m^2]
𝒻S(t=t₀) = 𝒻F(t)/4

#global temperature [K]
𝒻T(f, t=t₀) = T₀ + (1/a₀)*( 𝒻S(t)*(1 - α₀) - OLR₀ + b₀*log(f/f₀) )

#fraction of carbon in the atmosphere [-]
# Mills, Benjamin, et al. "Timing of Neoproterozoic glaciations linked to transport-limited global weathering." Nature geoscience 4.12 (2011): 861-864.
𝒻ϕ(C) = 0.78*C/(C + h₀)

#partial pressure of CO2 [Pa]
𝒻P(C) = 𝒻ϕ(C)*(C*1e12)*𝛍*𝐠/𝐒ₑ

#molar concentration of CO2 [ppmv]
function 𝒻f(C)
    @assert C > 0
    P = 𝒻P(C)
    return 1e6*P/(P + P₀)
end

#convert carbon reservior [Tmole] directly to temperature [K]
C2T(C, t=t₀) = 𝒻T(𝒻f(C), t)

#latent heat of vaporization of water [J/m^3]
𝒻L(T) = 1.918e9*(T/(T - 33.91))^2

#global precipitation [m/s]
𝒻p(T, t=t₀) = max(min(p₀*(1 + ϵ*(T - T₀)), (1 - α₀)*𝒻S(t)/𝒻L(T)), 1e-12)

#global runoff [m/s]
𝒻q(T, t=t₀) = Γ*𝒻p(T, t)

##

function 𝒻mac(C=C₀, t=t₀; Λ=6.1837709746872e-5, β=0.2)
    #CO₂ concentration, temperature, runoff
    f = 𝒻f(C)
    T = 𝒻T(f, t)
    q = 𝒻q(T, t)
    #weathering rate [mole/second/m²]
    w = mac(q, T, f, 11.1, T₀, f₀, Λ=Λ, β=β)
    #global weathering [teramole/year]
    w*(0.3*𝐒ₑ*yr/1e12)
end

function 𝒻whak(C=C₀, t=t₀; k=0.2287292550091995, β=0.0)
    #CO₂ concentration, temperature, runoff
    f = 𝒻f(C)
    T = 𝒻T(f, t)
    q = 𝒻q(T, t)
    #weathering rate [mole/second/m²]
    w = whak(q, T, f, k, 11.1, T₀, f₀, β)
    #global weathering [teramole/year]
    w*(0.3*𝐒ₑ*yr/1e12)
end

function 𝒻Wₑ(𝒻W::F, W=W₀, t=t₀) where {F<:Function}
    Roots.find_zero(
        x -> 𝒻W(exp10(x), t) - W,
        (-20, 20),
        Roots.Bisection()
    ) |> exp10
end

𝒻Tₑ(𝒻W, W=W₀, t=t₀) = 𝒻T(𝒻Wₑ(𝒻W, W, t) |> 𝒻f, t)

##

fig, axs = plt.subplots(1, 2, figsize=(12,5), sharey=true)

ax = axs[1]
fracs = LinRange(0.2, 5, 50)
ax.semilogx(
    fracs,
    T₀ .- 𝒻Tₑ.(𝒻whak, fracs*W₀),
    label="WHAK",
    linewidth=2.5,
    solid_capstyle="round"
)
fracs = LinRange(0.2, 2.5, 50)
ax.semilogx(
    fracs,
    T₀ .- 𝒻Tₑ.(𝒻mac, fracs*W₀),
    label="MAC",
    linewidth=2.5,
    solid_capstyle="round"
)

ax.legend(frameon=false, loc="lower left")
ax.set_xlabel("Weatherability Factor \$(\\mathcal{X})\$", labelpad=12)
ylabel("Temperature\nResponse\n[K]", ax)
ax.set_xticks([0.2, 0.5, 1, 2, 5])
ax.set_xticklabels(["0.2x", "0.5x", "no\nchange", "2x", "5x"])
ax.axhline(0, linestyle=":", linewidth=1, alpha=0.5, zorder=-1)
ax.minorticks_off()

ax = axs[2]
Δwhak = [T₀ - 𝒻Tₑ(𝒻whak, 2W₀), T₀ - 𝒻Tₑ(𝒻whak, W₀/2)]
Δmac = [T₀ - 𝒻Tₑ(𝒻mac, 2W₀), T₀ - 𝒻Tₑ(𝒻mac, W₀/2)]
ax.bar(
    [0, 1],
    [Δwhak[2] - Δwhak[1], Δmac[2] - Δmac[1]],
    0.5,
    [Δwhak[1], Δmac[1]],
    color=["C0", "C1"]
)
ax.set_xlim(-0.5, 1.5)
yl = max(maximum(Δwhak .|> abs), maximum(Δmac .|> abs))
ax.set_xticks([0, 1])
ax.set_xticklabels(["WHAK", "MAC"])
ax.axhline(0, linestyle=":", linewidth=1, alpha=0.5, zorder=-1)
for x ∈ ([-0.25, 0.25], [0.75, 1.25])
    ax.plot(
        x,
        [0,0],
        linestyle="-",
        color="gray",
        linewidth=1.5,
        zorder=1
    )
end
for (i,(top, bot)) ∈ enumerate([Δwhak, Δmac])
    ax.text(i-1, top, round(bot, digits=1), ha="center", va="bottom", color="k", fontsize=12)
    ax.text(i-1, bot*0.98, round(top, digits=1), ha="center", va="top", color="k", fontsize=12)
end
ylabel("Doubling/Halving\nTemperature\nSensitivity\n[K]", ax)

fig.tight_layout()
plt.subplots_adjust(bottom=0.15, top=0.85)
fig.savefig(joinpath("..", "plots", "temperature_response_curves"), dpi=500, pad_inches=1)
#plt.close(fig)

##

ens = combine(
    groupby(
        CSV.File(
            joinpath(
                "..",
                "data",
                "segmented_weathering.csv"
            )
        ) |> DataFrame,
        [:n, :p]
    ),
    #aggregate by summing and converting units
    [:whak, :mac] .=> (x -> sum(x)*365*24*3600/1e12) .=> [:whak, :mac]
)

ens[!,:WHAK] = 𝒻Tₑ.(𝒻whak, ens[!,:whak] .- mean(ens[!,:whak]) .+ 7)
ens[!,:MAC] = 𝒻Tₑ.(𝒻mac, ens[!,:mac] .- mean(ens[!,:mac]) .+ 7)
df = stack(ens, [:WHAK, :MAC])

fig, ax = plt.subplots(1, 1, figsize=(6,5))

Seaborn.kdeplot(
    x=df.value,
    hue=df.variable,
    cut=0,
    bw_adjust=0.75,
    linewidth=2.5,
    ax=ax
)
ax.set_xlabel("Equilibrium Temperature [K]")
ax.set_ylabel(nothing)
ax.set_yticks([])
ax.spines["left"].set_visible(false)
leg = ax.get_legend()
leg.set_frame_on(false)

fig.tight_layout()
fig.subplots_adjust(wspace=0.2)
fig.savefig(joinpath("..", "plots", "temperature_response_distributions"), dpi=500)