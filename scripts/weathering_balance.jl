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
const ð  = 9.8
#molar mass of CO2 [kg/mole]
const ð = 0.044
#the Earth's mean radius [m]
const ðâ = 6.371e6
#the Earth's surface area [m^2]
const ðâ = 4Ï*ðâ^2

##

#modern insolation [W/m^2]
const Fâ = 1366.0
#temporal scaling of insolation [Gyr]
const tâ = 4.5
#preindustrial temperature [K]
const Tâ = 288.0
#planetary albedo [dimensionless]
const Î±â = 0.3
#OLR temperature scaling [W/m^2/K]
const aâ = 2.0
#OLR CO2 scaling [W/m^2]
const bâ = 5.0
#equilibrium outgoing longwave radiation [W/m^2]
const OLRâ = 239.05
#preindustrial atmospheric CO2 concentration
const fâ = 285.0
#reference atmospheric pressure excluding CO2
const Pâ = 1e5 - 28.5
#default parameter for CO2 ocean-atmosphere partioning [teramole]
const hâ = 2.3269250670587494e20/1e12
#default weathering rate [teramole/yr]
const Wâ = 7.0
#reference total carbon in ocean-atmosphere [teramole]
const Câ = 3.193e18/1e12
#reference precipitation rate [m/s]
const pâ = 1.0/yr
#fractional change in precipitation per K change in temperature
const Ïµ = 0.03
#proportionality between precipitation and runoff
const Î = 0.2

##

#instellation over time [W/m^2]
ð»F(t) = Fâ/(1 + (2/5)*(1 - t/tâ))

#surface area averaged radiation [W/m^2]
ð»S(t=tâ) = ð»F(t)/4

#global temperature [K]
ð»T(f, t=tâ) = Tâ + (1/aâ)*( ð»S(t)*(1 - Î±â) - OLRâ + bâ*log(f/fâ) )

#fraction of carbon in the atmosphere [-]
# Mills, Benjamin, et al. "Timing of Neoproterozoic glaciations linked to transport-limited global weathering." Nature geoscience 4.12 (2011): 861-864.
ð»Ï(C) = 0.78*C/(C + hâ)

#partial pressure of CO2 [Pa]
ð»P(C) = ð»Ï(C)*(C*1e12)*ð*ð /ðâ

#molar concentration of CO2 [ppmv]
function ð»f(C)
    @assert C > 0
    P = ð»P(C)
    return 1e6*P/(P + Pâ)
end

#convert carbon reservior [Tmole] directly to temperature [K]
C2T(C, t=tâ) = ð»T(ð»f(C), t)

#latent heat of vaporization of water [J/m^3]
ð»L(T) = 1.918e9*(T/(T - 33.91))^2

#global precipitation [m/s]
ð»p(T, t=tâ) = max(min(pâ*(1 + Ïµ*(T - Tâ)), (1 - Î±â)*ð»S(t)/ð»L(T)), 1e-12)

#global runoff [m/s]
ð»q(T, t=tâ) = Î*ð»p(T, t)

##

function ð»mac(C=Câ, t=tâ; Î=6.1837709746872e-5, Î²=0.2)
    #COâ concentration, temperature, runoff
    f = ð»f(C)
    T = ð»T(f, t)
    q = ð»q(T, t)
    #weathering rate [mole/second/mÂ²]
    w = mac(q, T, f, 11.1, Tâ, fâ, Î=Î, Î²=Î²)
    #global weathering [teramole/year]
    w*(0.3*ðâ*yr/1e12)
end

function ð»whak(C=Câ, t=tâ; k=0.2287292550091995, Î²=0.0)
    #COâ concentration, temperature, runoff
    f = ð»f(C)
    T = ð»T(f, t)
    q = ð»q(T, t)
    #weathering rate [mole/second/mÂ²]
    w = whak(q, T, f, k, 11.1, Tâ, fâ, Î²)
    #global weathering [teramole/year]
    w*(0.3*ðâ*yr/1e12)
end

function ð»Wâ(ð»W::F, W=Wâ, t=tâ) where {F<:Function}
    Roots.find_zero(
        x -> ð»W(exp10(x), t) - W,
        (-20, 20),
        Roots.Bisection()
    ) |> exp10
end

ð»Tâ(ð»W, W=Wâ, t=tâ) = ð»T(ð»Wâ(ð»W, W, t) |> ð»f, t)

##

fig, axs = plt.subplots(1, 2, figsize=(12,5), sharey=true)

ax = axs[1]
fracs = LinRange(0.2, 5, 50)
ax.semilogx(
    fracs,
    Tâ .- ð»Tâ.(ð»whak, fracs*Wâ),
    label="WHAK",
    linewidth=2.5,
    solid_capstyle="round"
)
fracs = LinRange(0.2, 2.5, 50)
ax.semilogx(
    fracs,
    Tâ .- ð»Tâ.(ð»mac, fracs*Wâ),
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
Îwhak = [Tâ - ð»Tâ(ð»whak, 2Wâ), Tâ - ð»Tâ(ð»whak, Wâ/2)]
Îmac = [Tâ - ð»Tâ(ð»mac, 2Wâ), Tâ - ð»Tâ(ð»mac, Wâ/2)]
ax.bar(
    [0, 1],
    [Îwhak[2] - Îwhak[1], Îmac[2] - Îmac[1]],
    0.5,
    [Îwhak[1], Îmac[1]],
    color=["C0", "C1"]
)
ax.set_xlim(-0.5, 1.5)
yl = max(maximum(Îwhak .|> abs), maximum(Îmac .|> abs))
ax.set_xticks([0, 1])
ax.set_xticklabels(["WHAK", "MAC"])
ax.axhline(0, linestyle=":", linewidth=1, alpha=0.5, zorder=-1)
for x â ([-0.25, 0.25], [0.75, 1.25])
    ax.plot(
        x,
        [0,0],
        linestyle="-",
        color="gray",
        linewidth=1.5,
        zorder=1
    )
end
for (i,(top, bot)) â enumerate([Îwhak, Îmac])
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

ens[!,:WHAK] = ð»Tâ.(ð»whak, ens[!,:whak] .- mean(ens[!,:whak]) .+ 7)
ens[!,:MAC] = ð»Tâ.(ð»mac, ens[!,:mac] .- mean(ens[!,:mac]) .+ 7)
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