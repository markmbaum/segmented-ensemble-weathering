#===
This module can be used by taking the following steps
 1. navigate to the same directory as the module
 2. activate the module's environment in the package manager with "activate ."
 3. get all the dependencies in the package manager with the "instantiate" command
 4. import the module with
      include("EnsembleTools.jl")
      using .EnsembleTools
If GEOCLIM gets updated, we need to enter "up GEOCLIM" in the package manager
===#

module EnsembleTools

using CSV
using GEOCLIM
using DataFrames
import DataFrames

export Trial
export trialparams, readtrial, readtrials, ensembleframe, weatheringframe, writeframe, CO2groups

struct Trial
    n::Int64 #the trial index or number (same for the different colors of a configuration)
    p::Int64 #configuration color
    L::Int64 #global land fraction
    co2::Int64 #carbon dioxide concentration [ppm]
    clim::Climatology #Climatology object
end

function ncext(fn::AbstractString)::String
    if fn[end-2:end] != ".nc"
        return string(fn * ".nc")
    else
        return string(fn)
    end
end

#determines whether something can be parsed to an integer
function isint(x)::Bool
    try
        parse(Int, x)
        return true
    catch
        return false
    end
end

function trialparams(trial::AbstractString)::NTuple{4,Int64}
    #break trial names on underscores
    @assert occursin('_', trial) "no underscores found in trial name"
    s = split(trial, '_')
    #configuration index
    i = findfirst(isint, s)
    n = parse(Int64, s[i])
    #continential color
    i = findfirst(x->occursin("p-", x), s)
    p = parse(Int64, s[i][end-1:end])
    #total land fraction
    i = findfirst(x->occursin('L', x), s)
    L = parse(Int64, s[i][end-1:end])
    #CO2 concentration
    i = findfirst(x->occursin("ppm", x), s)
    co2 = parse(Int64, replace(s[i], "ppm"=>""))
    return n, p, L, co2
end

function readtrial(dir::String;
                   fnr::String="ROF_landfrac.nc",
                   vr::String="QRUNOFF",
                   nullr::Real=1e36,
                   convr::Real=1e-3,
                   fnT::String="ATM.nc",
                   vT::String="TS",
                   fnf::String="ROF_landfrac.nc",
                   vf::String="landfrac",
                   fnlat::String="",
                   latname::String="lat")::Trial
    #first interpret the trial parameters
    n, p, L, co2 = trialparams(basename(dir))
    #read the climatology fields
    c = Climatology(
        ncext(joinpath(dir, fnr)),
        vr,
        nullr,
        convr,
        ncext(joinpath(dir, fnT)),
        vT,
        ncext(joinpath(dir, fnf)),
        vf,
        fnlat=isempty(fnlat) ? "" : ncext(joinpath(dir, fnlat)),
        latname=latname
    )
    return Trial(n, p, L, co2, c)
end

function readtrials(topdir::String; kwargs...)::Vector{Trial}
    #get subdirectory names
    dirs = joinpath.(topdir, readdir(topdir))
    #read information from all of them
    readtrial.(dirs; kwargs...)
end

function DataFrames.DataFrame(trials::Vector{Trial})::DataFrame
    DataFrame(
        "n" => getfield.(trials, :n),
        "p" => getfield.(trials, :p),
        "L" => getfield.(trials, :L)/100.0,
        "co2" => getfield.(trials, :co2),
        "clim" => getfield.(trials, :clim)
    )
end

ensembleframe(topdir; kwargs...) = readtrials(topdir; kwargs...) |> DataFrame

function weatheringframe(direns::String; kwargs...)::DataFrame
    #read all the results into a dataframe
    df = ensembleframe(direns; kwargs...)
    #weathering rates
    df[!,:godderis] = godderis.(df.clim, 0.18, 48200, 288.15)
    df[!,:whak] = whak.(df.clim, df.co2*1e-6, 0.18, 11.1, 288.15, 285e-6)
    df[!,:mac] = mac.(df.clim, df.co2*1e-6, 11.1, 288.15, 285e-6, Λ=0.0084, L=1.5)
    #whole planet runoff, temperature, etc
    df[!,:q] = meanlandrunoff.(df.clim)
    df[!,:r] = totallandrunoff.(df.clim)
    df[!,:T] = meanlandtemperature.(df.clim)
    df[!,:O] = meanoceandistance.(df.clim)
    #various quantities for tropical slices of different sizes
    for L ∈ Int64[5, 10, 15, 20]
        #string version of the tropical boundary latitude
        s = string(L)
        #land fraction, runoff, and temperature
        df[!,'L'*s] = landfraction.(df.clim, L)
        df[!,'r'*s] = totallandrunoff.(df.clim, L)
        df[!,'q'*s] = meanlandrunoff.(df.clim, L)
        df[!,'T'*s] = meanlandtemperature.(df.clim, L)
        df[!,'O'*s] = meanoceandistance.(df.clim, L)
    end
    #a "center of latitude" value which may be useful for supercontinent configs
    df[!,:meanlat] = meanlandlatitude.(df.clim)
    df[!,:meanabslat] = meanabslandlatitude.(df.clim)
    #return the frame with climatologies
    return df
end

#write frame to file without the climatology column
function writeframe(df::DataFrame, path::String)::Nothing
    CSV.write(path, select(df, Not(:clim)))
    println("file written: $path")
    return nothing
end

#select only instances with multiple CO2 values
function CO2groups(df::DataFrame)::Vector{DataFrame}
    [sl for sl ∈ groupby(df, [:n, :L, :p]) if size(sl,1) > 1]
end

end