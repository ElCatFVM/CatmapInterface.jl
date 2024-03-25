__precompile__()
module CatmapInterface

using PyCall

function __init__()
    @pyinclude(joinpath(@__DIR__, "../data/parameter_data.py"))
end

using Artifacts
using Catalyst
using DelimitedFiles
using LessUnitful
using RuntimeGeneratedFunctions
using DocStringExtensions
using LinearAlgebra
using JSON
RuntimeGeneratedFunctions.init(@__MODULE__)
include("utils.jl")
include("ideal-gas-model.jl")
include("harmonic-model.jl")
include("species.jl")
export AbstractSpecies, GasSpecies, AdsorbateSpecies, SiteSpecies, TStateSpecies, FictiousSpecies
include("interface.jl")
export CatmapParams, parse_catmap_input
include("corrections.jl")
include("reaction_network.jl")
export create_reaction_network, generate_function, liquidize

end