__precompile__()
module CatmapInterface



using Artifacts: Artifacts, @artifact_str
using Catalyst: Catalyst, @parameters, @species, @variables, Equation,
ModelingToolkit, Num, ODESystem, Reaction, ReactionSystem,
SymbolicUtils, Symbolics,  build_function, complete,
equations, expand_derivatives, netstoichmat, numreactions,
parameters, reactionrates, species, speciesmap,
substitute, unknowns
using DelimitedFiles: DelimitedFiles, readdlm
using DocStringExtensions: DocStringExtensions, SIGNATURES, TYPEDEF, TYPEDFIELDS
using JSON: JSON
using LessUnitful: LessUnitful, @local_phconstants, @local_unitfactors, @ufac_str
using LinearAlgebra: LinearAlgebra, I, Symmetric, convert, eigvals
using PyCall: PyCall, @py_str, @pyinclude, keys
using RuntimeGeneratedFunctions: RuntimeGeneratedFunctions, @RuntimeGeneratedFunction, drop_expr
using Symbolics: setmetadata
using Ploynomials

function __init__()
    @pyinclude(joinpath(@__DIR__, "../data/parameter_data.py"))
end

RuntimeGeneratedFunctions.init(@__MODULE__)
include("utils.jl")
export conserve_pressures!
include("ideal-gas-model.jl")
include("harmonic-model.jl")
include("species.jl")
export AbstractSpecies, GasSpecies, AdsorbateSpecies, SiteSpecies, TStateSpecies, FictiousSpecies
include("interface.jl")
export CatmapParams, parse_catmap_input
include("corrections.jl")
include("reaction_network.jl")
export create_reaction_network, reversible_kinetics_electrochemical, generate_function, liquidize, paramsidx, compose_poly, generate_echem_TS_energies
end
