__precompile__()
"""
CatmapInterface

CatmapInterface implements the functionality of a subset of the Python package [CatMAP](https://catmap.readthedocs.io) that is used in the Python package [CatINT](https://catint.readthedocs.io).
Whereas [CatINT](https://catint.readthedocs.io) uses an iterative approach for combining the Poisson-Nernst-Planck transport model with a microkinetic model of the surface reactions at the electrode, the output of CatmapInterface can be directly plugged into the functionality of [LiquidElectrolytes](https://j-fu.github.io/LiquidElectrolytes.jl) to solve the coupled system.
"""
module CatmapInterface
    using Artifacts: Artifacts, @artifact_str
    using Catalyst: Catalyst, @parameters, @species, @variables, Equation,
        Num, ODESystem, Reaction, ReactionSystem, complete,
        equations, expand_derivatives, netstoichmat, numreactions,
        parameters, reactionrates, species,
        substitute, substitute_in_deriv,
        unknowns, tosymbol
    using DelimitedFiles: DelimitedFiles, readdlm
    using DocStringExtensions: DocStringExtensions, SIGNATURES, TYPEDEF, TYPEDFIELDS
    using JSON: JSON
    using LessUnitful: LessUnitful, @local_phconstants, @local_unitfactors, @ufac_str
    using LinearAlgebra: LinearAlgebra, I, Symmetric, eigvals
    using ModelingToolkitBase: varmap_to_vars
    using PreallocationTools: DiffCache, get_tmp
    using PyCall: PyCall, @py_str, @pyinclude, keys
    using RuntimeGeneratedFunctions: RuntimeGeneratedFunctions
    using SciMLBase: ODEProblem
    using Symbolics: Symbolics, SymbolicUtils
    using SymbolicIndexingInterface: getname

    function __init__()
        return @pyinclude(joinpath(@__DIR__, "../data/parameter_data.py"))
    end

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
    export create_reaction_network, generate_function, liquidize, paramsidx

end
