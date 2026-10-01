__precompile__()
"""
    CatmapInterface

$(read(joinpath(@__DIR__, "..", "README.md"), String))
"""
module CatmapInterface
    using Artifacts: Artifacts, @artifact_str
    using Catalyst: Catalyst, @parameters, @species, @variables, @independent_variables, Equation,
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
    using SciMLBase: ODEProblem
    using SciMLPublic: @public
    using Symbolics: Symbolics, SymbolicUtils, getmetadata, VariableDefaultValue, hasmetadata
    using SymbolicIndexingInterface: getname

    function __init__()
        return @pyinclude(joinpath(@__DIR__, "../data/parameter_data.py"))
    end

    include("utils.jl")
    include("ideal-gas-model.jl")
    include("harmonic-model.jl")
    include("species.jl")
    export AbstractSpecies, GasSpecies, LocalGasSpecies, AdsorbateSpecies, SiteSpecies, TStateSpecies, FictiousSpecies
    include("interface.jl")
    export CatmapParams, parse_catmap_input
    include("corrections.jl")
    include("reaction_network.jl")
    export create_reaction_network, generate_function, liquidize
    export parameter_cache, parameter_defaults
    @public ReactionTerm, ReactionTermParameterCache
    export unknown_indexes, parameter_indexes, parameter_dict
end
