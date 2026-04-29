module test_steady_state
using Catalyst: Catalyst, netstoichmat, numreactions, reactionrates, species,
speciesmap, symmap_to_varmap, ode_model
using CatmapInterface: CatmapInterface, CatmapParams, create_reaction_network,
parse_catmap_input
using ModelingToolkit: ModelingToolkit, ODESystem, Symbolics, substitute, get_bindings, complete
using OrdinaryDiffEqRosenbrock: OrdinaryDiffEqRosenbrock, Rodas5P,
SteadyStateProblem, solve
using PyCall: PyCall
using SteadyStateDiffEq: SteadyStateDiffEq, DynamicSS
using Test: Test, @test, @testset

using CatmapInterface

include("Utils.jl")
using .Utils

## Models

struct ModelInstance
    name::String
    path::String
    catmap_params::CatmapParams
    rn::Catalyst.ReactionSystem
    ss_params_iter
end
function ModelInstance(; name, path, test_params_path)
    catmap_params   = parse_catmap_input(path)
    rn              = create_reaction_network(catmap_params; conserve_pressures=true)
    ss_params_iter  = load_test_params(test_params_path, catmap_params)
    ModelInstance(name, path, catmap_params, rn, ss_params_iter)
end

## Model definitions

const model_instances = [
    ModelInstance(; name="Au-model-hbond"   , path=joinpath(@__DIR__, "..", "data", "Au-model-hbond"  , "catmap_CO2R_template.mkm"), test_params_path=joinpath(@__DIR__, "..", "data", "Au-model-hbond", "test_params.csv")),
    ModelInstance(; name="Au-model-simple"  , path=joinpath(@__DIR__, "..", "data", "Au-model-simple" , "catmap_CO2R_template.mkm"), test_params_path=joinpath(@__DIR__, "..", "data", "Au-model-simple", "test_params.csv")),
    #ModelInstance(; name="Liu-model-simple" , path=joinpath(@__DIR__, "..", "data", "Liu-model-simple", "catmap_CO2R_template.mkm"), test_params_path=joinpath(@__DIR__, "..", "data", "Liu-model-simple", "test_params.csv")),
]

## Steady Steate Solver

"""
ssolve(odesys::ModelingToolkit.ODESystem, params::SSParams; solver=DynamicSS(Rodas5P()), maxiters=1e6)

Solve the `odesys` for the steady state with the parameter set in `params`. 
"""
function ssolve(odesys::ModelingToolkit.ODESystem, params; solver=DynamicSS(Rodas5P()), maxiters=1e6)
    (; u0, ps) = params
    ssprob = SteadyStateProblem(
        odesys, 
        symmap_to_varmap(odesys, u0), 
        symmap_to_varmap(odesys, ps)
    )
    solve(ssprob, solver; maxiters, reltol=1.0e-10)
end

function test_steady_state_with_catmap(rn::Catalyst.ReactionSystem, odesys, catmap_params, template_file_path, ss_params; solver=DynamicSS(Rodas5P()), maxiters=1e6, rtol=1.0e-3)
    (; species_list) = catmap_params
    
    ssol = ssolve(odesys, ss_params; solver, maxiters)
    catmap_ssol = catmap_ssolve(template_file_path, ss_params) # note that the types of ssol and catmap_ssol are different!!!

    (; ps)  = ss_params
    nr          = numreactions(rn)
    rrs_sym     = reactionrates(rn)
    stoichmat   = netstoichmat(rn)

    #=
    rrs_num     = zeros(nr)

    for (ir, rr_sym) in enumerate(rrs_sym)
        rrs_num[ir] = substitute(substitute(rr_sym, merge(get_bindings(rn), Dict(symmap_to_varmap(rn, ps)))), Dict(sp => ssol[sp] for sp in species(rn)))
    end

    spmap = let 
        spmap = speciesmap(rn)
        Dict(Symbolics.tosymbol(k; escape=false) => v for (k, v) in spmap)
    end
    =#
    @testset "species=$(string(species))" for (species, θcatmap) in catmap_ssol
        sp = species_list[string(species)]
        θ  = ssol[species] #* sp.n_sites
        @test isapprox(θ, θcatmap; rtol)
        #@test isapprox(sum(stoichmat[spmap[species], :] .* rrs_num), 0.0; atol=1.0e-7)
    end
end

"""
runtests()

For each model instance run a testset over model parameters where the steady state solutions without transport are compared to CatMAP simulations. 
"""
function runtests()
    @testset "Stationary Solutions" begin
	@testset "Model=$(model_instance.name)" for model_instance in model_instances
            odesys = complete(ode_model(model_instance.rn))
            @testset "$(repr(ss_params.ps))" for ss_params in model_instance.ss_params_iter
                test_steady_state_with_catmap(model_instance.rn, odesys, model_instance.catmap_params, model_instance.path, ss_params)
	    end
        end
    end
end
end
