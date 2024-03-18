module test_steady_state
using CatmapInterface
using Catalyst
using ModelingToolkit
using DifferentialEquations
using PyCall
using ..Utils
using Test

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
    rn              = create_reaction_network(catmap_params)
    ss_params_iter  = load_test_params(test_params_path, catmap_params)
    ModelInstance(name, path, catmap_params, rn, ss_params_iter)
end

## Model definitions

const model_instances = [
    ModelInstance(; name="Au-model-hbond"   , path=joinpath("..", "data", "Au-model-hbond"  , "catmap_CO2R_template.mkm"), test_params_path=joinpath("..", "data", "Au-model-hbond", "test_params.csv")),
    ModelInstance(; name="Au-model-simple"  , path=joinpath("..", "data", "Au-model-simple" , "catmap_CO2R_template.mkm"), test_params_path=joinpath("..", "data", "Au-model-simple", "test_params.csv")),
    #ModelInstance(; name="Liu-model-simple" , path=joinpath("..", "data", "Liu-model-simple", "catmap_CO2R_template.mkm"), test_params_path=joinpath("..", "data", "Liu-model-simple", "test_params.csv")),
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
    solve(ssprob, solver; maxiters)
end

function test_steady_state_with_catmap(odesys, template_file_path, ss_params; solver=DynamicSS(Rodas5P()), maxiters=1e6, rtol=1.0e-4)
    ssol = ssolve(odesys, ss_params; solver, maxiters)
    catmap_ssol = catmap_ssolve(template_file_path, ss_params)
    @testset "species=$(string(species))" for (species, cov) in catmap_ssol
        @test isapprox(ssol[species], cov; rtol)
    end
end

function runtests()
    @testset "Stationary Solutions" begin
		@testset "Model=$(model_instance.name)" for model_instance in model_instances
            CatmapInterface.conserve_pressures!(model_instance.rn, model_instance.catmap_params)
            odesys = convert(ODESystem, model_instance.rn)
            @testset "$(repr(ss_params.ps))" for ss_params in model_instance.ss_params_iter
                test_steady_state_with_catmap(odesys, model_instance.path, ss_params)
			end
        end
    end
end

end