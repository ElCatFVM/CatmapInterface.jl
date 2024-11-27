module test_free_energies
using CatmapInterface
using Test
include("Utils.jl")
using .Utils

const eV = 1.602176634e-19

struct ModelInstance
    name::String
    path::String
    catmap_params::CatmapParams
    interface_params_iter
end
function ModelInstance(; name, path, test_params_path)
    @assert isfile(path) && isfile(test_params_path)

    catmap_params = parse_catmap_input(path)
    interface_params_iter = load_test_params(test_params_path, catmap_params; only_interface_params=true)
    ModelInstance(name, path, catmap_params, interface_params_iter)
end

## Model definitions

const model_instances = [
    ModelInstance(; name="Au-model-hbond"           , path=joinpath(@__DIR__, "..", "data", "Au-model-hbond"          , "catmap_CO2R_template.mkm"), test_params_path=joinpath(@__DIR__, "..", "data", "Au-model-hbond"         , "test_params.csv")),
    ModelInstance(; name="Au-model-simple"          , path=joinpath(@__DIR__, "..", "data", "Au-model-simple"         , "catmap_CO2R_template.mkm"), test_params_path=joinpath(@__DIR__, "..", "data", "Au-model-simple"        , "test_params.csv")),
    ModelInstance(; name="Liu-model-simple"         , path=joinpath(@__DIR__, "..", "data", "Liu-model-simple"        , "catmap_CO2R_template.mkm"), test_params_path=joinpath(@__DIR__, "..", "data", "Liu-model-simple"       , "test_params.csv")),
    ModelInstance(; name="Liu-model-first-order"    , path=joinpath(@__DIR__, "..", "data", "Liu-model-first-order"   , "catmap_CO2R_template.mkm"), test_params_path=joinpath(@__DIR__, "..", "data", "Liu-model-first-order"  , "test_params.csv")),
]

function test_free_energies_with_catmap(catmap_params, template_file_path, interface_params; rtol=1.0e-5)
    free_energies = Dict([ sp => 0.0	for sp in keys(catmap_params.species_list)	])
    CatmapInterface.compute_free_energies!(free_energies, catmap_params, interface_params)
    catmap_free_energies = Dict([ sp => 0.0	for sp in keys(catmap_params.species_list)	])
    compute_catmap_free_energies!(catmap_free_energies, template_file_path, interface_params)
    @testset "species=$sp" for sp in keys(free_energies)
        @test isapprox(free_energies[sp] / eV, catmap_free_energies[sp]; rtol)
    end
end

function runtests()
    @testset "Free Energies" begin
		@testset "Model=$(model_instance.name)" for model_instance in model_instances
            @testset "$(convert(String, interface_params))" for interface_params in model_instance.interface_params_iter
                test_free_energies_with_catmap(model_instance.catmap_params, model_instance.path, interface_params; rtol=1.0e-5)
			end
        end
    end
end

end
