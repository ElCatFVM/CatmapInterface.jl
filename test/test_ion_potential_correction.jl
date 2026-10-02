module test_ion_potential_correction

using CatmapInterface: CatmapInterface, AdsorbateSpecies, TStateSpecies, parse_catmap_input
using Test: Test, @test, @testset, @test_throws

const eV = 1.602176634e-19

const models = [
    ("Au-model-simple", joinpath(@__DIR__, "..", "data", "Au-model-simple", "catmap_CO2R_template.mkm")),
    ("Au-model-hbond", joinpath(@__DIR__, "..", "data", "Au-model-hbond", "catmap_CO2R_template.mkm")),
]

# Copy of the template with extra Python lines appended. The energy table path is made
# absolute so that the copy can live in a temporary directory.
function write_variant(dir, template_path, extra_lines)
    energies_file = joinpath(dirname(abspath(template_path)), "catmap_CO2R_energies.txt")
    path = joinpath(dir, "variant.mkm")
    write(path, read(template_path, String) * "\ninput_file = '$energies_file'\n" * extra_lines * "\n")
    return path
end

function free_energies(catmap_params, interface_params)
    energies = Dict(s => 0.0 for s in keys(catmap_params.species_list))
    CatmapInterface.compute_free_energies!(energies, catmap_params, interface_params)
    return energies
end

function interface_params(catmap_params; ϕ_we, ϕ, σ, local_pH)
    θ = Dict(s => 0.0 for (s, sp) in catmap_params.species_list if isa(sp, AdsorbateSpecies))
    return CatmapInterface.InterfaceParams(; θ, ϕ_we, ϕ, σ, local_pH)
end

# expected shift (in units of ϕ·eV) of the free energy when switching the flag on
function expected_shift(s, sp)
    if s == "H_g"
        return 1.0
    elseif s == "ele_g" || s == "OH_g" || (isa(sp, TStateSpecies) && occursin("ele", sp.species_name))
        return -1.0
    else
        return 0.0
    end
end

function test_shifts(name, template_path)
    return @testset "$name" begin
        mktempdir() do dir
            params_off = parse_catmap_input(write_variant(dir, template_path, "ion_potential_correction = False"))
            params_on = parse_catmap_input(write_variant(dir, template_path, "ion_potential_correction = True"))
            @test !params_off.ion_potential_correction
            @test params_on.ion_potential_correction
            # make sure the species affected by the flag are present in the model
            @test haskey(params_on.species_list, "ele_g")
            @test haskey(params_on.species_list, "OH_g")
            @test any(isa(sp, TStateSpecies) && occursin("ele", sp.species_name) for sp in values(params_on.species_list))

            @testset "ϕ_we=$ϕ_we, ϕ=$ϕ" for (ϕ_we, ϕ) in [(-1.0, -0.3), (-0.5, 0.2), (-1.2, 0.0)]
                G_off = free_energies(params_off, interface_params(params_off; ϕ_we, ϕ, σ = -0.05, local_pH = 7.0))
                G_on = free_energies(params_on, interface_params(params_on; ϕ_we, ϕ, σ = -0.05, local_pH = 7.0))
                @testset "species=$s" for (s, sp) in params_on.species_list
                    @test isapprox((G_on[s] - G_off[s]) / eV, expected_shift(s, sp) * ϕ; atol = 1.0e-10)
                end
            end
        end
    end
end

function runtests()
    return @testset "ion_potential_correction" begin
        for (name, path) in models
            test_shifts(name, path)
        end

        @testset "default and reset between parses" begin
            name, template_path = models[1]
            mktempdir() do dir
                @test !parse_catmap_input(template_path).ion_potential_correction
                @test parse_catmap_input(write_variant(dir, template_path, "ion_potential_correction = True")).ion_potential_correction
                # a flag set by a previously parsed file must not leak into the next one
                @test !parse_catmap_input(template_path).ion_potential_correction
            end
        end

        @testset "invalid value" begin
            name, template_path = models[1]
            mktempdir() do dir
                @test_throws ArgumentError parse_catmap_input(write_variant(dir, template_path, "ion_potential_correction = 'False'"))
            end
        end
    end
end

end
