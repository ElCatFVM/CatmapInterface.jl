module test_local_gas_species

using Test: Test, @test, @testset, @test_throws
using CatmapInterface: CatmapInterface, LocalGasSpecies, localgasspecies,
    GasSpecies, AdsorbateSpecies

function test_local_gas_species_struct()
    return @testset "LocalGasSpecies struct construction & validation" begin
        # 1. Default constructor
        sp = LocalGasSpecies(; species_name = "COlocal", formation_energy = 0.0)
        @test sp.species_name == "COlocal"
        @test sp.formation_energy == 0.0
        @test sp.pressure == 1.0
        @test sp.site == "b"
        @test sp.frequencies == Float64[]
        @test sp.parent_gas == "CO_g"
        @test ismissing(sp.henry_const)

        # 2. Custom values
        sp2 = LocalGasSpecies(; species_name = "CO2local", formation_energy = 1.5, pressure = 2.0, site = "diff", frequencies = [1000.0, 2000.0], parent_gas = "CO2_g", henry_const = 0.03)
        @test sp2.species_name == "CO2local"
        @test sp2.formation_energy == 1.5
        @test sp2.pressure == 2.0
        @test sp2.site == "diff"
        @test sp2.frequencies == [1000.0, 2000.0]
        @test sp2.parent_gas == "CO2_g"
        @test sp2.henry_const == 0.03

        # 3. Domain validation errors
        @test_throws DomainError LocalGasSpecies(; species_name = "COlocal", formation_energy = 0.0, pressure = -1.0)
        @test_throws DomainError LocalGasSpecies(; species_name = "COlocal", formation_energy = 0.0, frequencies = [-100.0])
        @test_throws DomainError LocalGasSpecies(; species_name = "COlocal", formation_energy = 0.0, henry_const = -0.5)
    end
end

function test_local_gas_species_helper()
    return @testset "localgasspecies helper function" begin
        sp1 = GasSpecies(; species_name = "CO2", formation_energy = 0.0, pressure = 1.0, frequencies = Float64[], henry_const = missing)
        sp2 = LocalGasSpecies(; species_name = "COlocal", formation_energy = 0.0)

        dict = Dict("CO2_g" => sp1, "COlocal_b" => sp2)
        filtered = localgasspecies(dict)
        @test haskey(filtered, "COlocal_b")
        @test !haskey(filtered, "CO2_g")
        @test length(filtered) == 1
    end
end

function runtests()
    return @testset "LocalGasSpecies Tests" begin
        test_local_gas_species_struct()
        test_local_gas_species_helper()
    end
end

end
