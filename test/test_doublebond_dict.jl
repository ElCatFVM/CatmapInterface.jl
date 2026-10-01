module test_doublebond_dict

using Test: Test, @test, @testset
using PyCall: PyCall, @pyinclude, @py_str
using CatmapInterface: CatmapInterface

function test_doublebond_dict_loading()
    return @testset "doublebond_dict loading & content" begin
        @pyinclude(joinpath(@__DIR__, "..", "data", "parameter_data.py"))
        doublebond_dict = py"doublebond_dict"
        @test !isnothing(doublebond_dict)
        @test haskey(doublebond_dict, "COOH")
        @test haskey(doublebond_dict, "OH")
        @test haskey(doublebond_dict, "OCCO")
        @test doublebond_dict["COOH"] == 0.25
        @test doublebond_dict["OCCO"] == 0.3
    end
end

function runtests()
    return @testset "doublebond_dict Tests" begin
        test_doublebond_dict_loading()
    end
end

end
