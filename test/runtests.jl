using Test
import ExampleJuggler
using ExampleJuggler: cleanexamples, @testmodules, @testscripts
import ExplicitImports
import CatmapInterface

@testset "ExplicitImports" begin
    @test ExplicitImports.check_no_implicit_imports(CatmapInterface) === nothing
    @test ExplicitImports.check_no_stale_explicit_imports(CatmapInterface, ignore=(:setmetadata,)) === nothing
end

ExampleJuggler.verbose!(true)

function run_tests_from_directory(testdir, prefix)
    @info "Directory $(testdir):"
    examples = filter(ex -> length(ex) >= length(prefix) && ex[1:length(prefix)] == prefix, basename.(readdir(testdir)))
    @info examples
    @testmodules(testdir, examples)
end

function run_all_tests(; run_notebooks = false, notebooksonly = false)
    if !notebooksonly
        @testset "basictest" begin
            run_tests_from_directory(@__DIR__, "test_")
        end
    end
end

run_all_tests()


# const Cgap = 0.2 # in F/m^2
# models = [
#     (;  
#         model                   = "CO₂-Reduction on Au with hbond corrections", 
#         catmap_template_path    = "../data/Au-model-hbond/catmap_CO2R_template.mkm", 
#         params_set              = map(
#             row -> (; zip([:θ                                                                                                                                                               ,:ϕ_we      ,:local_pH  ,:T     ,:ϕ_pzc     ,:ϕ     ,:σ     ], [row; Cgap * (row[2] - row[6] - row[5])])...),
#             eachrow([
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.80
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.72
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.64
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.56
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.48
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.40
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.32
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.16
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.08
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.24
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        +0.00
#             ])
#         )
#     ),
#     (;  
#         model                   = "CO₂-Reduction on Au with simple corrections", 
#         catmap_template_path    = "../data/Au-model-simple/catmap_CO2R_template.mkm", 
#         params_set              = map(
#             row -> (; zip([:θ                                                                                                                                                               ,:ϕ_we      ,:local_pH  ,:T     ,:ϕ_pzc     ,:ϕ     ,:σ     ], [row; Cgap * (row[2] - row[6] - row[5])])...),
#             eachrow([
#                             Dict("CO2_t" => 0.1, "COOH_t" => 0.2, "CO_t" => 0.05)   -0.8        6.8         298.0   0.16        -0.80
#                             Dict("CO2_t" => 0.1, "COOH_t" => 0.2, "CO_t" => 0.05)   -0.8        6.8         298.0   0.16        -0.72
#                             Dict("CO2_t" => 0.1, "COOH_t" => 0.2, "CO_t" => 0.05)   -0.8        6.8         298.0   0.16        -0.64
#                             Dict("CO2_t" => 0.1, "COOH_t" => 0.2, "CO_t" => 0.05)   -0.8        6.8         298.0   0.16        -0.56
#                             Dict("CO2_t" => 0.1, "COOH_t" => 0.2, "CO_t" => 0.05)   -0.8        6.8         298.0   0.16        -0.48
#                             Dict("CO2_t" => 0.1, "COOH_t" => 0.2, "CO_t" => 0.05)   -0.8        6.8         298.0   0.16        -0.40
#                             Dict("CO2_t" => 0.1, "COOH_t" => 0.2, "CO_t" => 0.05)   -0.8        6.8         298.0   0.16        -0.32
#                             Dict("CO2_t" => 0.1, "COOH_t" => 0.2, "CO_t" => 0.05)   -0.8        6.8         298.0   0.16        -0.16
#                             Dict("CO2_t" => 0.1, "COOH_t" => 0.2, "CO_t" => 0.05)   -0.8        6.8         298.0   0.16        -0.08
#                             Dict("CO2_t" => 0.1, "COOH_t" => 0.2, "CO_t" => 0.05)   -0.8        6.8         298.0   0.16        -0.24
#                             Dict("CO2_t" => 0.1, "COOH_t" => 0.2, "CO_t" => 0.05)   -0.8        6.8         298.0   0.16        +0.00
#             ])
#         )
#     ),
#     (;  
#         model                   = "CO₂-Reduction on Cu with simple corrections", 
#         catmap_template_path    = "../data/Liu-model-simple/catmap_CO2R_template.mkm", 
#         params_set              = map(
#             row -> (; zip([:θ                                                                                                                                                               ,:ϕ_we      ,:local_pH  ,:T     ,:ϕ_pzc     ,:ϕ     ,:σ     ], [row; Cgap * (row[2] - row[6] - row[5])])...),
#             eachrow([
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.80
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.72
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.64
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.56
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.48
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.40
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.32
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.16
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.08
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.24
#                             Dict("OCCO_t" => 0.01, "CO2_t" => 0.01, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        +0.00
#             ])
#         )
#     ),
#     (;  
#         model                   = "CO₂-Reduction on Cu with simple corrections and first order adsorbate corrections", 
#         catmap_template_path    = "../data/Liu-model-first-order/catmap_CO2R_template.mkm", 
#         params_set              = map(
#             row -> (; zip([:θ                                                                                                                                                               ,:ϕ_we      ,:local_pH  ,:T     ,:ϕ_pzc     ,:ϕ     ,:σ     ], [row; Cgap * (row[2] - row[6] - row[5])])...),
#             eachrow([
#                             Dict("OCCO_t" => 0.2, "CO2_t" => 0.1, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.80
#                             Dict("OCCO_t" => 0.2, "CO2_t" => 0.1, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.72
#                             Dict("OCCO_t" => 0.2, "CO2_t" => 0.1, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.64
#                             Dict("OCCO_t" => 0.2, "CO2_t" => 0.1, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.56
#                             Dict("OCCO_t" => 0.2, "CO2_t" => 0.1, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.48
#                             Dict("OCCO_t" => 0.2, "CO2_t" => 0.1, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.40
#                             Dict("OCCO_t" => 0.2, "CO2_t" => 0.1, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.32
#                             Dict("OCCO_t" => 0.2, "CO2_t" => 0.1, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.16
#                             Dict("OCCO_t" => 0.2, "CO2_t" => 0.1, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.08
#                             Dict("OCCO_t" => 0.2, "CO2_t" => 0.1, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        -0.24
#                             Dict("OCCO_t" => 0.2, "CO2_t" => 0.1, "H_t" => 0.01, "CHOH_t" => 0.01, "CHO_t" => 0.01, "COOH_t" => 0.01, "CO_t" => 0.01, "CH_t"    => 0.01, "OCCOH_t" => 0.01)   -0.8        6.8         298.0   0.16        +0.00
#             ])
#         )
#     ),
# ]

# @testset "CatmapInterface.jl" begin
#     @testset "$model" for (model, catmap_template_path, params_set) in models
#         @testset "$(params.ϕ_we)" for params in params_set
#             test_free_energies(catmap_template_path, params)
#         end
#     end
# end
