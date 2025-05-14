### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 7db4d4bd-bc8b-4949-a241-537bcea259f0
begin
    using Pkg
    Pkg.activate(@__DIR__)
    using Revise
	using PyCall ## 추가
    using LiquidElectrolytes
    using CatmapInterface
    using VoronoiFVM
    using LessUnitful
    using ExtendableGrids, GridVisualize
    using DelimitedFiles
    using PlutoUI, HypertextLiteral
    using PreallocationTools
    using Latexify
    using Catalyst
    using Printf
    using LinearAlgebra
    using Plots
    using DifferentialEquations
	using GlobalSensitivity## 추가
    if isdefined(Main, :PlutoRunner)
        using CairoMakie
        default_plotter!(CairoMakie)
        CairoMakie.activate!(type = "svg")
    end
end;

# ╔═╡ 1d12956c-9ffb-4f42-8551-903213ec8dda
begin
	bayopt = pyimport("bayes_opt")
	Logger= pyimport("bayes_opt.logger")
	JSONLogger=Logger.JSONLogger
	Eve=pyimport("bayes_opt.event")
	load=pyimport("bayes_opt.util")
end

# ╔═╡ 406d9080-c50f-446c-9d51-35e15eafa98e
math = pyimport("math")

# ╔═╡ 29313e7c-4d7f-448c-80f0-32f1beb4e4bd
md"""
## Setup
"""

# ╔═╡ 74f51da6-9d34-480d-b59a-b37e1cf57714
md"""
### Units
"""

# ╔═╡ 578bd9a2-6ba8-4448-8c9e-4ba458e135e1
@unitfactors mol dm m s K pm μm bar Pa eV μF V cm μA mA Å C;

# ╔═╡ 96d0009a-77bb-4e65-a87d-e85a3660ae26
md"""
### Data
"""

# ╔═╡ 1e188adb-f2fa-4de5-998f-6b8319b37d8f
begin
    @phconstants N_A c_0 k_B e h
    const F = N_A * e

    # geometrical constants
    const hmin = 1.0e-6 * μm
    const hmax = 1.0e-1 * μm
    const nref = 0
    const L = 80.0 * μm  ## 80이었음
    const Γ_we = 1
    const Γ_bulk = 2

    # bulk constants
    const pH = 6.8
    const T = 298.0 * K
    const Hcp_CO = 9.7e-6 * mol / (m^3 * Pa)
    const Hcp_CO2 = 3.3e-4 * mol / (m^3 * Pa)
    # species not involved in reactions
    const ikplus = 1
    # species involved in buffer reactions but not in surface reactions
    const ibufferstart = 2
    const ihplus = 2
    #const ihco3 	= 3
    #const ico3 		= 4
    # species involved in buffer reactions and surface reactions
    const isurfacestart = 3
    const ico2 = 3
    const iohminus = 4
    const ibufferend = 4
    # species involved in surface reactions but not in buffer reactions
    const ico = 5
    const ih2 = 6
    const ich4 = 7
    const ich3ch2oh = 8
    const nc = 8
    ## reaction rate constants for bulk reactions
    ### CO2 + OH- <=> HCO3-
    const kbe1 = 4.44e7 / (mol / dm^3)
    const kbf1 = 5.93e3 / (mol / dm^3) / s
    const kbr1 = kbf1 / kbe1
    ### HCO3- + OH- <=> CO3-- + H2O
    const kbe2 = 4.66e3 / (mol / dm^3)
    const kbf2 = 1.0e8 / (mol / dm^3) / s
    const kbr2 = kbf2 / kbe2
    ### CO2 + H20 <=> HCO3- + H+
    const kae1 = 4.44e-7 * (mol / dm^3)
    const kaf1 = 3.7e-2 / s
    const kar1 = kaf1 / kae1
    ### HCO3- <=> CO3-- + H+ 
    const kae2 = 4.66e-5 / (mol / dm^3)
    const kaf2 = 59.44e3 / (mol / dm^3) / s
    const kar2 = kaf2 / kae2
    ### autoprotolyse
    const kwe = 1.0e-14 * (mol / dm^3)^2
    const kwf = 2.4e-5 * (mol / dm^3) / s
    const kwr = kwf / kwe
    const aH₂O = 1.0 #* mol/dm^3

    const scheme = :μex

    # surface constants
    const S = 0.169061707523246e20 / m^2## 9.64e4 / N_A * 0.169061707523246e20 * C*mol/m^2
    const C_gap = 20 * μF / cm^2
    const ϕ_pzc = 0.16 * V
    const ico2_t = 9
    const icooh_t = 10
    const ico_t = 11
    const icho_t = 12
    const ichoh_t = 13
    const ich_t = 14
    const iocco_t = 15
    const ioccoh_t = 16
    const ih_t = 17
    const isurfaceend = 17
    const na = 9
    const M0 = 18.0153 * ufac"g/mol"
    const v0 = N_A * (8.2 * Å)^3 # 1 / (55.4 * ufac"M")

    const species_dict = Dict(
        "K⁺" => ikplus,
        #"HCO₃⁻" => ihco3,
        #"CO₃²⁻" => ico3,
        "CO₂" => ico2,
        "OH⁻" => iohminus,
        "H⁺" => ihplus,
        "CO" => ico,
        "CO_t" => ico_t,
        "COOH_t" => icooh_t,
        "CO2_t" => ico2_t,
        "CHO_t" => icho_t,
        "CHOH_t" => ichoh_t,
        "H_t" => ih_t,
        "H2" => ih2,
        "CH_t" => ich_t,
        "CH4" => ich4,
        "OCCO_t" => iocco_t,
        "OCCOH_t" => ioccoh_t,
        "CH3CH2OH" => ich3ch2oh,
    )

    const species_dict_catmap = Dict(
        "H_t" => ih_t,
        "H2_aq" => ih2,
        "OH_g" => iohminus,
        "CO2_aq" => ico2,
        "CO_aq" => ico,
        "CO_t" => ico_t,
        "COOH_t" => icooh_t,
        "CO2_t" => ico2_t,
        "CHO_t" => icho_t,
        "CHOH_t" => ichoh_t,
        "CH_t" => ich_t,
        "CH4_aq" => ich4,
        "OCCO_t" => iocco_t,
        "OCCOH_t" => ioccoh_t,
        "CH3CH2OH_aq" => ich3ch2oh,
    )
end;

# ╔═╡ 893faaec-7039-40ce-ac74-ab833726ac14
begin
    @kwdef struct BulkSpecies
        name::String
        z::Int
        D::Float64
        c_bulk::Union{Nothing,Float64}
        κ::Float64
        a::Float64
        v::Float64
        M::Float64
        color::Symbol
    end
    function BulkSpecies(; name, z, c_bulk = nothing, D, κ = 0.0, a = 0.0, color)
        D *= m^2 / s
        c_bulk = isnothing(c_bulk) ? nothing : c_bulk * mol / dm^3
        a *= Å
        v = N_A * a^3
        M = M0 * v
        BulkSpecies(name, z, D, c_bulk, κ, a, v, M, color)
    end
    function make_eneutral(
        bulk_species::Vector{BulkSpecies};
        name,
        z,
        D,
        κ = 0.0,
        a = 0.0,
        v = N_A * (a * Å)^3,
        M = M0 * v,
        color,
    )
        a *= Å
        c_bulk = -mapreduce(x -> x.c_bulk * x.z, +, bulk_species) / z
        BulkSpecies(name, z, D, c_bulk, κ, a, v, M, color)
    end
    function create_markdown(bulk::Vector{BulkSpecies})
        table = """
      | Name | z | D | c_bulk | a | v | M | κ | color |
      |------|---|---|--------|---|---|---|---|-------|
      """
        for sp in bulk
            (; name, z, D, c_bulk, a, v, M, κ, color) = sp
            table *= @sprintf(
                "| %s | %i | %1.3e | %1.3e | %1.3e | %1.3e | %1.3e | %1.2f | %s |\n",
                name,
                z,
                D,
                c_bulk,
                a,
                v,
                M,
                κ,
                color
            )
        end
        Markdown.parse(table)
    end
end;

# ╔═╡ 080a139a-bafc-4cac-ab16-a202905f179f
begin
    const bulk = let
        bulk = [
            #BulkSpecies(;name="HCO₃⁻", z=-1, D=1.185e-9, c_bulk=0.091, color=:brown),
            #BulkSpecies(;name="CO₃²⁻", z=-2, D=0.923e-9, c_bulk=2.68e-5, color=:violet),
            BulkSpecies(;
                name = "CO₂",
                z = 0,
                D = 1.91e-9,
                c_bulk = 0.033,
                a = 4.0,
                color = :red,
            ),
            BulkSpecies(;
                name = "OH⁻",
                z = -1,
                D = 5.273e-9,
                c_bulk = 10^(pH - 14),
                color = :green,
            ),
            BulkSpecies(;
                name = "H⁺",
                z = 1,
                D = 9.310e-9,
                c_bulk = 10^(-pH),
                color = :gray,
            ),
            BulkSpecies(;
                name = "CO",
                z = 0,
                D = 2.23e-9,
                c_bulk = 0.0,
                a = 4.0,
                color = :blue,
            ),
            BulkSpecies(;
                name = "H2",
                z = 0,
                D = 1.23e-9,
                c_bulk = 0.0,
                a = 4.0,
                color = :darkred,
            ),
            BulkSpecies(;
                name = "CH4",
                z = 0,
                D = 1.88e-9,
                c_bulk = 0.0,
                a = 4.0,
                color = :gold,
            ),
            BulkSpecies(;
                name = "CH3CH2OH",
                z = 0,
                D = 1.23e-9,
                c_bulk = 0.0,
                a = 4.0,
                color = :greenyellow,
            ),
        ]
        push!(
            bulk,
            make_eneutral(bulk; name = "K⁺", z = 1, D = 1.957e-9, a = 8.2, color = :orange),
        )
        sort(bulk, by = x -> species_dict[x.name])
    end
    create_markdown(bulk)
end

# ╔═╡ 92017c84-03a2-4fbb-bc70-4b43031d7f44
md"""
### Solver Control
"""

# ╔═╡ 1a4651c5-f895-4d77-a332-eca09cbb2a7e
md"""
Consider the bicarbonate buffer system in base and acid as well as autoprotolysis of water:
"""

# ╔═╡ 7d1dc93a-71a4-4164-bcf9-6c25db235fa9
md"""
#### Surface Reactions
"""

# ╔═╡ ecea209d-addc-4059-9920-fb9b4709c8f1
begin
    begin
        struct ModelInstance
            name::String
            path::String
            catmap_params::CatmapParams
            rn::Catalyst.ReactionSystem
        end
        function ModelInstance(; name, path)
            catmap_params = parse_catmap_input(path) ## change the path of this template # & template_interaction
            rn = create_reaction_network(catmap_params)
            CatmapInterface.conserve_pressures!(rn, catmap_params)
            ModelInstance(name, path, catmap_params, rn)
        end
    end
    const model_instances = Dict([
        "Cu-model-ideal" =>
            ModelInstance(; name = "Cu-model-ideal", path = "catmap_CO2R_template.mkm"),
        "Cu-model-interaction" => ModelInstance(;
            name = "Cu-model-interaction",
            path = "catmap_CO2R_template_interaction.mkm",
        ),
        #"Liu-model-simple" => ModelInstance(; name="Liu-model-simple" , path=joinpath("..", "data", "Liu-model-simple", "catmap_CO2R_template.mkm")),
    ])
end;

# ╔═╡ d4986ba9-295c-4170-bcda-a3db7d8eed06
md"""
Model: $(@bind model_name PlutoUI.Select(collect(keys(model_instances)), default="Cu-model-ideal")) 
"""

# ╔═╡ 986879c4-fe6b-4a5c-a19e-5dfcc1d1471a
model_instance = model_instances[model_name];

# ╔═╡ 2fd1776f-df92-406c-8dec-42b85ffd741f
standard_C2=[4.200367971642985e-6, 1.3519488975303392e-5, 4.348740865247605e-5, 0.00013982402428373673, 0.00044889453245214165, 0.001434754235272873, 0.0045232549793875355, 0.013680242869588577, 0.03693099954218972, 0.07825784799151869, 0.07910058289101625, 0.07750805494069611, 0.07493403286707202, 0.07089716788949194, 0.06477480428704059, 0.055928690603671116, 0.0440901267791845, 0.03008789586585555]

# ╔═╡ e55a0993-98ef-47f8-bda4-b7fbf25f44cb
standard_log_C2=[-5.376712661732128, -4.869039724030498, -4.361636470501147, -3.854418202571268, -3.3478556842849705, -2.843222484578618, -2.3445489300476243, -1.8639061923830185, -1.432608937954682, -1.106472098912449, -1.1018203161809683, -1.1106531615644522, -1.1253208935224892, -1.1493711131178095, -1.1885938905333644, -1.2523653482071881, -1.3556586523775538, -1.5216081826753467]

# ╔═╡ 50ae1a36-a89d-44b0-a132-1b1d6a37a976
interaction_log_C2=[
0.413935
0.837041
1.19707
1.5007
1.76485
2.00369
2.22593
2.43518
2.62865
2.79146
2.88936
2.77839
2.48463
1.93109
1.27044
0.593806
-0.0857046
-0.76458]

# ╔═╡ aa46b63b-a4e9-4ac0-89b6-5cbd46f4a66e
standard_methane =[0.03517379974439584, 0.1181608265717288, 0.31456423829930885, 0.8329455218334136, 2.2051283212750015, 5.8371372403280235, 15.448204852419577, 40.86631312972923, 108.00875522354842, 284.99655106115296, 749.9948200340624, 1965.3596003060495, 5115.738561374798, 13172.961283356291, 33332.831737103406, 81994.8667651421, 192725.95029055176, 421468.97782896005]

# ╔═╡ aaead4f8-9b22-4072-a0f2-3fb8b9e36b18
standard_log_methane=[-1.4537807133183585, -0.9275264796605184, -0.5022906523169238, -0.07938340236372113, 0.34343386709188883, 0.7661999043803112, 1.1888780198136448, 1.6113654583712522, 2.0334589609606635, 2.4548396043447456, 2.8750582638738442, 3.293442024502979, 3.7089083418555835, 4.119683415391238, 4.5228722100169225, 4.9137866645135375, 5.2849401957618705, 5.624765613946485]

# ╔═╡ f7655e56-9dd7-452c-9d62-3c4474a4c846
interaction_log_CH4=[
-1.82396
-1.2999
-0.878507
-0.461007
-0.0446877
0.370855
0.785861
1.20052
1.61508
2.03011
2.44647
2.86627
3.28105
3.68834
4.09193
4.49418
4.89537
5.29336
]

# ╔═╡ 026598b5-48bd-47f6-8cac-342a7ac5a11a
standard_CO=[21.349335654478413, 21.35763203137659, 21.360663713562012, 21.3635256117202, 21.364318331791324, 21.362637697310788, 21.356634089152084, 21.341260254294678, 21.306562152572912, 21.237032971780145, 21.111170707197008, 20.897575108214312, 20.547643816284104, 19.98650787512449, 19.10405370440713, 17.751693937344164, 15.761340427203553, 13.02022185029321]

# ╔═╡ 413d05d8-4529-4ac0-9416-d7da01a15446
standard_log_CO=[1.3293843652587878, 1.3295530998211895, 1.3296147428630416, 1.3296729256676474, 1.3296890404031851, 1.3296548750729484, 1.3295328067977668, 1.3292200620109118, 1.3285133812668644, 1.327093841317527, 1.3245123175131486, 1.3200958948086785, 1.312762028827505, 1.3007369190280693, 1.2811255303189824, 1.249239801477872, 1.1975931494088072, 1.1146183841966362]

# ╔═╡ c5397573-cb3c-4ad3-baa1-64cc8a41df0a
interaction_log_CO=[
0.959327
0.957348
0.953513
0.948116
0.941607
0.934333
0.92653
0.91838
0.910139
0.902363
0.895923
0.892928
0.884906
0.869398
0.850179
0.829629
0.808025
0.783214
]

# ╔═╡ fed77c8a-d35c-4e67-b112-903c159a8221
const energy=model_instance.catmap_params.species_list

# ╔═╡ 1bde0c68-15c0-4f44-a33d-4bca9c39e66f
U_eq=Dict()

# ╔═╡ f552af14-3a76-4a1c-9e7f-16fa0d34cd3c
for i in model_instance.catmap_params.reactions
	if isnothing(i.tstate)
		nothing
	else
	    j=[pair.first for pair in  i.educts] ## name 
	    k=[pair.second for pair in i.educts] ## coeff
	    l=[pair.first for pair in i.products]
	    m=[pair.second for pair in i.products]
	    U_eq["$(i.tstate.components[1][1])"]= 	(sum(energy["$(j)"].formation_energy*k for (j,k) in zip(j,k))-sum(energy["$(l)"].formation_energy*m for (l,m) in zip(l,m)))/e
		 
    end
end

# ╔═╡ d82b8194-4c18-4bd4-92c4-03fa3c7a0f1d
# ╠═╡ disabled = true
#=╠═╡
begin
	BayesianOpt= bayopt.BayesianOptimization
	domain_red= bayopt.SequentialDomainReductionTransformer ##https://bayesian-optimization.github.io/BayesianOptimization/2.0.3/domain_reduction.html
end
  ╠═╡ =#

# ╔═╡ 511aed8f-95a1-4b01-ae19-b85fb16ea888
#=╠═╡
bounds_transformer = domain_red(minimum_window=0.5)

  ╠═╡ =#

# ╔═╡ fd7bf3d1-0fe4-46ce-8024-e75a6c63dde8
md"""
#### Debugging
"""

# ╔═╡ fc23f850-cd11-45cf-a769-ad4347b2b6cd
@bind params_input PlutoUI.combine() do Child
    input_params = [
        (; name = "ϕ_we", range = -1.85:0.05:-1.0, default = -1.85, unit = "V"),
        (; name = "local_pH", range = 5.0:1.0:9.0, default = 7.0, unit = ""),
        (; name = "aH2O_g", range = 0.9:0.05:1.0, default = 1.0, unit = ""),
		#(; name = "CO2_g", range = 0.0:0.05:1.0, default=1.0, unit ="") ## added
    ]
    for (s, sp) in model_instance.catmap_params.species_list
        if (isa(sp, GasSpecies) && s ≠ "H2O_g")
            push!(
                input_params,
                (; name = "γ$s", range = 0.2:0.2:1.4, default = 1.0, unit = ""),
            )
        end
    end

    params_input = [
        md""" $(name) : $(Child(name, PlutoUI.Slider(range; default=default, show_value=true))) $unit
        """ for (; name, range, default, unit) in input_params
    ]
    md"""
    #### Input Parameters:
    $(params_input)
    """
end

# ╔═╡ 038e4bfb-c117-44ee-8ac3-b6c19d4696a7
@bind pressures_input PlutoUI.combine() do Child
    input_params =
        @NamedTuple{name::String, range::StepRangeLen, default::Float64, unit::String}[]
    for (s, sp) in model_instance.catmap_params.species_list
		if (isa(sp, GasSpecies) && s == "CO2_g") 
		    push!(input_params, (; name = s, range = 0.0:0.1:1, default=1.0, unit = "bar"))
		else
        if (isa(sp, GasSpecies) && s ≠ "H2O_g") || (isa(sp, FictiousSpecies) && s ≠ "ele_g")
            push!(input_params, (; name = s, range = 0.0:1e-4:1, default = 0, unit = "bar"))
        end
    end
	end

    pressures_input = [
        md""" $(name) : $(Child(name, PlutoUI.Slider(range; default=default, show_value=true))) $unit
        """ for (; name, range, default, unit) in input_params
    ] 
    md"""
    #### Input Pressures:
    $(pressures_input)
    """
end


# ╔═╡ 11432324-db74-443f-a027-f9d89a5da45f
const θinit = 1.0e-20 ## 1.0e-10 

# ╔═╡ 355b8c18-c6f0-4b68-beb8-5c0698dd1204
begin 
	const pressures = Dict(pairs(pressures_input))
	const (; rn, catmap_params) = model_instance
	const species_list = catmap_params.species_list
	const params = Dict(pairs(params_input))
	const θ0 = Dict(Symbol(s) => θinit for (s, sp) in species_list if isa(sp, 		          AdsorbateSpecies))
	const rxns=reactions(rn)
	init_points =10
	n_iter =70
	random_state=10
	maxiters=1e6  
	abstol=1e-6 #1e-8
	reltol=1e-4 #1e-6
	best_result = Ref{Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}}}((-Inf, [], [], []))
	CH4= best_result[][2]
	CO= best_result[][3]
	CH3CH2OH= best_result[][4]
end

# ╔═╡ 97266262-a38c-4a8f-ad47-709f04ea4a4d
begin
function loss_function(; E1, E2, E3, E6, β1, β2, β3, β6)
	result = nothing
	t = @elapsed begin
	rn.:ECHΔOHΔele_t= E1*e + energy["CHOH_t"].formation_energy-β1*e*U_eq["CHΔOHΔele_t"]# 0.85
	rn.:ECOOHΔH2OΔele_t= E2*e +energy["COOH_t"].formation_energy+energy["H2O_g"].formation_energy-β2*e*U_eq["COOHΔH2OΔele_t"] # -0.09
	rn.:EH2OΔCOΔele_t= E3*e + energy["CO_t"].formation_energy+energy["H2O_g"].formation_energy-β3*e*U_eq["H2OΔCOΔele_t"]#0.409
	rn.:EOCCOΔH2OΔele_t= E6*e+energy["OCCO_t"].formation_energy+energy["H2O_g"].formation_energy-β6*e*U_eq["OCCOΔH2OΔele_t"] #-0.65
	rn.:βCHΔOHΔele_t=β1 # 0.5
	rn.:βCOOHΔH2OΔele_t=β2 # 0.5
	rn.:βH2OΔCOΔele_t= β3  #0.5
	rn.:βOCCOΔH2OΔele_t = β6 #0.6	
 (Δϕs, currs_CH4, currs_CO, currs_CH3CH2OH) = currentvoltage(rn, catmap_params, pressures, θ0, params)
	result= -(norm(interaction_log_CH4-currs_CH4)+norm(interaction_log_C2-currs_CH3CH2OH)+norm(interaction_log_CO-currs_CO))
	 if result > best_result[][1]
        best_result[] = (result, currs_CH4, currs_CO, currs_CH3CH2OH)
    end
	end
		open("opt_inter_domain_reduction_with_iter=$(n_iter)_init=$(init_points).txt", "a") do io
        println(io, "elapsed: $(round(t, digits=4))s, E1: $(E1), E2: $(E2), E3: $(E3), E6: $(E6), β1: $(β1), β2: $(β2), β3: $(β3), β6: $(β6)")
    end
	return result
	end
	   pbounds= Dict("E1" => (0.0,2.0), "E2" => (0.0,2.0), "E3" => (0.0,2.0), "E6" => (0.0,2.0), "β1" => (0.1,0.9), "β2" => (0.1,0.9), "β3" => (0.1,0.9), "β6" => (0.1,0.9))
end

# ╔═╡ 23c7a8cb-3734-413c-8216-5dc6a12c645b
#=╠═╡
begin
	optimizer = BayesianOpt(
	    f=loss_function,
	    pbounds=pbounds, 
	    random_state=random_state,
		#bounds_transformer=bounds_transformer
	)
end
  ╠═╡ =#

# ╔═╡ d00de1d2-72c3-48fd-b1e1-6ebb41b73e4c
#=╠═╡
begin
	new_optimizer = BayesianOpt(
	    f=loss_function,
	    pbounds=pbounds, 
	    random_state=random_state,
		bounds_transformer=bounds_transformer
	)
end
  ╠═╡ =#

# ╔═╡ 65b2abd0-4365-45b9-b47d-9f3070384cc2
#=╠═╡
load.load_logs(new_optimizer, logs=["log/opt_inter_no_domain_reduction_with_iter=700_init=100_original.log.json"])
  ╠═╡ =#

# ╔═╡ d83489ca-587d-466e-a95c-579240771eb2
#=╠═╡
print(new_optimizer.max)
  ╠═╡ =#

# ╔═╡ 3725f3e0-4a7d-47c0-8130-f940a6b1ee7c
#=╠═╡
optimized=new_optimizer.max["params"]
  ╠═╡ =#

# ╔═╡ 751ce0f4-ff2b-4fea-a9a3-84a7d8b6f4ff
#=╠═╡
loss_function(; β1=optimized["β1"], β2=optimized["β2"], β3=optimized["β3"], β6=optimized["β6"], E1=optimized["E1"], E2=optimized["E2"], E3=optimized["E3"], E6=optimized["E6"]) 
  ╠═╡ =#

# ╔═╡ 12323fda-3cc3-4015-9d1e-4b7747650b28
#=╠═╡
best_params=new_optimizer.max
  ╠═╡ =#

# ╔═╡ 427af424-c026-4a30-95a6-37be0a59514b
#=╠═╡
begin
	logger = Logger.JSONLogger(path="log/opt_inter_domain_reduction_with_iter=$(n_iter)_init=$(init_points).log")
	new_optimizer.subscribe(Eve.Events.OPTIMIZATION_STEP, logger)
end
  ╠═╡ =#

# ╔═╡ 798eec1c-1383-45bb-ae55-b2620ce185fb
#=╠═╡
new_optimizer[:maximize](
		init_points=0,
		n_iter= n_iter
	)
  ╠═╡ =#

# ╔═╡ f816dd6d-2fb7-4c1c-8a73-1fb20ecc1359
#=╠═╡
optimizer[:maximize](
		init_points=init_points,
		n_iter= n_iter
	)
  ╠═╡ =#

# ╔═╡ ee158f3b-ae5d-4129-a1b9-4665c0215084
#=╠═╡
open("opt_inter_domain_reduction_with_iter=$(n_iter)_init=$(init_points).txt", "a") do io
        println(io, "reltol: $(reltol), abstol:$(abstol)")
		println(io, "iter:$(n_iter), init:$(init_points)")
		println(io, "random_state:$(random_state)")
		println(io, "best_result:$(best_params)" )
end
  ╠═╡ =#

# ╔═╡ a43e66b7-28d9-462b-9983-03a72623f3e7
#=╠═╡
open("log/opt_inter_domain_reduction_with_iter=$(n_iter)_init=$(init_points).log.json", "a") do io
        println(io, "reltol: $(reltol), abstol:$(abstol)")
		println(io, "iter:$(n_iter), init:$(init_points)")
		println(io, "random_state:$(random_state)")
		println(io, "best_result:$(best_params)" )
end
  ╠═╡ =#

# ╔═╡ ae47cd84-fff5-485e-8aba-81f8cd13a469
voltage = [
    -1.0,
    -1.05,
    -1.1,
    -1.15,
    -1.2,
    -1.25,
    -1.3,
    -1.35,
    -1.4,
    -1.45,
    -1.5,
    -1.55,
    -1.6,
    -1.65,
    -1.7,
    -1.75,
    -1.8,
    -1.85,
]

# ╔═╡ 548a2d56-77ff-4d75-981b-6d9cc81ea32d
let
    f1 = Figure(resolution = (1200, 400))
    ax4 = Axis(
        f1[1, 1];
        title = "Current-Voltage Curve_catmapInterface_methane",
        xlabel = "Δϕ [V]",
        ylabel = "I #[mA/cm^2]",
		titlesize = 15,
       
        limits = (-1.85, -1.0, -5, 5),
    )
    lines!(ax4, voltage,  best_result[][2], color=:red, linewidth =5)
	lines!(ax4, voltage, interaction_log_CH4, color=:blue, linewidth=3 )
	
	ax5 = Axis(
        f1[1, 2],
        title = "Current-Voltage Curve_catmapInterface_co",
        xlabel = "Δϕ [V]",
        ylabel = "I #[mA/cm^2]",
		titlesize = 15,
		
        limits = (-1.85, -1.0, -5, 5),  # Adjust limits as needed
    )
    lines!(ax5, voltage,  best_result[][3], color = :red, linewidth =5 )
	lines!(ax5, voltage, interaction_log_CO , color = :blue, linewidth =3)

	ax6 = Axis(
        f1[1, 3],
        title = "Current-Voltage Curve_catmapInterface_c2",
        xlabel = "Δϕ [V]",
        ylabel = "I #[mA/cm^2]",
		titlesize = 15,
		
        limits = (-1.85, -1.0, -10, 5),  # Adjust limits as needed
    )
    lines!(ax6, voltage,  best_result[][4], color = :red, linewidth =5 )
	lines!(ax6, voltage, interaction_log_C2 , color = :blue, linewidth =3)
	
	
	f1
	
end

# ╔═╡ bfec4e58-d8a9-40f9-a3e2-43b29b94fecb
begin
	x_value=[]
	y_value=[]
end

# ╔═╡ 36ab1db5-edea-4d28-8251-bbd47654acfd
function get_name(dic)
        a= dic[1]
	    b= keys(a["params"])
	return b
end

# ╔═╡ 8bc53188-6ba4-4c0d-b2be-544d29904822
#=╠═╡
key_params=collect(get_name(new_optimizer[:res])) ## vector of key
  ╠═╡ =#

# ╔═╡ 49a36f61-2036-4ea9-b508-612df810794b
#=╠═╡
begin
	for (i,res) in enumerate(new_optimizer[:res])
		b=i
		a=res["target"]
		x_value=push!(x_value,b)
		y_value=push!(y_value,a)
	end
end
  ╠═╡ =#

# ╔═╡ 2c1a42dd-d980-4a87-963d-367ea9b55f5c
dict_of_params=Dict()

# ╔═╡ fe2a8f6f-cd0a-4dd4-88f2-47f483374df2
#=╠═╡
for i in key_params
	dict_of_params[i]=[]
end
  ╠═╡ =#

# ╔═╡ ef62ab21-e48d-4f56-98bd-1eb5a578d5c0
dict_of_params

# ╔═╡ 32c79f65-5042-4a74-ac54-ee87223d87e3
#=╠═╡
for (j, res) in enumerate(new_optimizer[:res])
	for i in key_params
dict_of_params[i]=push!(dict_of_params[i],res["params"][i])
	end
end
  ╠═╡ =#

# ╔═╡ a71b0077-0f01-471c-a132-2c783126e7e7
dict_of_params

# ╔═╡ 4865c4aa-2db8-46f4-803c-edc51aa9d2be
set_theme!(Theme(
    fontsize = 24,             # sets base font size
    Axis = (
        titlesize = 40,
        xlabelsize = 30,
        ylabelsize = 30,
        xticklabelsize = 16,
        yticklabelsize = 25,
    )
))

# ╔═╡ 33a3b669-475f-48a8-aae5-168df8436bf5
#=╠═╡
let 
	  f_1 = Figure(resolution = (1000, 2000))
    ax = Axis(
        f_1[1, 1];
        title = "target",
        xlabel = "iteration",
        ylabel = "target",
        limits = (0, length(x_value), -1000, 0),
    )
    lines!(ax, x_value, y_value, color=:red)
for (k,i) in enumerate(sort(key_params)[1:4])
	ax = Axis(
		f_1[k+1,1];
       title = "$(i)",
       xlabel = "iteration",
       ylabel = "value", 
	limits = (0, length(x_value), 0, 2),
    )
    lines!(ax, x_value, dict_of_params[i], color=:blue)
	
end
	f_1
end
  ╠═╡ =#

# ╔═╡ 724eba9d-23ea-4681-a6e0-9a6eac94d5b6
#=╠═╡
let 
	  f_2 = Figure(resolution = (1000, 2000))
	for (k,i) in enumerate(sort(key_params)[5:end])
	ax = Axis(
		f_2[k,1];
       title = "$(i)",
       xlabel = "iteration",
       ylabel = "value", 
	limits = (0, length(x_value), 0, 1),
    )
    lines!(ax, x_value, dict_of_params[i], color=:blue)
	
end
	f_2
end
  ╠═╡ =#

# ╔═╡ bdec29c2-d876-4545-9e17-7d8e9c40d61b
function electrontransfer_CO(rn::Catalyst.ReactionSystem, ssol, params)
    rxns = reactions(rn)
    ir_forward = findfirst(rxns) do r
        (; products, substrates) = r
        !isempty(products) &&
            Symbolics.tosymbol.(products; escape = false) == [:CO_g]
    end
    ir_back = findfirst(rxns) do r
        (; products, substrates) = r
        !isempty(products) &&
            Symbolics.tosymbol.(products; escape = false) == [:CO_t]
    end
    curr = 0.0
    if !isnothing(ir_forward)
        curr += substitute(
            substitute(rxns[ir_forward].rate, Dict(sp => ssol[sp] for sp in species(rn))),
            merge(rn.defaults, symmap_to_varmap(rn, params)),
        )
    end
    if !isnothing(ir_back)
        curr -= substitute(
            substitute(rxns[ir_back].rate, Dict(sp => ssol[sp] for sp in species(rn))),
            merge(rn.defaults, symmap_to_varmap(rn, params)),
        )
    end
    
    #println(curr * ph"e" * S * 2 * 0.1) #current_density
    #println( curr / S)
    return log10(abs(curr * ph"e" * S * 2)) ## 10is for A/m^2 to mA/cm^2 but we don't have to there is ufac in the plotting
end

# ╔═╡ e37d77b4-37d6-46cf-bfde-42a673c83cbd
function electrontransfer_CH3CH2OH(rn::Catalyst.ReactionSystem, ssol, params)
    rxns = reactions(rn)
    ir_forward = findfirst(rxns) do r
        (; products, substrates) = r
        !isempty(products) &&
            Symbolics.tosymbol.(products; escape = false) == [:CH3CH2OH_g, :OH_g]
    end
    ir_back = findfirst(rxns) do r
        (; products, substrates) = r
        !isempty(products) &&
            Symbolics.tosymbol.(products; escape = false) == [:OCCOH_t]
    end
    curr = 0.0
    if !isnothing(ir_forward)
        curr += substitute(
            substitute(rxns[ir_forward].rate, Dict(sp => ssol[sp] for sp in species(rn))),
            merge(rn.defaults, symmap_to_varmap(rn, params)),
        )
    end
    if !isnothing(ir_back)
        curr -= substitute(
            substitute(rxns[ir_back].rate, Dict(sp => ssol[sp] for sp in species(rn))),
            merge(rn.defaults, symmap_to_varmap(rn, params)),
        )
    end
  
    #println(curr * ph"e" * S * 12 * 0.1) #current_density
    #println( curr / S)
    return log10(abs(curr * ph"e" * S * 12)) ## 10is for A/m^2 to mA/cm^2 but we don't have to there is ufac in the plotting ## log10(abs("")) 
end

# ╔═╡ d96a98d6-208d-42b1-bc06-5ed71e3f0a8e
function electrontransfer_CH4(rn::Catalyst.ReactionSystem, ssol, params)
    rxns = reactions(rn)
    ir_forward = findfirst(rxns) do r
        (; products, substrates) = r
        !isempty(products) &&
            Symbolics.tosymbol.(products; escape = false) == [:CH4_g, :OH_g]
    end
    ir_back = findfirst(rxns) do r
        (; products, substrates) = r
        !isempty(products) &&
            Symbolics.tosymbol.(products; escape = false) == [:CH_t]
    end
    curr = 0.0
    if !isnothing(ir_forward)
        curr += substitute(
            substitute(rxns[ir_forward].rate, Dict(sp => ssol[sp] for sp in species(rn))),
            merge(rn.defaults, symmap_to_varmap(rn, params)),
        )
    end
    if !isnothing(ir_back)
        curr -= substitute(
            substitute(rxns[ir_back].rate, Dict(sp => ssol[sp] for sp in species(rn))),
           merge(rn.defaults, symmap_to_varmap(rn, params)),
        )
    end
 
    #println(curr * ph"e" * S * 8 * 0.1) #current_density
    #println( curr / S)
    return log10(abs(curr * ph"e" * S * 8)) ## 10is for A/m^2 to mA/cm^2 but we don't have to there is ufac in the plotting
end

# ╔═╡ 053f971e-f217-49da-a515-477e1393ee54
begin
    function currentvoltage(
        rn,
        catmap_params,
        pressures,
        θ0,
        params;
        #solver = DynamicSS(Rodas5P()), solver & iteration is defined in EnsembleThreads
        #maxiters = 1e4, ##1e6에서 바꿈
    )
        nr = numreactions(rn)
        stoichmat = netstoichmat(rn)
        rrs_sym = reactionrates(rn)
        rrs_num = zeros(nr)

        Δϕs = -1.0:-0.05:params[:ϕ_we]
		params[:ϕ]=0 ## dummy variable for setting SteadyStateProblem
        u0 = symmap_to_varmap(rn, merge(pressures, θ0))
        currs_CH4 = zeros(length(Δϕs))
		currs_CO = zeros(length(Δϕs))
		currs_CH3CH2OH = zeros(length(Δϕs))
		ssprob = SteadyStateProblem(rn, u0, symmap_to_varmap(rn, params))
		function prob_func(prob, i, repeat)
			#println("Trajectory $i is running on thread $(Threads.threadid())")
			params_i = copy(params) ## need to checked(gpt)
			params_i[:ϕ]=params_i[:ϕ_we]-Δϕs[i]
        	remake(prob, p = symmap_to_varmap(rn, params_i))
		end
		ensemble_prob = EnsembleProblem(ssprob, prob_func = prob_func)
		ssol = solve(ensemble_prob , DynamicSS(Rodas5P()); maxiters=maxiters , abstol=abstol, reltol=reltol, ensemblealg=EnsembleThreads(), trajectories = length(Δϕs))
            currs_CH4 = [electrontransfer_CH4(rn, ssol[i], params) for i in 1:length(Δϕs)]
			currs_CO =  [electrontransfer_CO(rn, ssol[i], params) for i in 1:length(Δϕs)]
			currs_CH3CH2OH =  [electrontransfer_CH3CH2OH(rn, ssol[i], params) for i in 1:length(Δϕs)]
		
            #push!(coverage_from_interface, ssol[11])  ## 5 
        #println(collect(Δϕs))
        return collect(Δϕs), currs_CH4, currs_CO, currs_CH3CH2OH
      end## log랑 abs 붙임
end

# ╔═╡ 1393048e-4a39-44a4-92e9-27de6add2632
begin
function peak_cases(p)
	result = nothing
	#ps = [E1 => p[1], E2 => p[2], :E3 => p[3], :E6 => p[4], :β1 => p[5], :β2 => p[6], :β3 => p[7], :β6 => p[8]]
	#ps=[:E1, :E2, :E3, :E6, :β1, :β2, :β3, :β6]
	rn.:ECHΔOHΔele_t= p[1]*e + energy["CHOH_t"].formation_energy-p[5]*e*U_eq["CHΔOHΔele_t"]# 0.85
	rn.:ECOOHΔH2OΔele_t= p[2]*e +energy["COOH_t"].formation_energy+energy["H2O_g"].formation_energy-p[6]*e*U_eq["COOHΔH2OΔele_t"] # -0.09
	rn.:EH2OΔCOΔele_t= p[3]*e + energy["CO_t"].formation_energy+energy["H2O_g"].formation_energy-p[7]*e*U_eq["H2OΔCOΔele_t"]#0.409
	rn.:EOCCOΔH2OΔele_t= p[4]*e+energy["OCCO_t"].formation_energy+energy["H2O_g"].formation_energy-p[8]*e*U_eq["OCCOΔH2OΔele_t"] #-0.65
	rn.:ECOΔ_t= p[9]*e+ energy["CO_g"].formation_energy+energy["_t"].formation_energy
	rn.:EOCΔCO_t = p[10]*e + energy["CO_t"].formation_energy*2
	rn.:ECO2Δ_t = p[11]*e +energy["CO_g"].formation_energy + 2* energy["_t"].formation_energy
	rn.:βCHΔOHΔele_t=p[5] # 0.5
	rn.:βCOOHΔH2OΔele_t=p[6] # 0.5
	rn.:βH2OΔCOΔele_t= p[7] #0.5
	rn.:βOCCOΔH2OΔele_t = p[8] #0.6	
 (Δϕs, currs_CH4, currs_CO, currs_CH3CH2OH) = currentvoltage(rn, catmap_params, pressures, θ0, params)
	result= -(norm(interaction_log_CH4-currs_CH4)+norm(interaction_log_C2-currs_CH3CH2OH)+norm(interaction_log_CO-currs_CO))
	return result
	end
end

# ╔═╡ 764462d4-dd65-43e1-ad9f-bb7b79efbdc2
global_sens = gsa(peak_cases, Sobol(), [(0.0,2.0), (0.0,2.0), (0.0,2.0), (0.0,2.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.1, 0.9), (0.1, 0.9), (0.1, 0.9), (0.1, 0.9)]; samples = 100)

# ╔═╡ 7f7e6146-85f4-48c8-b4a5-f0e7d2acf077
begin
	tbl = (cat = [1,2,3,4,5,6,7,8,9,10,11],
	       height= global_sens.ST)
	
	barplot(tbl.cat, tbl.height,
	        axis = (xticks = (1:11, ["E1", "E2", "E3", "E6", "Beta1" , "Beta2", "Beta3", "Beta6", "E4", "E5", "E7"]),
	                title = "Total Order Indices"),
	        )
end

# ╔═╡ 6dc5664c-4a72-4584-a64c-47a1cb3b7517
begin
	tbl2 = (cat = [1,2,3,4,5,6,7,8],
	       height= global_sens.S1)
	
	barplot(tbl2.cat, tbl2.height,
	        axis = (xticks = (1:8, ["E1", "E2", "E3", "E6", "Beta1" , "Beta2", "Beta3", "Beta6"]),
	                title = "First Order Indices"),
	        )
end

# ╔═╡ Cell order:
# ╠═7db4d4bd-bc8b-4949-a241-537bcea259f0
# ╠═1d12956c-9ffb-4f42-8551-903213ec8dda
# ╠═406d9080-c50f-446c-9d51-35e15eafa98e
# ╟─29313e7c-4d7f-448c-80f0-32f1beb4e4bd
# ╟─74f51da6-9d34-480d-b59a-b37e1cf57714
# ╠═578bd9a2-6ba8-4448-8c9e-4ba458e135e1
# ╟─96d0009a-77bb-4e65-a87d-e85a3660ae26
# ╠═1e188adb-f2fa-4de5-998f-6b8319b37d8f
# ╠═893faaec-7039-40ce-ac74-ab833726ac14
# ╠═080a139a-bafc-4cac-ab16-a202905f179f
# ╟─92017c84-03a2-4fbb-bc70-4b43031d7f44
# ╟─1a4651c5-f895-4d77-a332-eca09cbb2a7e
# ╠═7d1dc93a-71a4-4164-bcf9-6c25db235fa9
# ╠═d4986ba9-295c-4170-bcda-a3db7d8eed06
# ╠═ecea209d-addc-4059-9920-fb9b4709c8f1
# ╠═986879c4-fe6b-4a5c-a19e-5dfcc1d1471a
# ╠═2fd1776f-df92-406c-8dec-42b85ffd741f
# ╠═e55a0993-98ef-47f8-bda4-b7fbf25f44cb
# ╠═50ae1a36-a89d-44b0-a132-1b1d6a37a976
# ╠═aa46b63b-a4e9-4ac0-89b6-5cbd46f4a66e
# ╠═aaead4f8-9b22-4072-a0f2-3fb8b9e36b18
# ╠═f7655e56-9dd7-452c-9d62-3c4474a4c846
# ╠═026598b5-48bd-47f6-8cac-342a7ac5a11a
# ╠═413d05d8-4529-4ac0-9416-d7da01a15446
# ╠═c5397573-cb3c-4ad3-baa1-64cc8a41df0a
# ╠═fed77c8a-d35c-4e67-b112-903c159a8221
# ╠═1bde0c68-15c0-4f44-a33d-4bca9c39e66f
# ╠═f552af14-3a76-4a1c-9e7f-16fa0d34cd3c
# ╠═355b8c18-c6f0-4b68-beb8-5c0698dd1204
# ╠═97266262-a38c-4a8f-ad47-709f04ea4a4d
# ╠═1393048e-4a39-44a4-92e9-27de6add2632
# ╠═764462d4-dd65-43e1-ad9f-bb7b79efbdc2
# ╠═7f7e6146-85f4-48c8-b4a5-f0e7d2acf077
# ╠═6dc5664c-4a72-4584-a64c-47a1cb3b7517
# ╠═d82b8194-4c18-4bd4-92c4-03fa3c7a0f1d
# ╠═511aed8f-95a1-4b01-ae19-b85fb16ea888
# ╠═23c7a8cb-3734-413c-8216-5dc6a12c645b
# ╠═d00de1d2-72c3-48fd-b1e1-6ebb41b73e4c
# ╠═65b2abd0-4365-45b9-b47d-9f3070384cc2
# ╠═427af424-c026-4a30-95a6-37be0a59514b
# ╠═798eec1c-1383-45bb-ae55-b2620ce185fb
# ╠═f816dd6d-2fb7-4c1c-8a73-1fb20ecc1359
# ╠═d83489ca-587d-466e-a95c-579240771eb2
# ╠═3725f3e0-4a7d-47c0-8130-f940a6b1ee7c
# ╠═12323fda-3cc3-4015-9d1e-4b7747650b28
# ╠═751ce0f4-ff2b-4fea-a9a3-84a7d8b6f4ff
# ╠═ee158f3b-ae5d-4129-a1b9-4665c0215084
# ╠═a43e66b7-28d9-462b-9983-03a72623f3e7
# ╠═fd7bf3d1-0fe4-46ce-8024-e75a6c63dde8
# ╠═fc23f850-cd11-45cf-a769-ad4347b2b6cd
# ╠═038e4bfb-c117-44ee-8ac3-b6c19d4696a7
# ╠═11432324-db74-443f-a027-f9d89a5da45f
# ╠═ae47cd84-fff5-485e-8aba-81f8cd13a469
# ╠═548a2d56-77ff-4d75-981b-6d9cc81ea32d
# ╠═bfec4e58-d8a9-40f9-a3e2-43b29b94fecb
# ╠═36ab1db5-edea-4d28-8251-bbd47654acfd
# ╠═8bc53188-6ba4-4c0d-b2be-544d29904822
# ╠═49a36f61-2036-4ea9-b508-612df810794b
# ╠═2c1a42dd-d980-4a87-963d-367ea9b55f5c
# ╠═fe2a8f6f-cd0a-4dd4-88f2-47f483374df2
# ╠═ef62ab21-e48d-4f56-98bd-1eb5a578d5c0
# ╠═32c79f65-5042-4a74-ac54-ee87223d87e3
# ╠═a71b0077-0f01-471c-a132-2c783126e7e7
# ╠═4865c4aa-2db8-46f4-803c-edc51aa9d2be
# ╠═33a3b669-475f-48a8-aae5-168df8436bf5
# ╠═724eba9d-23ea-4681-a6e0-9a6eac94d5b6
# ╠═bdec29c2-d876-4545-9e17-7d8e9c40d61b
# ╠═e37d77b4-37d6-46cf-bfde-42a673c83cbd
# ╠═d96a98d6-208d-42b1-bc06-5ed71e3f0a8e
# ╠═053f971e-f217-49da-a515-477e1393ee54
