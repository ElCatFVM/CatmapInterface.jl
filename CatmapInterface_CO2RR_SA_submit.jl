### A Pluto.jl notebook ###
# v0.20.10

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
	using CSV
    using DataFrames
    using Query
    using DifferentialEquations ## 추가
	using GlobalSensitivity
	using QuasiMonteCarlo
	using Statistics
    using CairoMakie
end

# ╔═╡ 83b85f0b-acf0-42a3-940b-f4a7f1978f57
pathof(CatmapInterface)

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

# ╔═╡ 217a122f-564e-4002-8678-979d83bebdfe
begin 
	  @phconstants N_A c_0 k_B e h
      const F = N_A * e
	  const pH = 6.8
      const T = 298.0 * K
	  const S = 0.169061707523246e20 / m^2## 9.64e4 / N_A * 0.169061707523246e20 *      C*mol/m^2
      const C_gap = 20 * μF / cm^2
      const ϕ_pzc = 0.16 * V
	  const M0 = 18.0153 * ufac"g/mol"
      const v0 = N_A * (8.2 * Å)^3 # 1 / (55.4 * ufac"M"
end

# ╔═╡ 96d0009a-77bb-4e65-a87d-e85a3660ae26
md"""
### Data
"""

# ╔═╡ 7d1dc93a-71a4-4164-bcf9-6c25db235fa9
md"""
#### Model
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

# ╔═╡ e55a0993-98ef-47f8-bda4-b7fbf25f44cb
standard_log_C2=[-1.52161
-1.35566
-1.25237
-1.18859
-1.14937
-1.12532
-1.11065
-1.10182
-1.10647
-1.43261
-1.86391
-2.34455
-2.84323
-3.34786
-3.85442
-4.36163
-4.86901
-5.37661]

# ╔═╡ 50ae1a36-a89d-44b0-a132-1b1d6a37a976
interaction_log_C2=[
-0.764591
-0.085715
0.593796
1.27043
1.93108
2.48462
2.77839
2.88935
2.79146
2.62865
2.43518
2.22593
2.00368
1.76484
1.5007
1.19706
0.837043
0.413982]

# ╔═╡ aaead4f8-9b22-4072-a0f2-3fb8b9e36b18
standard_log_CH4=[5.62477
5.28495
4.91379
4.52288
4.11969
3.70891
3.29345
2.87506
2.45485
2.03347
1.61137
1.18888
0.766206
0.34344
-0.0793755
-0.502248
-0.926041
-1.39315]

# ╔═╡ f7655e56-9dd7-452c-9d62-3c4474a4c846
interaction_log_CH4=[
5.29337
4.89538
4.49418
4.09193
3.68835
3.28106
2.86628
2.44647
2.03011
1.61509
1.20052
0.785868
0.370861
-0.0446811
-0.461
-0.878467
-1.29841
-1.7634]

# ╔═╡ 413d05d8-4529-4ac0-9416-d7da01a15446
standard_log_CO=[1.11462
1.1976
1.24924
1.28113
1.30074
1.31277
1.3201
1.32452
1.3271
1.32852
1.32922
1.32954
1.32966
1.32969
1.32968
1.32962
1.32957
1.32944]

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

# ╔═╡ 9906d038-7065-4ee0-afba-e0d4f2a3dfb1
interaction_log_H2=[2.53813
1.90478
1.2698
0.634859
0.00415143
-0.593421
-1.10918
-1.56713
-2.06638
-2.58304
-3.10722
-3.63227
-4.15497
-4.67263
-5.18483
-5.69041
-6.23563
-6.77025]

# ╔═╡ 75bf9c9a-d43b-4be2-b4c1-7f41bf57dbfb
standard_log_H2=[2.73001
2.13734
1.52908
0.911137
0.28756
-0.338478
-0.963421
-1.58127
-2.18067
-2.7445
-3.25979
-3.73131
-4.17586
-4.60816
-5.03825
-5.49128
-6.01887
-6.49772
]	

# ╔═╡ 7127dd37-2a18-418b-8eb0-41d876be57c9
begin
    file_path = "./data_png/combined_data_ideal_mode.csv" ## interaction:  combined_products_interaction_mode.csv
    df1 = CSV.read(file_path, DataFrame)

end

# ╔═╡ cb35926c-e456-4f64-9b04-cffd5b956750
kinetic_TOF = Dict()																

# ╔═╡ 00fa4e3c-db3f-468b-b11e-9e91a9f129ca
electron_transfer=Dict{String, Int}()

# ╔═╡ ca466e6a-b689-4ad6-9776-394f9a8f44c5
begin
    products = ["CH4_g", "CO_g", "H2_g", "CH3CH2OH_g"]
	elec_trans = [8, 2, 2, 12]
end

# ╔═╡ 1e51d925-dc2a-498a-a9b8-bcfc280f0707
for (i,j) in zip(products, elec_trans)
	electron_transfer[i]= j
end

# ╔═╡ ece46b21-a7f9-4af5-be64-71c2ca3263f9
begin
    for k in products
        result = @from i in df1 begin
            @where i.Product == k && i.pressure == 0.0
            @select (i.Voltage, i.TOF_Coverage)
            @collect DataFrame
        end
        kinetic_TOF[k] = result
    end
end

# ╔═╡ 771eb140-004d-4585-a8fa-85d452810185
current_density=Dict{String, Vector}()

# ╔═╡ cf2a7c0f-c6af-4cc7-929d-18b6063306e2
function TOF_to_current(products)
	for product in products
	current_density[product]=[]
	for i in kinetic_TOF[product][:,2]
		a = log10.(i * F / N_A * 0.169061707523246e20 * electron_transfer[product])
		current_density[product]=push!(current_density[product], a)
	end
	end
end

# ╔═╡ 3d36eb3b-712a-4465-ba01-677ad93ffc3d
TOF_to_current(products)

# ╔═╡ c60249a7-a4ba-4885-b53c-88fcff08c081
md"""
### Input
"""

energy= model_instance.catmap_params.species_list


# ╔═╡ 4433770f-45a5-4df5-bee5-cc47f9b9dd4f
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

@bind params_input PlutoUI.combine() do Child
    input_params = [
        (; name = "ϕ_we", range = -1.85:0.05:-1.0, default = -1.85, unit = "V"),
        (; name = "local_pH", range = 5.0:0.2:9.0, default = 6.8, unit = ""),
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

	for (s ,sp) in model_instance.catmap_params.species_list
		if (isa(sp, TStateSpecies)) && occursin("Δele", s)
			push!(
				input_params,
				(; name = "β$s", range = 0.0:0.01:1.0, default = sp.β, unit = ""),
			)
		end
	end

	for (s, sp) in model_instance.catmap_params.species_list
		if (isa(sp, TStateSpecies))
			push!(
				input_params,
				(; name = "E$s", range = -1e10:1e5:1e10, default = energy[s].formation_energy , unit = "J"),
			)
		end
	end

	 pressures_input = [
        md""" $(name) : $(Child(name, PlutoUI.Slider(range; default=default, show_value=true))) $unit
        """ for (; name, range, default, unit) in input_params
    ] 
    md"""
    ##### Input Params:
    $(pressures_input)
    """
end

# ╔═╡ 9e765a39-95ab-497b-8443-22714db21cb5
md"""
### Functions for kinetic
"""

# ╔═╡ b1c1a32b-686b-4cc7-ab36-3b063edced77
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

# ╔═╡ 65a43473-fd54-4dd7-a02e-11ed66824d76
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

# ╔═╡ cec44348-9032-4a96-ab5b-1e8373b1287f
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

# ╔═╡ 4bd582fd-fe32-4a2f-a71a-9dc74887d1d7
function electrontransfer_H2(rn::Catalyst.ReactionSystem, ssol, params) ## look more
    rxns = reactions(rn)
	 ir_forward_tafel = findfirst(rxns) do r
        (; products, substrates) = r
      !isempty(products) &&
           Symbolics.tosymbol.(products; escape = false) == [:H2_g]
    end
	
    ir_back_tafel = findfirst(rxns) do r
        (; substrates, products) = r
        !isempty(products) && Symbolics.tosymbol.(substrates; escape=false) == [:H2_g]
    end
	
    ir_forward_hey = findfirst(rxns) do r
        (; products, substrates) = r
      !isempty(products) &&
           Symbolics.tosymbol.(products; escape = false) == [:H2_g, :OH_g]
    end
	
    ir_back_hey = findfirst(rxns) do r
        (; products, substrates) = r
    !isempty(products) &&
          Symbolics.tosymbol.(products; escape = false) == [:H_t]
    end
    curr = 0.0
	curr_1=0.0
	 if !isnothing(ir_forward_tafel)
        curr += substitute(
            substitute(rxns[ir_forward_tafel].rate, Dict(sp => ssol[sp] for sp in species(rn))),
            merge(rn.defaults, symmap_to_varmap(rn, params)),
        )
    end
    if !isnothing(ir_forward_hey)
        curr_1 += substitute(
            substitute(rxns[ir_forward_hey].rate, Dict(sp => ssol[sp] for sp in species(rn))),
            merge(rn.defaults, symmap_to_varmap(rn, params)),
        )
    end
	if !isnothing(ir_back_tafel)
       curr -= substitute(
            substitute(rxns[ir_back_tafel].rate, Dict(sp => ssol[sp] for sp in species(rn))),
           merge(rn.defaults, symmap_to_varmap(rn, params)),
       )
    end
  if !isnothing(ir_back_hey)
       curr_1 -= substitute(
            substitute(rxns[ir_back_hey].rate, Dict(sp => ssol[sp] for sp in species(rn))),
           merge(rn.defaults, symmap_to_varmap(rn, params)),
       )
    end
	#println(ir_forward_hey, ir_back_hey)
    #println(curr * ph"e" * S * 8 * 0.1) #current_density
    #println( curr / S)
    return log10(abs(curr_1 * ph"e" * S * 2+ curr * ph"e" * S * 2)) ## 10is for A/m^2 to mA/cm^2 but we don't have to there is ufac in the plotting
end

# ╔═╡ 1aeb6857-85a9-447c-9b9a-6fa880f97b39
md"""
### Current_voltage function
"""

function currentvoltage(rn, catmap_params, pressures, θ0, params; solver=DynamicSS(Rodas5P()), maxiters=maxiters, abstol = abstol, reltol = reltol)
	nr = numreactions(rn)
	stoichmat = netstoichmat(rn)
	rrs_sym = reactionrates(rn)
	rrs_num = zeros(nr)

	Δϕs = params[:ϕ_we]:0.05:-1.0
	u0 = symmap_to_varmap(rn, merge(pressures, θ0)) 
		currs_CH4 = zeros(length(Δϕs))
		currs_CO = zeros(length(Δϕs))
		currs_CH3CH2OH = zeros(length(Δϕs))
		currs_H2 = zeros(length(Δϕs))
	for (iΔϕ, Δϕ) in enumerate(Δϕs)
		 params[:ϕ] = params[:ϕ_we] - Δϕ
		if catmap_params.electrochemical_thermo_mode == :hbond_surface_charge_density
			params[:σ] = surface_charge_relation(Δϕ)
		end
		ssprob = SteadyStateProblem(
			rn, 
			u0, 
			symmap_to_varmap(rn, params)
		)
		ssol = solve(ssprob, solver; maxiters, abstol, reltol)
		currs_CH4[iΔϕ] = electrontransfer_CH4(rn, ssol, params)
		currs_CO[iΔϕ] =  electrontransfer_CO(rn, ssol, params) 
		currs_CH3CH2OH[iΔϕ] =  electrontransfer_CH3CH2OH(rn, ssol, params) 
		currs_H2[iΔϕ] =  electrontransfer_H2(rn, ssol, params) 
		u0 = Dict(sp => ssol[sp] for sp in species(rn))
	end
        return collect(Δϕs), currs_CH4, currs_CO, currs_CH3CH2OH, currs_H2
	end

# ╔═╡ a405cf23-07fe-41df-83b0-2fbf6f46959b
md"""
#### chunk_parallel
"""

# ╔═╡ 43217a9b-dd21-4dff-a55e-b15d71bc6864
# ╠═╡ disabled = true
#=╠═╡
function merge_currentvoltage_results(results)
    Δϕs      = Float64[]
    currs_CH4    = Float64[]
    currs_CO     = Float64[]
    currs_CH3CH2OH = Float64[]
    currs_H2     = Float64[]

    for (ϕs, ch4, co, etoh, h2) in results
        append!(Δϕs, ϕs)
        append!(currs_CH4, ch4)
        append!(currs_CO, co)
        append!(currs_CH3CH2OH, etoh)
        append!(currs_H2, h2)
    end

    return Δϕs, currs_CH4, currs_CO, currs_CH3CH2OH, currs_H2
end

  ╠═╡ =#

# ╔═╡ 04f74348-a536-4f11-ac47-cbf323344c8d
md"""
#### serial
"""

# ╔═╡ 2a22ec17-9bfa-49f2-ac63-bde821706f63
md"""
### Classify parameters
"""

# ╔═╡ c600733a-5f3f-44ed-b6f5-876704103231
function parameters_to_string(parameters)
	params_name=[]
	for p in parameters
		params=string(p)
		push!(params_name, params)
	end
	return params_name
end				

# ╔═╡ 6522d4cc-0644-4e4f-9bb3-e304a5c42eff
function pcet_name(name)
    pcet_barrier = []
    chem_barrier = []
    pcet_beta = []
    for i in name
        i_replaced = replace(i, "E" => "G") ## replace E to G
        if contains(i_replaced, "ele") && contains(i_replaced, "G") && contains(i_replaced, "Δ")
            push!(pcet_barrier, i_replaced)
        elseif contains(i_replaced, "G") && contains(i_replaced, "Δ")
            push!(chem_barrier, i_replaced)
        elseif contains(i_replaced, "ele") && contains(i_replaced, "β") && contains(i_replaced, "Δ")
            push!(pcet_beta, i_replaced)
        end
    end
    return pcet_barrier, chem_barrier, pcet_beta
end


# ╔═╡ 1bde0c68-15c0-4f44-a33d-4bca9c39e66f
U_eq=Dict()

# ╔═╡ 1a25c56d-cde5-40ca-a102-4becd2a54d10




# ╔═╡ f552af14-3a76-4a1c-9e7f-16fa0d34cd3c
function get_U_eq(reactions)
	for i in reactions
		if isnothing(i.tstate)
		nothing
		else
	    j=[pair.first for pair in  i.educts] ## name 
	    k=[pair.second for pair in i.educts] ## coeff
	    l=[pair.first for pair in i.products]
	    m=[pair.second for pair in i.products]
	    	U_eq["$(i.tstate.components[1][1])"]=(sum(energy["$(j)"].formation_energy*k for (j,k) in zip(j,k))-sum(energy["$(l)"].formation_energy*m for (l,m) in zip(l,m)))/e
		 
    	end
	end
end

# ╔═╡ d2236c4b-5c72-4b9b-9533-22c767c44b84
get_U_eq(model_instance.catmap_params.reactions)

# ╔═╡ c9d38446-019f-4f40-afda-544f40ff19f0
md"""
### loss function
"""

# ╔═╡ a127638f-1bd0-43d4-bf42-1074e967e064
const θinit = 1.0e-10 ## 1.0e-10 

# ╔═╡ 355b8c18-c6f0-4b68-beb8-5c0698dd1204
begin 
	const pressures = Dict(pairs(pressures_input))
	const (; rn, catmap_params) = model_instance
	const species_list = catmap_params.species_list
 	const	params = Dict(pairs(params_input))
	const θ0 = Dict(Symbol(s) => θinit for (s, sp) in species_list if isa(sp, AdsorbateSpecies))
	const rxns=reactions(rn)
	maxiters=1e6  
	abstol=1e-6 #1e-8 default, for making prior use 1e-6
	reltol=1e-4 #1e-6 default  for making prior use 1e-4
	best_result = Ref{Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}}((-Inf, [], [], [], []))
	CH4= best_result[][2]
	CO= best_result[][3]
	CH3CH2OH= best_result[][4]
	H2= best_result[][5] ## 1.0e-10 
	const warning_count = Ref(0)
	const warning_params = Ref(Vector{NamedTuple}())
end

# ╔═╡ 0eb97e69-b983-4abc-af7f-0551c6913fa6
# ╠═╡ disabled = true
#=╠═╡
function currentvoltage_serial_chunk(
    rn, catmap_params, pressures, θ0, params, Δϕs_chunk;
    solver = DynamicSS(TRBDF2()), maxiters = maxiters, abstol = abstol, reltol = reltol
) ##Rodas5p()9
    nr = numreactions(rn)
    stoichmat = netstoichmat(rn)
    rrs_sym = reactionrates(rn)
    rrs_num = zeros(nr)

    u0 = symmap_to_varmap(rn, merge(pressures, θ0))
    currs_CH4 = zeros(length(Δϕs_chunk))
    currs_CO = zeros(length(Δϕs_chunk))
    currs_CH3CH2OH = zeros(length(Δϕs_chunk))
    currs_H2 = zeros(length(Δϕs_chunk))

    for (iΔϕ, Δϕ) in enumerate(Δϕs_chunk)
        local_params = deepcopy(params)
        local_params[:ϕ] = local_params[:ϕ_we] - Δϕ
        if catmap_params.electrochemical_thermo_mode == :hbond_surface_charge_density
            local_params[:σ] = surface_charge_relation(Δϕ)
        end

        ssprob = SteadyStateProblem(rn, u0, symmap_to_varmap(rn, local_params))
            ssol = solve(ssprob, solver; maxiters, abstol, reltol)
			if ssol.retcode != :Success
    			@warn "Solver failed: retcode = $(ssol.retcode)"
    		return (collect(Δϕs_chunk), zeros(length(Δϕs_chunk)), zeros(length(Δϕs_chunk)), zeros(length(Δϕs_chunk)), zeros(length(Δϕs_chunk)))
			end
            currs_CH4[iΔϕ]       = electrontransfer_CH4(rn, ssol, local_params)
            currs_CO[iΔϕ]        = electrontransfer_CO(rn, ssol, local_params)
            currs_CH3CH2OH[iΔϕ]  = electrontransfer_CH3CH2OH(rn, ssol, local_params)
            currs_H2[iΔϕ]        = electrontransfer_H2(rn, ssol, local_params)

            u0 = Dict(sp => ssol[sp] for sp in species(rn))  # update for next loop
    end

    return collect(Δϕs_chunk), currs_CH4, currs_CO, currs_CH3CH2OH, currs_H2
end

  ╠═╡ =#

# ╔═╡ 37a2514d-b546-492b-90cb-38d316d91972
params_name = parameters_to_string(Catalyst.parameters(rn))

# ╔═╡ d5fa068a-0713-4292-b62b-eca37f10706a
pcet_barrier, chem_barrier, pcet_beta = pcet_name(params_name)

# ╔═╡ 5eb78291-1fad-41bd-af32-2d89e9ceb98b
Catalyst.parameters(rn)

# ╔═╡ 8bd2dc6f-57fb-40fe-8efd-b59eb4b7bce9
corrections = Dict(zip(keys(catmap_params.species_list), zeros(Float64, length(catmap_params.species_list))))

# ╔═╡ 6b9cf7b2-6362-4b86-95f7-f056999e1c75
function compute_thermo_corrections!(corrections, catmap_params::CatmapParams)
    (; gas_thermo_mode, adsorbate_thermo_mode) = catmap_params
    gas_thermo_correction!             = getfield(CatmapInterface, Symbol(gas_thermo_mode))
    adsorbate_thermo_correction!       = getfield(CatmapInterface, Symbol(adsorbate_thermo_mode))
    gas_thermo_correction!(corrections, catmap_params)
    adsorbate_thermo_correction!(corrections, catmap_params)
	return corrections
end

# ╔═╡ 9fdf0ea0-d488-4996-92b8-181160f03a54
thermo_corrections = compute_thermo_corrections!(corrections, catmap_params::CatmapParams)

# ╔═╡ a123d077-2b8d-410a-a875-f8b5e0e850b0
md"""
### Checking the loss function
"""

# ╔═╡ 1f43f628-7c35-4f11-9b12-27e73214a35b
begin
	function loss_function(parameters::NamedTuple)
		result = nothing
		t= @elapsed begin
		β_names= pcet_beta[begin:end] ## to exclude HER
		E_pcet_names= pcet_barrier[begin:end] ## to exclude HER
		E_chem_names= chem_barrier[begin:end] ## to exclude HER
		β_values = Dict(name => parameters[Symbol(name)] for name in β_names)
	    barrier_pcet_values = Dict(name => parameters[Symbol(name)] for name in E_pcet_names)
		barrier_chem_values = Dict(name => parameters[Symbol(name)] for name in E_chem_names)
		local_params=copy(params)
		for (i,rxn) in enumerate(model_instance.catmap_params.reactions)
			if isnothing(rxn.tstate)
				continue
			end
			tstate_name = rxn.tstate.components[1][1]
			#if tstate_name in  ["H2OΔele_t", "HΔH2OΔele_t", "HΔH_t"] ## remove HER
			#	continue
			#end
			
			if "ele_g" in  [pair.first for pair in  rxn.educts] ## pcet_step
				j=[pair.first for pair in  rxn.educts] ## name 
		    	k=[pair.second for pair in rxn.educts] ## coeff
				val = barrier_pcet_values["G"*tstate_name]*e +sum(energy["$(j)"].formation_energy*k+thermo_corrections["$(j)"]*k for (j,k) in zip(j,k))-thermo_corrections[tstate_name]-β_values["β"*tstate_name]*e*U_eq[tstate_name]
		    local_params[Symbol("E"*tstate_name)] = val
			else ## chem_step
				j=[pair.first for pair in  rxn.educts] ## name 
		    	k=[pair.second for pair in rxn.educts] ## coeff
				val = barrier_chem_values["G"*tstate_name]*e +sum(energy["$(j)"].formation_energy*k+thermo_corrections["$(j)"]*k for (j,k) in zip(j,k))-thermo_corrections[tstate_name]
			local_params[Symbol("E"*tstate_name)] = val	
			end
		end
		for (i,j) in β_values
	  		local_params[Symbol(i)]= j
	    end
		
	 (Δϕs, currs_CH4, currs_CO, currs_CH3CH2OH, currs_H2) = currentvoltage(rn, catmap_params, pressures, θ0, local_params)
		result= -(norm(standard_log_CH4-currs_CH4)+norm(standard_log_C2-currs_CH3CH2OH)+norm(standard_log_CO-currs_CO)+norm(standard_log_H2-currs_H2))
		end
		 if !isfinite(result)
			@warn "Non-finite or zero current, returning fallback"
			 warning_count[] += 1
			 push!(warning_params[], parameters)
			return -1e6
		end
		 if result > best_result[][1]
	       
			 best_result[] = (result, currs_CH4, currs_CO, currs_CH3CH2OH, currs_H2)
	    end
	#		open("log_txt/HER_serial_low_tol_inter_abstol=$(abstol)_reltol=$(reltol)_with_iter=$(n_iter)_init=$(init_points)_randstate=$(random_state).txt", "a") do io
    #    println(io, "elapsed: $(round(t, digits=4))s, E1: $(parameters[1]), E2: $(parameters[2]), E3: $(parameters[3]), E4: $(parameters[4]),E5: $(parameters[5]), E6: $(parameters[6]),E7: $(parameters[7]), E8: $(parameters[8]), E9: $(parameters[9]), E10:$(parameters[10]) β1: $(parameters[11]), β2: $(parameters[12]), β3: $(parameters[13]), β4: $(parameters[14]), β5: $(parameters[15]), β6: $(parameters[16])")
	#end
		return result
end
end

# ╔═╡ 312fa03e-a913-4a97-a3db-0cf436150a41
function loss_function(; kwargs...)
    loss_function(NamedTuple(kwargs))
end


# ╔═╡ 7ec99492-254e-4bcf-ad07-075eb4b6f2eb
begin
	pbounds= Dict("GH2OΔele_t" => (0.0   ,1.8), "GHΔH2OΔele_t" => (0.0,1.8), "GH2OΔCOΔele_t" => (0.0,1.8), "GCHΔOHΔele_t" => (0.0,1.8), "GOCCOΔH2OΔele_t" => (0.0,1.8), "GCOOHΔH2OΔele_t" => (0.0,1.8), "GHΔH_t" => (0.0, 1.8), "GCOΔ_t" => (0.0,1.8), "GOCΔCO_t" => (0.0,1.8), "GCO2Δ_t" => (0.0,1.8), "βH2OΔele_t" => (0.1,0.9), "βHΔH2OΔele_t" => (0.1, 0.9), "βH2OΔCOΔele_t" => (0.1,0.9), "βCHΔOHΔele_t" => (0.1,0.9), "βOCCOΔH2OΔele_t" => (0.1,0.9), "βCOOHΔH2OΔele_t" => (0.1,0.9))
end

# ╔═╡ c1bf4c12-5292-4188-8839-d4fc27348f12
parameters = (
    GH2OΔele_t=1.351,
GHΔH2OΔele_t=1.52281, GH2OΔCOΔele_t=0.966143, GCHΔOHΔele_t=0.803231, GOCCOΔH2OΔele_t=0.508165, GCOOHΔH2OΔele_t=0.699272, GHΔH_t= 0.786548, GCOΔ_t= 0.328544, GOCΔCO_t= 0.3998, GCO2Δ_t=0.103707, 
   βH2OΔele_t=0.5,
βHΔH2OΔele_t=0.89, βH2OΔCOΔele_t=0.5, βCHΔOHΔele_t=0.5, βOCCOΔH2OΔele_t=0.6, βCOOHΔH2OΔele_t=0.5
) ## theortical value

# ╔═╡ f76bc64f-244b-4658-87db-238274ca8093
begin
	a= loss_function(parameters)
	println(a)
end

# ╔═╡ a2692d58-37df-4e2e-8e88-22240caf555f
md"""
### Sensitivity Analysis
"""

# ╔═╡ 73d7d8a2-e044-4caa-ad8c-3de4e557bd04


begin
	function peak_cases(p)
			result = nothing
		    parameters = (
   GH2OΔele_t=p[1],
GHΔH2OΔele_t=p[2], GH2OΔCOΔele_t = p[3], GCHΔOHΔele_t = p[4], GOCCOΔH2OΔele_t = p[5], GCOOHΔH2OΔele_t = p[6], GHΔH_t= p[7],
    GCOΔ_t = p[8], GOCΔCO_t= p[9], GCO2Δ_t = p[10],
   βH2OΔele_t=p[11],
βHΔH2OΔele_t=p[12], βH2OΔCOΔele_t = p[13], βCHΔOHΔele_t = p[14], βOCCOΔH2OΔele_t = p[15], βCOOHΔH2OΔele_t = p[16]
)
			β_names= pcet_beta[begin:end] ## to exclude HER
			E_pcet_names= pcet_barrier[begin:end] ## to exclude HER
			E_chem_names= chem_barrier[begin:end] ## to exclude HER
			β_values = Dict(name => parameters[Symbol(name)] for name in β_names)
		    barrier_pcet_values = Dict(name => parameters[Symbol(name)] for name in E_pcet_names)
			barrier_chem_values = Dict(name => parameters[Symbol(name)] for name in E_chem_names)
			local_params=copy(params)
			for (i,rxn) in enumerate(model_instance.catmap_params.reactions)
				if isnothing(rxn.tstate)
					continue
				end
				tstate_name = rxn.tstate.components[1][1]
				
				if "ele_g" in  [pair.first for pair in  rxn.educts] ## pcet_step
					j=[pair.first for pair in  rxn.educts] ## name 
			    	k=[pair.second for pair in rxn.educts] ## coeff
					val = barrier_pcet_values["G"*tstate_name]*e +sum(energy["$(j)"].formation_energy*k+thermo_corrections["$(j)"]*k for (j,k) in zip(j,k))-thermo_corrections[tstate_name]-β_values["β"*tstate_name]*e*U_eq[tstate_name]
		    local_params[Symbol("E"*tstate_name)] = val
				else ## chem_step
					j=[pair.first for pair in  rxn.educts] ## name 
			    	k=[pair.second for pair in rxn.educts] ## coeff
					val = barrier_chem_values["G"*tstate_name]*e +sum(energy["$(j)"].formation_energy*k+thermo_corrections["$(j)"]*k for (j,k) in zip(j,k))-thermo_corrections[tstate_name]
			local_params[Symbol("E"*tstate_name)] = val	
			end
			end
			for (i,j) in β_values
		  		local_params[Symbol(i)]= j
		    end 	
			 (Δϕs, currs_CH4, currs_CO, currs_CH3CH2OH, currs_H2)= currentvoltage(rn, catmap_params, pressures, θ0, local_params)  # adapt as needed

             result= -(norm(standard_log_CH4-currs_CH4)+norm(standard_log_C2-currs_CH3CH2OH)+norm(standard_log_CO-currs_CO)+norm(standard_log_H2-currs_H2))

			return (isfinite(result) && result ≠ -Inf) ? result : NaN
    		end
end
  ╠═╡ =#

# ╔═╡ 78ea7115-6d0d-4336-9cff-7b83c20675eb
bounds= [(0.01, 1.8), (0.01, 1.8), (0.01, 1.8), (0.01, 1.8), (0.01, 1.8), (0.01, 1.8), (0.01, 1.8), (0.01, 1.8), (0.01, 1.8), (0.01, 1.8), (0.1, 0.9), (0.1, 0.9), (0.1, 0.9), (0.1, 0.9), (0.1, 0.9), (0.1, 0.9)] 

# ╔═╡ 8fb9de64-0286-41db-9094-403f59e1b00f
begin
	samples = 5000
	lb = first.(bounds)
	ub = last.(bounds)
end

# ╔═╡ bb527498-120a-4c77-9116-e1cff829bf0c
begin
	sampler = SobolSample()
	A, B = QuasiMonteCarlo.generate_design_matrices(samples, lb, ub, sampler)
end

# ╔═╡ c1ce4aaa-1612-43b4-9f17-1cd83272d1e4
function peak_cases_batch(X::AbstractMatrix)
    return [peak_cases(vec(col)) for col in eachcol(X)]
end

  ╠═╡ =#

# ╔═╡ f90f10b9-0f03-4e01-9671-996b2be34be8

begin
		Y_A = peak_cases_batch(A)
		Y_B = peak_cases_batch(B)
		
		# 2. Keep only valid indices (where both A and B results are good)
		valid_idx = findall(i -> isfinite(Y_A[i]) && isfinite(Y_B[i]), 1:length(Y_A))
		
		# 3. Filter inputs and outputs
		A_valid = A[:, valid_idx]
		B_valid = B[:, valid_idx]
	res = gsa(peak_cases_batch, Sobol(), A_valid, B_valid, batch=true, Ei_estimator = :Sobol2007)
end
  ╠═╡ =#

# ╔═╡ 66d26e21-8762-4a4d-95b3-8e185a7cee45
println(res)
  ╠═╡ =#

# ╔═╡ 3b6ca314-cbff-4fbf-998f-5f0d5a7846f1
# ╠═╡ disabled = true
#=╠═╡
@show var(filter(!isnan, Y_A)), var(filter(!isnan, Y_B))
  ╠═╡ =#

# ╔═╡ 0f9aa6b7-744d-475a-a998-da03d7b4a2e7
set_theme!(Theme(
    fontsize = 16,             # sets base font size
    Axis = (
        titlesize = 20,
        xlabelsize = 10,
        ylabelsize = 30,
        xticklabelsize = 15,
        yticklabelsize = 25,
    )
))

# ╔═╡ c3c91dd2-d01d-44f8-ab5e-597b6bf32e72

begin
f = Figure()
	tbl = (cat = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],
	       height= res.ST)
	
    barplot(f[1,1], tbl.cat, tbl.height,
	        axis = (xticks = (1:16, ["G1", "G2", "G3", "G4",
             "G5", "G6", "G7", "G8", "G9", "G10",
             "β1", "β2", "β3", "β4", "β5", "β6"]),
	                title = "Total Order Indices"),
	        )
end
  ╠═╡ =#

# ╔═╡ f06d046f-ab0a-4949-8c7b-864b42e736b7
begin
	tbl1 = (cat = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],
	       height= res.S1)

barplot(f[1,2], tbl1.cat, tbl1.height,
	        axis = (xticks = (1:16, ["G1", "G2", "G3", "G4",
             "G5", "G6", "G7", "G8", "G9", "G10",
             "β1", "β2", "β3", "β4", "β5", "β6"]),
	                title = "First Order Indices"),
	        )
save("SA_fig.svg", f)
end
  ╠═╡ =#

# ╔═╡ ae47cd84-fff5-485e-8aba-81f8cd13a469
voltage = [
  -1.85
-1.8
-1.75
-1.7
-1.65
-1.6
-1.55
-1.5
-1.45
-1.4
-1.35
-1.3
-1.25
-1.2
-1.15
-1.1
-1.05
-1.0
]

# ╔═╡ 0e8788fb-4835-4b21-adc5-65a90b9426cc
md"""
### Plots(from Makie)
"""

# ╔═╡ 548a2d56-77ff-4d75-981b-6d9cc81ea32d
#let 
#    f1 = Figure(resolution = (1200, 1200))
#    ax4 = Axis(
#        f1[1, 1];
#        title = "Current-Voltage Curve_catmapInterface_CH4",
#        xlabel = "Δϕ [V]",
#        ylabel = "I #[mA/cm^2]",
#		titlesize = 20,
#       
#        limits = (-1.85, -1.0, -5, 5),
#    )
#    lines!(ax4, voltage,  best_result[][2], color=:red, linewidth =5)
#	lines!(ax4, kinetic_TOF["H2_g"][:,1], current_density["CH4_g"][:,1] , color=:blue, linewidth=3 )
#	
#	ax5 = Axis(
#        f1[1, 2],
#        title = "Current-Voltage Curve_catmapInterface_CO",
#        xlabel = "Δϕ [V]",
#        ylabel = "I #[mA/cm^2]",
#		titlesize = 20,
#		
#        limits = (-1.85, -1.0, -5, 7),  # Adjust limits as needed
#    )
#    lines!(ax5, voltage,   best_result[][3], color = :red, linewidth =5 )
#	lines!(ax5, kinetic_TOF["H2_g"][:,1] ,  current_density["CO_g"][:,1] , color = :blue, linewidth =3)
#
#	ax6 = Axis(
#        f1[2, 1],
#        title = "Current-Voltage Curve_catmapInterface_C2",
#        xlabel = "Δϕ [V]",
#        ylabel = "I #[mA/cm^2]",
#		titlesize = 20,
#		
#        limits = (-1.85, -1.0, -10, 5),  # Adjust limits as needed
#    )
#    lines!(ax6, voltage,   best_result[][4], color = :red, linewidth =5 )
#	lines!(ax6, kinetic_TOF["H2_g"][:,1] , current_density["CH3CH2OH_g"][:,1]  , color = :blue, linewidth =3)
#	
#	ax7 = Axis(
#        f1[2, 2],
#        title = "Current-Voltage Curve_catmapInterface_H2",
#        xlabel = "Δϕ [V]",
#        ylabel = "I #[mA/cm^2]",
#		titlesize = 20,
#		
#        limits = (-1.85, -1.0, -20, 6),  # Adjust limits as needed
#    )
#    lines!(ax7, voltage,  best_result[][5], color = :red, linewidth =5 )
#	lines!(ax7, kinetic_TOF["H2_g"][:,1]  , current_density["H2_g"] , color = :blue, linewidth =3)
#	
#	f1
#	
#end


# ╔═╡ 4bc1093a-950a-4ca2-8294-105b4b31bedd
#best_result[][5]

# ╔═╡ 416ccdf0-eae1-41e0-b60a-aaa8cbcade12
# ╠═╡ disabled = true
#=╠═╡
function currentvoltage(rn, catmap_params, pressures, θ0, params;)
    Δϕs = params[:ϕ_we]:0.05:-1.0
    N = Threads.nthreads()
    
    chunks = Iterators.partition(Δϕs, cld(length(Δϕs), N))

    tasks = map(chunks) do chunk
        Threads.@spawn begin
            try
             return currentvoltage_serial_chunk(rn, catmap_params, pressures, θ0, params, collect(chunk))
            catch e
				@warn "Thread crash: $e"
                n = length(chunk)
                return (collect(chunk), zeros(n), zeros(n), zeros(n), zeros(n))
    		end
			end
	end

	local results
	
    try
        results = fetch.(tasks)
    catch e
        @warn "Outer fetch error in currentvoltage: $e"
    	return (Δϕs, zeros(length(Δϕs)), zeros(length(Δϕs)), zeros(length(Δϕs)), zeros(length(Δϕs)))
	else
		return merge_currentvoltage_results(results)
    end
end
  ╠═╡ =#

# ╔═╡ b1d76daa-452e-46a5-96a8-b2cf335deaab



# ╔═╡ Cell order:
# ╠═7db4d4bd-bc8b-4949-a241-537bcea259f0
# ╠═83b85f0b-acf0-42a3-940b-f4a7f1978f57
# ╟─29313e7c-4d7f-448c-80f0-32f1beb4e4bd
# ╟─74f51da6-9d34-480d-b59a-b37e1cf57714
# ╠═578bd9a2-6ba8-4448-8c9e-4ba458e135e1
# ╠═217a122f-564e-4002-8678-979d83bebdfe
# ╟─96d0009a-77bb-4e65-a87d-e85a3660ae26
# ╟─7d1dc93a-71a4-4164-bcf9-6c25db235fa9
# ╠═d4986ba9-295c-4170-bcda-a3db7d8eed06
# ╠═ecea209d-addc-4059-9920-fb9b4709c8f1
# ╠═986879c4-fe6b-4a5c-a19e-5dfcc1d1471a
# ╠═e55a0993-98ef-47f8-bda4-b7fbf25f44cb
# ╠═50ae1a36-a89d-44b0-a132-1b1d6a37a976
# ╠═aaead4f8-9b22-4072-a0f2-3fb8b9e36b18
# ╠═f7655e56-9dd7-452c-9d62-3c4474a4c846
# ╠═413d05d8-4529-4ac0-9416-d7da01a15446
# ╠═c5397573-cb3c-4ad3-baa1-64cc8a41df0a
# ╠═9906d038-7065-4ee0-afba-e0d4f2a3dfb1
# ╠═75bf9c9a-d43b-4be2-b4c1-7f41bf57dbfb
# ╠═7127dd37-2a18-418b-8eb0-41d876be57c9
# ╠═cb35926c-e456-4f64-9b04-cffd5b956750
# ╠═00fa4e3c-db3f-468b-b11e-9e91a9f129ca
# ╠═ca466e6a-b689-4ad6-9776-394f9a8f44c5
# ╠═1e51d925-dc2a-498a-a9b8-bcfc280f0707
# ╠═6daa6e42-e868-4d3f-bc40-598796580068
# ╠═ece46b21-a7f9-4af5-be64-71c2ca3263f9
# ╠═771eb140-004d-4585-a8fa-85d452810185
# ╠═cf2a7c0f-c6af-4cc7-929d-18b6063306e2
# ╠═3d36eb3b-712a-4465-ba01-677ad93ffc3d
# ╠═c60249a7-a4ba-4885-b53c-88fcff08c081
# ╠═6ad2e68e-37a2-4155-9766-26478e7c07a2
# ╠═4433770f-45a5-4df5-bee5-cc47f9b9dd4f
# ╠═9e765a39-95ab-497b-8443-22714db21cb5
# ╠═b1c1a32b-686b-4cc7-ab36-3b063edced77
# ╠═65a43473-fd54-4dd7-a02e-11ed66824d76
# ╠═cec44348-9032-4a96-ab5b-1e8373b1287f
# ╠═4bd582fd-fe32-4a2f-a71a-9dc74887d1d7
# ╠═1aeb6857-85a9-447c-9b9a-6fa880f97b39
# ╠═a405cf23-07fe-41df-83b0-2fbf6f46959b
# ╠═416ccdf0-eae1-41e0-b60a-aaa8cbcade12
# ╠═0eb97e69-b983-4abc-af7f-0551c6913fa6
# ╠═43217a9b-dd21-4dff-a55e-b15d71bc6864
# ╠═04f74348-a536-4f11-ac47-cbf323344c8d
# ╠═b1d76daa-452e-46a5-96a8-b2cf335deaab
# ╠═2a22ec17-9bfa-49f2-ac63-bde821706f63
# ╠═c600733a-5f3f-44ed-b6f5-876704103231
# ╠═37a2514d-b546-492b-90cb-38d316d91972
# ╠═5eb78291-1fad-41bd-af32-2d89e9ceb98b
# ╠═6522d4cc-0644-4e4f-9bb3-e304a5c42eff
# ╠═d5fa068a-0713-4292-b62b-eca37f10706a
# ╠═1bde0c68-15c0-4f44-a33d-4bca9c39e66f
# ╠═1a25c56d-cde5-40ca-a102-4becd2a54d10
# ╠═f552af14-3a76-4a1c-9e7f-16fa0d34cd3c
# ╠═d2236c4b-5c72-4b9b-9533-22c767c44b84
# ╟─c9d38446-019f-4f40-afda-544f40ff19f0
# ╠═a127638f-1bd0-43d4-bf42-1074e967e064
# ╠═355b8c18-c6f0-4b68-beb8-5c0698dd1204
# ╠═8bd2dc6f-57fb-40fe-8efd-b59eb4b7bce9
# ╠═6b9cf7b2-6362-4b86-95f7-f056999e1c75
# ╠═9fdf0ea0-d488-4996-92b8-181160f03a54
# ╟─a123d077-2b8d-410a-a875-f8b5e0e850b0
# ╠═1f43f628-7c35-4f11-9b12-27e73214a35b
# ╠═312fa03e-a913-4a97-a3db-0cf436150a41
# ╠═7ec99492-254e-4bcf-ad07-075eb4b6f2eb
# ╠═c1bf4c12-5292-4188-8839-d4fc27348f12
# ╠═f76bc64f-244b-4658-87db-238274ca8093
# ╟─a2692d58-37df-4e2e-8e88-22240caf555f
# ╠═73d7d8a2-e044-4caa-ad8c-3de4e557bd04
# ╠═78ea7115-6d0d-4336-9cff-7b83c20675eb
# ╠═8fb9de64-0286-41db-9094-403f59e1b00f
# ╠═bb527498-120a-4c77-9116-e1cff829bf0c
# ╠═c1ce4aaa-1612-43b4-9f17-1cd83272d1e4
# ╠═f90f10b9-0f03-4e01-9671-996b2be34be8
# ╠═66d26e21-8762-4a4d-95b3-8e185a7cee45
# ╠═3b6ca314-cbff-4fbf-998f-5f0d5a7846f1
# ╠═0f9aa6b7-744d-475a-a998-da03d7b4a2e7
# ╠═c3c91dd2-d01d-44f8-ab5e-597b6bf32e72
# ╠═f06d046f-ab0a-4949-8c7b-864b42e736b7
# ╠═ae47cd84-fff5-485e-8aba-81f8cd13a469
# ╟─0e8788fb-4835-4b21-adc5-65a90b9426cc
# ╠═548a2d56-77ff-4d75-981b-6d9cc81ea32d
# ╠═4bc1093a-950a-4ca2-8294-105b4b31bedd
