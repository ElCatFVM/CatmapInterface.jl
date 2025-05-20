### A Pluto.jl notebook ###
# v0.20.8

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

# ╔═╡ f5a97e3e-e116-11ee-210f-cfc211ee8cbc
begin
    import Pkg as _Pkg
    haskey(ENV, "PLUTO_PROJECT") && _Pkg.activate(ENV["PLUTO_PROJECT"])
	#using Revise
    using Test
	using PlutoUI
	using CatmapInterface
	using Catalyst
	using ModelingToolkit
	using DifferentialEquations
	using LessUnitful
    using CairoMakie
    CairoMakie.activate!(; type = "svg", visible = false)
end;

# ╔═╡ 11ab36a0-88b5-408a-8405-88c995254081
md"""
#### Models
"""

# ╔═╡ e088298a-03b4-4678-a552-ca8b40cfc818
begin
	struct ModelInstance
	    name::String
	    path::String
	    catmap_params::CatmapParams
	    rn::Catalyst.ReactionSystem
	end
	function ModelInstance(; name, path)
	    catmap_params = parse_catmap_input(path)
	    rn            = create_reaction_network(catmap_params, conserve_pressures=true)
	    ModelInstance(name, path, catmap_params, rn)
	end
end;

# ╔═╡ bc7f352d-a5e1-4a08-a473-e0abf38410e2
const model_instances = Dict([
	    "Au-model-hbond" => ModelInstance(; name="Au-model-hbond"   , path=joinpath("..", "data", "Au-model-hbond"  , "catmap_CO2R_template.mkm")),
	    "Au-model-simple" => ModelInstance(; name="Au-model-simple"  , path=joinpath("..", "data", "Au-model-simple" , "catmap_CO2R_template.mkm")),
	    #"Liu-model-simple" => ModelInstance(; name="Liu-model-simple" , path=joinpath("..", "data", "Liu-model-simple", "catmap_CO2R_template.mkm")),
]);

# ╔═╡ fcd74eb7-329a-4c80-9fd8-1d3ce1d2060c
md"""
Model: $(@bind model_name PlutoUI.Select(collect(keys(model_instances)), default="Au-model-hbond")) 
"""

# ╔═╡ ebbe29d4-c470-4fc4-8d20-127fbdb3d09a
model_instance = model_instances[model_name];

# ╔═╡ d8bc0043-3960-420f-8237-024cb6c0e2c4
md"""
#### Structural Parameters
"""

# ╔═╡ 3942627a-621d-4e3e-ab4d-6dee080557c9
begin
	const θinit          = 1.0e-10
	const ϕ_pzc          = 0.16 * ufac"V"
	const C_gap          = 0.2 * ufac"F"
	const temp 	         = 298.0 * ufac"K"
	const S              = 9.61e-5 / ph"N_A" * (1.0e10)^2 * ufac"mol/m^2"
	surface_charge_relation(Δϕ) = round(C_gap * (Δϕ - ϕ_pzc); digits=6)
end;

# ╔═╡ f7088bbe-7553-48b4-8bca-6590222dc013
@bind params_input PlutoUI.combine() do Child
	input_params = [
		(; name="ϕ_we", range=-1.5:0.1:0.0, default=-1.5, unit="V"),
		(; name="local_pH", range=5.0:1.0:9.0, default=7.0, unit=""),
		(; name="aH2O_g", range=0.9:0.05:1.0, default=1.0, unit="")
	]
	for (s, sp) in model_instance.catmap_params.species_list
		if (isa(sp, GasSpecies) && s ≠ "H2O_g")
			push!(input_params, (; name="γ$s", range=0.2:0.2:1.4, default=1.0, unit=""))
		end
	end
	
	params_input = [
		md""" $(name) : $(Child(name, PlutoUI.Slider(range; default=default, show_value=true))) $unit
		"""
		for (; name, range, default, unit) in input_params
	]
	md"""
	#### Input Parameters:
	$(params_input)
	"""
end

# ╔═╡ 1bdb40fa-e405-4264-9f0c-d1ddca4fd8c8
@bind pressures_input PlutoUI.combine() do Child
	input_params = @NamedTuple{name::String, range::StepRangeLen, default::Float64, unit::String}[]
	for (s, sp) in model_instance.catmap_params.species_list
		if (isa(sp, GasSpecies) && s ≠ "H2O_g") || (isa(sp, FictiousSpecies) && s ≠ "ele_g")
			push!(input_params, (; name=s, range=0.0:1.0e-7:1.0e-5, default=1.0e-7, unit="bar"))
		end
	end

	pressures_input = [
		md""" $(name) : $(Child(name, PlutoUI.Slider(range; default=default, show_value=true))) $unit
		"""
		for (; name, range, default, unit) in input_params
	]
	md"""
	#### Input Pressures:
	$(pressures_input)
	"""
end

# ╔═╡ 66fbc742-4a5f-4c22-bb35-35acdcdb88a0
md"""
#### Current-Voltage Curve
"""

# ╔═╡ 6ff0e95b-ed42-4520-a29b-2fd419f1cc8e
function electrontransfer(rn::Catalyst.ReactionSystem, ssol, params)
	rxns = reactions(rn)
	irOH_g = findfirst(rxns) do r
		(; products, substrates) = r
		isempty(products) && Symbolics.tosymbol.(substrates; escape=false) == [:OH_g]
	end
	irH_g = findfirst(rxns) do r
		(; products, substrates) = r
		isempty(products) && Symbolics.tosymbol.(substrates; escape=false) == [:H_g]
	end
	curr = 0.0
	if !isnothing(irOH_g)
		curr += substitute(
			substitute(rxns[irOH_g].rate, Dict(sp => ssol[sp] for sp in species(rn))), symmap_to_varmap(rn, params)
		)
	end
	if !isnothing(irH_g)
		curr -= substitute(
			substitute(rxns[irOH_g].rate, Dict(sp => ssol[sp] for sp in species(rn))), symmap_to_varmap(rn, params)
		)
	end
	return curr * ph"N_A" * ph"e" * S
end

# ╔═╡ 84b7244f-dc0a-46d9-9f1b-947189d10111
function currentvoltage(rn, catmap_params, pressures, θ0, params; solver=DynamicSS(Rodas5P()), maxiters=1e6)
	nr = numreactions(rn)
	stoichmat = netstoichmat(rn)
	rrs_sym = reactionrates(rn)
	rrs_num = zeros(nr)

	Δϕs = 0.0:-0.1:params[:ϕ_we]
	u0 = symmap_to_varmap(rn, merge(pressures, θ0)) 
	currs = zeros(length(Δϕs))
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
		ssol = solve(ssprob, solver; maxiters)
		currs[iΔϕ] = electrontransfer(rn, ssol, params)
		u0 = Dict(sp => ssol[sp] for sp in species(rn))
	end
	return collect(Δϕs), currs
end

# ╔═╡ e284e23f-c1ab-4022-b929-d11634e71fe6
let
	(; rn, catmap_params) = model_instance
	(; species_list)      = catmap_params
	pressures = Dict(pairs(pressures_input))
	θ0 = Dict(Symbol(s) => θinit for (s, sp) in species_list if isa(sp, AdsorbateSpecies))
	params = Dict(pairs(params_input))
	(Δϕs, currs) = currentvoltage(rn, catmap_params, pressures, θ0, params)
	f = Figure()
	ax = Axis(f[1, 1]; title="Current-Voltage Curve", xlabel="Δϕ [V]", ylabel="I [mA/cm^2]", yscale=log10, limits=(-1.5, -0.5, 1.0e-10, 1.0e4))
	lines!(ax, Δϕs[currs .> 0.0], currs[currs .> 0.0] ./ (ufac"mA/cm^2"))
	f
end

# ╔═╡ Cell order:
# ╠═f5a97e3e-e116-11ee-210f-cfc211ee8cbc
# ╟─11ab36a0-88b5-408a-8405-88c995254081
# ╠═e088298a-03b4-4678-a552-ca8b40cfc818
# ╠═bc7f352d-a5e1-4a08-a473-e0abf38410e2
# ╟─fcd74eb7-329a-4c80-9fd8-1d3ce1d2060c
# ╠═ebbe29d4-c470-4fc4-8d20-127fbdb3d09a
# ╟─d8bc0043-3960-420f-8237-024cb6c0e2c4
# ╠═3942627a-621d-4e3e-ab4d-6dee080557c9
# ╟─f7088bbe-7553-48b4-8bca-6590222dc013
# ╟─1bdb40fa-e405-4264-9f0c-d1ddca4fd8c8
# ╟─66fbc742-4a5f-4c22-bb35-35acdcdb88a0
# ╠═e284e23f-c1ab-4022-b929-d11634e71fe6
# ╠═84b7244f-dc0a-46d9-9f1b-947189d10111
# ╠═6ff0e95b-ed42-4520-a29b-2fd419f1cc8e
