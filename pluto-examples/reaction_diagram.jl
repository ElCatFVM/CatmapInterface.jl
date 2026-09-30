### A Pluto.jl notebook ###
# v0.20.24

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

# ╔═╡ 252a1d28-d7eb-11ee-1f6a-990ae6a3184f
begin
    # using Revise
    using Pkg
    Pkg.activate(@__DIR__)
	using Revise
    using Test
	using PlutoUI
	using Format
	using CatmapInterface
	using LessUnitful
    using CairoMakie
    CairoMakie.activate!(; type = "svg", visible = false)
end;

# ╔═╡ af65d97c-dae1-4ea8-8a95-58bbd2ffd80a
begin
	include("Utils.jl")
	using .Utils
end

# ╔═╡ b591ebbf-e98a-4e79-abb3-77b523c6e652
md"""
#### Improve Printing of Parsed Reactions
"""

# ╔═╡ 667684f7-3835-4a4d-80ab-d9229b9c1898
begin
	function Base.show(io::Base.IO, reaction::CatmapInterface.ParsedReaction)
		educt_string = ""
		for (educt, c) in reaction.educts
			educt_string *= "$(c > 1 ? c : "") $educt + "
		end
	
		product_string = ""
		for (product, c) in reaction.products
			product_string *= "$(c > 1 ? c : "") $product + "
		end
	
		if !isnothing(reaction.tstate)
			tstate_string = ""
			for (s, c) in reaction.tstate.components
				tstate_string *= "$(c > 1 ? c : "") $s + "
			end
			print(io, educt_string[1:end-3] * " <-> " * tstate_string[1:end-3] * " <-> " * product_string[1:end-3])
		else
			print(io, educt_string[1:end-3] * " <-> " * product_string[1:end-3])
		end
		nothing
	end
	CatmapInterface.ParsedReaction(["H2O_g" => 1, "ele_g" => 1, "_t" => 1], ["H_t" => 1, "OH_g" => 1], CatmapInterface.TState(["H2O-ele_t" => 1],nothing, 0.5))
end

# ╔═╡ 617f2253-a80c-4feb-b292-47337bd36752
md"""
#### Models
"""

# ╔═╡ 1cfbb9ef-8e32-4482-bf70-c5a63427d518
function listmodels(datadir)
	ismodeldir(f) = isdir(joinpath(datadir, f)) && occursin("model", f)
	ismodeltemplate(f) = splitext(f)[end] == ".mkm"
	modeldirs = filter(ismodeldir, readdir(datadir))
	modeltemplates = map(modeldirs) do d
		filter(ismodeltemplate, readdir(joinpath(datadir, d); join=true))[end]
	end
	return Dict(zip(modeldirs, modeltemplates))
end

# ╔═╡ 2273a044-adb2-45b1-b166-88b47f30ca68
const models = listmodels("../data")

# ╔═╡ 5b294ca8-d63f-4277-ad20-5327e418a219
md"""
#### Structural Parameters
"""

# ╔═╡ bb3b204b-5b55-4a70-ab14-ac4de457e563
begin
	const θ 		= 0.2:0.2:0.2
	const ϕ_we 		= -1.0:0.5:0.0
	const ϕ 		= -1.0:0.5:0.0
	const local_pH 	= 6.0:2.0:8.0
	const ϕ_pzc 	= 0.16 * ufac"V"
	const C_gap 	= 0.2 * ufac"F"
	const temp 		= 298 * ufac"K"
end;

# ╔═╡ 66ed4181-1393-4ed0-9555-b3608b01f223
surface_charge_relation(Δϕ) = round(C_gap * (Δϕ - ϕ_pzc); digits=6)

# ╔═╡ bd93f2e2-9368-4ea1-a58b-80ac76e15f59
md"""
#### Model
$(@bind model_name PlutoUI.Select(collect(keys(models))))
"""

# ╔═╡ 29d5bcf8-6332-4ab6-91ed-9d6cf30e5421
function energyplot(catmap_params, free_energies, ireaction)
	function get_free_energy(state)
		Gf = 0.0
		for (species, c) in state
			Gf += c * free_energies[species]
		end
		return Gf
	end
	function get_label(state)
		label = ""
		for (species, c) in state
			label *= (c > 0) ? "$(c > 1 ? c : "") $species + " : "" 
		end
		return label[1:end-3]
	end

	(; educts, products, tstate) = catmap_params.reactions[ireaction]

	Gfs = Float64[]
	labels = String[]
	hlines = Tuple{Float64, Float64}[]
	
	# educts
	push!(Gfs, get_free_energy(educts))
	push!(labels, get_label(educts))
	push!(hlines, (0.0, 0.33))
	# tstate
	if !isnothing(tstate)
		push!(Gfs, get_free_energy(tstate.components))
		push!(labels, get_label(tstate.components))
		push!(hlines, (0.33, 0.67))
	end
	# products
	push!(Gfs, get_free_energy(products))
	push!(labels, get_label(products))
	push!(hlines, (0.66, 1.0))

	f = Figure(; size=(800, 300))
	ax = Axis(f[1, 1],
	    title = "Energy Diagram",
	    xlabel = "Reaction Path",
	    ylabel = "Free Energy [eV]",
		limits = (0, 1, -1, 5)
	)
	Gfs ./= ufac"eV"
	hlines!(ax, Gfs, xmin=first.(hlines), xmax=last.(hlines))
	text!(ax,first.(hlines).+0.02, Gfs, text=labels)
	f
end

# ╔═╡ 17b3ef9b-fa14-4bd8-841a-ad3d55bfeacf
function reformulate(params_input)
	(; ϕ_we, ϕ, local_pH) = params_input
	θdict = Dict{String, Float64}()
	for (k, v) in pairs(params_input)
		k = String(k)
		if k[1] == 'θ'
			θdict[k[3:end]] = v
		end
	end
	CatmapInterface.InterfaceParams(; θ=θdict, ϕ_we=ϕ_we, ϕ=ϕ, σ=surface_charge_relation(ϕ_we-ϕ), local_pH=local_pH)
end

# ╔═╡ 4e20776d-261a-40c8-9720-dc151de9599f
const catmap_params_dict = Dict([
	model_name => parse_catmap_input(model_template_path)
	for (model_name, model_template_path) in models
])

# ╔═╡ b983297c-697d-4058-a2ba-0003bd55d8dd
md"""
#### Reaction:
$(@bind ireaction PlutoUI.Select([i => r for (i, r) in enumerate(catmap_params_dict[model_name].reactions)]))
"""

# ╔═╡ c765b5a4-e706-4151-83b7-a27ee59ea068
@bind params_input PlutoUI.combine() do Child
	input_params = [
		(; name="ϕ_we", range=ϕ_we, unit="V"),
		(; name="ϕ", range=ϕ, unit="V"),
		(; name="local_pH", range=local_pH, unit="")
	]
	catmap_params = catmap_params_dict[model_name]
	ads_species = CatmapInterface.adsorbatespecies(catmap_params.species_list)
	append!(input_params, [
		(; name="θ$(first(p))", range=θ, unit="")
		for p in ads_species
	])

	params_input = [
		md""" $(name) : $(Child(name, PlutoUI.Slider(range, show_value=true))) $unit
		"""
		for (; name, range, unit) in input_params
	]
	md"""
	#### Input Parameters:
	$(params_input)
	"""
end

# ╔═╡ 01e29eca-6a92-4b32-a6a4-79248402a87a
const params_dict = let
	ps = []
	for (model_name, catmap_params) in catmap_params_dict
		ads_species = CatmapInterface.adsorbatespecies(catmap_params.species_list)
		θiter = Utils.θProductIterator(
			ads_name => copy(θ) for ads_name in keys(ads_species)
		)
		push!(ps, model_name => 
			Utils.InterfaceParamsProductIterator(; θ=θiter, ϕ_we=copy(ϕ_we), ϕ=copy(ϕ), 						local_pH=copy(local_pH), surface_charge_relation)
		)
	end
	Dict(ps)
end

# ╔═╡ 04dab72e-a847-4470-816b-4c94931019fb
md"""
### Free Energy Calculations
"""

# ╔═╡ 9b06dd0d-721a-4947-b749-f59d722d3420
begin 
	const free_energies_dict = Dict(
		model_name => Dict{CatmapInterface.InterfaceParams, Dict{String, Float64}}()
		for model_name in keys(models)
	)
	for model_name in keys(models)
		catmap_params 	= catmap_params_dict[model_name]
		params 			= params_dict[model_name]
		CatmapInterface.compute_free_energies!(
			free_energies_dict[model_name], catmap_params, params
		)
	end
end

# ╔═╡ e52df405-05f3-447d-907a-d1106aa47b15
let
	catmap_params = catmap_params_dict[model_name]
	free_energies = free_energies_dict[model_name][reformulate(params_input)]
	energyplot(catmap_params, free_energies, ireaction)
end

# ╔═╡ 0b857ada-a26a-4637-8f57-8262a0c13ca2
md"""
### Tests
"""

# ╔═╡ 6fe83143-ae95-426b-8242-bde4f830fb07
begin
	const catmap_free_energies_dict = Dict(
		model_name => Dict{CatmapInterface.InterfaceParams, Dict{String, Float64}}()
		for model_name in keys(models)
	)
	for (model_name, catmap_template_path) in models
		params = params_dict[model_name]
		Utils.compute_catmap_free_energies!(
			catmap_free_energies_dict[model_name], catmap_template_path, params
		)
	end
end

# ╔═╡ 79ba095b-dc84-44c5-863e-d7d70322254d
md"""
__Relative tolerance__: $(@bind rtol PlutoUI.Select([1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4,1.0e-5, 1.0e-6]))
"""

# ╔═╡ a656e35a-9184-4241-8aa9-301379bcf728
function runtests(models, free_energies_dict, catmap_free_energies_dict; rtol=1.0e-5)
	@testset "model=$model_name" for model_name in keys(models)
		free_energies_ps = free_energies_dict[model_name]
		catmap_free_energies_ps = catmap_free_energies_dict[model_name]
		params = params_dict[model_name]
		@testset "$(convert(String, intparams))" for intparams in params
			free_energies = free_energies_ps[intparams]
			catmap_free_energies = catmap_free_energies_ps[intparams]
			@testset "species=$species" for (species, free_energy) in free_energies
		        @test isapprox(free_energy/ufac"eV", catmap_free_energies[species]; rtol)
		    end 
		end
	end
	nothing
end

# ╔═╡ f6a9addf-c580-4ae9-b9c2-0b978d6216b9
runtests(models, free_energies_dict, catmap_free_energies_dict; rtol)

# ╔═╡ Cell order:
# ╠═252a1d28-d7eb-11ee-1f6a-990ae6a3184f
# ╠═af65d97c-dae1-4ea8-8a95-58bbd2ffd80a
# ╟─b591ebbf-e98a-4e79-abb3-77b523c6e652
# ╟─667684f7-3835-4a4d-80ab-d9229b9c1898
# ╟─617f2253-a80c-4feb-b292-47337bd36752
# ╟─1cfbb9ef-8e32-4482-bf70-c5a63427d518
# ╠═2273a044-adb2-45b1-b166-88b47f30ca68
# ╟─5b294ca8-d63f-4277-ad20-5327e418a219
# ╠═bb3b204b-5b55-4a70-ab14-ac4de457e563
# ╠═66ed4181-1393-4ed0-9555-b3608b01f223
# ╟─bd93f2e2-9368-4ea1-a58b-80ac76e15f59
# ╟─29d5bcf8-6332-4ab6-91ed-9d6cf30e5421
# ╟─e52df405-05f3-447d-907a-d1106aa47b15
# ╟─b983297c-697d-4058-a2ba-0003bd55d8dd
# ╟─c765b5a4-e706-4151-83b7-a27ee59ea068
# ╟─17b3ef9b-fa14-4bd8-841a-ad3d55bfeacf
# ╠═4e20776d-261a-40c8-9720-dc151de9599f
# ╟─01e29eca-6a92-4b32-a6a4-79248402a87a
# ╟─04dab72e-a847-4470-816b-4c94931019fb
# ╠═9b06dd0d-721a-4947-b749-f59d722d3420
# ╟─0b857ada-a26a-4637-8f57-8262a0c13ca2
# ╠═6fe83143-ae95-426b-8242-bde4f830fb07
# ╟─79ba095b-dc84-44c5-863e-d7d70322254d
# ╟─a656e35a-9184-4241-8aa9-301379bcf728
# ╠═f6a9addf-c580-4ae9-b9c2-0b978d6216b9
