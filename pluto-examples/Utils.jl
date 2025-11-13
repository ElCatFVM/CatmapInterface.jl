module Utils
using PyCall
using CatmapInterface
using ModelingToolkit
using Catalyst
using OrdinaryDiffEqRosenbrock
using SteadyStateDiffEq
using Format

const SymmapType    = Vector{Pair{Symbol, Float64}}
const SSParamsType  = @NamedTuple{u0::SymmapType, ps::SymmapType}
const temp = 298.0

if haskey(ENV, "CATMAP")
    pushfirst!(pyimport("sys")."path", ENV["CATMAP"])
end

"""
instantiate_catmap_template!(instance_file_path::String, template_file_path::String, params, T::Float64)

Instantiate a template file by inserting the parameters in the `params`.
"""
function instantiate_catmap_template!(instance_file_path::String, template_file_path::String, params::SSParamsType, T)
    (; u0, ps) = params
    ps = Dict(ps)
    ϕ_we, ϕ, local_pH = (ps[:ϕ_we], ps[:ϕ], ps[:local_pH])
    σ = get(ps, :σ, nothing)
    instance_string = open(template_file_path, "r") do template_file
        read(template_file, String)
    end

	replacements = [
		r"descriptor_ranges.?=.*" =>SubstitutionString("descriptor_ranges = [[$ϕ_we, $ϕ_we], [$T, $T]]"),
		r"voltage_diff_drop.?=.*" => SubstitutionString("voltage_diff_drop = $ϕ"),
		r"pH.?=.*" => SubstitutionString("pH = $local_pH"),
	]
    if !isnothing(σ)
        push!(replacements, r"\nsigma_input.?=.*" => SubstitutionString("\\nsigma_input = $(σ/0.01)")) # in μF/cm^2
    end

    for (sym, val) in u0
        s = string(sym)
        push!(replacements, Regex("^species_definitions\\[['|\"]\\Q$s\\E['|\"]\\].?=.*?{(?<before>.*?)['|\"]pressure['|\"].*?:.*?[0-9e\\-\\.]*(?<after>.*?)}", "m") => SubstitutionString("species_definitions['$s'] = {\\g<before>'pressure':$(val)\\g<after>}"))
    end

    instance_string = replace(instance_string, replacements...)

    open(instance_file_path, "w") do instance_file
        write(instance_file, instance_string)
    end

    return instance_file_path
end



py"""
from catmap import ReactionModel
import catmap

def catmap_kinetic_model(setup_file, theta):

    # ReactionModel is the main class that is initialized with a setup-file
    model = ReactionModel(setup_file=setup_file)

    # some solver parameters have to be set manually (?!)
    import mpmath as mp
    model.solver._mpfloat = mp.mpf
    model.solver._math = mp
    model.solver._matrix = mp.matrix

    model.solver.compile() # compiles all templates, here (rate_constants) are needed

    # Set up interaction model.
    if model.adsorbate_interaction_model == 'first_order':
        interaction_model = \
            catmap.thermodynamics.FirstOrderInteractions(model)
        interaction_model.get_interaction_info()
        response_func = interaction_model.interaction_response_function
        if not callable(response_func):
            int_function = getattr(interaction_model,
                                    response_func+'_response')
            interaction_model.interaction_response_function = int_function
        model.thermodynamics.__dict__['adsorbate_interactions'] = interaction_model
    elif model.adsorbate_interaction_model in ['ideal',None]:
        model.thermodynamics.adsorbate_interactions = None
    else:
        raise AttributeError(
                'Invalid adsorbate_interaction_model specified.')

    descriptor_values = [descriptor_range[0] for descriptor_range in model.descriptor_ranges]
    rxn_parameters = model.scaler.get_rxn_parameters(descriptor_values)
    n_tot = len(model.adsorbate_names) + len(model.transition_state_names)
    energies = rxn_parameters[:n_tot]
    if len(rxn_parameters) == n_tot + n_tot**2:
        interaction_vector = rxn_parameters[-n_tot**2:]
    elif len(rxn_parameters) == n_tot:
        interaction_vector = [0]*n_tot**2
    F = model.interaction_response_function if hasattr(model, 'interaction_response_function') else None
    theta = [theta[species] for species in model.adsorbate_names] + len(model.transition_state_names) * [0.0]
    _, Gf, _ = model.interaction_function(theta, energies, interaction_vector, F)

    energies = {
        k:v for k, v in zip(
            model.adsorbate_names + model.transition_state_names,
            Gf
        )
    }
    for k, v in zip(model.gas_names, model.solver._gas_energies):
        energies[k] = v
    for k, v in zip(model.site_names, model.solver._site_energies):
        energies[k] = v

    return energies
"""
catmap_kinetic_model = py"catmap_kinetic_model"

function compute_catmap_free_energies!(free_energies, catmap_instance_path, params)
    (; θ) = params
    starting_dir = pwd()
    catmap_instance_path = abspath(catmap_instance_path)
    cd(dirname(catmap_instance_path))
    local tmp
	try
        tmp = catmap_kinetic_model(catmap_instance_path, θ)
	finally
        cd(starting_dir)
    end
	for (k, v) in tmp
		if length(k) == 1 && k != "g"
			free_energies["_$k"] = convert(Float64, v)
		else
			free_energies[k] = convert(Float64, v)
		end
	end
    nothing
end

function compute_catmap_free_energies!(
	free_energies::Dict{CatmapInterface.InterfaceParams, Dict{String, T}}, 
	catmap_template_path, 
	params
) where {T <: Real}
	(dname, fname) = splitdir(catmap_template_path)
	catmap_instance_path = joinpath(dname,
		replace(fname, 
			r"(?<pre>.*)_template(?<post>.*)" => s"\g<pre>_instance\g<post>"
		)
	)
	for intparams in params
		CatmapInterface.instantiate_catmap_template!(catmap_instance_path, catmap_template_path, intparams, temp)
		free_energies[intparams] = Dict{String, Float64}()
		compute_catmap_free_energies!(free_energies[intparams], catmap_instance_path, intparams)
		rm(catmap_instance_path)
	end
end

# Steady State

## Catmap

py"""
from catmap import ReactionModel

def runcatmap(setup_file):
	model = ReactionModel(setup_file=setup_file)
	return model.run()
"""
	
"""
runcatmap(instance_file_path::String)

Run CatMAP on the microkinetic model at `instance_file_path`.
"""
function runcatmap(instance_file_path)
    @assert isfile(instance_file_path)
    currdir = pwd()
    cd(dirname(instance_file_path))
    try
        py"runcatmap"(splitdir(instance_file_path)[end])
    finally
        cd(currdir)
    end
end

"""
get_coverage_map(logfile_path::String)

Get the coverage map from CatMAP's logfile at `logfile_path`
"""
function get_coverage_map(logfile_path)
    @assert isfile(logfile_path)
	currdir = pwd()
	newdir = dirname(logfile_path)
	cd(newdir)
	local cmap
	try
		@pyinclude(splitdir(logfile_path)[end])
		labels = Symbol.(py"output_labels"["coverage"])
		coverages = py"coverage_map"[2]
		coverages = py"float".(coverages)
		cmap = Dict(zip(labels, coverages))
	finally
		cd(currdir)
	end
	cmap
end

"""
ssolve!(ssols::Dict{String, Dict{Utils.SSParamsType, SciMLBase.NonlinearSolution}}, template_file_path::String, params_iter)

Solve the microkinetic model specified in `template_file_path` for the steady state using CatMAP for each parameter set in `params_iter` and add it to `ssols`. 
"""
function ssolve!(ssols, template_file_path::String, params_iter)
    @assert isfile(template_file_path)
	for params in params_iter
		instance_file_path = joinpath(dirname(template_file_path), "test.mkm")
		instantiate_catmap_template!(
			instance_file_path, 
			template_file_path, 
			params, 
			temp
		)
		try
			runcatmap(instance_file_path)
			logfile = splitext(instance_file_path)[1] * ".log"
			@pyinclude(instance_file_path)
			datafile = joinpath(dirname(instance_file_path), py"data_file")
			try
				ssols[params] = get_coverage_map(logfile)
			finally
				rm(logfile)
				rm(datafile)
			end
		finally
			rm(instance_file_path)
		end
	end
end

## CatmapInterface



"""
ssolve!(ssols::Dict{String, Dict{Utils.SSParamsType, SciMLBase.NonlinearSolution}}, odesys::ModelingToolkit.ODESystem, params_iter)

Solve the `odesys` for the steady state for each parameter set in `params_iter` and add it to `ssols`. 
"""
function ssolve!(ssols, odesys::ModelingToolkit.ODESystem, params_iter)
	for params in params_iter
		(; u0, ps) = params
		ssprob = SteadyStateProblem(
			odesys, 
			symmap_to_varmap(odesys, u0), 
			symmap_to_varmap(odesys, ps)
		)
		ssols[params] = solve(ssprob, DynamicSS(Rodas5P()); maxiters=1e6)
	end
end


# Product Iterators

## InterfaceParamsProductIterator
struct InterfaceParamsProductIterator
    prod
    surface_charge_relation
    function InterfaceParamsProductIterator(; θ, ϕ_we, ϕ, local_pH, surface_charge_relation)
        new(Iterators.product(θ, ϕ_we, ϕ, local_pH), surface_charge_relation)
    end
end
function Iterators.iterate(iter::InterfaceParamsProductIterator)
    o = Iterators.iterate(iter.prod)
    if isnothing(o)
        return nothing
    end
    (first_item, initial_state) = o
    (θ, ϕ_we, ϕ, local_pH) = first_item
    σ = iter.surface_charge_relation(ϕ_we - ϕ)
    return CatmapInterface.InterfaceParams(; θ, ϕ_we, ϕ, σ, local_pH), initial_state
end
function Iterators.iterate(iter::InterfaceParamsProductIterator, state)
    next = Iterators.iterate(iter.prod, state)
    if isnothing(next)
        return nothing
    end
    (next_item, next_state) = next
    (θ, ϕ_we, ϕ, local_pH) = next_item
    σ = iter.surface_charge_relation(ϕ_we - ϕ)
    return CatmapInterface.InterfaceParams(θ, ϕ_we, ϕ, σ, local_pH), next_state
end
Base.eltype(iter::InterfaceParamsProductIterator) = CatmapInterface.InterfaceParams
Base.length(iter::InterfaceParamsProductIterator) = length(iter.prod)

## θProductIterator
struct θProductIterator
    prod
    θnames
    θProductIterator(θs) = new(Iterators.product(last.(θs)...), first.(θs))
end
function Base.iterate(iter::θProductIterator)
    o = Iterators.iterate(iter.prod)
    if isnothing(o)
        return nothing
    end
    (first_item, initial_state) = o
    return Dict(zip(iter.θnames, first_item)), initial_state
end
function Base.iterate(iter::θProductIterator, state)
    next = Iterators.iterate(iter.prod, state)
    if isnothing(next)
        return nothing
    end
    (next_item, next_state) = next
    return Dict(zip(iter.θnames, next_item)), next_state
end
Base.eltype(iter::θProductIterator) = Dict{String, eltype(iter.prod)}
Base.length(iter::θProductIterator) = length(iter.prod)


# Helpful Stuff

function Base.convert(::Type{String}, intparams::CatmapInterface.InterfaceParams)
	(; θ, ϕ_we, ϕ, σ, local_pH) = intparams
	fmt = "%.5f"
	θstr = mapreduce(*, θ) do kv
		(k, v) = first(kv), last(kv)
		"θ$k=$(cfmt(fmt, v)),"
	end
	otherstr = mapreduce(*,["ϕ_we"=>ϕ_we, "ϕ"=>ϕ, "σ"=>σ, "local_pH"=>local_pH]) do kv
		(k, v) = first(kv), last(kv)
		"$k=$(cfmt(fmt, v)),"
	end
	return θstr * otherstr[1:end-1]
end

function Base.convert(::Type{CatmapInterface.InterfaceParams}, s)
	ps = split.(split(s, ","), "=")
	ps = Dict(pname => parse(Float64, pvalue) for (pname, pvalue) in ps)
	θdict = Dict{String, Float64}()
	for (k, v) in ps
		k = String(k)
		if k[1] == 'θ'
			θdict[k[3:end]] = v
		end
	end
	CatmapInterface.InterfaceParams(; θ=θdict, ϕ_we=ps["ϕ_we"], ϕ=ps["ϕ"], σ=ps["σ"], local_pH=ps["local_pH"])
end

begin
	function Base.isequal(a::CatmapInterface.InterfaceParams{T}, b::CatmapInterface.InterfaceParams{T}) where {T <: Real}
		a.θ == b.θ && a.ϕ_we == b.ϕ_we && a.ϕ == b.ϕ && a.local_pH == b.local_pH 
	end
	function Base.hash(a::CatmapInterface.InterfaceParams{T}) where {T <: Real}
		hash((; θ=a.θ, ϕ_we=a.ϕ_we, ϕ=a.ϕ, σ=a.σ, local_pH=a.local_pH))
	end
end

# Data Loading
function listmodels(datadir)
	ismodeldir(f) = isdir(joinpath(datadir, f)) && occursin("model", f)
	ismodeltemplate(f) = splitext(f)[end] == ".mkm"
	modeldirs = filter(ismodeldir, readdir(datadir))
	modeltemplates = map(modeldirs) do d
		filter(ismodeltemplate, readdir(joinpath(datadir, d); join=true))[end]
	end
	return Dict(zip(modeldirs, modeltemplates))
end


# Exports 
export compute_catmap_free_energies!
export InterfaceParamsProductIterator
export θProductIterator
export listmodels
export conserve_pressures!, ssolve!

end;
