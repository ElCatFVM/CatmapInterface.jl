module Utils
using CatmapInterface: CatmapInterface, AdsorbateSpecies, FictiousSpecies, GasSpecies
using DelimitedFiles: DelimitedFiles, readdlm, writedlm
using Format: Format, cfmt
using PyCall: PyCall, @py_str, @pyinclude, keys, pyimport

const Symmap    = Vector{Pair{Symbol, Float64}}
const SSParams  = @NamedTuple{u0::Symmap, ps::Symmap}
const temp      = 298.0

if haskey(ENV, "CATMAP")
    pushfirst!(pyimport("sys")."path", ENV["CATMAP"])
end

begin ## helper functions for CatMAP testing
    """
    instantiate_catmap_template!(instance_file_path::String, template_file_path::String, params, T::Float64)
    
    Instantiate a template file by inserting the parameters in the `params`.
    """
    function instantiate_catmap_template!(instance_file_path::String, template_file_path::String, params::SSParams; T=temp)
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
    
    begin ## CatMAP energy calculations
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
        
        function compute_catmap_free_energies!(free_energies::Dict{String, T}, catmap_template_path::String, params::CatmapInterface.InterfaceParams) where {T <: Real}
            # instantiate template
            (dname, fname) = splitdir(catmap_template_path)
            catmap_instance_path = joinpath(dname,
                replace(fname, 
                    r"(?<pre>.*)_template(?<post>.*)" => s"\g<pre>_instance\g<post>"
                )
            )
            CatmapInterface.instantiate_catmap_template!(catmap_instance_path, catmap_template_path, params, temp)

            # compute free energies
            (; θ) = params
            starting_dir = pwd()
            catmap_instance_path = abspath(catmap_instance_path)
            cd(dirname(catmap_instance_path))
            local tmp
            try
                tmp = catmap_kinetic_model(catmap_instance_path, θ)
            finally
                cd(starting_dir)
                rm(catmap_instance_path)
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
            for intparams in params
                free_energies[intparams] = Dict{String, Float64}()
                compute_catmap_free_energies!(free_energies[intparams], catmap_template_path, intparams)
            end
        end 
    end

    begin ## CatMAP steady state calculations and helpers
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
		py"""
		import numpy as np
		$$(read(splitdir(logfile_path)[end], String))
		"""
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
        catmap_ssolve(template_file_path::String, params::SSParams)
    
        Solve the microkinetic model specified in `template_file_path` for the steady state using CatMAP for each parameter set in `params_iter` and add it to `ssols`. 
        """
        function catmap_ssolve(template_file_path::String, params::SSParams)
            @assert isfile(template_file_path)
            instance_file_path = joinpath(dirname(template_file_path), "test.mkm")
            instantiate_catmap_template!(
                instance_file_path, 
                template_file_path, 
                params
            )
            local ssol
            try
                runcatmap(instance_file_path)
                logfile = splitext(instance_file_path)[1] * ".log"
                @pyinclude(instance_file_path)
                datafile = joinpath(dirname(instance_file_path), py"data_file")
                try
                    ssol = get_coverage_map(logfile)
                finally
                    rm(logfile)
                    rm(datafile)
                end
            finally
                rm(instance_file_path)
            end
            return ssol
        end 
    end
end


begin ## some helper functions for writing and loading the test parameters
    surface_charge_relation = let
        ϕ_pzc = 0.16
        C_gap = 0.2
        Δϕ->round(C_gap * (Δϕ - ϕ_pzc); digits=6)
    end
    
    
    function write_test_params(test_params_path, catmap_params; surface_charge_relation, ϕ_we_iter=-1.0:0.5:0.0, ϕ_iter=-1.0:0.5:0.0, local_pH_iter=6.0:1.0:8.0, u0gas=1.0, θ0=1.0e-10)
        (; species_list) = catmap_params
        
        params_list = Dict{Symbol, Float64}[]
        for (ϕ_we, ϕ, local_pH) in Iterators.product(ϕ_we_iter, ϕ_iter, local_pH_iter)
            params = [
                :ϕ_we     => ϕ_we,
                :ϕ        => ϕ,
                :local_pH => local_pH,
                :aH2O_g => 1.0,
                :σ => surface_charge_relation(ϕ_we - ϕ)
            ]
            for (s, sp) in species_list
                if isa(sp, GasSpecies) && s ≠ "H2O_g"
                    push!(params, Symbol("γ$s") => 1.0)
                    push!(params, Symbol(s) => u0gas)
                elseif isa(sp, FictiousSpecies) && s ≠ "ele_g"
                    push!(params, Symbol(s) => u0gas)
                elseif isa(sp, AdsorbateSpecies)
                    push!(params, Symbol(s) => θ0)
                end
            end
            push!(params_list, Dict(params))
        end
        first_params = first(params_list)
        header = [p_sym for p_sym in keys(first_params)]
        open(test_params_path, "w") do file
            writedlm(file, [[header]; [[params[sym] for sym in header] for params in params_list]])
        end
    end
    
    function load_test_params(test_params_path, catmap_params; only_interface_params=false)
        @assert isfile(test_params_path)
        
        (data_cells, header_cells) = readdlm(test_params_path, Float64; header=true)
        if any(isnan, data_cells)
            throw(ArgumentError("Not all values in the file at $(last(splitdir(test_params_path))) are numeric"))        
        end
        header_cells = Symbol.(header_cells)
        (; species_list) = catmap_params
        p_syms  = [:ϕ_we, :ϕ, :local_pH, :aH2O_g]
        if catmap_params.electrochemical_thermo_mode == :hbond_surface_charge_density
            push!(p_syms, :σ)
        end
        u0_syms = Symbol[]
        θ_syms  = Symbol[]
        for (s, sp) in species_list
            if isa(sp, GasSpecies) && s ≠ "H2O_g"
                push!(p_syms, Symbol("γ$s"))
                push!(u0_syms, Symbol(s))
            elseif isa(sp, FictiousSpecies) && s ≠ "ele_g"
                push!(u0_syms, Symbol(s))
            elseif isa(sp, AdsorbateSpecies)
                push!(u0_syms, Symbol(s))
                push!(θ_syms, Symbol(s))
            end
        end
        for p_sym in p_syms
            if p_sym ∉ header_cells
                throw(ArgumentError("The column for the parameter $(string(p_sym)) is missing."))
            end
        end
        for u0_sym in u0_syms
            if u0_sym ∉ header_cells
                throw(ArgumentError("The column for the initial valye of $(string(u0_sym)) is missing."))
            end
        end
    
        local test_params
        if only_interface_params == true
            (ϕ_we, ϕ, σ, local_pH) = zeros(4)
            test_params = CatmapInterface.InterfaceParams[]
            for row in eachrow(data_cells)
                θ = Pair{String, Float64}[]
                for (sym, val) in zip(header_cells, row)
                    if sym in θ_syms
                        push!(θ, string(sym) => val)
                    elseif sym == :ϕ_we
                        ϕ_we = val
                    elseif sym == :ϕ
                        ϕ = val
                    elseif sym == :σ
                        σ = val
                    elseif sym == :local_pH
                        local_pH = val
                    end
                end
                push!(test_params, CatmapInterface.InterfaceParams(; ϕ_we, ϕ, σ, local_pH, θ=Dict(θ)))
            end
        else
            test_params = SSParams[]
            for row in eachrow(data_cells)
                u0 = Symmap()
                ps = Symmap()
                for (sym, val) in zip(header_cells, row)
                    if sym in p_syms
                        push!(ps, sym => val)
                    elseif sym in u0_syms
                        push!(u0, sym => val)
                    end
                end
                push!(test_params, SSParams((u0, ps)))
            end 
        end
        test_params
    end    
end 

begin # product iterators
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
end


begin # helpful Base extensions
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

    function Base.convert(::Type{CatmapInterface.InterfaceParams}, s::String)
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
export Symmap, SSParams
export instantiate_catmap_template!, catmap_ssolve
export compute_catmap_free_energies!
export InterfaceParamsProductIterator
export θProductIterator
export load_test_params
export listmodels

end;
