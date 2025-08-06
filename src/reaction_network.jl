
"""
$(TYPEDEF)

Interface parameters in the free energy model

$(TYPEDFIELDS)
"""

@kwdef struct InterfaceParams{T <: Real}
    """
    Coverages of the adsorbates
    """
    θ::Dict{String, T}
    """
    Electric potential at the working electrode in V
    """
    ϕ_we::T
    """
    Electric potential at the reaction plane in V
    """
    ϕ::T
    """
    Surface charge density on the electrode in Cm⁻²
    """
    σ::T
    """
    pH-value at the reaction plane
    """
    local_pH::T
end

"""
$(SIGNATURES)

Compute the Gibbs free energies of all species specified in the `catmap_params` by applying the specified correction modes.
"""
function compute_free_energies!(free_energies, catmap_params::CatmapParams, formation_energies, θ, σ, ϕ_we, ϕ, local_pH, β = Dict([s => sp.β for (s, sp) in catmap_params.species_list if isa(sp, TStateSpecies)]))
    (; adsorbate_interaction_params, gas_thermo_mode, adsorbate_thermo_mode, electrochemical_thermo_mode) = catmap_params
    (; adsorbate_interaction_model) = adsorbate_interaction_params

    for (s,formation_energy) in formation_energies
        free_energies[s] += formation_energy
    end

    adsorbate_interaction_correction!  = getfield(@__MODULE__, Symbol(adsorbate_interaction_model, "_adsorbate_interaction"))
    gas_thermo_correction!             = getfield(@__MODULE__, gas_thermo_mode)
    adsorbate_thermo_correction!       = getfield(@__MODULE__, adsorbate_thermo_mode)
    electrochemical_thermo_correction! = getfield(@__MODULE__, electrochemical_thermo_mode)

    adsorbate_interaction_correction!(free_energies, catmap_params, θ)

    thermo_corrections = Dict(zip(keys(free_energies), zeros(valtype(free_energies), length(free_energies))))
    gas_thermo_correction!(thermo_corrections, catmap_params)
    adsorbate_thermo_correction!(thermo_corrections, catmap_params)
    # electrochemical corrections
    electrochemical_thermo_correction!(thermo_corrections, catmap_params, σ, ϕ_we, ϕ, local_pH, β)
    for (species, thermo_correction) in thermo_corrections
        free_energies[species] += thermo_correction
    end
    nothing
end

"""
$(SIGNATURES)

Compute the Gibbs free energies of all species specified in the `catmap_params` by applying the specified correction modes, but without ad-ad interaction
"""
function compute_free_energies_only_ad_TS!(free_energies, catmap_params::CatmapParams, barriers, θ, σ, ϕ_we, ϕ, local_pH, β = Dict([s => sp.β for (s, sp) in catmap_params.species_list if isa(sp, TStateSpecies)]))
    (;adsorbate_interaction_params, species_list) = catmap_params
    (;adsorbate_interaction_model) = adsorbate_interaction_params
    
    for (s, barrier) in barriers
        free_energies[s] += barrier ## formation_energy of TS shoulbe 0 or maybe it can be the value of barrier.
    end

    adsorbate_interaction_correction!  = getfield(@__MODULE__, Symbol(adsorbate_interaction_model, "_adsorbate_interaction"))
    adsorbate_interaction_correction!(free_energies, catmap_params, θ)

    nothing
end



"""
$(SIGNATURES)

Compute the Gibbs free energies of all species specified in the `catmap_params` by applying the specified correction modes.
"""
function compute_free_energies!(free_energies::Dict{String, T}, catmap_params::CatmapParams, interface_params::InterfaceParams) where {T <: Real}
    (; θ, σ, ϕ_we, ϕ, local_pH) = interface_params
    (; adsorbate_interaction_params, gas_thermo_mode, adsorbate_thermo_mode, electrochemical_thermo_mode) = catmap_params
    (; adsorbate_interaction_model) = adsorbate_interaction_params

    for (s, (; formation_energy)) in catmap_params.species_list
        free_energies[s] += formation_energy
    end
    β = Dict([s => sp.β for (s, sp) in catmap_params.species_list if isa(sp, TStateSpecies)])

    adsorbate_interaction_correction!  = getfield(@__MODULE__, Symbol(adsorbate_interaction_model, "_adsorbate_interaction"))
    gas_thermo_correction!             = getfield(@__MODULE__, gas_thermo_mode)
    adsorbate_thermo_correction!       = getfield(@__MODULE__, adsorbate_thermo_mode)
    electrochemical_thermo_correction! = getfield(@__MODULE__, electrochemical_thermo_mode)

    adsorbate_interaction_correction!(free_energies, catmap_params, θ)

    thermo_corrections = Dict(zip(keys(free_energies), zeros(valtype(free_energies), length(free_energies))))
    gas_thermo_correction!(thermo_corrections, catmap_params)
    adsorbate_thermo_correction!(thermo_corrections, catmap_params)
    # electrochemical corrections
    electrochemical_thermo_correction!(thermo_corrections, catmap_params, σ, ϕ_we, ϕ, local_pH, β)
    for (species, thermo_correction) in thermo_corrections
        free_energies[species] += thermo_correction
    end
    nothing
end

"""
$(SIGNATURES)

Compute the Gibbs free energies of all species specified in the `catmap_params` by applying the specified correction modes.
"""
function CatmapInterface.compute_free_energies!(free_energies::Dict{InterfaceParams, Dict{String, T}}, catmap_params, params) where {T <: Real}
	for intparams in params
		free_energies[intparams] = Dict([ sp => 0.0	for sp in keys(catmap_params.species_list)	])
		CatmapInterface.compute_free_energies!(free_energies[intparams], catmap_params, intparams)
	end
end

"""
$(SIGNATURES)

Computes the rate of an elementary reaction that passes through a transition state.

The rate law is based on the Arrhenius relation.
The change in Gibbs free energy between the initial and transition state is the needed activation energy.
A `prefactor` and the product of the activities `activprod` complete the rate law.
"""
function ratelaw_TS(prefactor, Gf_IS, Gf_TS, T, activprod)
    @local_phconstants k_B
    prefactor * exp(- Gf_TS/ (k_B * T)) * exp(Gf_IS / (k_B * T)) * activprod
end

"""
$(SIGNATURES)

Create a [ReactionSystem](https://docs.sciml.ai/Catalyst/stable/api/catalyst_api/#Catalyst.ReactionSystem) for the specified microkinetic model.

For each elementary reaction the [`CatmapInterface.ratelaw_TS`](@ref) is used.
No separate rate equations for the solvent (= H₂O) and the active sites are created but their activities can be specified as parameters.
For ficitious gases (OH⁻ and H⁺) and adsorbates the activity coefficients are assumed to be 1.
The activity coefficients of the gaseous species can specified as parameters.
The thermodynamical corrections to the DFT-data of the formation energies are applied according to the specified modes.
New modes can be added by the user by adding a function with the same name to the module. 

if `conserve_pressures==true`,  conserve the pressures of the gaseous and fictious species involved in the heterogeneous reaction network.
The pressures of the gaseous and fictious species are conserved by adding an additional (production/elimination) reaction for each species.
"""
function create_reaction_network(catmap_params::CatmapParams; conserve_pressures = false)
    (; species_list, T) = catmap_params

    @parameters σ ϕ_we ϕ local_pH
    @variables t
    vars        = Dict{String, Num}() # converages and concentrations
    θ           = Dict{String, Num}() # coverages
    activ_coefs = Dict{String, Num}()
    β           = Dict{String, Num}() # transition state beta 
    formation_energies= Dict{String, Num}()# I changed this 8/1
    barriers = Dict{String, Num}() ## barrier
    for (s, sp) in species_list
        if isa(sp, SiteSpecies) # the coverage of the free sites of site type is 1 - sum(coverages of adsorbates on site)
            vars[s] = Num(1)
        end
    end
    for (s, sp) in species_list
        Es= Symbol("E$s")
        formation_energies[s] = first(@parameters $Es = sp.formation_energy)
        if s =="H2O_g" # the solvent is assumed to have constant activity
            as      = Symbol("a$s")
            vars[s] = first(@parameters $as)
        elseif (isa(sp, FictiousSpecies) && s ≠ "ele_g") # fictious species and adsorbates have no activity coeff
            ss          = Symbol(s)
            vars[s]     = first(@species $ss(t))

        elseif isa(sp, AdsorbateSpecies)
            ss                  = Symbol(s)
            vars[s]             = first(@species $ss(t))
            θ[s]                = vars[s] #* Num(sp.n_sites)
            vars["_$(sp.site)"]-= vars[s] #* Num(sp.n_sites)
        elseif (isa(sp, GasSpecies) && s ≠  "H2O_g")
            ss              = Symbol(s)
            vars[s]         = first(@species $ss(t))
            gs              = Symbol("γ$s")
            activ_coefs[s]  = first(@parameters $gs)
        elseif isa(sp, TStateSpecies)
            Gas   = Symbol("Ga$s")
            barriers[s] = first(@parameters $Gas = sp.barrier) # should be checked.
            βs   = Symbol("β$s")
            β[s] = first(@parameters $βs = sp.β) # note: default value not included when generate_function is used!
        end
    end
    
    free_energies = Dict(zip(keys(species_list), fill(Num(0.0), length(species_list))))
    TS_list = Dict(k => v for (k,v) in species_list if v isa TStateSpeciesdd)
    free_energies_only_ad_TS= Dict(zip(keys(TS_list), fill(Num(0.0), length(TS_list))))
    compute_free_energies!(free_energies, catmap_params::CatmapParams, formation_energies, θ, σ, ϕ_we, ϕ, local_pH, β)
    compute_free_energies_only_ad_TS!(free_energies_only_ad, catmap_params::CatmapParams, Ga, θ, σ, ϕ_we, ϕ, local_pH, β)


    function process_reaction_side(reactants)
        @local_unitfactors mol dm
        Gf = Num(0.0)
        rs = Num[]
        γs = Int[]
        a  = Num(1.0)
        for (reactant, factor) in reactants
            sp = species_list[reactant]
            if (reactant =="H2O_g" || isa(sp, SiteSpecies))
                a *= vars[reactant]^factor
            elseif (isa(sp, FictiousSpecies) && reactant ≠ "ele_g") # activity is assumed 1 b/c their influence is in rate constant
                push!(rs, vars[reactant])
                push!(γs, factor)
            elseif isa(sp, AdsorbateSpecies) # activity coefficients are assumed to be 1
                push!(rs, vars[reactant])
                push!(γs, factor)
                a *= (vars[reactant])^factor #  (vars[reactant] * Num(sp.n_sites))^factor
            elseif (isa(sp, GasSpecies) && reactant ≠ "H2O_g")
                push!(rs, vars[reactant])
                push!(γs, factor)
                a *= (activ_coefs[reactant] * vars[reactant])^factor
            end
            Gf += factor * free_energies[reactant]
        end
        Gf, rs, γs, a
    end

    function process_reaction_side_only_ad_TS(reactants) # process_reaction without ad-ad interaction
        barrier_ad_corr = Num(0.0)
        for (reactant, factor) in reactants
            barrier_ad_corr += factor * free_energies_only_TS[reactant]
        end
        barrier_ad_corr
    end

    

    function compute_reversiblepotential(Gf_IS, Gf_FS, surface_charge_relation, ϕ_we, ϕ) ## While calculating revpot, energies[OH_g], energies[H_g] should be replaced by the pH-indepedent value(it's in _get_echem_corrections in catmap) & we have to thinks about is it okay to inlclude ad-ad interaction in Gf_FS, Gf_IS in this funciton.
        ΔGf_r = substitute(Gf_FS - Gf_IS, Dict(surface_charge_relation))  
        symbolic_solve(ΔGf_r ~ 0, ϕ_we - ϕ)
    end


            
            




    rxs = Reaction[]
    for ((; educts, products, tstate), prefactor) in zip(catmap_params.reactions, catmap_params.prefactors)
        (Gf_IS, es, αs, af) = process_reaction_side(educts)
        (Gf_FS, ps, βs, ar) = process_reaction_side(products)
        (barrier_ad_corr, es, αs, af) = process_reaction_side_for_only_TS(tstate.components) ## adding ad-ad correction for TS
        
        Gf_TS= if isnothing(tstate)
            max(Gf_IS, Gf_FS)
        elseif !isnothing(tstate.barrier)
            surface_charge_relation = σ => C_gap*(ϕ_we - ϕ - ϕ_pzc)
            ϕ_rev = compute_reversiblepotential(Gf_IS, Gf_FS, σ => C_gap * (ϕ_we - ϕ - ϕ_pzc), ϕ_we - ϕ)
            if beta_mode == :simple
                ΔGf_r = substitute(Gf_FS - Gf_IS, Dict(surface_charge_relation)) ## do not need to subtrac ΔGf_r at revpot, since it is just 0.
                Gf_IS + barrier_ad_corr + tstate.β*ΔGf_r
            elseif beta_mode == :effectie_surface_charging
                Gf_IS + barrier_ad_corr + tstate.β*e*(ϕ_we - ϕ - ϕ_rev) ## e should be defined
            else
                throw(ArgumentError("$beta_mode is not a valid beta-mode"))
            end
        else
            max(Gf_IS, Gf_FS, mapreproduce(x-> free_energies[first(x)]^last(x), +, tstate.components))
        rxn_f = Reaction(ratelaw_TS(prefactor, Gf_IS, Gf_TS, T, af), es, ps, αs, βs; only_use_rate=true)
        rxn_r = Reaction(ratelaw_TS(prefactor, Gf_FS, Gf_TS, T, ar), ps, es, βs, αs; only_use_rate=true)
        push!(rxs, rxn_f)
        push!(rxs, rxn_r)
    end

    rn=ReactionSystem(rxs, t, name = :microkinetics, combinatoric_ratelaws=false)

    if conserve_pressures
	stoichmat = netstoichmat(rn)
	rr = reactionrates(rn)
	nr = numreactions(rn)
 	for (isp, s) in enumerate(species(rn))
 	    sp = species_list[string(Symbolics.operation(Symbolics.value(s)))]
            if isa(sp, GasSpecies) || isa(sp, FictiousSpecies)
		R = sum([stoichmat[isp ,i] * rr[i] for i in 1:nr])
	        r = Reaction(R, [s], nothing; only_use_rate=true)
                push!(rxs, r)
	    end
        end
        rn1=ReactionSystem(rxs, t, name = :microkinetics, combinatoric_ratelaws=false)
        @assert all(map(enumerate(species(rn1))) do (isp, s)
                    sp = species_list[string(Symbolics.operation(Symbolics.value(s)))]
                    new_stoichmat = netstoichmat(rn1)
                    new_rr = reactionrates(rn1)
                    new_nr = numreactions(rn1)
                    isequal(sum([new_stoichmat[isp ,i] * new_rr[i] for i in 1:new_nr]), isa(sp, GasSpecies) || isa(sp, FictiousSpecies) ? Num(0.0) : sum([stoichmat[isp ,i] * rr[i] for i in 1:nr]))
                    end)
        return complete(rn1)
    else
        return complete(rn)
    end
end

"""
$(SIGNATURES)

Transform the (micro-)kinetic model from surface/gas-reactions to surface/electrolyte-reactions using Henry's law.
"""
function liquidize(odesys::ODESystem, catmap_params::CatmapParams)
    @local_unitfactors bar
    (; species_list) = catmap_params

    sts     = unknowns(odesys)
    ps      = parameters(odesys)

    usubs = Pair{SymbolicUtils.BasicSymbolic{Real}, SymbolicUtils.BasicSymbolic{Real}}[]
    csubs = Pair{Num, Num}[]
    psubs = Pair{SymbolicUtils.BasicSymbolic{Real}, SymbolicUtils.BasicSymbolic{Real}}[]
    @variables t
    for st in sts
        sp = species_list[string(Symbolics.operation(Symbolics.value(st)))]
        if isa(sp, GasSpecies)
            (; henry_const) = sp
            ss          = Symbol("$(sp.species_name)_aq")
            var         = first(@variables $ss(t))
            push!(usubs, st => Symbolics.value(var))
            push!(csubs, st => var / henry_const / bar)

            gs          = Symbol("γ$(sp.species_name)_aq")
            activ_coef  = first(@parameters $gs)
            p = getproperty(odesys, Symbol("γ$(sp.species_name)_g"); namespace=false)
            p = Symbolics.value(p)
            push!(psubs, p => Symbolics.value(activ_coef))
        end
    end

    new_eqs = Equation[]
    for eq in equations(odesys)
        lhs = expand_derivatives(substitute(eq.lhs, Dict(usubs)))
        rhs = substitute(eq.rhs, Dict(csubs..., psubs...))
        push!(new_eqs, Equation(lhs, rhs))
    end

    # JF: this runs into the fact that the symbolic tools now require to work with `ifelse()` instead of
    # `if ... then ... else ... end`. The later seems to be used deep down in some packages.
    # structural_simplify(ODESystem(new_eqs, t, replace(sts, usubs...), replace(ps, psubs...); name=odesys.name))
    
    ODESystem(new_eqs, t, replace(sts, usubs...), replace(ps, psubs...); name=odesys.name)
end

"""
    $(SIGNATURES)

Create index map of odesys parameters as a `Dict{Symbol,Int}`.
E.g. with `pidx=paramsidx(odesys)`, the index of `odesys.σ` can be accessed via `pidx[:σ]`.
"""
function paramsidx(odesys)
    pidx=Dict{Symbol, Int}()
    px=Catalyst.parameters(odesys)
    for i=1:length(px)
	pidx[Catalyst.getname(px[i])]=i
    end
    return pidx
end


"""
$(SIGNATURES)

Generate a mutating function from a `ReactionSystem` that computes the concentration fluxes due to the reaction.
"""
function generate_function(rn::ReactionSystem; dvs::Vector{Tval}=species(rn), ps::Vector{Tval}=parameters(rn)) where {Tval <: Union{SymbolicUtils.BasicSymbolic{Real}, Num}}
    @assert Set(dvs) == Set(species(rn))
    @assert Set(ps)  == Set(parameters(rn))

    species_map = speciesmap(rn)
    sys = convert(ODESystem, rn; combinatoric_ratelaws=false)
    eqs = equations(sys)
    rhss = [-1 * eqs[species_map[dv]].rhs for dv in dvs] # multiply by -1 because the orientation assumed in VoronoiFVM physics functions

    u = map(x -> ModelingToolkit.time_varying_as_func(Symbolics.value(x), sys), dvs)
    p = map(x -> ModelingToolkit.time_varying_as_func(Symbolics.value(x), sys), ps)
    t = ModelingToolkit.get_iv(sys)

    pre, sol_states = ModelingToolkit.get_substitutions_and_solved_unknowns(sys, no_postprocess = false)

    f_expr = build_function(rhss, u, p, t; postprocess_fbody = pre, states = sol_states)[2]
    drop_expr(@RuntimeGeneratedFunction(@__MODULE__, f_expr))
end

"""
$(SIGNATURES)

Generate a mutating function from a `ODESystem` that computes the concentration fluxes due to the reaction.
"""
function generate_function(sys::ODESystem; dvs=unknowns(sys), ps=parameters(sys))
    @assert Set(dvs) == Set(unknowns(sys))
    @assert Set(ps)  == Set(parameters(sys))


    #state_map = Dict(zip(unknowns(sys), length(unknowns(sys))))
    state_map = Dict([st => i for (i, st) in enumerate(unknowns(sys))])
    eqs = equations(sys)
    rhss = [-1 * eqs[state_map[dv]].rhs for dv in dvs] # multiply by -1 because the orientation assumed in VoronoiFVM physics functions

    u = map(x -> ModelingToolkit.time_varying_as_func(Symbolics.value(x), sys), dvs)
    p = map(x -> ModelingToolkit.time_varying_as_func(Symbolics.value(x), sys), ps)
    t = ModelingToolkit.get_iv(sys)

    pre, sol_states = ModelingToolkit.get_substitutions_and_solved_unknowns(sys, no_postprocess = false)

    f_expr = build_function(rhss, u, p, t; postprocess_fbody = pre, states = sol_states)[2]
    drop_expr(@RuntimeGeneratedFunction(@__MODULE__, f_expr))
end





