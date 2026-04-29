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
function compute_free_energies!(
        free_energies, Ga, catmap_params::CatmapParams, formation_energies,
        θ, σ, ϕ_we, ϕ, local_pH,
        β = Dict([s => sp.β for (s, sp) in catmap_params.species_list if isa(sp, TStateSpecies)]);
        symbolic_formation_energies::Bool = false
    )
    (; adsorbate_interaction_params, gas_thermo_mode, adsorbate_thermo_mode, electrochemical_thermo_mode, beta_mode) = catmap_params
    (; adsorbate_interaction_model) = adsorbate_interaction_params

    if symbolic_formation_energies && beta_mode == :none
        for (s, formation_energy) in formation_energies
            free_energies[s] += formation_energy
        end
    elseif symbolic_formation_energies && (beta_mode == :effective_surface_charging || beta_mode == :simple)
        for (s, formation_energy) in formation_energies
            free_energies[s] += formation_energy
        end
        merge!(free_energies, Ga)
    else
        for (s, (; formation_energy)) in catmap_params.species_list
            free_energies[s] += formation_energy
        end
    end


    adsorbate_interaction_correction! = getfield(@__MODULE__, Symbol(adsorbate_interaction_model, "_adsorbate_interaction"))
    gas_thermo_correction! = getfield(@__MODULE__, gas_thermo_mode)
    adsorbate_thermo_correction! = getfield(@__MODULE__, adsorbate_thermo_mode)
    electrochemical_thermo_correction! = getfield(@__MODULE__, electrochemical_thermo_mode)

    adsorbate_interaction_correction!(free_energies, catmap_params, θ)

    thermo_corrections = Dict(zip(keys(free_energies), zeros(valtype(free_energies), length(free_energies))))
    gas_thermo_correction!(thermo_corrections, catmap_params)
    adsorbate_thermo_correction!(thermo_corrections, catmap_params)
    # electrochemical corrections
    electrochemical_thermo_correction!(thermo_corrections, catmap_params, σ, ϕ_we, ϕ, local_pH, β)
    if symbolic_formation_energies && (beta_mode == :effective_surface_charging || beta_mode == :simple)
        for (species, thermo_correction) in thermo_corrections
            if !haskey(Ga, species)
                free_energies[species] += thermo_correction
            end
        end
    else
        for (species, thermo_correction) in thermo_corrections
            free_energies[species] += thermo_correction
        end
    end
    #    @show free_energies
    return nothing
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

    adsorbate_interaction_correction! = getfield(@__MODULE__, Symbol(adsorbate_interaction_model, "_adsorbate_interaction"))
    gas_thermo_correction! = getfield(@__MODULE__, gas_thermo_mode)
    adsorbate_thermo_correction! = getfield(@__MODULE__, adsorbate_thermo_mode)
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
    return nothing
end

"""
$(SIGNATURES)

Compute the Gibbs free energies of all species specified in the `catmap_params` by applying the specified correction modes.
"""
function CatmapInterface.compute_free_energies!(free_energies::Dict{InterfaceParams, Dict{String, T}}, catmap_params, params) where {T <: Real}
    for intparams in params
        free_energies[intparams] = Dict([ sp => 0.0    for sp in keys(catmap_params.species_list)    ])
        CatmapInterface.compute_free_energies!(free_energies[intparams], catmap_params, intparams)
    end
    return
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
    return prefactor * exp(- Gf_TS / (k_B * T)) * exp(Gf_IS / (k_B * T)) * activprod
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
function create_reaction_network(catmap_params::CatmapParams; conserve_pressures = false, symbolic_formation_energies = true)
    (; species_list, T) = catmap_params

    @parameters σ ϕ_we ϕ local_pH C_gap ϕ_pzc
    @variables t
    vars = Dict{String, Num}() # converages and concentrations
    θ = Dict{String, Num}() # coverages
    activ_coefs = Dict{String, Num}()
    β = Dict{String, Num}() # transition state beta
    formation_energies = Dict{String, Num}() # I changed this 8/1
    Ga = Dict{String, Num}() ## barrier
    for (s, sp) in species_list
        if isa(sp, SiteSpecies) # the coverage of the free sites of site type is 1 - sum(coverages of adsorbates on site)
            vars[s] = Num(1)
        end
    end
    for (s, sp) in species_list
        Es = Symbol("E$s")
        formation_energies[s] = sp.formation_energy
        if s == "H2O_g" # the solvent is assumed to have constant activity
            as = Symbol("a$s")
            vars[s] = first(@parameters $as)
            formation_energies[s] = sp.formation_energy
        elseif (isa(sp, FictiousSpecies) && s ≠ "ele_g") # fictious species and adsorbates have no activity coeff
            ss = Symbol(s)
            vars[s] = first(@species $ss(t))
            formation_energies[s] = sp.formation_energy
        elseif isa(sp, AdsorbateSpecies)
            ss = Symbol(s)
            vars[s] = first(@species $ss(t))
            θ[s] = vars[s] #* Num(sp.n_sites)
            vars["_$(sp.site)"] -= vars[s] #* Num(sp.n_sites)
            formation_energies[s] = sp.formation_energy
        elseif (isa(sp, GasSpecies) && s ≠ "H2O_g")
            ss = Symbol(s)
            vars[s] = first(@species $ss(t))
            gs = Symbol("γ$s")
            activ_coefs[s] = first(@parameters $gs)
            formation_energies[s] = sp.formation_energy
        elseif isa(sp, TStateSpecies)
            Es = Symbol("E$s")
            formation_energies[s] = first(@parameters $Es = sp.formation_energy)
            Gas = Symbol("Ga$s")
            Ga[s] = first(@parameters $Gas = sp.barrier) # should be checked.
            βs = Symbol("β$s")
            β[s] = first(@parameters $βs = sp.β) # note: default value not included when generate_function is used!
        end
    end

    free_energies = Dict(zip(keys(species_list), fill(Num(0.0), length(species_list))))
    compute_free_energies!(free_energies, Ga, catmap_params::CatmapParams, formation_energies, θ, σ, ϕ_we, ϕ, local_pH, β; symbolic_formation_energies)
    #    @show free_energies

    function process_reaction_side(reactants)
        @local_unitfactors mol dm
        Gf = Num(0.0)
        rs = Num[]
        γs = Int[]
        a = Num(1.0)
        for (reactant, factor) in reactants
            sp = species_list[reactant]
            if (reactant == "H2O_g" || isa(sp, SiteSpecies))
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
        return Gf, rs, γs, a
    end

    function compute_reversiblepotential(Gf_IS, Gf_FS, surface_charge_relation, ϕ_we, θ, local_pH) ## While calculating revpot, energies[OH_g], energies[H_g] should be replaced by the pH-indepedent value(it's in _get_echem_corrections in catmap) & we have to thinks about is it okay to inlclude ad-ad interaction in Gf_FS, Gf_IS in this funciton.
        @local_unitfactors μF cm
        ΔGf_r = substitute(Gf_FS - Gf_IS, Dict(surface_charge_relation))
        ΔGf_r = substitute(ΔGf_r, Dict(collect(values(θ)) .=> 0)) # exclude ad-ad interaction
        ΔGf_r = substitute(ΔGf_r, Dict(local_pH => 0)) # exclude ph_dependece
        ΔGf_r = substitute(ΔGf_r, Dict(ϕ => 0)) # assume that potential at reaction_plane is 0
        ΔGf_r = Symbolics.expand(ΔGf_r)
        variable = Symbolics.get_variables(ΔGf_r)
        if !any(v -> isequal(v, ϕ_we), variable) ## for no_surface charge & no electron transfer
            return 0
        else
            sol = Symbolics.symbolic_solve(ΔGf_r ~ 0, ϕ_we)
        end
        return sol
    end


    rxs = Reaction[]
    for ((; educts, products, tstate), prefactor) in zip(catmap_params.reactions, catmap_params.prefactors)
        number_electron = get(Dict(educts), "ele_g", 0.0)
        (Gf_IS, es, αs, af) = process_reaction_side(educts)
        (Gf_FS, ps, βs, ar) = process_reaction_side(products)
        @local_phconstants e
        @local_unitfactors eV cm μF
        Gf_TS = if isnothing(tstate)
            max(Gf_IS, Gf_FS)
        elseif isnothing(tstate.barrier) || catmap_params.beta_mode == :none
            max(Gf_IS, Gf_FS, mapreduce(x -> free_energies[first(x)]^last(x), +, tstate.components))
        elseif !isnothing(tstate.barrier)
            surface_charge_relation = σ => C_gap * (ϕ_we - ϕ - ϕ_pzc)
            ϕ_rev = compute_reversiblepotential(Gf_IS, Gf_FS, surface_charge_relation, ϕ_we, θ, local_pH)
            tstate_name = first(tstate.components)[1]
            tstate_factor = first(tstate.components)[2]
            if catmap_params.beta_mode == :simple
                ΔGf_r = substitute(Gf_FS - Gf_IS, Dict(surface_charge_relation))
                ΔGf_r = substitute(ΔGf_r, Dict(local_pH => 0))
                ΔGf_r = substitute(ΔGf_r, Dict(collect(values(θ)) .=> 0))
                IS_no_int = substitute(Gf_IS, Dict(collect(values(θ)) .=> 0))
                max(Gf_IS, Gf_FS, IS_no_int + free_energies[tstate_name] * tstate_factor + β[tstate_name] * ΔGf_r)

            elseif catmap_params.beta_mode == :effective_surface_charging
                IS_no_int = substitute(Gf_IS, Dict(collect(values(θ)) .=> 0))
                max(Gf_IS, Gf_FS, IS_no_int + free_energies[tstate_name] * tstate_factor + β[tstate_name] * e * (ϕ_we - ϕ - ϕ_rev[1]))
            else
                throw(ArgumentError("$beta_mode is not a valid beta-mode"))
            end
        else
            throw(ArgumentError("$beta_mode is not defined in mkm file"))
        end
        rxn_f = Reaction(ratelaw_TS(prefactor, Gf_IS, Gf_TS, T, af), es, ps, αs, βs; only_use_rate = true)
        rxn_r = Reaction(ratelaw_TS(prefactor, Gf_FS, Gf_TS, T, ar), ps, es, βs, αs; only_use_rate = true)
        push!(rxs, rxn_f)
        push!(rxs, rxn_r)
    end
    rn = ReactionSystem(rxs, t, name = :microkinetics, combinatoric_ratelaws = false)
    if conserve_pressures
        stoichmat = netstoichmat(rn)
        rr = reactionrates(rn)
        nr = numreactions(rn)
        for (isp, s) in enumerate(species(rn))
            sp = species_list[string(Symbolics.operation(Symbolics.value(s)))]
            if isa(sp, GasSpecies) || isa(sp, FictiousSpecies)
                R = sum([stoichmat[isp, i] * rr[i] for i in 1:nr])
                r = Reaction(R, [s], nothing; only_use_rate = true)
                push!(rxs, r)
            end
        end
        rn1 = ReactionSystem(rxs, t, name = :microkinetics, combinatoric_ratelaws = false)
        @assert all(
            map(enumerate(species(rn1))) do (isp, s)
                sp = species_list[string(Symbolics.operation(Symbolics.value(s)))]
                new_stoichmat = netstoichmat(rn1)
                new_rr = reactionrates(rn1)
                new_nr = numreactions(rn1)
                isequal(sum([new_stoichmat[isp, i] * new_rr[i] for i in 1:new_nr]), isa(sp, GasSpecies) || isa(sp, FictiousSpecies) ? Num(0.0) : sum([stoichmat[isp, i] * rr[i] for i in 1:nr]))
            end
        )
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

    sts = unknowns(odesys)
    ps = parameters(odesys)

    usubs = Pair{SymbolicUtils.BasicSymbolic, SymbolicUtils.BasicSymbolic}[]
    csubs = Pair{Num, Num}[]
    psubs = Pair{SymbolicUtils.BasicSymbolic, SymbolicUtils.BasicSymbolic}[]
    @variables t
    for st in sts
        sp = species_list[string(Symbolics.operation(Symbolics.value(st)))]
        if isa(sp, GasSpecies)
            (; henry_const) = sp
            ss = Symbol("$(sp.species_name)_aq")
            var = first(@variables $ss(t))
            push!(usubs, st => Symbolics.value(var))
            push!(csubs, st => var / henry_const / bar)

            gs = Symbol("γ$(sp.species_name)_aq")
            activ_coef = first(@parameters $gs)
            p = getproperty(odesys, Symbol("γ$(sp.species_name)_g"); namespace = false)
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

    return ODESystem(new_eqs, t, replace(sts, usubs...), replace(ps, psubs...); name = nameof(odesys))
end

"""
    $(SIGNATURES)

Create index map of odesys parameters as a `Dict{Symbol,Int}`.
E.g. with `pidx=paramsidx(odesys)`, the index of `odesys.σ` can be accessed via `pidx[:σ]`.
"""
function paramsidx(odesys)
    pidx = Dict{Symbol, Int}()
    px = Catalyst.parameters(odesys)
    for i in 1:length(px)
        pidx[Catalyst.getname(px[i])] = i
    end
    return pidx
end


"""
$(SIGNATURES)

Generate a mutating function from a `ReactionSystem` that computes the concentration fluxes due to the reaction.
"""
function generate_function(
        rn::ReactionSystem;
        dvs = species(rn),
        ps = parameters(rn)
    )
    @assert Set(dvs) == Set(species(rn))
    @assert Set(ps) == Set(parameters(rn))

    #    species_map = speciesmap(rn)
    sys = ode_model(rn; combinatoric_ratelaws = false)
    prob = ODEProblem(sys, zeros(length(dvs)), (0, 1.0), ps)
    return function (f, u, p, t)
        prob.f(f, u, p, t)
        f .*= -1
        return nothing
    end
    #=
    eqs = equations(sys)
    rhss = [-1 * eqs[species_map[dv]].rhs for dv in dvs] # multiply by -1 because the orientation assumed in VoronoiFVM physics functions

    u = dvs # map(x -> ModelingToolkit.time_varying_as_func(Symbolics.value(x), sys), dvs)
    p = ps #map(x -> ModelingToolkit.time_varying_as_func(Symbolics.value(x), sys), ps)
    t = ModelingToolkit.get_iv(sys)

    # pre, sol_states = ModelingToolkit.get_substitutions_and_solved_unknowns(sys, no_postprocess = false)

    f_expr = build_function(rhss, u, p, t; postprocess_fbody = pre, states = sol_states)[2]
    return drop_expr(@RuntimeGeneratedFunction(@__MODULE__, f_expr))
    =#
end

"""
$(SIGNATURES)

Generate a mutating function from a `ODESystem` that computes the concentration fluxes due to the reaction.
"""
function generate_function(sys::ODESystem; dvs = unknowns(sys), ps = parameters(sys))
    @assert Set(dvs) == Set(unknowns(sys))
    @assert Set(ps) == Set(parameters(sys))


    #state_map = Dict(zip(unknowns(sys), length(unknowns(sys))))
    state_map = Dict([st => i for (i, st) in enumerate(unknowns(sys))])
    eqs = equations(sys)
    rhss = [-1 * eqs[state_map[dv]].rhs for dv in dvs] # multiply by -1 because the orientation assumed in VoronoiFVM physics functions

    u = map(x -> ModelingToolkit.time_varying_as_func(Symbolics.value(x), sys), dvs)
    p = map(x -> ModelingToolkit.time_varying_as_func(Symbolics.value(x), sys), ps)
    t = ModelingToolkit.get_iv(sys)

    pre, sol_states = ModelingToolkit.get_substitutions_and_solved_unknowns(sys, no_postprocess = false)

    f_expr = build_function(rhss, u, p, t; postprocess_fbody = pre, states = sol_states)[2]
    return drop_expr(@RuntimeGeneratedFunction(@__MODULE__, f_expr))
end
