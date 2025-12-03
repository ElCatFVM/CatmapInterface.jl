"""
$(SIGNATURES)

Correct the relative Gibbs free energies of formation due to ideal adsorbate interactions.

Model for the excess chemical potential that assumes ideal adsorbate interactions so that the corrections are zero.
"""
function ideal_adsorbate_interaction(energies, catmap_params, coverages)
    nothing
end

    
# c_0 = 0.0
# if θ_tot <= x0
#     c_0 = 0.0
# elseif θ_tot <= x1
#     α = slope/(2*(x1-x0))
#     c_0 = (α * (θ_tot - x0)^2)/θ_tot
#     #dC = alpha*(1-(x0/theta_tot)**2)
#     #d2C = (2*alpha*x0**2)/(theta_tot**3)
# else
#     c_0 = slope*(θ_tot - cutoff)/θ_tot
#     # dC = slope*(cutoff/(theta_tot**2))
#     # d2C = (-2*slope*cutoff)/(theta_tot**3)
# end
# c_0

function smooth_piecewise_linear_interaction_response(θ_tot, interaction_response_params::InteractionResponseParams)
    (; slope, cutoff, smoothing) = interaction_response_params
    θ_1 = min(max(0.0, θ_tot - (cutoff - smoothing)), 2 * smoothing)
    θ_2 = max(0.0, θ_tot - (cutoff + smoothing))

    α   = smoothing > 0 ? slope / (4 * smoothing) : 0.0
    (α * θ_1^2 + slope * θ_2)
end

function linear_interaction_response(θ_tot, interaction_response_params::InteractionResponseParams)
    (; slope) = interaction_response_params
    smooth_piecewise_linear_interaction_response(θ_tot, InteractionResponseParams(; slope, cutoff=0.0, smoothing=0.0))
end

function piecewise_linear_interaction_response(θ_tot, interaction_response_params::InteractionResponseParams)
    (; slope, cutoff) = interaction_response_params
    smooth_piecewise_linear_interaction_response(θ_tot, InteractionResponseParams(; slope, cutoff, smoothing=0.0))
end

geometric_mean_cross_interaction(ϵ_s, ϵ_os) = √(ϵ_s * ϵ_os)
arithemtic_mean_cross_interaction(ϵ_s, ϵ_os) = (ϵ_s + ϵ_os) / 2
neglect_cross_interaction(ϵ_s, ϵ_os) = 0.0

intermediate_state_transition_cross_interaction(ϵ_es, ϵ_ps) = 0.5 * (ϵ_es + ϵ_ps)
final_state_transition_cross_interaction(ϵ_es, ϵ_ps) = ϵ_ps
initial_state_transition_cross_interaction(ϵ_es, ϵ_ps) = ϵ_es
neglect_transition_cross_interaction(ϵ_rs) = 0.0

"""
    get_interaction_term(s::String, sp::AdsorbateSpecies, os::String, osp::AdsorbateSpecies, cross_interaction_function)

Extract the interaction term if needed by applying the `cross_interaction_function`
"""
function _get_interaction_term(s::String, sp::AdsorbateSpecies, os::String, osp::AdsorbateSpecies, cross_interaction_function, species_list)
    ϵ = nothing
    if s == os
        ϵ = sp.self_interaction_param
    else
        ϵ = get(sp.cross_interaction_params, os, nothing)
        if isnothing(ϵ)
            ϵ = get(osp.cross_interaction_params, s, nothing)
        end
        if isnothing(ϵ)
            ϵ_s  = sp.self_interaction_param
            ϵ_os = osp.self_interaction_param
            ϵ    = cross_interaction_function(ϵ_s, ϵ_os)
        end
    end 
    return ϵ
end

"""
    get_interaction_term(s::String, sp::TStateSpecies, os::String, osp::AdsorbateSpecies, cross_interaction_function)

Extract the interaction term if needed by applying the `transition_state_cross_interaction_function`
"""
function _get_interaction_term(s::String, sp::TStateSpecies, os::String, osp::AdsorbateSpecies, cross_interaction_function, species_list)
    ϵ = nothing
    ϵ = get(sp.cross_interaction_params, os, nothing)
    if isnothing(ϵ)
        ϵ = get(osp.cross_interaction_params, s, nothing)
    end
    if isnothing(ϵ)
        ϵ_es = 0.0
        ϵ_ps = 0.0
        for (reactant, factor) in sp.between_species
            if isa(species_list[reactant], AdsorbateSpecies)
                if factor > 0
                    ϵ_ps += factor * species_list[reactant].cross_interaction_params[os]
                else
                    ϵ_es += -factor * species_list[reactant].cross_interaction_params[os]
                end
            end
        end
        ϵ = cross_interaction_function(ϵ_es, ϵ_ps)
    end
    return ϵ
end

function _coverage_of_site(site, species_list, θ)
    θ_tot = 0.0
    for (s, sp) in species_list
        if isa(sp, AdsorbateSpecies) && (site == sp.site)
            θ_tot += θ[s]
        end
    end
    θ_tot
end

"""
$(SIGNATURES)

Add correction terms to the adsorbation energies based on first order interactions between the adsorbates. 

Model for the excess chemical potential that correct the relative Gibbs free energies of formation due to first-order adsorbate interactions.
This adsorbation interaction model is expained in ![CatMAP's documentation](https://catmap.readthedocs.io/en/latest/topics/including_adsorbate_adsorbate_interactions.html#coverage-dependent-adsorption-eneriges) and in this ![issue](https://github.com/smaasz/CatmapInterface.jl/issues/10).
"""
function first_order_adsorbate_interaction(energies, catmap_params::CatmapParams, θ)
    @local_unitfactors eV
    (; species_list, adsorbate_interaction_params) = catmap_params
    (; interaction_response_function, cross_interaction_mode, transition_state_cross_interaction_mode) = adsorbate_interaction_params

    response_function                           = getfield(@__MODULE__, Symbol(interaction_response_function, "_interaction_response"))
    cross_interaction_function                  = getfield(@__MODULE__, Symbol(cross_interaction_mode, "_cross_interaction"))
    transition_state_cross_interaction_function = getfield(@__MODULE__, Symbol(transition_state_cross_interaction_mode, "_transition_cross_interaction"))

    for (Species, interaction_function) in [(AdsorbateSpecies, cross_interaction_function), (TStateSpecies, transition_state_cross_interaction_function)]
        for (s, sp) in species_list
            if isa(sp, Species)
                (; site)                        = sp
                θ_tot                           = _coverage_of_site(site, species_list, θ)
                (; interaction_response_params) = species_list["_$site"]
                response_value                  = response_function(θ_tot, interaction_response_params)
                for (os, osp) in species_list
                    if isa(osp, AdsorbateSpecies)
                        ϵ = _get_interaction_term(s, sp, os, osp, interaction_function, species_list)
                        sp.cross_interaction_params[os] = ϵ
                        osp.cross_interaction_params[s] = ϵ
                        energies[s] += response_value * ϵ * ((θ[os])/(θ_tot + 1.0e-50)) * eV
                    end
                end
            end
        end
    end
end

"""
$(SIGNATURES)

Add thermodynamic correction terms for all gas species using the ideal gas approximation.

An ab-initio statistical model to correct the DFT-energies in order to obtain relative Gibbs free energies of formation. A reference for the approach is given in the book 'Essentials of Computational Chemistry: Theories and Models', 2nd edition, by Cramer, C.J. and published by Wiley.
"""
function ideal_gas(energies::Dict{String, Tval}, catmap_params::CatmapParams) where Tval <: Real 
    
    ideal_gas_params = py"ideal_gas_params"
    (; species_list, T) = catmap_params
    for (s, sp) in species_list
        if isa(sp, GasSpecies)
            (; species_name, frequencies) = sp
            (; symmetrynumber, geometry, spin) = get_ideal_gas_params(species_name)
            (; numbers, positions, masses) = get_molecule_spec(species_name)
            ideal_gas = IdealGas(; elements=numbers, masses, positions, symmetrynumber, frequencies, spin, geometry, temperature=T)
            energies[s] += enthalpy(ideal_gas) - T * entropy(ideal_gas)
        end
    end
end

"""
$(SIGNATURES) 

Add thermodynamic correction terms for all adsorbed species using the harmonic adsorbate approximation.

An ab-initio statistical model to correct the DFT-energies of adsorbed species in order to obtain relative Gibbs free energies of formation.
"""
function harmonic_adsorbate(energies, catmap_params::CatmapParams)
    (; species_list, T) = catmap_params
    for (s, sp) in species_list
        if isa(sp, AdsorbateSpecies)
            frequencies = sp.frequencies
            harmonic_phase = HarmonicPhase(; frequencies, temperature = T)
            energies[s] += enthalpy(harmonic_phase) - T * entropy(harmonic_phase)
        end
    end
    for (s, sp) in species_list
        if isa(sp, TStateSpecies)
            if !isempty(sp.frequencies)
                frequencies = sp.frequencies
                harmonic_phase = HarmonicPhase(; frequencies, temperature = T)
                energies[s] += enthalpy(harmonic_phase) - T * entropy(harmonic_phase)
            else
                (; between_species) = sp
                for bs in first.(between_species)
                    if isa(species_list[bs], AdsorbateSpecies)
                        thermo_correction = energies[bs]
                        energies[s] += 0.5 * thermo_correction
                    end
                end
            end
        end

    end
    nothing
end

"""
$(SIGNATURES) 

Add electrochemical pH-correction terms for all influenced species.
"""
function _get_echem_corrections(energies, catmap_params::CatmapParams, σ, ϕ_we, ϕ, local_pH, β)
    @local_unitfactors eV
    (; species_list, potential_reference_scale, T) = catmap_params
    if haskey(species_list, "H_g") || haskey(species_list, "OH_g")
        G_H2 = energies["H2_g"]
        if potential_reference_scale == "SHE"
            G_H = G_H2 * 0.5 - (0.0592 * local_pH / 298.14 * T) * eV
        elseif potential_reference_scale == "RHE"
            G_H = G_H2 * 0.5
        end
        G_H2O = energies["H2O_g"]
        G_OH = G_H2O - G_H
        if haskey(species_list, "H_g")
            energies["H_g"] += G_H
        end
        if haskey(species_list, "OH_g")
            energies["OH_g"] += G_OH
        end
    end
    nothing
end

"""
$(SIGNATURES) 

Add electrochemical correction terms to the relative Gibbs free energies of formation of 'free' electrons and transition states including 'free' electrons.

The model assumes a linear capacitor model for the double layer between the electrode surface and the inner Helmholtz plane. In this model the corrections accomodate the Frumkin effects of proton-coupled electron transfers. For a reference see 'Double layer charging driven carbon dioxide adsorption limits the rate of electrochemical carbon dioxide reduction on Gold' by Ringe, S. et al. and published in Nature Communications.
"""
function simple_electrochemical(energies, catmap_params::CatmapParams, σ, ϕ_we, ϕ, local_pH, β)
    @local_unitfactors eV
    (; species_list, Uref) = catmap_params
    # simple_electrochem_corrections
    if haskey(energies, "ele_g")
        energies["ele_g"] += - (ϕ_we-ϕ) * eV
    end

    #if haskey(energies, "H_g")
    #   energies["H_g"] += - ϕ * eV
    #end 

    #if haskey(energies, "OH_g") ## this can be expanded to all ionspecies.
    #   energies["OH_g"] +=  ϕ * eV
    #end 

    for (s, sp) in species_list
        if isa(sp, TStateSpecies) && occursin("ele", sp.species_name)
            energies[s] += (-(ϕ_we-ϕ) + β[s] * (ϕ_we -ϕ - Uref)) * eV
        end
    end
    # pH_correction
    _get_echem_corrections(energies, catmap_params, σ, ϕ_we, ϕ, local_pH, β)
    nothing
end

"""
$(SIGNATURES) 

Add electrochemical correction terms to the relative Gibbs free energies of formation of (polarized) adsorbed species.

The model assumes a linear capacitor model for the double layer between the electrode surface and the inner Helmholtz plane. The surface charge dependence is fitted to a quadratic model. For a reference see 'Double layer charging driven carbon dioxide adsorption limits the rate of electrochemical carbon dioxide reduction on Gold' by Ringe, S. et al. and published in Nature Communications.
"""
function hbond_surface_charge_density(energies, catmap_params::CatmapParams, σ, ϕ_we, ϕ, local_pH, β)
    @local_unitfactors eV
    # simple_electrochem
    simple_electrochemical(energies, catmap_params, σ , ϕ_we, ϕ, local_pH, β)
    hbond_dict = py"hbond_dict"
    (; species_list, beta_mode) = catmap_params
    #hbond_electrochemical
    #for (s, sp) in species_list
    #   if isa(sp, AdsorbateSpecies)
    #      (; species_name) = sp
    #      if haskey(hbond_dict, species_name)
    #         energies[s] += hbond_dict[species_name] * eV
    #      end
    #   end
    #end
#    @show energies
    # hbond_surface_charge_density
    for (s, sp) in species_list
        if (isa(sp, AdsorbateSpecies) || isa(sp, TStateSpecies))
            (; a, b) = sp.sigma_params
            if beta_mode==:simple||beta_mode==:effective_surface_charging
               energies[s] += (a * σ) * eV
            else
               energies[s] += (a * σ + b * σ^2) * eV
            end
        end
    end
    nothing
end

