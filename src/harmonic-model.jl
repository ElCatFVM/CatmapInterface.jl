"""
$(TYPEDEF)

Data struct for a harmonic phase

$(TYPEDFIELDS)
"""
@kwdef struct HarmonicPhase
    """
    List of normal modes of vibration in m⁻¹ (wavenumber)
    """
    frequencies::Vector{Float64}
    """
    Temperature of a canonical ensemble
    """
    temperature::Float64
end


"""
$(SIGNATURES)

Compute the entropy of an harmonic phase based on ab-initio data.

In this (ab-initio) statiscal model all degrees of freedom are approximated harmonically.
"""
function entropy(harmonicphase::HarmonicPhase)
    @local_phconstants k_B h c_0

    T = harmonicphase.temperature
    ω = harmonicphase.frequencies .* c_0

    vibrational = k_B * sum(@. h * ω / (k_B * T * (exp(h * ω / (k_B * T)) - 1)) - log(1 - exp(-h * ω / (k_B * T))))

    vibrational
end

"""
$(SIGNATURES)

Compute the enthalpy of an ideal gas with 0 eV electronic ground-state energy

In this (ab-initio) statiscal model all degrees of freedom are approximated harmonically.
"""
function enthalpy(harmonicphase::HarmonicPhase)
    @local_phconstants k_B h R c_0

    T = harmonicphase.temperature
    ω = harmonicphase.frequencies .* c_0
    
    zpe = h * sum(ω) / 2
    vibrational = h * sum(@. ω / (exp(h * ω / (k_B * T)) - 1))

    zpe + vibrational
end