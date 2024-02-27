"""
$(TYPEDEF)

Data struct for a harmonic phase

$(TYPEDFIELDS)
"""
@kwdef struct HarmonicPhase
    """
    Molecule geometry, one of monoatomic/linear/nonlinear
    """
    geometry::MoleculeGeometry
    """
    List of normal modes of vibration
    """
    frequencies::Vector{Float64}
    """
    Temperature of a canonical ensemble
    """
    temperature::Float64
end


"""
Compute the entropy of an harmonic phase
"""
function entropy(harmonicphase::HarmonicPhase)
    @local_phconstants k_B h

    T = harmonicphase.temperature
    geometry = harmonicphase.geometry
    if geometry == monoatomic
        ω = []
    elseif geometry == linear
        ω = sort(harmonicphase.frequencies, rev=true)[1:(3 * N - 5)]
    else
        ω = sort(harmonicphase.frequencies, rev=true)[1:(3 * N - 6)]
    end

    vibrational = k_B * sum(@. h * ω / (k_B * T * (exp(h * ω / (k_B * T)) - 1)) - log(1 - exp(-h * ω / (k_B * T))))

    vibrational
end

"""
Compute the enthalpy of an ideal gas with 0 eV electronic ground-state energy
"""
function enthalpy(harmonicphase::HarmonicPhase)
    @local_phconstants k_B h R

    T = harmonicphase.temperature
    geometry = harmonicphase.geometry
    if geometry == monoatomic
        ω = []
    elseif geometry == linear
        ω = sort(harmonicphase.frequencies, rev=true)[1:(3 * N - 5)]
    else
        ω = sort(harmonicphase.frequencies, rev=true)[1:(3 * N - 6)]
    end

    zpe = h * sum(ω) / 2
    vibrational = h * sum(@. ω / (exp(h * ω / (k_B * T)) - 1))

    zpe + vibrational
end