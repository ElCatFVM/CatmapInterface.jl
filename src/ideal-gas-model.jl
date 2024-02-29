# py"""
# from ase.collections import g2

# def moleculedata(formula):
#     a = g2[formula]
#     return a.numbers, a.get_masses(), a.get_positions()
# """

# function molecule_geometry(geo)
#     if geo == monoatomic
#         return "monoatomic"    
#     elseif geo == linear
#         return "linear"
#     elseif geo == nonlinear
#         return "nonlinear"
#     end
# end


"""
"""
@enum MoleculeGeometry monoatomic linear nonlinear

"""
$(TYPEDEF)

Data struct for an ideal gas

$(TYPEDFIELDS)
"""
@kwdef struct IdealGas
    """
    List of atomic numbers of the molecule's atoms
    """
    elements::Vector{Int}
    """
    List of the masses of the molecule's atoms
    """
    masses::Vector{Float64}
    """
    Matrix of 3D-positions with one row per atom in the molecule
    """
    positions::Matrix{Float64}
    """
    """
    symmetrynumber::Int32
    """
    List of normal modes of vibration
    """
    frequencies::Vector{Float64}
    """

    """
    spin::Float64
    """
    Molecule geometry, one of monoatomic/linear/nonlinear
    """
    geometry::MoleculeGeometry
    """
    Temperature of a canonical ensemble
    """
    temperature::Float64
end

function totalmass(idealgas::IdealGas)
    sum(idealgas.masses)
end

function momentsofinertia(idealgas::IdealGas)
    m = idealgas.masses
    x = idealgas.positions
    M = sum(m)
    
    cm = sum(x .* m, dims=1) / M
    x .-= cm

    IM = zeros(3,3)
    for i in 1:length(m)
        IM += m[i] * ((x[i, :])' * x[i, :] * I - x[i, :] * (x[i, :])') 
    end
    eigvals(IM)
end

"""
Compute the entropy of an ideal gas
"""
function entropy(idealgas::IdealGas)
    @local_phconstants k_B h R ħ
    @local_unitfactors bar eV atm

    N = length(idealgas.elements)
    M = totalmass(idealgas)
    T = idealgas.temperature
    σ = idealgas.symmetrynumber
    geometry = idealgas.geometry
    if geometry == monoatomic
        ω = []
    elseif geometry == linear
        ω = sort(idealgas.frequencies, rev=true)[1:(3 * N - 5)]
    else
        ω = sort(idealgas.frequencies, rev=true)[1:(3 * N - 6)]
    end
    S = idealgas.spin
    (IA, IB, IC) = momentsofinertia(idealgas)
    if geometry == linear
        @assert isapprox(IA, 0.0, atol = 1e-6) && isapprox(IB, IC, rtol = 1e-6)
    end
    P° = 1 * bar

    translational = k_B * (log((2 * π * M * k_B * T / h^2)^(3/2) * k_B * T / P°) + 5 / 2)
    rotational = 0.0
    if geometry == linear
        rotational +=  k_B * (log(2 * IC * k_B * T / (ħ^2 * σ)) + 1)
    elseif geometry == nonlinear
        rotational += k_B * (log(sqrt(π * IA * IB * IC) / σ * (2 * k_B * T / ħ^2)^(3/2)) + 3/2)
    end
    vibrational = k_B * sum(@. h * ω / (k_B * T * (exp(h * ω / (k_B * T)) - 1)) - log(1 - exp(-h * ω / (k_B * T))))
    electronic = k_B * log(2 * S + 1)

    translational + rotational + vibrational + electronic
end

"""
Compute the enthalpy of an ideal gas with 0 eV electronic ground-state energy
"""
function enthalpy(idealgas::IdealGas)
    @local_phconstants k_B h R
    @local_unitfactors eV

    N = length(idealgas.elements)
    T = idealgas.temperature
    geometry = idealgas.geometry
    if geometry == monoatomic
        ω = []
    elseif geometry == linear
        ω = sort(idealgas.frequencies, rev=true)[1:(3 * N - 5)]
    else
        ω = sort(idealgas.frequencies, rev=true)[1:(3 * N - 6)]
    end

    zpe = h * sum(ω) / 2
    translational = 3 / 2 * k_B * T
    rotational = 0
    if geometry == linear
        rotational += k_B * T
    elseif  geometry == nonlinear
        rotational += 3 / 2 * k_B * T
    end
    vibrational = h * sum(@. ω / (exp(h * ω / (k_B * T)) - 1))
    corr = k_B * T

    zpe + translational + rotational + vibrational + corr
end
