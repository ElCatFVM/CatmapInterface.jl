# from IUPAC 2016
const atomic_masses = [
    1.008,          # H
    4.002602,       # He
    6.94,           # Li
    9.0121831,      # Be
    10.81,          # B
    12.011,         # C
    14.007,         # N
    15.999,         # O
    18.998403163,   # F
    20.1797,        # Ne
    22.98976928,    # Na
    24.305,         # Mg
    26.9815385,     # Al
    28.085,         # Si 
    30.973761998,   # P
    32.06,          # S 
    35.45,          # Cl 
    39.948,         # Ar
    39.0983,        # K
    40.078,         # Ca
    44.955908,      # Sc
    47.867,         # Ti
    50.9415,        # V
    51.9961,        # Cr
    54.938044,      # Mn
    55.845,         # Fe
    58.933194,      # Co
    58.6934,        # Ni
    63.546,         # Cu
    65.38,          # Zn
    69.723,         # Ga
    72.630,         # Ge
    74.921595,      # As
    78.971,         # Se
    79.904,         # Br 
    83.798,         # Kr
    85.4678,        # Rb
    87.62,          # Sr
    88.90584,       # Y
    91.224,         # Zr
    92.90637,       # Nb
    95.95,          # Mo
    97.90721,       # 98Tc
    101.07,         # Ru
    102.90550,      # Rh
    106.42,         # Pd
    107.8682,       # Ag
    112.414,        # Cd
    114.818,        # In
    118.710,        # Sn
    121.760,        # Sb
    127.60,         # Te
    126.90447,      # I
    131.293,        # Xe
    132.90545196,   # Cs
    137.327,        # Ba
    138.90547,      # La
    140.116,        # Ce
    140.90766,      # Pr
    144.242,        # Nd
    144.91276,      # 145Pm
    150.36,         # Sm
    151.964,        # Eu
    157.25,         # Gd
    158.92535,      # Tb
    162.500,        # Dy
    164.93033,      # Ho
    167.259,        # Er
    168.93422,      # Tm
    173.054,        # Yb
    174.9668,       # Lu
    178.49,         # Hf
    180.94788,      # Ta
    183.84,         # W
    186.207,        # Re
    190.23,         # Os
    192.217,        # Ir
    195.084,        # Pt
    196.966569,     # Au
    200.592,        # Hg
    204.38,         # Tl
    207.2,          # Pb
    208.98040,      # Bi
    208.98243,      # 209Po
    209.98715,      # 210At
    222.01758,      # 222Rn
    223.01974,      # 223Fr
    226.02541,      # 226Ra
    227.02775,      # 227Ac
    232.0377,       # Th
    231.03588,      # Pa
    238.02891,      # U
    237.04817,      # 237Np
    244.06421,      # 244Pu
    243.06138,      # 243Am
    247.07035,      # 247Cm
    247.07031,      # 247Bk
    251.07959,      # 251Cf
    252.0830,       # 252Es
    257.09511,      # 257Fm
    258.09843,      # 258Md
    259.1010,       # 259No
    262.110,        # 262Lr
    267.122,        # 267Rf
    268.126,        # 268Db
    271.134,        # 271Sg
    270.133,        # 270Bh
    269.1338,       # 269Hs
    278.156,        # 278Mt
    281.165,        # 281Ds
    281.166,        # 281Rg
    285.177,        # 285Cn
    286.182,        # 286Nh
    289.190,        # 289Fl
    289.194,        # 289Mc
    293.204,        # 293Lv
    293.208,        # 293Ts
    294.214,        # 294Og
] .* ufac"u"

@kwdef struct MoleculeSpec
    name::String
    numbers::Vector{Int}
    masses::Vector{Float64} 
    positions::Matrix{Float64}
end

function parse_molecule_spec(js)

    name = try
        js["key_value_pairs"]["name"]
    catch e
        if isa(e, KeyError)
            throw("Molecule entry misses 'name' entry")
        else
            rethrow(e)
        end
    end

    numbers = try
        convert(Vector{Int}, js["numbers"])
    catch e
        if isa(e, KeyError)
            throw("Molecule entry misses 'numbers' entry")
        elseif isa(e, MethodError)
            throw("The 'numbers' entry in for the molecule is invalid")
        else
            rethrow(e)
        end
    end
    @assert all(map(number -> 1 ≤ number ≤ length(atomic_masses), numbers)) "Invalid atoms in the molecule"
    masses = map(number->atomic_masses[number], numbers)

    positions = try
        convert.(Vector{Float64}, js["positions"])
    catch
        if isa(e, KeyError)
            throw("Molecule entry with id=$id misses 'positions' entry")
        elseif isa(e, MethodError)
            throw("The 'positions' entry in for the molecule with id=$id is invalid")
        else
            rethrow(e)
        end
    end
    @assert length(positions) == length(numbers) "Too many/few positions"
    @assert all(length.(positions) .== 3) "Every positions must be a coordinate in 3D"
    positions = stack(positions; dims=1) * ufac"Å"
 
    return MoleculeSpec(; name, numbers, masses, positions)
end

function parse_molecule_data(filename)
    js = open(filename) do file
        JSON.parse(file)
    end

    # check ids entry
    ids = Int64[]
    if haskey(js, "ids")
        try
            ids = convert(Vector{Int}, js["ids"])
        catch e
            if isa(e, MethodError)
                throw(ArgumentError("Entry 'ids' in the molecule database is invalid."))
            end
        end
    else
        throw(ArgumentError("Missing 'ids' entry in the molecule database."))
    end
    delete!(js, "ids")
    haskey(js, "nextid") && delete!(js, "nextid") # also delete nextid from json

    @assert all(map(id->haskey(js, "$id"), ids)) "Invalid molecule database."

    molecule_specs = Dict{String, MoleculeSpec}()
    for id in ids
        molecule_spec = parse_molecule_spec(js["$id"])
        molecule_specs[molecule_spec.name] = molecule_spec
    end
    return molecule_specs
end

const molecule_specs = parse_molecule_data(joinpath(readdir(artifact"ase_collections", join=true)[1], "ase", "collections", "g2.json"))


"""
$(SIGNATURES)

Get the specification of the molecule of type 'name'
"""
function get_molecule_spec(name)
    if !haskey(molecule_specs, name)
        throw(ArgumentError("The molecule $name is not included in the molecule database"))
    end
    return  molecule_specs[name]
end

"""
$(SIGNATURES)

Get the symmetrynumber, geometry, and spin of the ideal gas' molecules of type 'name'
"""
function get_ideal_gas_params(name)
    (symmetrynumber, geometry, spin) = try
        py"ideal_gas_params"[name * "_g"]
    catch e
        if isa(e, KeyError)
            throw(ArgumentError("$name has no specified ideal gas parameters"))
        else
            rethrow(e)
        end
    end
    if geometry == "monoatomic"
        geometry = monoatomic    
    elseif geometry == "linear"
        geometry = linear
    elseif geometry == "nonlinear"
        geometry = nonlinear
    else
        throw(ArgumentError("The molecular geometry $geometry is not valid"))
    end

    return (; symmetrynumber = symmetrynumber, geometry = geometry, spin = spin)
end

"""
    instantiate_catmap_template!(instance_file_path, template_file_path, params)

Instantiate a template file by inserting the parameters in the `params`.
"""
function instantiate_catmap_template!(instance_file_path, template_file_path, params, T)
    (; σ, ϕ_we, ϕ, local_pH) = params
    instance_string = open(template_file_path, "r") do template_file
        read(template_file, String)
    end

	replacements = [
		r"descriptor_ranges.?=.*" =>"descriptor_ranges = [[$ϕ_we, $ϕ_we], [$T, $T]]",
		r"voltage_diff_drop.?=.*" => "voltage_diff_drop = $ϕ",
		r"pH.?=.*" => "pH = $local_pH",
		r"\nsigma_input.?=.*" => "\nsigma_input = $σ/0.01", # in μF/cm^2
	]
    instance_string = replace(instance_string, replacements...)

    open(instance_file_path, "w") do instance_file
        write(instance_file, instance_string)
    end

    return instance_file_path
end

"""
conserve_pressures!(rn::Catalyst.ReactionSystem, catmap_params::CatmapInterface.CatmapParams)

Conserve the pressures of the gaseous and fictious species involved in the heterogeneous reaction network `rn`.

The pressures of the gaseous and fictious species are conserved by adding an additional (production/elimination) reaction for each species.
"""
function conserve_pressures!(rn, catmap_params)
	(; species_list) = catmap_params
	stoichmat = netstoichmat(rn)
	rr = reactionrates(rn)
	nr = numreactions(rn)
	for (isp, s) in enumerate(species(rn))
		sp = species_list[string(Symbolics.operation(Symbolics.value(s)))]
		if isa(sp, GasSpecies) || isa(sp, FictiousSpecies)
			R = sum([stoichmat[isp ,i] * rr[i] for i in 1:nr])
			addreaction!(rn, Reaction(R, [s], nothing; only_use_rate=true))
		end
	end
    @assert all(map(enumerate(species(rn))) do (isp, s)
        sp = species_list[string(Symbolics.operation(Symbolics.value(s)))]
        new_stoichmat = netstoichmat(rn)
        new_rr = reactionrates(rn)
        new_nr = numreactions(rn)
        isequal(sum([new_stoichmat[isp ,i] * new_rr[i] for i in 1:new_nr]), isa(sp, GasSpecies) || isa(sp, FictiousSpecies) ? Num(0.0) : sum([stoichmat[isp ,i] * rr[i] for i in 1:nr]))
    end)
end


function rename_tstate(text::AbstractString; without_site=false)
    if without_site
        re = r"(?<before>(([A-Z]+[1-9]?)+|ele))-(?<after>(([A-Z]+[1-9]?)+|ele|\w|$))"
    else
        re = r"(?<before>(([A-Z]+[1-9]?)+|ele))-(?<after>(([A-Z]+[1-9]?)+|ele|\*?_[a-z]))"
    end
    old_text = text
    text = replace(text, re => s"\g<before>Δ\g<after>") 
    while old_text ≠ text
        old_text = text
        text = replace(text, re => s"\g<before>Δ\g<after>") 
    end
    return text
end
