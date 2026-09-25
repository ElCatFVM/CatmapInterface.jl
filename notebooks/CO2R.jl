### A Pluto.jl notebook ###
# v1.0.3

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

# ╔═╡ 8bd435fe-89c6-4b96-99fe-0ae42b8bb6bd
if isdefined(Main, :PlutoRunner)
    using Pkg
    docsdir = joinpath(@__DIR__, "..", "docs")
    if isdir(docsdir)
        Pkg.activate(docsdir)
    end
end


# ╔═╡ 91ac9e35-71eb-4570-bef7-f63c67ce3881
begin
    using Revise
    using LiquidElectrolytes
    using CatmapInterface
    using Catalyst: unknowns
    using VoronoiFVM
    using LessUnitful
    using ExtendableGrids, GridVisualize
    using DelimitedFiles
    using PlutoUI, HypertextLiteral
    using PreallocationTools
    using Latexify
    using Catalyst
    using Printf
    using Test
    using FileIO
    if isdefined(Main, :PlutoRunner)
        using CairoMakie
        default_plotter!(CairoMakie)
        CairoMakie.activate!(type = "svg")
    end
end;

# ╔═╡ beae1479-1c0f-4a55-86e1-ad2b50174c83
md"""
## Setup
"""

# ╔═╡ ab2184fc-0279-46d9-9ee4-88fe3e732789
md"""
### Units
"""

# ╔═╡ 7316901c-d85d-48e9-87dc-3614ab3d81a5
@unitfactors mol dm m s K μm bar Pa eV μF V cm μA mA Å;

# ╔═╡ 6b7cfe87-8190-40a5-8d25-e39ef8d55db5
md"""
### Data
"""

# ╔═╡ 5a146a44-03dc-45f3-ae15-993d11c2edac
begin
    @phconstants N_A c_0 k_B e h
    const F = N_A * e

    const voltages = (-1.5:0.1:-0.0) * V

    # geometrical constants
    const hmin = 1.0e-6 * μm
    const hmax = 1.0 * μm
    const nref = 0
    const L = 80.0 * μm
    const Γ_we = 1
    const Γ_bulk = 2

    # bulk constants
    const pH = 6.8
    const T = 298.0 * K
    const Hcp_CO = 9.7e-6 * mol / (m^3 * Pa)
    const Hcp_CO2 = 3.3e-4 * mol / (m^3 * Pa)
    # species not involved in reactions
    const ikplus = 1
    # species involved in buffer reactions but not in surface reactions
    const ibufferstart = 2
    const ihplus = 2
    const ihco3 = 3
    const ico3 = 4
    # species involved in buffer reactions and surface reactions
    const isurfacestart = 5
    const ico2 = 5
    const iohminus = 6
    const ibufferend = 6
    # species involved in surface reactions but not in buffer reactions
    const ico = 7
    const nc = 7
    ## reaction rate constants for bulk reactions
    ### CO2 + OH- <=> HCO3-
    const kbe1 = 4.44e7 / (mol / dm^3)
    const kbf1 = 5.93e3 / (mol / dm^3) / s
    const kbr1 = kbf1 / kbe1
    ### HCO3- + OH- <=> CO3-- + H2O
    const kbe2 = 4.66e3 / (mol / dm^3)
    const kbf2 = 1.0e8 / (mol / dm^3) / s
    const kbr2 = kbf2 / kbe2
    ### CO2 + H20 <=> HCO3- + H+
    const kae1 = 4.44e-7 * (mol / dm^3)
    const kaf1 = 3.7e-2 / s
    const kar1 = kaf1 / kae1
    ### HCO3- <=> CO3-- + H+
    const kae2 = 4.66e-5 / (mol / dm^3)
    const kaf2 = 59.44e3 / (mol / dm^3) / s
    const kar2 = kaf2 / kae2
    ### autoprotolyse
    const kwe = 1.0e-14 * (mol / dm^3)^2
    const kwf = 2.4e-5 * (mol / dm^3) / s
    const kwr = kwf / kwe
    const aH₂O = 1.0 #* mol/dm^3

    # surface constants
    const S = 9.61e-5 / N_A * (1.0e10)^2 * mol / m^2
    const C_gap = 20 * μF / cm^2
    const ϕ_pzc = 0.16 * V
    const ico_t = 8
    const icooh_t = 9
    const ico2_t = 10
    const isurfaceend = 10
    const na = 3 # CO_t, CO2_t, COOH_t
    const M0 = 18.0153 * ufac"g/mol"
    const v0 = N_A * (8.2 * Å)^3 # 1 / (55.4 * ufac"M") #

    const species_dict = Dict(
        "K⁺" => ikplus,
        "HCO₃⁻" => ihco3,
        "CO₃²⁻" => ico3,
        "CO₂" => ico2,
        "OH⁻" => iohminus,
        "H⁺" => ihplus,
        "CO" => ico,
        "CO_t" => ico_t,
        "COOH_t" => icooh_t,
        "CO2_t" => ico2_t,
    )

    const species_dict_catmap = Dict(
        "OH_g" => iohminus,
        "CO2_aq" => ico2,
        "CO_aq" => ico,
        "CO_t" => ico_t,
        "COOH_t" => icooh_t,
        "CO2_t" => ico2_t,
    )
end;

# ╔═╡ 00947475-c96e-4ecc-a1ef-5be5e3e3c864
begin
    @kwdef struct BulkSpecies
        name::String
        z::Int
        D::Float64
        c_bulk::Union{Nothing, Float64}
        κ::Float64
        a::Float64
        v::Float64
        M::Float64
        color::Symbol
    end
    function BulkSpecies(; name, z, c_bulk = nothing, D, κ = 0.0, a = 0.0, color)
        D *= m^2 / s
        c_bulk = isnothing(c_bulk) ? nothing : c_bulk * mol / dm^3
        a *= Å
        v = N_A * a^3
        M = M0 * v
        return BulkSpecies(name, z, D, c_bulk, κ, a, v, M, color)
    end
    function make_eneutral(
            bulk_species::Vector{BulkSpecies}
            ; name, z, D, κ = 0.0, a = 0.0, v = N_A * (a * Å)^3, M = M0 * v, color
        )
        a *= Å
        c_bulk = -mapreduce(x -> x.c_bulk * x.z, +, bulk_species) / z
        return BulkSpecies(name, z, D, c_bulk, κ, a, v, M, color)
    end
end;

# ╔═╡ ed1812f4-fdab-4fb5-88e1-0ece3c1e26b1
begin
    const bulk = let
        bulk = [
            BulkSpecies(; name = "HCO₃⁻", z = -1, D = 1.185e-9, c_bulk = 0.091, color = :brown),
            BulkSpecies(; name = "CO₃²⁻", z = -2, D = 0.923e-9, c_bulk = 2.68e-5, color = :violet),
            BulkSpecies(; name = "CO₂", z = 0, D = 1.91e-9, c_bulk = 0.033, a = 0.0, color = :red),
            BulkSpecies(; name = "OH⁻", z = -1, D = 5.273e-9, c_bulk = 10^(pH - 14), color = :green),
            BulkSpecies(; name = "H⁺", z = 1, D = 9.31e-9, c_bulk = 10^(-pH), a = 0.0, color = :gray),
            BulkSpecies(; name = "CO", z = 0, D = 2.23e-9, c_bulk = 0.0, a = 0.0, color = :blue),
        ]
        push!(bulk, make_eneutral(bulk; name = "K⁺", z = 1, D = 1.957e-9, a = 8.2, color = :orange))
        sort(bulk, by = x -> species_dict[x.name])
    end
end;

# ╔═╡ 06f52599-7006-4a5c-ba86-0b668b6952c9
begin
    function create_markdown(bulk::Vector{BulkSpecies})
        table = """
        | Name | z | D | c_bulk | a | v | M | κ | color |
        |------|---|---|--------|---|---|---|---|-------|
        """
        for sp in bulk
            (; name, z, D, c_bulk, a, v, M, κ, color) = sp
            table *= @sprintf("| %s | %i | %1.3e | %1.3e | %1.3e | %1.3e | %1.3e | %1.2f | %s |\n", name, z, D, c_bulk, a, v, M, κ, color)
        end
        return Markdown.parse(table)
    end
    create_markdown(bulk)
end

# ╔═╡ ca22e3fe-5cb7-4910-b9fa-890fd2d20e4b
md"""
### Solver Control
"""

# ╔═╡ 4c95d645-f909-492b-a425-927c093ae31a
solver_control = (;
    max_round = 4,
    maxiters = 20,
    tol_round = 1.0e-9,
    verbose = "a",
    reltol = 1.0e-8,
    tol_mono = 1.0e-10,
)

# ╔═╡ 4b64e168-5fe9-4202-9657-0d4afc237ddc
md"""
### Reaction Description
"""

# ╔═╡ de2c826d-6c05-47cf-b5f5-44a00ea9889c
md"""
#### Buffer System
"""

# ╔═╡ d8f00649-e2ed-4bdd-853f-05268f0d5353
md"""
Consider the bicarbonate buffer system in base and acid as well as autoprotolysis of water:
"""

# ╔═╡ 47b36c81-b57e-4dd0-a22f-999e4fd3ac9f
begin
    @variables t
    @species HCO₃⁻(t), CO₃²⁻(t), CO₂(t), OH⁻(t), H⁺(t)
    @parameters γHCO₃⁻ γCO₃²⁻ γCO₂ γOH⁻ γH⁺
    buffer_rn = @reaction_network buffer begin
        ($kbf1 * γCO₂ * γOH⁻, $kbr1 * γHCO₃⁻), CO₂ + OH⁻ <--> HCO₃⁻
        ($kbf2 * γHCO₃⁻ * γOH⁻, $kbr2 * $aH₂O * γCO₃²⁻), HCO₃⁻ + OH⁻ <--> CO₃²⁻
        ($kaf1 * $aH₂O * γCO₂, $kar1 * γHCO₃⁻ * γH⁺), CO₂ <--> HCO₃⁻ + H⁺
        ($kaf2 * γHCO₃⁻, $kar2 * γCO₃²⁻ * γH⁺), HCO₃⁻ <--> CO₃²⁻ + H⁺
        ($kwf * $aH₂O, $kwr * γH⁺ * γOH⁻), ∅ <--> H⁺ + OH⁻
    end
    odesys_buffer = convert(ODESystem, buffer_rn; combinatoric_ratelaws = false)
    const f_buffer! = CatmapInterface.generate_function(
        buffer_rn;
        dvs = [H⁺, HCO₃⁻, CO₃²⁻, CO₂, OH⁻],
        ps = [γH⁺, γHCO₃⁻, γCO₃²⁻, γCO₂, γOH⁻]
    )
    latexify(buffer_rn; env = :chemical)
end

# ╔═╡ 1e877f17-0219-45f1-b640-3a25ae085dbd
latexify(odesys_buffer) # hide

# ╔═╡ 8a1047fa-e483-40d9-8904-7576f30acfb4
begin
    const γ_cache = DiffCache(zeros(nc), 12)

    function reaction(
            f,
            u,
            node,
            data
        )
        (; ip, iϕ, v0, v, M0, M, κ, ε_0, ε, RT, nc, pscale, p_bulk) = data

        # compute activity coefficients according to the approach in Ringe et al.
        γ = get_tmp(γ_cache, u[ico2])
        γ .= 1.0 / (1 - v[ikplus] * u[ikplus] / (mol / dm^3))

        # compute activity coefficients according to the approach in Dreyer et al.
        # p = u[ip] * pscale-p_bulk
        # c0, bar_c = c0_barc(u, data)
        # for ic in 1:nc
        #	Mrel = M[ic] / M0
        #	barv = v[ic] + κ[ic] * v0
        #	tildev = barv - Mrel * v0
        # 	γ[ic] = exp(tildev * p / (RT)) * (bar_c / c0)^Mrel*(1/bar_c) /v0
        # end


        @views f_buffer!(
            f[ibufferstart:ibufferend],
            u[ibufferstart:ibufferend],
            γ[ibufferstart:ibufferend],
            nothing
        )
        return nothing
    end
end;

# ╔═╡ 8912f990-6b02-467a-bd11-92f94818b1c7
md"""
#### Surface Reactions
"""

# ╔═╡ a8157cc1-1761-4b11-a37c-9e12a9ca695e
md"""
A microkinetic modeling approach is taken:

The reaction mechanism for the $CO_2$ reduction is divided into four elementary reactions at the electrode surface:

1. Adsorption of $CO_2$ molecules at the oxygen atoms
${CO_2}_{(aq)} + * \rightleftharpoons {CO_2 *}_{(ad)}$

2. First proton-coupled electron transfer
${CO_2*}_{(ad)} + H_2O_{(l)} + e^- \rightleftharpoons COOH*_{(ad)} + OH^-_{(aq)}$

3. Second proton-coupled electron transfer
$COOH*_{(ad)} + e^- \rightleftharpoons CO*_{(ad)} + OH^{-}_{(aq)}$
with the transition state: $*CO-OH^{TS}$

4. Desorption of $CO$
$*CO_{(ad)} \rightleftharpoons CO_{(aq)} + *$
"""

# ╔═╡ 489ead3b-04b8-44bb-9d73-7b1d13cf5346
const symbolic_formation_energies = true


# ╔═╡ 6b5cf93c-0df3-4a18-8786-502361736838
begin
    catmap_params = parse_catmap_input(joinpath(@__DIR__, "..", "data", "models", "Au", "catmap_CO2R_template.mkm"))
    rn, _ = create_reaction_network(catmap_params; symbolic_formation_energies)
    odesys0 = convert(ODESystem, rn; combinatoric_ratelaws = false)
    odesys = CatmapInterface.liquidize(odesys0, catmap_params)
    vars = unknowns(odesys)
    const f_microkinetics! = CatmapInterface.generate_function(
        odesys;
        dvs = sort(vars, by = x -> species_dict_catmap[string(operation(x))])
    )
    const pidx = paramsidx(odesys)
    latexify(odesys)
end

# ╔═╡ d2c0642d-dfa5-4a76-bd36-ac4a735a3299
md"""
##### Reaction Rates
"""

# ╔═╡ 06d45088-ab8b-4e5d-931d-b58701bf8464
md"""
For the calculation of the reaction rates a __mean field approach__ based is applied.

The acitivities of H⁺, OH⁻ and e⁻ are set to be zero. The dependence on the rates on the pH-value, applied voltage and surface charges are contained in the reaction rate constants.

The surface charging relation $σ = σ(U)$ is given by the Robin boundary condition

$σ(U) = C_{gap} (ϕ_{we} - ϕ_{pzc} - ϕ^\ddagger)$

where the gap capacitance between the working electrode and the reaction plane ($\ddagger$) is given by $C_{gap} = 20~μF/cm^2$. The potential of zero current is measured to be $ϕ_{pzc} = 0.16~V$.

__Question is the pH-dependence only in the reaction rate constants (i.e. activity of OH⁻ must be set to 0)?__
"""

# ╔═╡ 91113083-d80e-4528-be41-82d10f6860fc
begin
    const nparams = maximum([p.second for p in pairs(pidx)])
    const ps_cache = DiffCache(zeros(nparams), 13)
    const us_cache = DiffCache(zeros(isurfaceend - isurfacestart + 1), 13)

    function we_breactions(
            f,
            u,
            bnode,
            data
        )
        (; ip, iϕ, v0, v, M0, M, κ, RT, nc, pscale, p_bulk, ϕ_we) = data

        γ_co2 = 1.0 / (1 - v[ikplus] * u[ikplus] / (mol / dm^3))
        γ_co = 1.0 / (1 - v[ikplus] * u[ikplus] / (mol / dm^3))
        σ = C_gap * (ϕ_we - u[iϕ] - ϕ_pzc)
        local_pH = -log10(u[ihplus] / (mol / dm^3))

        ps = get_tmp(ps_cache, u[iϕ])
        ps[pidx[:σ]] = σ
        ps[pidx[:γCO2_aq]] = γ_co2
        ps[pidx[:aH2O_g]] = aH₂O
        ps[pidx[:ϕ]] = u[iϕ]
        ps[pidx[:ϕ_we]] = ϕ_we
        ps[pidx[:local_pH]] = local_pH
        ps[pidx[:γCO_aq]] = γ_co
        ps[pidx[:βCOOHΔH2OΔele_t]] = 0.59
        ps[pidx[:prefactor_CO2_t]] = 1.0e13
        ps[pidx[:prefactor_COOH_t]] = 1.0e13
        ps[pidx[:prefactor_COOHΔH2OΔele_t]] = 1.0e13
        ps[pidx[:prefactor_CO_g]] = 1.0e8
        if symbolic_formation_energies
            ps[pidx[:ECO2_t]] = 0.657600203 * e
            ps[pidx[:ECOOH_t]] = 0.128214079 * e
            ps[pidx[:ECO_t]] = -0.02145440850567823 * e
            ps[pidx[:ECOOHΔH2OΔele_t]] = 0.95 * e
        end
        @views f_microkinetics!(
            f[isurfacestart:isurfaceend],
            u[isurfacestart:isurfaceend],
            ps,
            nothing
        )

        # conversion from turnover frequency (appropriate for change in coverage) to
        # production rate (per unit area) (approprite for change in concentration)
        # by S = number of free catalyst sites in mole per unit area
        f[ico2] *= S
        f[iohminus] *= S
        return f[ico] *= S
    end
end

# ╔═╡ d0093605-0e35-4888-a93c-8456c698e6f0
md"""
##### Reaction Rate Constants
"""

# ╔═╡ f0b5d356-6b97-4878-98de-bee5f380d41a
md"""
The details of the rate constant calculations can be found in the documentation of CatMAP and the CatMAP specification for this model: "./models/Au/catmap_CO2R_template.mkm". Here only a few import points are highlighted.

According to the transition state theory the kinetic model assumes that the dynamics of a reaction from an educt complex to a product complex can be projected onto a path parametrized by a reaction coordinate. On the reaction path an activated complex can be identified. The rate constant is interpreted according to the Eyring equation
```math
k = κ \frac{k_B T}{h} exp\left(-ΔG^\ddagger\right)
```
where the $κ$ is a transmission coefficient and $ΔG^\ddagger$ is the change in the Gibbs free energy between the educt and the transition state for canonical ensembles in the respective states. For the transition state only degrees of freedom perpendicular to the reaction coordinate are considered and treated as harmonic degrees of freedom.

Only the second concerted proton-electron transfer (CPET) is assumed to pass through a distinct activated complex. For the other reactions no kinetic barrier is assumed and $ΔG^\ddagger$ is interpreted as the change in the relative Gibbs free energy of formation between the educt and product state $Δᵣμ^{⦵}$ in case that the reaction does not occur spontaneously $Δᵣμ^{⦵} > 0$.
The (entropy-free) pre-exponential factors $κ k_B T/ h$ are approximated as $10^{13}$ for all but the adsorption of CO where it set to $10^8$.
"""

# ╔═╡ e3eda42f-e2f3-4c10-81c4-610246ca528d
md"""
The Gibbs free energies are based on potential energies calculated by applying DFT-methods. In order to obtain thermodynamic values ab-initio statistical models are used.
"""

# ╔═╡ 7b87aa2a-dbaf-441c-9ad7-444abf15f664
md"""
The free energies $ΔG_f$ of the surface species are corrected according to the (excess) surface charge density $σ$ by a fitted quadratic model:

$ΔG_f(σ) = a_σ~σ + b_σ~σ^2$
"""

# ╔═╡ 6e4c792e-e169-4b49-89d0-9cf8d5ac8c04
md"""
### Nernst-Planck Half-Cell
"""

# ╔═╡ e7e0eb0d-fe3e-4f1d-876f-cc13a9aaf84c
grid = let
    X = geomspace(0, L, hmin, hmax)
    simplexgrid(X)
end

# ╔═╡ e510bce3-d33f-47bb-98d6-121eee8f2252
celldata = ElectrolyteData(;
    nc = nc,
    na = na,
    z = getproperty.(bulk, :z),
    D = getproperty.(bulk, :D),
    T = T,
    eneutral = false,
    κ = getproperty.(bulk, :κ),
    c_bulk = getproperty.(bulk, :c_bulk),
    v0 = v0,
    v = getproperty.(bulk, :v),
    M0 = M0,
    M = getproperty.(bulk, :M),
    Γ_we = Γ_we,
    Γ_bulk = Γ_bulk,
);

# ╔═╡ dc203e95-7763-4b13-8408-038b933c5c9c
function halfcellbc(
        f,
        u,
        bnode,
        data
    )

    (; Γ_we, Γ_bulk, ϕ_we, iϕ) = data

    bulkbcondition(f, u, bnode, data; region = Γ_bulk)

    # Robin b.c. for the Poisson equation
    boundary_robin!(f, u, bnode, iϕ, Γ_we, C_gap, C_gap * (ϕ_we - ϕ_pzc))

    if bnode.region == Γ_we
        we_breactions(f, u, bnode, data)
    end
    return nothing
end;

# ╔═╡ 842b074b-f808-48d8-8dc5-110ddd907f90
md"""
## Results
"""

# ╔═╡ 84d1270b-8df5-4d5d-a153-da4ffdb1d283
function simulate_CO2R(grid, celldata; voltages = (-1.5:0.1:0.0) * V, kwargs...)
    kwargs = merge(solver_control, kwargs)
    cell = PNPSystem(grid; bcondition = halfcellbc, reaction = reaction, celldata)
    ivresult = ivsweep(cell; voltages, store_solutions = true, kwargs...)

    return cell, ivresult
end;

# ╔═╡ 11b12556-5b61-42c2-a911-4ea98a0a1e85
cell, result = simulate_CO2R(grid, celldata; voltages);

# ╔═╡ 114d2324-5289-4e44-8d77-736a9bdec365
md"""
Show only pH: $(@bind useonly_pH PlutoUI.CheckBox(default=false))
"""

# ╔═╡ 659091d3-60b2-4158-80e2-cd28a492e870
(~, default_index) = findmin(abs, result.voltages .+ 0.9 * ufac"V");

# ╔═╡ 3bcb8261-5b98-4f4d-a9fe-fb71d5c5b476
md"""
$(@bind vindex PlutoUI.Slider(1:5:length(result.voltages), default=default_index))
"""

# ╔═╡ c4876d26-e841-4e28-8303-131d4635fc23
md"""
Potential at the working electrode 
$(vshow = result.voltages[vindex]; @sprintf("%+1.4f", vshow))
"""

# ╔═╡ c1d2305e-fb8b-4845-a414-08fff84aa9b0
md"""
### Plotting Functions
"""

# ╔═╡ 2ce5aa45-4aa5-4c2a-a608-f581266e55f0
begin
    function addplot(vis, sol, vshow)
        species = getproperty.(bulk, :name)
        colors = getproperty.(bulk, :color)

        scale = 1.0 / (mol / dm^3)
        title = @sprintf("Φ_we=%+1.2f [V vs. SHE]", vshow)

        return if useonly_pH
            i = findfirst(isequal("H⁺"), species)
            scalarplot!(
                vis,
                grid.components[XCoordinates] .+ 1.0e-14,
                log10.(sol[ihplus, :] * scale),
                color = colors[i],
                label = species[i],
                clear = true,
                title = title
            )
        else
            scalarplot!(
                vis,
                grid.components[XCoordinates] .+ 1.0e-14,
                log10.(sol[1, :] * scale),
                color = colors[1],
                label = species[1],
                clear = true,
                title = title
            )
            for ia in 2:nc
                scalarplot!(
                    vis,
                    grid.components[XCoordinates] .+ 1.0e-14,
                    log10.(sol[ia, :] * scale),
                    color = colors[ia],
                    label = species[ia],
                    clear = false,
                )
            end
        end
    end

    # function addplot(vis, df)

    # 	function extract_interpolation(df, i)
    # 		X = collect(skipmissing(df[!, 2*i-1]))
    # 		I = sortperm(X)
    # 		X .= X[I]
    # 		Y = collect(skipmissing(df[!, 2*i]))[I]
    # 		linear_interpolation(X, Y, extrapolation_bc=Line())
    # 	end

    # 	species = getproperty(bulk, :name)
    # 	colors = [:orange, :brown, :violet, :red, :blue, :green, :gray]

    # 	knots = grid.components[XCoordinates] .+ 1.0e-14
    # 	sol = [extract_interpolation(df, i) for i in 1:nc]

    # 	if useonly_pH
    # 		scalarplot!(vis,
    # 				    knots,
    # 				    log10.(sol[ihplus].(knots)),
    # 				    color = colors[ihplus],
    # 				    clear = false,
    # 					linewidth = 0,
    # 					label = "",
    # 					markershape = :cross,
    # 					markersize = 8,
    # 					markevery = 20)
    # 	else
    # 		for ia = 1:nc
    # 			scalarplot!(vis,
    # 					    knots,
    # 					    log10.(sol[ia].(knots)),
    # 					    color = colors[ia],
    # 					    clear = false,
    # 						linewidth = 0,
    # 						label = "",
    # 						markershape = :cross,
    # 						markersize = 8,
    # 						markevery = 20)
    # 		end
    # 	end
    # end

    function plot1d(result, celldata, vshow; df_compare = nothing)
        tsol = LiquidElectrolytes.voltages_solutions(result)
        vis = GridVisualizer(;
            size = (600, 300),
            clear = true,
            legend = :rt,
            limits = (-14, 2),
            xlimits = (10.0e-12, 80 * μm),
            xlabel = "Distance from electrode [m]",
            ylabel = "log c(aᵢ)",
            xscale = :log,
        )
        addplot(vis, tsol(vshow), vshow)
        if !isnothing(df_compare)
            addplot(vis, df_compare)
        end
        return reveal(vis)
    end

    function plot1d(result, celldata)
        tsol = LiquidElectrolytes.voltages_solutions(result)
        vis = GridVisualizer(;
            size = (600, 300),
            clear = true,
            legend = :rt,
            limits = (-14, 2),
            xlimits = (10.0e-12, 80 * μm),
            xlabel = "Distance from electrode [m]",
            ylabel = "log c(aᵢ)",
            xscale = :log,
        )

        vrange = result.voltages[end:-5:1]
        movie(vis, file = "concentrations.gif", framerate = 3) do vis
            for vshow_it in vrange
                addplot(vis, tsol(vshow_it), vshow_it)
                reveal(vis)
            end
        end
        return isdefined(Main, :PlutoRunner) && LocalResource("concentrations.gif")
    end
end

# ╔═╡ 5dd1a1e6-7db1-479e-a684-accec53ce06a
plot1d(result, celldata, vshow)

# ╔═╡ 15fadfc2-3cf8-4fda-9aed-a79c602b1d51
plot1d(result, celldata)

# ╔═╡ d5ab1a28-3a60-49d9-bb3e-ca589b1c79fd
begin
    curr(J, ix) = [F * abs(j[ix]) for j in J]

    function plotcurr(result; df = nothing)
        scale = 1 / (mol / dm^3)
        volts = result.voltages[result.voltages .< -0.4]
        vis = GridVisualizer(;
            size = (600, 400),
            tilte = "IV Curve",
            xlabel = "Φ_WE/(V vs. SHE)",
            ylabel = "I/(mA/cm²)",
            legend = :lb,
            yscale = :log,
        )

        scalarplot!(
            vis,
            volts,
            curr(result.j_we, iohminus)[result.voltages .< -0.4] .* cm^2 / mA;
            color = :green,
            clear = false,
            linestyle = :solid,
            label = "e⁻, we"
        )
        if !isnothing(df)
            scalarplot!(
                vis,
                df[:voltage],
                df[:current],
                clear = false,
                linewidth = 0,
                markershape = :cross,
                markersize = 8,
                markevery = 1,
                color = :red,
                label = "Ringe et. al"
            )
        end

        if isdefined(Main, :PlutoRunner)
            GridVisualize.save("CO2Rcompare.png", vis)
        end
        return reveal(vis)
    end
end

# ╔═╡ 1cd669ac-05eb-48b2-b457-8c395cd5807d
let
    table = readdlm(joinpath(@__DIR__, "..", "data", "IV-Ringe-digitized.csv"), ',', Float64, '\n')
    df = Dict(:voltage => table[:, 1], :current => table[:, 2])
    plotcurr(result; df = df)
end

# ╔═╡ 686ac3dc-c191-4575-ba0c-d4c2551474b5
md"""
### Regression Test
"""

# ╔═╡ 32eb1122-5013-4a8e-be54-18a30c151515
begin
    sresult = load(joinpath(@__DIR__, "..", "data", "regressionresults-CMI-v0.3.0.jld2"))["regressionresults"]

    vidxs_result = [findfirst(isequal(v), result.voltages) for v in voltages[1:(end - 1)]]
    vidxs_sresult = [findfirst(isequal(v), sresult.voltages) for v in voltages[1:(end - 1)]]

    if any(isnothing.(vidxs_sresult))
        throw(ArgumentError("For the full regression test use the applied voltages  -1.5:0.1:0.0"))
    end

    @testset begin
        @testset "Concentrations" begin
            @testset "$(bulk[ia].name)" for ia in 1:nc
                @testset "U=$(result.voltages[vidx_result])" for (vidx_result, vidx_sresult) in zip(vidxs_result, vidxs_sresult)
                    @test all(
                        isapprox(
                            result.solutions[vidx_result][ia, :], sresult.solutions[vidx_sresult][ia, :],
                            rtol = 1.0e-5
                        )
                    )
                end
            end
        end

        @testset "Currents" begin
            for (vidx_result, vidx_sresult) in zip(vidxs_result, vidxs_sresult)
                for (j_result, j_sresult) in zip(
                        result.j_we[vidx_result][iohminus],
                        sresult.j_we[vidx_sresult][iohminus]
                    )
                    @test isapprox(j_result, j_sresult, atol = 1.0e-13)
                end
            end
        end
    end
end;

# ╔═╡ d0985ca6-fef5-4b67-9ad6-f51d84b595b4
TableOfContents(title = "📚 Table of Contents", indent = true, depth = 4, aside = true)

# ╔═╡ 8ae53b8a-0fb3-4c1c-8e5f-a3782a85141c
begin
    hrule() = html"""<hr>"""
    function highlight(mdstring, color)
        return htl"""<blockquote style="padding: 10px; background-color: $(color);">$(mdstring)</blockquote>"""
    end

    macro important_str(s)
        return :(highlight(Markdown.parse($s), "#ffcccc"))
    end
    macro definition_str(s)
        return :(highlight(Markdown.parse($s), "#ccccff"))
    end
    macro statement_str(s)
        return :(highlight(Markdown.parse($s), "#ccffcc"))
    end

    html"""
            <style>
             h1{background-color:#dddddd;  padding: 10px;}
             h2{background-color:#e7e7e7;  padding: 10px;}
             h3{background-color:#eeeeee;  padding: 10px;}
             h4{background-color:#f7f7f7;  padding: 10px;}
            
    	     pluto-log-dot-sizer  { max-width: 655px;}
             pluto-log-dot.Stdout { background: #002000;
    	                            color: #10f080;
                                    border: 6px solid #b7b7b7;
                                    min-width: 18em;
                                    max-height: 300px;
                                    width: 675px;
                                        overflow: auto;
     	                           }
    	
        </style>
    """
end

# ╔═╡ Cell order:
# ╠═8bd435fe-89c6-4b96-99fe-0ae42b8bb6bd
# ╠═91ac9e35-71eb-4570-bef7-f63c67ce3881
# ╟─beae1479-1c0f-4a55-86e1-ad2b50174c83
# ╟─ab2184fc-0279-46d9-9ee4-88fe3e732789
# ╠═7316901c-d85d-48e9-87dc-3614ab3d81a5
# ╟─6b7cfe87-8190-40a5-8d25-e39ef8d55db5
# ╠═5a146a44-03dc-45f3-ae15-993d11c2edac
# ╠═00947475-c96e-4ecc-a1ef-5be5e3e3c864
# ╠═ed1812f4-fdab-4fb5-88e1-0ece3c1e26b1
# ╟─06f52599-7006-4a5c-ba86-0b668b6952c9
# ╟─ca22e3fe-5cb7-4910-b9fa-890fd2d20e4b
# ╟─4c95d645-f909-492b-a425-927c093ae31a
# ╟─4b64e168-5fe9-4202-9657-0d4afc237ddc
# ╟─de2c826d-6c05-47cf-b5f5-44a00ea9889c
# ╟─d8f00649-e2ed-4bdd-853f-05268f0d5353
# ╠═47b36c81-b57e-4dd0-a22f-999e4fd3ac9f
# ╟─1e877f17-0219-45f1-b640-3a25ae085dbd
# ╠═8a1047fa-e483-40d9-8904-7576f30acfb4
# ╟─8912f990-6b02-467a-bd11-92f94818b1c7
# ╟─a8157cc1-1761-4b11-a37c-9e12a9ca695e
# ╠═489ead3b-04b8-44bb-9d73-7b1d13cf5346
# ╠═6b5cf93c-0df3-4a18-8786-502361736838
# ╟─d2c0642d-dfa5-4a76-bd36-ac4a735a3299
# ╟─06d45088-ab8b-4e5d-931d-b58701bf8464
# ╠═91113083-d80e-4528-be41-82d10f6860fc
# ╟─d0093605-0e35-4888-a93c-8456c698e6f0
# ╟─f0b5d356-6b97-4878-98de-bee5f380d41a
# ╟─e3eda42f-e2f3-4c10-81c4-610246ca528d
# ╟─7b87aa2a-dbaf-441c-9ad7-444abf15f664
# ╟─6e4c792e-e169-4b49-89d0-9cf8d5ac8c04
# ╠═e7e0eb0d-fe3e-4f1d-876f-cc13a9aaf84c
# ╠═e510bce3-d33f-47bb-98d6-121eee8f2252
# ╠═dc203e95-7763-4b13-8408-038b933c5c9c
# ╟─842b074b-f808-48d8-8dc5-110ddd907f90
# ╠═84d1270b-8df5-4d5d-a153-da4ffdb1d283
# ╠═11b12556-5b61-42c2-a911-4ea98a0a1e85
# ╟─114d2324-5289-4e44-8d77-736a9bdec365
# ╟─659091d3-60b2-4158-80e2-cd28a492e870
# ╟─c4876d26-e841-4e28-8303-131d4635fc23
# ╟─3bcb8261-5b98-4f4d-a9fe-fb71d5c5b476
# ╠═5dd1a1e6-7db1-479e-a684-accec53ce06a
# ╠═15fadfc2-3cf8-4fda-9aed-a79c602b1d51
# ╠═1cd669ac-05eb-48b2-b457-8c395cd5807d
# ╟─c1d2305e-fb8b-4845-a414-08fff84aa9b0
# ╟─2ce5aa45-4aa5-4c2a-a608-f581266e55f0
# ╠═d5ab1a28-3a60-49d9-bb3e-ca589b1c79fd
# ╟─686ac3dc-c191-4575-ba0c-d4c2551474b5
# ╠═32eb1122-5013-4a8e-be54-18a30c151515
# ╟─d0985ca6-fef5-4b67-9ad6-f51d84b595b4
# ╟─8ae53b8a-0fb3-4c1c-8e5f-a3782a85141c
