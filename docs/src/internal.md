# Internal Documentation

## Corrections

```@docs
CatmapInterface.first_order_adsorbate_interaction
CatmapInterface.ideal_adsorbate_interaction
CatmapInterface.ideal_gas
CatmapInterface.harmonic_adsorbate
CatmapInterface.simple_electrochemical
CatmapInterface.hbond_surface_charge_density
CatmapInterface._get_echem_corrections
CatmapInterface._get_interaction_term
```

## Interface

```@docs
CatmapInterface.AdsorbateInteractionParams
CatmapInterface.TState
CatmapInterface.ParsedReaction
CatmapInterface.parse_reaction
CatmapInterface.parse_reactant_sum
CatmapInterface.parse_energy_table
CatmapInterface.specieslist
CatmapInterface.findspecies
CatmapInterface._get_adsorbate_interaction_params
```

## Reaction Network

```@docs
CatmapInterface.ratelaw_TS
CatmapInterface.compute_free_energies!
```

## Statistical Models

```@docs
CatmapInterface.MoleculeGeometry
CatmapInterface.IdealGas
CatmapInterface.HarmonicPhase
CatmapInterface.enthalpy
CatmapInterface.entropy
```

## Utilities

```@docs
CatmapInterface.get_molecule_spec
CatmapInterface.get_ideal_gas_params
CatmapInterface.instantiate_catmap_template!
```