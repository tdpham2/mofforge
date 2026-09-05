# Adsorption sites and guest placement

These lessons use core dependencies and the bundled IRMOF-1 structure.
They construct initial host–guest geometries. They do not calculate adsorption
energies, equilibrium loadings, or molecular dynamics trajectories.

## 1. Discover candidate sites

```bash
python examples/adsorbate/site_discovery.py
```

[site_discovery.py](site_discovery.py) demonstrates both site strategies:

1. Infer framework bonds and sample void space on a **1.0 Å grid**.
2. Return up to five void sites, including fractional/Cartesian coordinates and
   the nearest framework-atom distance.
3. Query the framework for under-coordinated metal sites.
4. Apply that heuristic to an explicitly constructed CuO3 motif as a positive
   control.
5. Set an impossible 100 Å clearance to demonstrate an empty site list.

The default framework has zero detected open-metal sites; the illustrative Cu
motif has one because Cu has three inferred neighbors, below the implementation's
minimum of four. `illustrative_CuO3.cif` is a teaching motif, not an experimental
MOF or an optimized complex.

The detector also returns zero sites for the existing MOF-74 test fixture.
It tests under-coordination relative to generic element ranges, not every
chemically meaningful definition of an open metal site. Empty results must be
interpreted alongside coordination and input preparation.

**Try:** reduce grid spacing to 0.5 Å and compare sites and runtime.
The reported nearest-atom distance is not a probe-accessible pore radius from
a full pore-characterization calculation. Coarse grids and clustering change
which candidate centers are returned.

## 2. Place built-in and custom molecules

```bash
python examples/adsorbate/placement.py --count 2 --seed 42
```

[placement.py](placement.py) starts from a clean 424-atom parent for each case:

1. Load CO2 from the built-in molecule menu.
2. Detect sites, filter their centers using a **3.0 Å** minimum intermolecular
   distance, and orient molecules with a local seeded RNG.
3. Save `CO2_loaded.cif` and inspect actual placement count and validation.
4. Construct a custom water `Crystal` and place it at an explicitly selected
   site with fixed orientation, producing `water_loaded.cif`.

The default CO2 case places two molecules and produces **430 atoms with zero
reported clashes** in the checked environment. `adsorbate_indices` identifies
each guest's atom indices in the combined crystal. These indices belong to
that in-memory ordering.

`requested` and `placed` may differ: filtering can remove candidate sites.
The spacing option constrains site centers, not every atom–atom separation.
Validation can therefore report clashes after a placement completes.
You can supply one `site` or a `sites` list; those arguments are mutually
exclusive. An automatically detected empty site set raises an error.

**Try:** compare a fixed and random orientation at the same site, or increase
`--count` and inspect actual count and contacts. Preserve the host and guest
geometry for a subsequent relaxation or adsorption simulation using the
appropriate external software.
