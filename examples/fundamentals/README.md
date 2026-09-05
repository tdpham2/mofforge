# Structure fundamentals

Start with a core installation; no optional extras or downloads are needed.
All distances below are in ångströms. Follow these lessons before selecting
fragments or trusting a coordination-based analysis.

## 1. Read, inspect, and export

```bash
python examples/fundamentals/structure_io.py
```

Read [structure_io.py](structure_io.py) alongside the output:

1. Load the bundled IRMOF-1 CIF into a `Crystal`.
2. Call `infer_bonds(periodic=True)` to create its graph, including cell-boundary
   contacts. Loading a CIF alone does not perform this step.
3. Load the ten-atom phenylene XYZ fragment, write both formats, and reparse them.
4. Inspect a Zn coordination number through `coordination_number`, which accounts
   for distinct periodic images rather than just counting graph edges.

The reference has **424 atoms, 512 graph bonds, and four neighbors at the first
Zn site**. `framework.cif` preserves the cell; `phenylene.xyz` contains molecular
Cartesian coordinates. Neither interchange format is a lossless serialization
of all in-memory metadata. The accompanying JSON manifests record export history.

**Try:** compare `frac_coords`, `cart_coords`, and `lattice.abc` for UiO-66.
If you replace the CIF with a partially occupied structure, resolve its disorder
before bond inference; do not silently discard occupancy information.

## 2. Understand periodic geometry

```bash
python examples/fundamentals/periodic_geometry.py
```

[periodic_geometry.py](periodic_geometry.py) constructs an explicitly illustrative
two-atom structure in a skewed cell. Fractional x coordinates 0.98 and 0.02
describe atoms on opposite sides of the conventional cell boundary.

1. Compare direct Cartesian separation (**7.68 Å**) with the minimum-image
   separation (**0.32 Å**).
2. Add an integer lattice translation to both atoms.
3. Wrap the translated structure and inspect `wrapped.cif`.

Periodic distances use the full lattice metric. Component-wise rounding in
fractional coordinates is not generally the nearest-image solution for skewed
cells. The exported motif is a coordinate demonstration, not a proposed material.

**Try:** change the off-diagonal lattice entries and compare the two distances.
When moving atoms through the API, use coordinate setters so associated geometry
is refreshed.

## 3. Inspect bond inference assumptions

```bash
python examples/fundamentals/bonding_rules.py
```

[bonding_rules.py](bonding_rules.py) compares default periodic inference,
nonperiodic inference, removal of all graph bonds, and an explicit Zn–O rule.

The custom rule deliberately lowers only the Zn–O cutoff to **1.0 Å**. Other
default rules are retained so the neighbor-search cutoff still covers other
species. This is a diagnostic perturbation, not a recommended bonding model.
The default graph has 512 bonds; the perturbed graph has 384 and loses Zn–O
coordination. `remove_bonds` changes the graph, not the atom coordinates.

**Try:** use a chemically motivated cutoff and inspect the affected neighbors.
A change in inferred connectivity changes search matches, solvent classification,
and open-metal-site detection downstream.

## 4. Compare symmetry and build a supercell

```bash
python examples/search/symmetry_analysis.py --seed 42
```

The existing [symmetry lesson](../search/symmetry_analysis.py) uses the small
NiPyC fixture, finds a PyC linker, modifies one location, and analyzes the
space group at `symprec=0.1` Å. It also checks that wrapping preserves periodic
positions, then replicates the modified structure with `[2, 1, 1]`.

Expect **27 parent atoms → 30 modified atoms → 60 supercell atoms**. Output files
are `NiPyC_CH3.cif` and `NiPyC_CH3_supercell.cif`. Replication uses a copy of the
pymatgen structure followed by fresh bond inference. Reported space-group labels
depend on the tolerance and numerical geometry.

**Try:** compare `symprec=0.01` and `0.1`. Interpret reduced symmetry after
selective substitution as a property of the chosen model, not a failed build.
