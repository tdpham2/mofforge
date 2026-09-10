# Polymer box construction

Mofforge constructs reproducible **initial geometries** for amorphous polymers.
It owns molecular inputs, Packmol execution, explicit topology edits, geometric
checks, and saved construction states. MatKit owns potential selection,
force-field assignment, MLIP evaluation, LAMMPS setup, relaxation, and MD.
A calling workflow can alternate construction and external relaxation.
Construction does not establish equilibration, permanent porosity, or a reaction
mechanism. SMARTS detection is an optional annotation aid; no named reaction
recipes are built in.

## Installation and engines

```bash
pip install 'mofforge[pop]'
pip install 'packmol==21.2.3'  # pinned native integration environment
mofforge pop-doctor --as-json
```

Packmol must support orthorhombic PBC (version 20.15.0 or newer). Resolve its
executable through `PopBuilder(packmol_bin=...)`, `[backends.pop].packmol_bin`
in `mofforge.toml`, `MOFFORGE_PACKMOL_BIN`, or PATH. Explicitly configured invalid
paths fail without falling back to another binary. The current Python
installation's `bin/packmol` is also searched. Packmol input is supplied through
a seekable stdin file because some builds rewind it.

RDKit is needed for molecular preparation and connection valence checks.
Native box loading and inspection require neither RDKit nor external engines.
Mofforge does not run LAMMPS or load MLIP checkpoints. Those settings belong to
MatKit. The previous pysimm backend has been removed.

## Inputs, labels, and density

Each registered component requires a positive integer count. Components may be
SMILES, a single-molecule MOL/SDF file, a finite connected `Crystal` with explicit
bond orders, or XYZ accompanied by an explicitly mapped molecular graph.
Coordinate inputs must include their hydrogens as atoms with coordinates.
Unresolved implicit hydrogens are rejected. SMILES preparation expands hydrogens
and performs seeded, bounded conformer preparation; supplied coordinates are
preserved. A wrapped finite `Crystal` is unwrapped using its bond images.

Exactly one sizing input is required:

- `box_lengths=[Lx, Ly, Lz]`, in angstroms, or
- `initial_packing_density`, in g/cm³, from which a cubic box is calculated using
  exact component counts and the masses of represented atoms.

There is no universal density default. Initial packing density remains in the
history; current mass, volume, and density change when atoms are deleted or an
external geometry update changes the cell. Neither is automatically an
experimental or equilibrated density. Optional `experimental_density` metadata
requires `value`, `unit`, and `measurement` fields. `material_state` records a
caller declaration (default `unspecified`); solvent-containing models must
include solvent as counted components.

Initial concentration, connectivity, and external preparation can influence
porosity. Compare independent seeds, initial densities, and cell sizes in the
external workflow; do not infer permanent porosity from a loosely packed box.
Only fully periodic orthorhombic construction boxes are supported.

## Packing with Python

```python
from mofforge.polymerize import PopBuilder

builder = PopBuilder()
builder.add_monomer("O", count=3, name="water")
builder.add_monomer("CCO", count=2, name="ethanol")
result = builder.pack(
    output_dir="packed",
    box_lengths=[20, 20, 20],  # illustrative geometry, not a material-density recipe
    minimum_distance=2.0,
    random_seed=42,
    material_state="solvent-containing",
)
print(result.to_dict())
```

Packing preserves intramolecular topology and independently checks periodic
separations, including contacts with a molecule's own images. The default
separation is 2.0 Å, with at most 0.01 Å numerical violation; the default timeout
is 300 seconds. Output geometry must preserve atom counts, element ordering, and
molecular geometry. A fresh run directory retains generated template files,
input, version probe, stdout, stderr, and a result report even on engine failure.
Display names never determine writable filenames.

`builder.prepare(random_seed=42)` returns templates without running Packmol.
Each template exposes `atoms`, `coordinates`, `bonds`, `molar_mass`, and
`to_dict()`. Unmapped input atoms have labels `a0`, `a1`, etc.; mapped atoms use
`map:NUMBER`. Hydrogens added from SMILES use `h:PARENT_LABEL:ORDINAL`.
Inspect these labels before defining connectors. Set them using
`builder.set_connectors(monomer, connectors)` or supply known labels to
`add_monomer(..., connectors=...)`.

An XYZ graph uses the coordinate atom order:

```json
{
  "atoms": [{"label": "h0", "species": "H"}, {"label": "h1", "species": "H"}],
  "bonds": [{"atom1": "h0", "atom2": "h1", "order": 1}]
}
```

Pass this as the component's `graph`. Atom records may include `formal_charge`
and `isotope`. XYZ tags alone do not supply topology.

## Explicit connections

A connector declares an attachment atom, orientation anchor, role, unique site
label, and replaceable atoms. The vector **attachment → anchor** is the outgoing
bond direction; the anchor must be a uniquely bonded neighbor. A replaceable
hydrogen is useful for the aromatic construction fixture below.

```python
from mofforge.polymerize import ConnectionRule, Connector, PopBuilder

builder = PopBuilder()
builder.add_monomer(
    "c1ccccc1", count=3,
    connectors=[
        Connector("left", "a0", "h:a0:1", "aryl", ("h:a0:1",)),
        Connector("right", "a3", "h:a3:1", "aryl", ("h:a3:1",)),
    ],
)
rule = ConnectionRule(
    name="aryl_construction_fixture", roles=("aryl", "aryl"),
    bond_order=1, bond_length=(1.45, 1.60),
    delete_atoms={"a": ["$replaceable"], "b": ["$replaceable"]},
)
result = builder.build(
    output_dir="connected", rules=[rule], box_lengths=[25, 25, 25],
    random_seed=42, target_conversion=2 / 3,
    candidate_attempt_budget=30, min_nonbonded_distance=1.5,
)
print(result.status, result.metadata)
```

This fixture constructs C18H14 with two new bonds; it is not a claimed chemical
synthesis mechanism. `build()` requires rules. Use `pack()` for packing alone,
or `connect(state, rules, ...)` to operate on an existing state.

Rules explicitly declare both deletion lists, even when empty. `$replaceable`
expands only the selected connector's declared atoms. Literal labels refer to
the selected molecular instance; `$atom` and `$anchor` reference its connector.
Optional `bond_changes` records contain `side` (`a` or `b`), `atom1`, `atom2`,
and `order`. Optional `charge_changes` records contain `side`, `atom`, and integer
`charge`. These edits must reference retained atoms and existing internal bonds.
No implicit capping, bond-order changes, charge assignment, or deletions occur.

One candidate attempt tests a compatible site pair and rule with up to twelve
seeded torsional orientations. Disconnected finite components may move rigidly;
percolating components cannot be moved as finite molecules. The rule's
`angle_tolerance` defaults to 30 degrees. Trials validate retained atoms after
edits, including local valence, bond lengths, severe overlaps, and the explicit
construction nonbonded separation. Bonded pairs and pairs separated by two bonds
are exempt from that nonbonded cutoff in their specific periodic images.

Each event consumes two original connectors. Conversion is cumulative consumed
connectors divided by the initial connector count, including after reloads.
Exhaustion or the attempt budget returns a valid `partial` state with
`success=False`. An invalid final geometry returns `failed`. The input state is
not mutated. Event history retains rejected trials and committed atom, bond,
and charge changes.

Same-component closures require `allow_cycles=True` and `minimum_cycle_size`
on every rule. Endpoints must already satisfy geometry limits. Finite-ring size
is evaluated with atom images; periodic winding is not a finite ring. Reported
periodic connectivity ranks range from zero (finite) to three.

## Saved boxes and external geometry

A native bundle contains `state.json`, `box.cif`, `box.xyz`, and `manifest.json`.
The manifest is published last and verifies every required file plus a versioned
state hash covering geometry, topology, sites, and metadata. CIF/XYZ alone do not
preserve topology. Partial construction states can be saved and reloaded.

```python
from mofforge.polymerize import ConstructionState, connect

state = ConstructionState.load(result.state_path)
# External code runs MatKit and maps its final coordinates back to these IDs.
records = [
    {"id": atom.id, "species": atom.species, "coordinates": xyz.tolist()}
    for atom, xyz in zip(state.atoms, state.coordinates, strict=True)
]
updated = state.update_geometry(
    parent_hash=state.state_hash,
    atoms=records,
    unwrapped=True,
    # box_lengths=[new_Lx, new_Ly, new_Lz],  # optional orthorhombic cell update
)
path = updated.save("updated")
```

Unwrapped positions must use the parent's per-atom image reference. Alternatively,
provide each atom's integer `wrapping_offset`; its unwrapped position is the
supplied wrapped coordinate plus that offset times the new box lengths, in the
same parent reference. Do not guess offsets from nearest distances. Atom order
may change in the update records, but IDs and species must match exactly once.

Geometry updates invalidate earlier geometric validation and derived results.
Subsequent connections validate retained geometry before committing. Updates
preserve declared topology and do not detect or import reactive MD chemistry.
A `state.crystal` is a detached representation; editing its graph does not edit
the authoritative native state. Legacy MOF graph-edit behavior is unchanged.

## CLI and MCP

Python `run_config`, CLI, and MCP use the same JSON shape:

```json
{
  "components": [{"source": "CC", "count": 3}],
  "packing": {"box_lengths": [20, 20, 20], "random_seed": 42}
}
```

```bash
mofforge pack --config packing.json --output packed --as-json
mofforge polymerize --config construction.json --output connected --as-json
```

For construction, add `connectors` to component definitions and a `connection`
object containing `rules`, `target_conversion`, `candidate_attempt_budget`, and
`min_nonbonded_distance`. For resumption, provide `state_path` and `connection`,
omitting `components` and `packing`. Relative paths are resolved from the calling
working directory. See the [executable example](../examples/polymerize/README.md).

MCP provides `mofforge_pack(config, output_dir)` and
`mofforge_polymerize(config, output_dir)`. The `pop` capability checks RDKit;
Packmol is required only when packing executes. `mofforge_list_reactions` lists
annotation patterns and reports an empty recipe list. CLI exits zero only on
completion. Partial and failed results retain status, diagnostics, and saved
state references where available across all interfaces.

Migration: `target_density` is replaced by `initial_packing_density`, scalar
`box_length` by `box_lengths`, and `n_monomers` by explicit component counts.
`forcefield`, `equilibrate`, `md_settings`, `lammps_bin`, and
`MOFFORGE_LAMMPS_BIN` are removed and produce migration errors when explicitly
used. Configure simulations in MatKit.
