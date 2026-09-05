# Construct frameworks and building blocks

Use the `build` extra; Python 3.11 is the builder CI baseline. SMILES conversion
also needs `chem`. Run the scripts from the repository root after installation.
Each construction report separates successful output generation from geometric
validation; a parseable CIF is not necessarily an optimized structure.

## 1. Build with Pormake

```bash
python examples/build/pormake_build.py
python examples/build/pormake_build.py --describe-topology pcu
python examples/build/pormake_build.py --list-bbs
```

[pormake_build.py](pormake_build.py) defaults to the integration-tested combination
**pcu + N108 + E1**, using accuracy 3:

1. Inspect the topology and register a node and edge with `MOFBuilder`.
2. Build the structure into the selected output directory.
3. Read `success`, `errors`, `outputs`, and the separate validation result.

The locally tested catalog produces a 123-atom structure. Consume returned
output paths, since catalog names affect filenames. Node N108 has six connection
points, matching pcu. Building-block connection points describe net connectivity;
they are not the coordination number of one metal atom.

`--node` and `--edge` are repeatable and accept catalog names or compatible
files. When you explicitly supply nodes, also supply desired edges; omitting
edges in that mode retains the original node-only behavior.
`--bb-dir`, `--status`, `--list-topologies`, and `--accuracy` remain available.

**Try:** register a compatible XYZ block using its absolute path. If construction
fails, inspect errors for connectivity or geometry mismatch before changing
accuracy. Raising accuracy cannot fix a four-connected node on a six-connected net.

## 2. Build with TOBACCO

```bash
python examples/build/tobacco_build.py
python examples/build/tobacco_build.py --tobacco-data /path/to/tobacco-data --status
python examples/build/tobacco_build.py --topology pcu --node 6c_Zn_1_Ch.cif --edge 1B_4H_Ch.cif
```

[tobacco_build.py](tobacco_build.py) defaults to the same **dmc** fixture used in
the real backend tests: nodes `4c_1Zn_Ch` and `triazole`, and edges `ntn_edge`,
`oxalic_edge`, and `squOxa_ch`. Their CIFs live in the external data bundle's
`tests/fixtures` directories. Explicit node/edge arguments support other recipes.

TOBACCO's package supplies code; runtime templates and blocks are separate.
The library resolves an explicit `--tobacco-data` / legacy `--tobacco-path`,
environment/TOML configuration, nearby installed data, or the existing verified
download-and-cache mechanism. For offline execution, supply prepared data
containing `template_database`, `nodes_database`, and `edges_database`.
The default lesson also requires the fixture subdirectories. Do not substitute
an incomplete custom catalog for that fixture silently.

One topology can produce several CIFs. Inspect each manifest, including files
whose names indicate a backend bond-check failure; the top-level validation
field describes the first parseable output. The local dmc run produced nine
outputs. `--parallel` requests the backend's parallel mode.

**Try:** list available nodes, edges, and topologies before adapting a recipe.
The [construction manual](../../docs/build.md) documents the data configuration
and upstream runtime-data provenance.

## 3. Convert a SMILES linker

```bash
python examples/build/smiles_building_blocks.py --seed 42
python examples/build/smiles_building_blocks.py --mode carboxylic --output-dir examples/_outputs/bdc_stripped
```

[smiles_building_blocks.py](smiles_building_blocks.py) detects BDC's two
carboxylate connection groups, embeds a molecular geometry, and writes
`linker.cif` for TOBACCO and `linker.xyz` for Pormake. Conversion needs RDKit
but does not require either construction backend.

| Mode | Geometry convention |
| --- | --- |
| auto | Detect carboxylates first, otherwise select direct endpoints |
| carboxylate | Keep carboxylate atoms and place X connection markers at oxygen centroids |
| direct | Mark connection atoms as X and remove their hydrogens |
| carboxylic | Remove COOH groups and mark the adjoining scaffold atoms as X |

The summary's `auto_detection` describes automatic detection; `requested_mode`
records any override. Choose a convention that matches the node building block,
avoiding duplicated or missing binding groups.

**Try:** convert an amino-substituted BDC linker using `--smiles`.
An invalid SMILES or wrong number of detected groups fails explicitly.

## 4. Select a reproducible random linker

```bash
python examples/build/pcu_zn_random_linker.py --seed 42
python examples/build/pcu_zn_random_linker.py --list-nodes
python examples/build/pcu_zn_random_linker.py --node N108 --edge E1 --seed 42
```

[pcu_zn_random_linker.py](pcu_zn_random_linker.py) reads the topology's required
connectivity, filters Zn-containing nodes accordingly, and samples from a sorted
catalog of two-connected edges with a local RNG. For pcu it selects six-connected
nodes, preferring N108 when available. It does not label every Zn node a
paddlewheel.

`--retries` limits distinct attempted edges. Every attempt retains its error,
output paths, and validation. A successful build stops the loop; exhausting all
candidates exits nonzero. The tested catalog selected E44 at seed 42, but catalog
changes can alter the selected name. The seed controls selection, not all
numerical operations inside the backend.

**Try:** explicitly request N109 on pcu and inspect the connectivity error.
For an end-to-end application, continue to [build and modify](../workflows/README.md).
