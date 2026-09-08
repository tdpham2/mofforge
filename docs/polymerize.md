# Amorphous Porous Organic Polymer (POP) Generation

[Polymerization lessons](../examples/polymerize/README.md) provide a runnable
walkthrough with bundled inputs and expected results.

mofforge generates **amorphous** porous organic polymers — CMPs, PIMs, HCPs —
which, unlike MOFs and crystalline COFs, have no periodic net. They are built by
**simulated polymerization**: pack monomers into a box, form covalent bonds
between nearby reactive sites under MD relaxation, then equilibrate to a
realistic density.

mofforge does not reimplement packing or the bond-formation loop. It orchestrates
[pysimm](https://pysimm.org), which drives **Packmol** (packing) and **LAMMPS**
(MD + Polymatic simulated polymerization). mofforge owns the thin layer around
that engine: monomer preparation, box sizing, conversion of the result into a
validated `Crystal`, provenance, and the MCP/CLI surface.

## Requirements

```bash
pip install 'mofforge[pop]'   # pysimm + rdkit
```

Plus the **Packmol** and **LAMMPS** binaries. mofforge resolves them at run time,
in order:

1. an explicit value (`packmol_bin` / `lammps_bin` kwargs, or `mofforge.toml`
   `[backends.pop]`),
2. the `MOFFORGE_PACKMOL_BIN` / `MOFFORGE_LAMMPS_BIN` environment variables,
3. a `PATH` search (`packmol`; `lmp`, `lmp_serial`, `lmp_mpi`, `lammps`).

A missing binary raises an actionable error only when a step needs it — the same
idiom used for the TOBACCO data directory. Check availability with:

```bash
mofforge pop-doctor
```

## Reactive sites and reactions

A monomer bonds through **reactive sites** — the POP analogue of a MOF
connection point / `!` anchor. mofforge ships a curated menu (in
`mofforge.polymerize.reactions`) of reactive groups and the reactions they form,
mirroring the curated functional-group menu used for linker functionalization:

| reactive site | example reaction |
|---|---|
| amine + aldehyde | imine (Schiff base) condensation |
| aryl halide + aryl halide | aryl-aryl coupling (C–C) |
| boronic acid + boronic acid | boroxine condensation |
| boronic acid + hydroxyl | boronate ester condensation |
| vinyl + vinyl | vinyl addition (C–C) |

Sites are detected from the monomer SMILES by SMARTS matching; supply them
explicitly (`ReactiveSite`) or via `!` anchor tags for file-based monomers.

## Python API

```python
from mofforge.polymerize import PopBuilder

builder = PopBuilder()
builder.add_monomer("NCCN", name="diamine", functionality=2)
builder.add_monomer("O=Cc1ccc(C=O)cc1", name="dialdehyde", functionality=2)

result = builder.build(
    output_dir="pop_out",
    target_density=0.8,     # g/cm^3, sizes the packing box
    forcefield="gaff2",     # gaff2 | dreiding | pcff
    target_conversion=0.95,
    equilibrate=True,
    random_seed=42,
)

if result.success:
    print(result.output_paths[0])          # P1 CIF
    print(result.metadata)                 # conversion, density, n_bonds, ...
    print(result.validation.summary())     # steric clashes, bond lengths
```

`POPResult` is field-for-field the same as the MOF builders' `BuildResult`, so a
POP is validated, manifested, and consumed like any other built structure.

## CLI

```bash
mofforge polymerize -m "NCCN" -m "O=Cc1ccc(C=O)cc1" \
    --target-density 0.8 --forcefield gaff2 -o pop_out --random-seed 42
mofforge pop-doctor
```

## MCP tools

With the `pop` capability (pysimm installed):

- `mofforge_polymerize` — generate a POP from monomer SMILES.
- `mofforge_list_reactions` — list curated reactive site types and reactions.

Both are hidden when pysimm is not installed, and report a clean error (rather
than raising) when Packmol or LAMMPS is missing.

## What happens under the hood

```
monomer SMILES ─► RDKit 3-D + reactive sites ─► pysimm System + force field
   ─► polymatic.pack (Packmol)  ─► polymatic.polymatic (Polymatic + LAMMPS)
   ─► System → Crystal (P1 box) ─► validate ─► CIF + provenance manifest
```
