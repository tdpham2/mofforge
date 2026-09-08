# Generate amorphous porous organic polymers (POPs)

Unlike MOFs and crystalline COFs, porous organic polymers (CMPs, PIMs, HCPs) are
*amorphous* covalent networks with no periodic net. mofforge builds them by
**simulated polymerization**: pack monomers into a box (Packmol), form covalent
bonds between nearby reactive sites and relax under MD (Polymatic + LAMMPS), then
equilibrate to a realistic density. The heavy lifting is delegated to
[pysimm](https://pysimm.org); mofforge prepares the monomers, orchestrates the
run, and converts the result into a validated P1 CIF.

Use the `pop` extra (`pip install 'mofforge[pop]'`) plus the **Packmol** and
**LAMMPS** binaries. mofforge resolves those binaries at run time from
`MOFFORGE_PACKMOL_BIN` / `MOFFORGE_LAMMPS_BIN`, a `mofforge.toml`
`[backends.pop]` section, or a `PATH` search.

## 1. Prepare a build (rdkit only)

```bash
python examples/polymerize/amorphous_pop.py
python examples/polymerize/amorphous_pop.py -m "NCCN" -m "O=Cc1ccc(C=O)cc1"
```

[amorphous_pop.py](amorphous_pop.py) runs the parts that need only rdkit:

1. Detect reactive sites on each monomer SMILES (a diamine and a dialdehyde by
   default) and confirm they form a compatible pair — here an **imine (Schiff
   base)** condensation.
2. Generate 3-D monomer geometry (RDKit embed + UFF) and compute each monomer's
   molar mass.
3. Size the cubic packing box for the requested `--target-density`.
4. Report tool availability via `doctor()` so you can see what is still needed
   for the full run.

The reactive-site menu and compatible reactions come from
`mofforge.polymerize.reactions` — the POP analogue of the curated functional
group menu used for linker functionalization.

## 2. Run the full polymerization (Packmol + LAMMPS)

Once the binaries are installed:

```bash
mofforge polymerize -m "NCCN" -m "O=Cc1ccc(C=O)cc1" \
    --target-density 0.8 --forcefield gaff2 -o pop_out --random-seed 42
mofforge pop-doctor          # report pysimm / Packmol / LAMMPS availability
```

or through the MCP server with the `mofforge_polymerize` and
`mofforge_list_reactions` tools (capability `pop`). Both write a P1 CIF and a
provenance manifest, and attach a validation report (steric clashes, bond
lengths) exactly like the MOF builders.

**Try:** swap in a trifunctional monomer (e.g. a trialdehyde,
`O=Cc1cc(C=O)cc(C=O)c1`) to form a cross-linked network instead of a linear
chain, and lower `--target-density` to open up the pore structure.
