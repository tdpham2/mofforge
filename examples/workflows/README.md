# Complete scientific workflows

These two lessons combine the focused tutorials into reviewable applications.
They run independently and write intermediate artifacts, summaries, and final
structures.

## 1. Build, modify, validate, and render

```bash
python examples/workflows/build_and_modify.py --seed 42
python examples/workflows/build_and_modify.py --seed 42 --render
```

[build_and_modify.py](build_and_modify.py) requires Pormake. The optional
`--render` step additionally requires Playwright, Chromium, and CDN access.

1. Build pcu from N108 and E1, using the tested construction recipe.
2. Infer the resulting framework graph and find a masked phenylene motif.
3. Replace one selected C-H site with a nitro group.
4. Write `modified.cif`, validate it, and record modification provenance.
5. Optionally render `modified.png` with the unit cell visible.

The tested catalog produces **123 parent atoms, nine matching locations, and
125 atoms after modification**. The selected fragment is a phenylene motif;
the complete E1 linker is not assumed to be BDC. The construction CIFs remain
in the `build` subdirectory.

A build failure or missing motif stops the workflow with an error instead of
silently skipping the modification. Geometric validation is retained even when
it finds issues; inspect it before choosing a relaxation calculation. The seed
controls fragment selection and does not make every external builder operation
deterministic.

**Try:** increase the number of modified locations in a copy of the script and
compare contacts. To change the linker, first confirm the new building block's
connection convention and query motif.

## 2. Screen, resolve, and place CO2

```bash
python examples/workflows/screen_and_place.py --seed 42
```

[screen_and_place.py](screen_and_place.py) requires only core dependencies in
demo mode:

1. Stage fictional CoRE metadata and its explicit geometry mappings.
2. Screen Zn records for PLD ≥ 3.8 Å, water-stability score ≥ 0.7, and the
   all-solvent-removed extension.
3. Resolve each candidate's local CIF.
4. Place up to two CO2 molecules at sampled void sites.
5. Save per-candidate outputs and validation, retaining missing/failure outcomes.

The demo selects **two records**: one produces a loaded IRMOF-1 geometry and one
intentionally has no CIF. Expect one `placed` and one `missing_structure`
entry. Fictional screening values demonstrate software behavior and are not
performance data for IRMOF-1.

Real mode uses your metadata and CIF collection together:

```bash
python examples/workflows/screen_and_place.py --data-path /data/coremof/metadata.csv --structures-dir /data/coremof/structures --metal Cu --limit 3 --seed 42
```

The script retains partial outcomes and fails if no candidate produces a loaded
structure. Read requested versus actual placement counts and validation for
every result. Empty screens, missing CIFs, and placement errors require different
responses; none should be labeled a successful adsorption simulation.

**Try:** narrow the pore-size filter and compare candidate counts before running
a larger campaign. The [MCP guide](../mcp/README.md) describes how the same broad
workflow can be submitted through ChemGraph's configured HPC fan-out interface.
