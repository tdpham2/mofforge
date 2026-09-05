# SMILES-driven linker functionalization

Install the `chem` extra. These two lessons use BDC SMILES
`O=C(O)c1ccc(C(=O)O)cc1` with the bundled IRMOF-1 CIF. The linker must describe
the structure being modified, including its explicit aromatic hydrogens.

## 1. Enumerate sites and choose a group

```bash
python examples/modify/agent_functionalization.py --group NH2 --site 0 --coverage 0.5 --seed 42
```

[agent_functionalization.py](agent_functionalization.py) demonstrates the same
discrete choices an agent can make without writing coordinates:

1. Enumerate aromatic C-H sites and their symmetry classes.
2. Read the curated functional-group menu.
3. Choose a site index, group, and fraction of matched rings.
4. Generate query/replacement fragments, install them, and inspect the result.

BDC provides four aromatic C-H sites in one symmetry class. Site indices describe
positions within the linker, while coverage selects matched framework locations.
At coverage 0.5, expect **12 of 24 matches** to be functionalized. Fractional
coverage uses rounded counts with at least one location for a positive fraction;
zero coverage is a no-change control.

`summary.json` contains site descriptions, available groups, match counts,
actual substitutions, validation, and any errors. `--output` overrides the CIF
path. The old `--campaign` option remains available and runs a small sweep.

**Try:** compare NH2, F, and NO2 at one site with the same seed. If the result
reports no query matches, first check linker identity, hydrogens, and bonding.
A missing optional dependency exits with an installation instruction.

## 2. Run a small campaign

```bash
python examples/modify/functionalization_campaign.py
python examples/modify/functionalization_campaign.py --groups NH2 F NO2 --coverages 0.25 0.5 1.0
```

[functionalization_campaign.py](functionalization_campaign.py) defaults to
**two groups × two coverages (0 and 0.5)**. The zero-coverage members are controls.
Each member receives a derived seed and its own identity-bearing CIF filename.

Read `results` in the returned order: successful valid structures are preferred,
then lower clash counts. This ranking does not calculate adsorption performance
or chemical stability. Retain error members for diagnosis; the script exits
nonzero if any combination fails.

The filenames encode run identity, so consume `output_cif` from results rather
than reconstructing names. Equal seeds reproduce choices within the same
scientific environment, not necessarily across versions of embedding libraries.

**Try:** add a coverage of 1.0 and compare clashes to the half-coverage structures.
The [provenance lesson](../pipeline/README.md) explains recorded recipes and
the distinction between seeded replay and durable workflow resumption.
