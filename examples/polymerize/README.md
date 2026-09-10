# Polymer box construction

Mofforge prepares explicit molecular topology, invokes Packmol, and constructs
connections using complete graph-edit rules. Simulation runs externally in
MatKit. The [construction guide](../../docs/polymerize.md) describes the API,
JSON configuration, density semantics, and geometry handoff.

```bash
pip install 'mofforge[pop]'
python examples/polymerize/amorphous_pop.py --output-dir pop_prepared
```

Default mode prepares three counted benzene templates with explicit para
connectors and writes their atom-label mapping plus `packing.json` and
`construction.json`. It needs only RDKit. This is a generic geometry fixture,
not a claimed chemical polymerization mechanism or material-density recipe.
The example uses a 25 Å cubic box; the library has no default volume or density.
You may supply `--initial-packing-density` in g/cm³ instead.

To execute the native stages:

```bash
pip install 'packmol==21.2.3'
mofforge pop-doctor --as-json
python examples/polymerize/amorphous_pop.py --mode pack --output-dir pop_packed
python examples/polymerize/amorphous_pop.py --mode connect --output-dir pop_connected
```

The connection example deliberately stops after one candidate attempt, reloads
its partial native state, imports an explicitly mapped translation, and resumes
to produce a finite C18H14 chain with two new bonds. The geometry update exercises
the handoff contract without claiming to run MD. Each saved bundle includes
native state JSON, P1 CIF, XYZ, and a verification manifest; Packmol attempts
also retain their input, version, stdout, stderr, and result report.

The generated configurations are shared by Python, CLI, and MCP:

```bash
mofforge pack --config pop_prepared/packing.json --output pop_cli_packed --as-json
mofforge polymerize --config pop_prepared/construction.json --output pop_cli_connected --as-json
```

The calling workflow owns physical preparation, comparisons across seeds and
initial densities, and subsequent adsorption or MD calculations in MatKit.
