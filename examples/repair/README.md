# Repair and solvent removal

These lessons require only core dependencies. The bundled structures are teaching
fixtures, so each script checks that its expected starting defect or guest exists.

## 1. Restore missing aromatic hydrogen atoms

```bash
python examples/repair/structure_repair.py --seed 42
```

[structure_repair.py](structure_repair.py) matches bare carbon rings in
IRMOF-1_noH and replaces them with complete phenylene fragments. Unlike a masked
substitution, the replacement here is a superset of the query scaffold.

Expect **328 → 424 atoms**, adding **96 H atoms**. Inspect `repaired_IRMOF-1.cif`
and the geometry report before using it further. `--fragment-path` permits a
custom query/replacement library, and `--output` selects the primary CIF.

**Try:** compare the repaired atom count and species with the bundled IRMOF-1
reference. This operation restores the particular aromatic hydrogens represented
by the fragment; it is not a general protonation or oxidation-state assignment.

## 2. Correct disorder and remove a named guest

```bash
python examples/repair/cleanup_and_repair.py --seed 42
```

[cleanup_and_repair.py](cleanup_and_repair.py) uses SIFSIX-2-Cu-i:

1. Match eight disordered ring locations using the tagged disorder fragment.
2. Replace them with ordered pyridyl rings: **200 → 136 atoms**.
3. Re-infer bonds and search for isolated acetylene molecules.
4. Delete eight guests: **136 → 104 atoms**.

`ordered_with_guests.cif` records the intermediate; `cleaned_SIFSIX.cif`
records the final structure. `disconnected_component=True` ensures the guest
query matches a whole isolated component instead of a similar framework motif.

**Try:** compare guest search with and without the disconnected-component
constraint before deleting anything. A changed fixture with no expected ring or
guest matches fails explicitly.

## 3. Classify solvent automatically

```bash
python examples/repair/solvent_removal.py --seed 42
```

[solvent_removal.py](solvent_removal.py) first places two CO2 molecules into clean
IRMOF-1 to provide a known loaded structure. Default desolvation removes two
three-atom components and restores **430 → 424 atoms**. Compare the automatically
cleaned structure with two alternatives:

- `min_atoms=3` retains every component containing at least three atoms,
  including these guests. It is a **retention** threshold.
- `n_framework_components=1` explicitly keeps the largest component;
  `keep_metal_containing=True` also protects metal-containing components.

The automatic heuristic keeps the largest component and others at least half
its size. It infers classification from connectivity and size, not solvent
identity or chemical charge. Disconnected framework units and counterions can
therefore be misclassified. In the SIFSIX example, a default cleanup can remove
two SiF6 units when inferred coordination leaves them disconnected; the
[cleanup pipeline](../pipeline/README.md) preserves those explicitly.

**Try:** replace CO2 with water and inspect removed formulas before accepting
the deletion. Use named guest removal when solvent identity is known.

Equivalent CLI for the generated loaded input:

```bash
mofforge desolvate --parent examples/_outputs/solvent_removal/loaded.cif --output examples/_outputs/solvent_removal/cli_cleaned.cif
```
