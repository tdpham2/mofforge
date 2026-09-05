# Fragment modification

These lessons use core dependencies and bundled XYZ fragments. Read the search
guide first: location and orientation indices refer to a specific match result.

The `!` marker identifies atoms removed during replacement. Untagged query
atoms form the alignment scaffold; they must be identifiable in the replacement.
Each lesson writes the child CIF, its manifest, and a summary containing a
separate validation report.

## 1. Install acetylamido groups

```bash
python examples/modify/linker_functionalization.py --nb-loc 6 --seed 42
```

[linker_functionalization.py](linker_functionalization.py) loads IRMOF-1, matches
`2-!-p-phenylene.xyz`, and installs
`2-acetylamido-p-phenylene.xyz` at six selected locations.

1. Inspect the tagged H in the query.
2. Find candidate rings in the parent graph.
3. Select locations with the local seeded RNG and align each replacement.
4. Write the child and inspect its validation.

The default changes **424 → 466 atoms**: each selected substitution adds seven
net atoms. `locations` in the summary records the actual selection. `--nb-loc 0`
means all available locations; a positive number is capped by available matches.
Use `--output` to override the primary CIF filename or `--output-dir` to collect
the artifacts together.

**Try:** compare seeds 42 and 7 at the same coverage. A count match alone does
not prove identical geometry; inspect selected indices and coordinates.
No matches usually means the input hydrogens or bonding graph differ from the
query.

## 2. Compare the five replacement modes

```bash
python examples/modify/selective_modification.py --seed 42
```

[selective_modification.py](selective_modification.py) applies every mode to the
same original parent, producing five CIFs:

| Mode | Selection | Default child atoms |
| --- | --- | ---: |
| all_optimal | Every location; lowest alignment error | 472 |
| random_locations | Eight seeded locations; optimal orientations | 440 |
| specific_locations | Locations 0, 5, 10, 15 | 432 |
| specific_orientations | Locations 0–3 with mappings 0–3 | 432 |
| random_orientations | Every location; seeded mapping choices | 472 |

Both explicit location and orientation indices are **zero-based**.
`ori=[0]` selects the first mapping; **omit `ori` to optimize**. The summaries
retain chosen locations, orientations, and validation for comparison.
Optimization here is rigid-fragment alignment, not a force-field relaxation.

**Try:** choose two explicit locations, inspect their possible mappings, and
compare the resulting steric contacts. Out-of-range indices should fail rather
than silently select another location.

## 3. Introduce missing-linker defects

```bash
python examples/modify/defect_engineering.py --loc 2 8 --seed 42
```

[defect_engineering.py](defect_engineering.py) uses UiO-66,
the masked `BDC.xyz` query, and `formate_caps.xyz`. Removing the linker core
while retaining/installing terminal caps changes **912 → 896 atoms** for the
two default locations. The output is `defected_UiO-66.cif`.

Review the remaining coordination and cap placement in the validation report.
The geometric construction does not determine a physically balanced defect
chemistry, charge state, or relaxation protocol.

**Try:** create a single defect with `--loc 0`, then compare cap geometry with the
two-defect model. For SMILES-based site selection and coverage campaigns, continue
to [functionalization](functionalization.md).
