# Interpret structure validation

```bash
python examples/validation/structure_validation.py
python examples/validation/structure_validation.py examples/data/crystals/UiO-66.cif
```

[structure_validation.py](structure_validation.py) requires core dependencies.
With no positional CIF, it inspects all bundled crystals. With a path, it
diagnoses that structure. Every run also checks deliberately overlapping atoms
and an empty/unchecked structure as negative controls.

## Read the result in order

1. Check `checks_performed` and `checks_skipped`. An unperformed check is not a
   passing check. Bond inference is used for ordered inputs; partial occupancy
   is diagnosed without forcing an ordered graph.
2. Read mandatory geometry/overlap errors and `steric_clashes`. These can make
   `is_valid` false.
3. Review `unusual_bonds`, `coordination_issues`, close contacts, and warnings,
   even if mandatory geometry checks passed.
4. Inspect the charge check. The example requests it, but ordinary CIFs often
   lack usable oxidation-state information; unknown charge is not neutrality.

The summary no longer announces that a structure “looks good” while suppressing
coordination caveats. For example, the UiO-66 fixture can pass mandatory geometry
checks while its inferred Zr coordination lies outside the generic range. The
SIFSIX disorder fixture contains real geometric overlaps detected by this
validator. Both results remain visible in `summary.json`.

The two controls must have `is_valid=false`. The script exits zero when it has
successfully diagnosed inputs and identified those controls; diagnostic findings
are expected lesson output. Invalid inputs that cannot be processed still fail.

## Adjust checks deliberately

Defaults use a **0.5 Å clash tolerance** and **0.3 relative bond-length tolerance**.
The current validator distinguishes severe overlaps from close contacts and
accounts for bonds/molecular membership where available. Do not interpret a
single tolerance as a universal chemistry criterion. Skipping a check should
be recorded and explained, rather than treated as evidence of validity.

**Try:** compare a raw structure with a modified one, keeping the check settings
fixed. Then inspect what changes when bonds or oxidation states are unavailable.

Equivalent CLI:

```bash
mofforge validate examples/data/crystals/IRMOF-1.cif --json
```

For batch outputs, run [artifact verification](../pipeline/README.md) separately.
A reproducibly generated or hash-verified artifact may still contain invalid
geometry. See the [reliability notes](../../docs/reliability.md) for the current
validation contract.
