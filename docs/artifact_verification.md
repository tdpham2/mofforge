# Artifact verification

[Executable verification lesson](../examples/pipeline/README.md) provides runnable workflows with bundled inputs and expected results.

This unreleased follow-up adds complete batch input records and a read-only
artifact verifier. It does not yet provide checkpoint resumption, lossless
reload, or atomic publication of a structure/manifest pair.

## Using the verifier

```bash
mofforge verify results/example.cif
mofforge verify results/example.cif --json
mofforge verify copied.cif --manifest original.cif.json --json
```

The default manifest location is `<artifact filename>.json`. The artifact path
you supply is the file checked; the original `output_path` recorded in the
manifest is informational. Verification never rewrites either file.

The CLI exits with 0 only when `is_verified` is true, and 1 for invalid,
unavailable, or incomplete verification. JSON output includes:

| Field | Meaning |
| --- | --- |
| `file_integrity` | The nonempty artifact's exact bytes match its recorded SHA-256 |
| `inputs_verified` | The descriptor, source files, scientific settings, and environment match |
| `is_verified` | Both checks passed and there are no errors |
| `identity_sha256` | The verified input descriptor's content identity, when available |
| `errors`, `warnings` | Reasons for failures or incomplete verification |

These fields describe artifact integrity and input compatibility. They do not
assert geometric validity, optimization convergence, or physical performance.
An invalid structure can still be an intact, reproducible workflow output.
Inspect its separate validation report before using it scientifically.

Python callers can check a specific proposed reuse with an expected descriptor:

```python
from mofforge import verify_artifact

report = verify_artifact("result.cif", expected_inputs=expected_descriptor)
if report.is_verified:
    print(report.identity_sha256)
else:
    print(report.errors, report.warnings)
```

`expected_descriptor` is the independently prepared descriptor for the requested
workflow. Without it, the verifier checks the recorded recipe and its current
inputs; it cannot infer a different workflow that the caller intends to run.

## Batch input records

New batch manifests contain an additive `input_descriptor` with its own schema
version, currently 1. The artifact and provenance schemas remain version 2.
The descriptor includes:

- Content hashes and resolved locations for the parent, query, replacement,
  and guest files actually referenced by the operations.
- Effective operation defaults and overrides, root and per-operation seeds,
  seed strategy, output format, and final validation settings.
- Bond padding, anchor tag, covalent and van der Waals radii, coordination
  ranges, and validation defaults captured before worker submission.
- Installed scientific package versions, platform, Python implementation,
  and a package-source fingerprint that detects changes in editable installs.

Workers execute with the captured scientific settings and resolved fragment
directory. Source files are checked again before outputs are written; changes
detected during execution fail that job. This check does not lock input files
or create an immutable snapshot of them.

File locations are excluded from the descriptor hash. Batch filenames still
include a source-instance identity so identical inputs from different locations
cannot overwrite one another. The existing `batch-v1` seed derivation is retained
for compatibility: derived seeds can depend on original paths, so relocating a
workflow can still change its effective seeds and therefore its input identity.
Explicit per-replacement seeds avoid that particular source of variation.

## Compatibility and limits

Batch output names change because their run identity now includes all consumed
file contents, effective settings, and the execution environment. Consume the
returned output paths instead of reconstructing filenames. Editing even an XYZ
comment changes its file identity; this is intentional for provenance integrity.

Existing schema-2 manifests without input descriptors remain readable. Their
file hash can pass, but the report warns that only file integrity was checked,
and `is_verified` remains false. Standalone provenance JSON is not an artifact
manifest. Other workflows that have not adopted complete input descriptors
likewise receive integrity-only verification. No missing inputs are invented and
no migration is performed automatically.

Full verification requires the recorded source files at their recorded locations
and matching settings and software in the current process. Missing inputs,
modified source code, changed package versions, and incompatible requested
parameters are explicit failures. File integrity remains separately inspectable.
Version strings and source hashes support reproducibility checks, not a promise
of bitwise identity across hardware or an authentication signature.

`structure_sha256` hashes the original in-memory serialization. The verifier
does not compare it against a CIF reparse: interchange formats can round
coordinates or omit metadata. Exact byte integrity, lossless reload, and
tolerance-based structural equivalence are distinct checks.

See the [roadmap](roadmap.md) for verified reload, durable artifact completion,
environment preflight, and resumable campaigns.
