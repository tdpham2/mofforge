# Pipelines, batches, provenance, and verification

These five core-only lessons connect individual operations into reproducible
workflows. Each runs independently with bundled inputs and writes its own
artifacts. Use seed 42 initially so results can be compared with the references.

## 1. Chain modifications and inspect intermediates

```bash
python examples/pipeline/multi_step_pipeline.py --seed 42
```

[multi_step_pipeline.py](multi_step_pipeline.py) queues two acetylamido
substitutions, two nitro substitutions, and validation on IRMOF-1. The root seed
is passed into `Pipeline`, which derives operation seeds.

`build_all` executes the recipe and returns intermediate crystals; the lesson
writes them to numbered CIFs and saves final `provenance.json`. Expect the first
replacement to produce 438 atoms and the second 442. Validation may add an
unchanged-geometry intermediate, so inspect operation history as well as counts.
`pipeline.validation_reports` exposes completed validation reports.

**Try:** place a validation step between substitutions. Repeated calls execute
the recipe from its parent; they are not checkpoint resumption.

## 2. Build a conservative cleanup pipeline

```bash
python examples/pipeline/cleanup_pipeline.py --seed 42
```

[cleanup_pipeline.py](cleanup_pipeline.py) chains disorder repair, named acetylene
removal, automatic desolvation, and validation for SIFSIX-2-Cu-i. The targeted
operations reduce 200 atoms to 104, as in the standalone repair lesson.

The desolvation step sets `min_atoms=7` and `keep_metal_containing=True`.
This intentionally preserves seven-atom SiF6 components that the inferred graph
can leave disconnected. Without that retention rule, automatic cleanup removes
14 additional atoms. The rule is specific to this teaching structure; it is
not a universal solvent cutoff.

**Try:** inspect the standalone solvent-removal metadata before choosing retention
rules for a different structure. Review the final 104-atom `cleaned.cif` and
the recorded operation sequence.

## 3. Execute a YAML batch

```bash
python examples/pipeline/batch_processing.py
python examples/pipeline/batch_processing.py --compare-parallel
mofforge batch --config examples/pipeline/batch.yaml
```

[batch_processing.py](batch_processing.py) copies the same parent into two
different input directories with the same filename. It loads the committed
[batch.yaml](batch.yaml), resolves all input/output paths, selects two nitro
substitutions with `mode: nb_loc_2`, and calls `run_batch`.

The resulting filenames include run identities, preventing same-name input
collisions. Read `output_path` and `success` from each result rather than
constructing output names. Execution success and final validation remain
separate. A failed or empty batch is an error for this lesson.

`--parallel 2` requests workers directly. `--compare-parallel` runs a serial
reference and a two-worker batch using the same inputs, then checks run IDs,
species, and fractional coordinates. The script has a main guard for platforms
that spawn workers.

The committed YAML command assumes the repository root because the current
batch API resolves relative paths against the working directory. The Python
lesson emits absolute paths and can run from anywhere. The supported operation
types are `replace`, `remove`, `desolvate`, and `validate`.

**Try:** change `mode` to `all_optimal` in a copied recipe. Use the YAML's actual
`mode` field; arbitrary replacement kwargs are not automatically forwarded.
Explicit per-operation seeds also avoid path-dependent seed changes for those
replacement steps.

## 4. Inspect and replay provenance

```bash
python examples/pipeline/provenance_replay.py --seed 42
```

[provenance_replay.py](provenance_replay.py) replaces three locations, saves the
operation record, reloads it with `Provenance.from_json`, and uses the recorded
seed to repeat the replacement from the original parent and fragments.

Expect identical species and coordinates within the same environment. The files
`original.cif` and `replayed.cif` need not have identical timestamps or manifest
bytes. `provenance.json` records selected locations, history, parameters, and
software versions; it is not itself a structure artifact manifest.

**Try:** change the replacement fragment and observe why a seed alone is not a
complete recipe. CIF reload can round coordinates and lose metadata; retain the
original inputs instead of treating exported CIF as a lossless checkpoint.

## 5. Verify outputs independently

```bash
python examples/pipeline/artifact_verification.py
```

[artifact_verification.py](artifact_verification.py) contrasts three cases:

| Case | File integrity | Complete verification |
| --- | --- | --- |
| Current batch output and its complete descriptor | Pass | Pass |
| Plain CIF export with no batch input descriptor | Pass | Incomplete |
| Deliberately modified copy checked against the original manifest | Fail | Fail |

The script changes only an output copy. Read the actual batch artifact path in
`summary.json`, then use:

```bash
mofforge verify /path/from/summary/result.cif --json
```

Full verification needs the recorded input files at their recorded locations
and compatible settings and software. For checking a proposed reuse, the API
also accepts an independently prepared `expected_inputs` descriptor.
Verification does not establish geometric validity or simulation convergence.
The CLI returns nonzero for both invalid and incomplete verification.

**Try:** preserve a batch output, then edit a copy of its source input and
compare reports. See [artifact verification](../../docs/artifact_verification.md)
for the precise compatibility limits.
