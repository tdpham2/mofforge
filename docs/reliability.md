# Research reliability in 0.2

[Validation and reliability lessons](../examples/validation/README.md) provides runnable workflows with bundled inputs and expected results.

Version 0.2 corrects scientific behavior while retaining the main Python entry
points. A validation pass means the required geometric checks ran without hard
errors. It does not establish energetic stability or experimental feasibility.

## Validation and migration

`ValidationReport.to_dict()` provides the structured report used by Python,
the CLI, and MCP. It includes `is_valid`, `errors`, `warnings`,
`checks_performed`, and `checks_skipped`, along with the detailed findings.

- Empty structures, nonfinite geometry, and singular cells cannot pass.
- An interatomic distance below `(1 - bond_tolerance)` times the sum of
  covalent radii is a severe overlap. The default `bond_tolerance` is 0.3.
  This mandatory check includes bonded atoms and periodic images, even when
  `check_clashes=False`.
- Van der Waals proximity is reported in `close_contacts` as an advisory
  warning. These radius-based estimates are not proof of an invalid structure.
- Unusual bond lengths and metal coordination are advisory. Their existing
  detailed fields remain available for applications with stricter policies.
- Missing oxidation states produce unknown charge (`None`) and an explanation
  under `checks_skipped`; they no longer silently produce zero charge.
- Partial occupancies survive transformations, but bond inference requires an
  explicitly ordered model. Validation reports the required overlap check as
  unavailable until occupancies are resolved; it does not pass the structure.
- A newly constructed `ValidationReport()` no longer represents a pass until
  the required checks have actually run.

`mofforge validate input.cif --json` emits the structured report and returns 1
for invalid or incomplete geometry validation. Batch CLI runs return 1 if no
inputs match, an operation fails, or final geometry validation fails. Workflow
artifacts and their reports remain available for inspection.

MCP validation retains the existing top-level count fields and adds the full
report under `validation`. Its `success` field describes whether the operation
executed; `is_valid` describes geometry. `BuildResult` uses the same distinction:
success requires existing, parseable outputs, while `validation` and
`metadata['validation_by_output']` describe their scientific checks.

## Periodic geometry and metadata

Geometric minimum images use the complete lattice metric, including skewed
cells. `nearest_image(displacement, lattice)` is the lattice-aware form; the
one-argument call retains component-wise fractional wrapping for compatibility.

`Crystal.bonds` remains the simple graph used for pattern matching.
`Crystal.periodic_bonds` stores each physical bond once as a `PeriodicBond`
containing atom indices, an integer lattice image, and a distance.
`coordination_number(index)` counts image-resolved neighbors, including two
neighbors per self-image bond. `bond_vectors(index)` returns their Cartesian
vectors. Graph degree alone is insufficient for primitive-cell coordination.

Use `set_frac_coords` or `set_cart_coords` to update coordinates and refresh bond
geometry. Slicing, wrapping, addition, and replacement preserve species
compositions, oxidation states, occupancies, site properties, and provenance.
Replacing atoms does not invent oxidation states for new fragments.

Adsorption spacing and clustering respect periodic boundaries. Internal
adsorbate contacts are separate from host/guest contacts; severe overlap checks
still apply to every atom pair.

## Reproducing a workflow

Pass `random_seed` to replacement, functionalization, fragment generation,
adsorption, campaigns, or `Pipeline`. The CLI replacement option is
`--random-seed`. Batch YAML accepts a top-level `random_seed` and per-replacement
overrides. If a seed is omitted, a seed is generated and recorded in the manifest.
Random operations use local generators; fragment generation also
sets the RDKit embedding seed. UFF optimization failure or nonconvergence is
reported instead of silently accepting an unoptimized embedding.

Seeded reproduction is guaranteed within a fixed software environment, not
across dependency versions. Use the committed `uv.lock` and `uv sync --locked`
with the extras needed by the workflow. Backend assembly algorithms retain
their own convergence behavior; a seed does not establish physical validity.

Workflow-generated CIF/XYZ artifacts with provenance have a neighboring
`<filename>.json` manifest with file and structure SHA-256 hashes, effective
operation parameters, software versions, provenance, and any validation report.
Provenance schema 2 uses flat history, avoiding recursive growth. Existing
provenance JSON files remain readable through `Provenance.from_json`.

The unreleased [artifact verification follow-up](artifact_verification.md) adds
complete batch input descriptors and `mofforge verify`. Older manifests receive
explicit integrity-only verification; verified reload remains future work.

Batch and campaign filenames now include a hash of their input identity and
parameters. Update scripts that assumed the old exact output names to consume
the returned output paths. Batch ordering is stable across serial and parallel
execution. Campaign coverages must be in `[0, 1]`; zero means no modification.

TOBACCO code and default data are pinned to commit
`effe43059cfd78015db616890a6d20c090633378`. Downloads verify the recorded archive
SHA-256 before extraction. Cache keys include repository, revision, and digest,
and completed directories are published atomically. Custom repository/tag
downloads additionally require `MOFFORGE_TOBACCO_DATA_SHA256`. Explicit local
datasets remain supported; consumed TOBACCO input files are hashed in manifests.

## Release checks and reference cases

CI separates core-only, chemistry, MCP, and builder installations. Core checks
cover Python 3.10 and 3.14. Builder checks use Python 3.11 because the pinned
Pormake `jaxlib` dependency has no Python 3.14 wheel. CI also installs the built
wheel into a clean environment and verifies its bundled fragments.

The regression suite checks the skewed-cell distance counterexample, primitive
versus supercell coordination, metadata preservation, exact duplicate atoms,
severe overlaps, unknown charges, isolated CO2 placement, flat provenance,
seeded geometry, serial/parallel batch equivalence, and CLI/MCP report agreement.

| Existing reference fixture | Atoms | Required geometry outcome |
| --- | ---: | --- |
| IRMOF-1 | 424 | Pass; no severe overlaps or unusual bond lengths |
| UiO-66 | 912 | Pass; no severe overlaps or unusual bond lengths |
| MOF-74 | 162 | Pass; no severe overlaps or unusual bond lengths |

These fixtures retain advisory close contacts; UiO-66 also demonstrates why
element-only coordination ranges are warnings. The explicit reference outcomes
are enforced in `tests/test_reliability.py`.

Run real backend builds with `MOFFORGE_RUN_BUILD_TESTS=1 pytest
tests/test_build_integration.py`. Optionally set `MOFFORGE_TEST_TOBACCO_DATA` to
a staged pinned source archive to avoid downloading data during that test.

For adsorption timing, run `python scripts/benchmark_adsorption.py`. It checks
clustering against exhaustive lattice-aware selection before reporting timings.

The reliability implementation is merged in [PR #11](https://github.com/tdpham2/mofforge/pull/11).
The [roadmap](roadmap.md) expands the follow-up work into artifact verification,
dataset resolution, environment diagnostics, resumable campaigns, and scientific
evaluation, with dependencies and acceptance criteria.
