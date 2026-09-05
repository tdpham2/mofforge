# mofforge roadmap

Updated 2026-09-05 against the merged 0.2 reliability implementation.

The next objective is to make screening studies restartable, easy to inspect,
and useful for scientific comparisons. Start with complete artifact records and
safe execution, then add better candidate selection and physical evaluation.
The commands and APIs below are proposals unless listed as completed.
The [next-session handoff](next-session-plan.md) records the working-tree state
and a concrete starting point for implementation.

## Completed foundation

[PR #11](https://github.com/tdpham2/mofforge/pull/11) is merged. All six CI jobs
passed: core on Python 3.10 and 3.14, chemistry, MCP, builders, and wheel
installation. This marks implementation and CI completion; package publication
is a separate release step.

- Lattice-aware minimum images and periodic coordination.
- Metadata preservation during in-memory transformations.
- Mandatory geometry checks and explicit errors, warnings, and skipped checks.
- Local recorded seeds, flat provenance, and generated artifact manifests.
- Stable batch ordering and parameter-based output names.
- Verified TOBACCO downloads and shared builder artifact validation.
- Scientific regressions, real backend checks, and a clustering benchmark.

See [0.2 behavior and migration notes](reliability.md) for the current contract.

## Priority and delivery order

Effort describes scope, not a delivery promise: small is one focused PR, medium
is two or three, and large should be split into independently useful milestones.
Release groupings are provisional.

| Order | Improvement | Why it matters | Target | Effort |
| --- | --- | --- | --- | --- |
| 1 | Complete artifact identity and verified reload | Prevent stale reuse and preserve history across sessions | Next 0.2 follow-up | Medium |
| 2 | Unambiguous dataset resolution | Ensure a screening result uses the intended CIF variant | Next 0.2 follow-up | Small |
| 3 | Environment diagnostics and workflow preflight | Find setup and configuration errors before expensive work | 0.3 | Medium |
| 4 | Resumable batches and campaigns | Recover interrupted studies without repeating completed work | 0.3 | Large |
| 5 | Broader scientific reference suite | Measure correctness and known limitations beyond smoke tests | Start now; ongoing | Medium |
| 6 | Better adsorption candidates and duplicate detection | Improve usable candidate yield and reduce repeated calculations | After items 1 and 5 | Medium |
| 7 | Property evaluation and optional relaxation | Rank candidates by the research objective | 0.4 proposal | Large |
| 8 | Measured scaling and study reports | Make larger studies affordable and results easy to compare | Alongside 0.3–0.4 | Medium |

## 1. Complete artifact identity and verified reload

This follow-up implements the first slice: complete batch input descriptors,
captured worker settings, and read-only artifact verification. See
[artifact verification](artifact_verification.md) for the API and compatibility
details. It is not yet a published release.

Remaining gaps: [CIF loading](../src/mofforge/core/crystal.py) creates a new load record
without restoring the adjacent manifest. Structure files and manifests are
written separately, and [provenance](../src/mofforge/provenance.py) has no verified
reload or replay entry point. Other workflow types still need complete input
descriptors before they can support verified reuse.

Build a common input descriptor containing consumed file hashes, effective
bonding and validation settings, seeds, and an environment fingerprint. Keep
file location separate from scientific identity. Define a versioned artifact
bundle with lossless structure metadata, species labels, periodic bonds, and
provenance; retain CIF/XYZ for interchange. Build explicit verified loading
on the new `mofforge verify` functionality, with migrations for existing records.
Stage outputs and mark a result complete only after its files and hashes verify.

Acceptance criteria:

- Editing a replacement XYZ at the same path changes the job identity.
- Write, reload, and continue a workflow without losing history, molecule IDs,
  occupancies, or periodic bond information.
- A changed CIF, missing sidecar, unsupported schema, or interrupted write is
  reported explicitly and cannot be reused as a completed checkpoint.
- Keep exact file integrity hashes separate from tolerance-based geometric
  comparisons; CIF rounding must not be mistaken for identical serialization.

## 2. Unambiguous dataset resolution

The [CoRE MOF resolver](../src/mofforge/coremof/structures.py) tries known filename
variants and can fall back to the first lexical prefix match. A request with
multiple ASR/FSR or other processing variants needs an explicit selection policy.

Return the exact identifier, dataset release, processing variant, resolved path,
and file hash. List candidates when the request is ambiguous. Keep any fuzzy
lookup as an explicit option and expose the same decision through Python and MCP.
An indexed local catalog can later replace repeated recursive scans.

Acceptance criteria: ambiguous fixtures never silently select a structure;
explicit variant requests resolve deterministically; manifests record the
selected dataset identity. Use synthetic or redistributable fixtures in CI.

## 3. Environment diagnostics and workflow preflight

[Capability detection](../src/mofforge/mcp/tool_selection.py) currently checks
module discoverability. That does not establish that a builder, dataset, browser,
or remote execution environment is usable. [Batch configuration](../src/mofforge/batch.py)
also defers several operation errors until individual inputs are processed.

Add a proposed `mofforge doctor --json` with per-capability readiness, versions,
configuration sources, dataset/cache status, and actionable error messages.
Default diagnostics should be offline and avoid downloads or job submission.
Share typed operation validation across batch YAML, Python, CLI, and MCP. A
preflight mode should resolve inputs, reject unknown keys, validate options,
and estimate job count and grid size before execution.

Acceptance criteria: malformed configuration fails before creating outputs;
missing browser binaries are distinguishable from a missing Python package;
offline diagnostics produce useful results without network access. Existing
Python entry points and documented configuration aliases remain supported.

## 4. Resumable batches and campaigns

[Batches](../src/mofforge/batch.py) use a process pool, while
[functionalization campaigns](../src/mofforge/functionalize/campaign.py) run a
serial sweep. Both return in-memory result lists and lack a durable run ledger.
Deterministic filenames alone are insufficient for safe resumption.

Start with a local SQLite ledger and one coordinator that records job identity,
root and child seeds, attempts, artifacts, and states such as pending, running,
succeeded, failed, and cancelled. Persist the generated root seed before the
first task starts. Add proposed `--resume` and `--retry-failed` options, bounded
submission, per-job logs, and progress counts. Reuse the same execution contract
for campaigns and the existing ChemGraph adapter.

Acceptance criteria:

- Interrupt a seeded study, restart with a different worker count, and obtain
  the same geometries and stable result ordering as an uninterrupted run.
- Reuse only completed artifacts whose input identity and hashes still match.
- Recover abandoned running jobs and record failed attempts without losing
  completed work; distinguish execution failure from invalid geometry.
- Native jobs have enforceable time limits in isolated worker processes;
  timing out a future alone is not treated as terminating its calculation.

Deliver local resume first, then parallel campaigns and cancellation, then
remote execution integration. Items 1 and 3 are prerequisites.

## 5. Broader scientific reference suite

The existing [reference regressions](../tests/test_reliability.py) cover three
frameworks and important geometry counterexamples. The
[real builder tests](../tests/test_build_integration.py) verify usable artifacts
and that checks ran; they do not yet enforce reference composition, connectivity,
coordination, or a reviewed expected validation verdict for each generated MOF.

Expand fixtures across metal families, skewed and primitive cells, missing
hydrogens, partial occupancy, guest loading, and disconnected fragments. Add
invariance checks under atom reordering, cell translations, and supercell
expansion. Record expected findings and their scientific rationale, including
cases where the supported outcome is an explicit inability to validate.

Acceptance criteria: each backend has reviewed reference outcomes with numeric
tolerances; intentionally broken structures trigger the expected finding;
geometry changes cannot silently reduce coverage by skipping a required check.
Add separate browser and ChemGraph integration jobs, and document any external
environment requirements rather than counting skipped tests as integration passes.

## 6. Better adsorption candidates and duplicate detection

[Placement](../src/mofforge/adsorbate/placement.py) currently chooses one
orientation per selected site, filters site centers, and validates afterward.
It can place fewer molecules than requested. Campaign ranking uses execution
status, geometry validity, and severe-clash count.

Generate a bounded number of seeded orientation/site trials, screen full
molecular geometry, and return explicit requested and achieved loading with
rejection reasons. Retain the present fast placement mode. Remove equivalent
candidates before expensive evaluation, with explicit tolerances and species,
occupancy, and host/guest comparison policies.

Evaluate the existing pymatgen dependency's
[StructureMatcher](https://pymatgen.org/pymatgen.core.html#pymatgen.core.structure_matcher.StructureMatcher)
for periodic equivalence; its comparison, scaling, and supercell options need a
deliberate policy. Keep this equivalence test separate from exact artifact hashes.

Acceptance criteria: report a loading shortfall explicitly; reproduce accepted
poses with a fixed seed; improve valid placement yield at a fixed trial budget
on the reference set; detect reordered equivalent structures without merging
chemically distinct functionalizations.

## 7. Property evaluation and optional relaxation

The current campaign ranking is a geometry filter. Add a pluggable evaluation
stage that records named metrics, units, methods, and unavailable values, then
ranks by a user-selected research objective. Start with an offline table of
candidate composition, density, requested/achieved modification, and validation.
Add pore descriptors through a tested adapter, and recompute properties after
modification instead of inheriting the parent's database values.

For physical evaluation, use an optional calculator adapter with explicit
element, charge, spin, and periodicity support. ASE provides
[optimizers with force convergence criteria, trajectories, and restart support](https://docs.ase-lib.org/ase/optimize.html),
which makes it a candidate for this adapter. Preserve both generated and relaxed
structures, calculator/model revision, constraints, convergence, and forces;
validate again after relaxation. Use the existing ChemGraph integration when
remote execution is needed.

Acceptance criteria: an unsupported or unconverged calculation cannot receive
a successful physical score; comparisons use the same documented protocol;
adsorption-energy results state their host and guest reference conventions;
benchmark results agree with reviewed references within declared tolerances.
Keep geometry validity, optimization convergence, and application performance
as separate results.

## 8. Measured scaling and study reports

The clustering benchmark is useful, but the
[void finder](../src/mofforge/adsorbate/sites.py) still materializes the full grid,
and [pattern matching](../src/mofforge/search/isomorphism.py) collects all
isomorphisms. Profile complete workflows before selecting further optimizations.

Stream grid generation, bound queued work, cache prepared fragments using their
full input identity, and add explicit limits or iterators for combinatorial
search. A limited search must report truncation. Measure peak memory, wall time,
throughput, and valid/unique candidate yield on fixed small and large cases.

Export the run ledger to CSV/JSON and a shareable report containing structures,
warnings, scores, failure reasons, and links to manifests. Users should be able
to compare one study without writing a custom script or requiring a hosted UI.

Acceptance criteria: documented resource budgets on the benchmark hardware;
optimized paths preserve reference results; reports account for every attempted
job and distinguish failed, invalid, incomplete, and successfully scored results.

## First implementation sequence

1. Hash all consumed batch inputs and add artifact verification — implemented in this follow-up.
2. Add verified reload and durable artifact completion.
3. Make dataset variant selection explicit.
4. Add doctor and shared configuration preflight.
5. Implement local resume with interruption regression tests.

Expand scientific references alongside these PRs. Start physical scoring after
the reference and artifact contracts are established. Defer a hosted dashboard,
additional builder backends, automatic disorder resolution, and learned ranking
until measured study outcomes justify their additional scope.
