# Next-session plan

Updated 2026-09-05 after implementing the first input identity and artifact
verification slice. Read the [roadmap](roadmap.md) for the full priorities and
[artifact verification](artifact_verification.md) for the new contract.

## Current state

- The 0.2 reliability work is merged in
  [PR #11](https://github.com/tdpham2/mofforge/pull/11). All six CI jobs passed.
- The current branch is `artifact-verification`, based on merge commit
  `1764de3d47fba03ed1c736fb160150bc4297c707`. Refresh remote state before further
  branching or PR creation.
- This follow-up groups input records, artifact verification, regression tests,
  and documentation into separate commits. Inspect `git status` before resuming
  and preserve any additional local changes.
- `dac_mof_design/` is pre-existing, untracked user work. Preserve it and exclude
  it from unrelated commits.
- User preference: group changes by topic and use clear PR titles, descriptions,
  and commit messages without conventional prefixes such as `feat:` or `fix:`.

## Implemented first slice

- [inputs.py](../src/mofforge/inputs.py) provides versioned descriptors covering
  consumed file hashes, effective options, scientific settings, seeds, and the
  execution environment, including editable package source changes.
- [batch.py](../src/mofforge/batch.py) records complete descriptors, captures
  settings before worker submission, and detects changed inputs before writing
  outputs. Existing seeded replacement geometry is preserved through the
  recorded `batch-v1` seed strategy. Output identities and names change.
- [artifacts.py](../src/mofforge/artifacts.py) exposes `verify_artifact` and
  `VerificationReport`; the CLI exposes `mofforge verify` and `--json`.
- File integrity and verified workflow inputs are separate. Legacy manifests
  without descriptors receive integrity-only verification, not a full pass.
- [test_artifacts.py](../tests/test_artifacts.py) covers changed fragments,
  partial/missing outputs, malformed and legacy manifests, changed scientific
  settings, read-only verification, and serial/parallel execution agreement.

The first slice does not implement resumption, lossless reload, or atomic
publication. Other workflow types still need complete descriptors before they
can support verified reuse.

## Next implementation slice: verified reload and durable completion

The objective is to continue a saved workflow without losing scientific metadata
or trusting a partially written result.

1. Define a versioned lossless bundle alongside CIF/XYZ interchange files. Store
   the native structure representation, species labels, periodic bonds, molecule
   IDs, and provenance. Keep its hash distinct from the interchange-file hash.
2. Add an explicit verified loading API. Validate schema, required files, and
   hashes before restoring metadata. A normal CIF load should not silently
   invent missing history or mark a legacy artifact fully verified.
3. Extend the verifier to cover every required bundle member and define the
   supported migration outcomes for old records. Retain read-only verification.
4. Stage outputs and publish a completion record only after all required members
   exist and verify. A pair of separate file renames is not a single atomic
   transaction; use a publication protocol that survives interruption.
5. Add interruption and round-trip regressions and update migration guidance.

Acceptance criteria:

- Write, reload, and continue a workflow with its history, occupancies, molecule
  IDs, species labels, site properties, and periodic bonds intact.
- Missing or modified bundle members and unsupported schemas prevent a complete
  verification result; incomplete outputs remain inspectable.
- An interrupted write cannot appear as a reusable completed result.
- CIF rounding does not cause false comparisons with the original native
  serialization; exact integrity and geometric equivalence stay separate.
- Generated and reloaded seeded workflows agree in a fixed software environment.

## Following priorities

1. Explicit CoRE MOF dataset and processing-variant selection.
2. Offline `mofforge doctor` and shared workflow configuration preflight.
3. Local SQLite run ledger, persistent root seeds, `--resume`, and
   `--retry-failed`, then bounded parallel campaigns and cancellation.
4. Broader scientific references, better adsorption trials, periodic duplicate
   detection, property evaluation, and optional relaxation.
5. Measured scaling and shareable study reports.

Expand scientific references alongside the infrastructure work. No particular
calculator model or hosted application has been selected.

## Verification and completion

The first slice passed 27 new artifact regressions. Full local suites passed:
Python 3.10 core (418 passed, 19 skipped), Python 3.11 with chemistry/MCP/builders
(541 passed, 9 skipped), and Python 3.14 with chemistry/MCP (492 passed,
13 skipped). Ruff, whitespace checks, and relative documentation links passed.
These are local results; consult the PR checks for the current remote CI status.

Use the committed lockfile and the extras required by changed paths. Python
3.11 was the tested builder environment; temporary environments from this
session may not persist. Start with artifact, provenance, batch, workflow, and
CLI regressions, then run the affected broader suites, Ruff, and `git diff --check`.
Exercise native builders when their behavior changes. Browser and ChemGraph
integration need separately configured environments.

Keep the roadmap/documentation changes identifiable separately from runtime
changes when preparing commits. Record actual validation and compatibility
changes in the eventual PR. Do not include the user's untracked work.

## Suggested prompt to resume

> Read docs/next-session-plan.md and docs/roadmap.md. Continue with verified
> reload and durable artifact completion, preserving the existing local work
> and dac_mof_design/. Group changes by topic and avoid prefixes in commit and
> PR messages.
