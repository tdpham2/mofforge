# Generic polymer-box construction

Updated: 2026-09-10. Scope: rework the POP functionality introduced by PR #17.

## Objective and boundary

Make mofforge build reproducible initial polymer boxes from user-supplied
molecules and explicit connection definitions. No experimental material recipe
is required. This replaces the earlier plans tied to styrene/divinylbenzene or
the Nd adsorption study.

Support two distinct operations in this PR:

1. **Pack:** place specified numbers of monomers, oligomers, or prebuilt finite
   polymer chains into a periodic box, preserving their existing topology.
2. **Connect:** create polymer connectivity through declared attachment sites
   and explicit graph edits, using bounded geometric placement and validation.

A packed monomer mixture is a valid packing result. Only the connect operation
reports newly formed polymer connectivity. Neither operation establishes
equilibration, permanent porosity, synthetic accessibility, or reaction kinetics.

Ownership is explicit:

| Owner | Responsibilities |
| --- | --- |
| Mofforge | Molecular geometry and topology, Packmol resolution and feature checks, input generation, execution, diagnostics, packing validation, explicit connection edits, and resumable structure artifacts. |
| MatKit | Potential selection and configuration, force-field assignment where needed, MLIP evaluation, LAMMPS setup, minimization, annealing, equilibration, production MD, DFT, charge prediction, and adsorption calculations. |
| Calling workflow | Coordination of construction and simulation stages, material preparation protocols, and comparisons across independent runs. |

Packmol performs geometric packing without an energy/force model. MLIP-based
relaxation is a subsequent MatKit stage. LAMMPS can run dynamics using forces
from a compatible MLIP; its required build packages and model configuration
belong to MatKit. See the [Packmol manual](https://m3g.github.io/packmol/userguide.shtml)
and [LAMMPS ML-IAP documentation](https://docs.lammps.org/pair_mliap.html).
MLIP checkpoints, potential parameters, and simulation settings are not packing
options in mofforge.

The basic workflow is mofforge packing followed by external MatKit preparation
and production simulation. Network construction can alternate bounded mofforge
connection steps with MatKit relaxation and explicit geometry imports. The
calling workflow owns this loop; `PopBuilder.build()` only composes mofforge's
packing and connection stages.

No MatKit code or dependency is added here. Automatic solvent removal, a workflow
engine, Nd placement, donor-site discovery, DFT cluster extraction, experimental
isotherms, and study orchestration are outside this PR.

## POP representation and density

Represent an amorphous POP as a periodic atomistic cell with explicit covalent
connectivity, including bonds across cell boundaries. A P1 CIF is an interchange
representation of this cell, not a claim of crystalline order. Packing supplied
finite chains is independently useful; constructing a crosslinked network from
precursors additionally requires connection rules.

The initial volume must be chosen, but the final experimental density need not
be known. Accept explicit box lengths or an initial packing density, deriving
`V_initial = total_packed_mass / initial_packing_density` with consistent units.
Do not provide a universal POP density default.

| Quantity | Meaning |
| --- | --- |
| `initial_packing_density` | Input used to size the initial box, or the equivalent value calculated from explicit box lengths. Retained as part of the packing history. |
| Current density | Current represented mass divided by current cell volume; recomputed after atom deletions and geometry/cell updates. |
| Final simulated density | A property assessed from the externally prepared structure or trajectory, with its preparation conditions recorded by the simulation workflow. |
| Experimental reference density | Optional comparison metadata with units and the measurement definition; not a required construction input or packing constraint. |

Record the caller's declared material state, such as dry/activated or
solvent-containing, without treating that declaration as proof of preparation.
Solvent must be supplied as explicit counted components and contributes to the
total packed mass. Selecting a solvent and defining a drying or activation
protocol belong to the external scientific workflow.

Initial density and connection history can influence the resulting pore
structure; low initial density alone does not establish permanent porosity.
Crosslinking degree and concentration during construction affect microporosity
in [published network-construction work](https://pubs.acs.org/doi/10.1021/ma200303p).
External preparation may use NPT simulations to evolve the volume or fixed-volume
simulations to retain an imposed density. The result depends on the chosen
potential and preparation procedure.

Recommend external comparisons across independent seeds, plausible initial
densities, and cell sizes. Assess final density, connectivity, accessible pore
volume, pore-size distribution, and the intended simulated property after
material-appropriate preparation. These comparisons and pore-property
calculations do not become a new mofforge workflow or prediction capability.

## Decisions from the review

- Keep `Crystal`, periodic geometry utilities, provenance, validation, and the
  `PopBuilder` facade. Replace the pysimm execution path in
  `src/mofforge/polymerize/engine.py`.
- Invoke Packmol directly for packing. Implement the connection selection and
  graph edits natively, following the staged construction idea used by Polymatic.
  Do not embed the Polymatic driver or port its force-field machinery. Algorithm
  inspiration does not imply numerical equivalence to Polymatic.
- Remove pysimm from the `pop` extra and regenerate the lockfile during
  implementation. Remove POP-specific LAMMPS configuration and the `forcefield`,
  `equilibrate`, and `md_settings` build options. Explicit use of removed options
  must raise an actionable error rather than silently change behavior.
- Treat SMARTS site detection as an optional annotation aid. The existing menu
  of five named reactions is not evidence that those reactions are implemented;
  do not advertise them as supported construction recipes.
- Limit initial packing and connection construction to fully periodic,
  orthorhombic boxes. Reject unsupported cell shapes explicitly. Existing MOF
  operations retain their current cell support.

## Implementation

### 1. Molecular inputs and explicit topology

- Extend `PopBuilder.add_monomer()` with an explicit positive integer `count`;
  the same operation accepts a prebuilt finite chain. Remove the implicit equal
  split and default twenty copies. Preserve registration order when assigning
  instance and atom identities.
- Support SMILES, single-molecule MOL/SDF files, and an in-memory `Crystal` with
  explicit connectivity. An XYZ input requires a companion molecular graph or
  an explicitly mapped molecular specification; `!` tags alone do not establish
  bonds or bond orders. Reject ambiguous inputs before invoking Packmol.
- Preserve supplied coordinates. For SMILES, retain RDKit conformer embedding
  and bounded local geometry cleanup; check its completion and record its seed.
  This preparation is not a user-facing scientific optimization capability.
- Require coordinate-bearing molecular inputs to represent hydrogens explicitly
  in both their graph and coordinates. Reject unresolved implicit hydrogens with
  guidance to prepare a complete molecular input. Expand hydrogens during SMILES
  preparation. Expose the stable prepared-template label mapping, including
  generated hydrogens, before connection edits are specified against it. Derive
  packed mass from the represented atoms so graph, coordinates, composition,
  and density calculations agree.
- Store immutable template atom IDs, instance IDs, bond orders, formal charges
  when supplied, and optional connector annotations. Derive functionality from
  the declared connectors; reject a contradictory explicit value.
- Stage inputs under generated internal names in a fresh run directory. Never
  write over a source file or derive writable paths directly from display names.
  Duplicate display names must not collide.

### 2. Working periodic packing

- Add `PopBuilder.pack(...)`. Require exactly one sizing input: `box_lengths`
  containing three positive lengths in angstroms, or a positive
  `initial_packing_density` in g/cm³. Density produces a cubic box using the
  exact component counts and represented masses. It is an initial packing
  choice, not a prediction of material density.
- Replace `target_density` with `initial_packing_density` across Python, CLI,
  MCP, configuration, and examples. Explicit use of `target_density` must raise
  an actionable rename error. Neither sizing input has a default, and providing
  both is an error. Record the equivalent initial density when lengths are used.
- Default the minimum distance between different packed molecules to 2.0 Å,
  expose it as a geometric parameter, and record the effective value. Validate
  counts, dimensions, density, distances, and seeds before running an engine.
- Run the resolved Packmol executable with an argument list, input on stdin,
  and a per-run working directory. Preserve stdout, stderr, generated input,
  version information, and the effective seed. Use a configurable timeout with
  a 300-second default; a timeout is a failed attempt with saved diagnostics.
- Use Packmol's orthorhombic `pbc` support, available from version 20.15.0.
  Establish and pin a real Packmol build in the integration-test environment;
  reject binaries without the required feature. The
  [Packmol manual](https://m3g.github.io/packmol/userguide.shtml) documents the
  periodic box and packing behavior.
- Replicate topology from the supplied templates; Packmol supplies coordinates.
  Check atom counts, element sequence, instance boundaries, and intramolecular
  geometry against the expected copies. Never infer new intermolecular bonds
  from proximity in the packed output.
- Independently validate distances across periodic boundaries, including
  nonbonded contacts with a molecule's own periodic images. Allow at most
  0.01 Å numerical violation of the requested packing distance. Reject malformed,
  incomplete, or geometrically invalid output even if the process exits zero.
  Do not silently enlarge the box, reduce separation, or drop molecules.

### 3. Native geometric connection steps

- Add `ConnectionRule` and `ConstructionState` types. A connector identifies
  its attachment atom, orientation anchor, role, and any replaceable atoms.
  A rule specifies compatible roles, the new bond order and permitted length
  interval, all atom deletions, internal bond-order changes, and formal-charge
  changes. Every edit references stable template labels. Reject incomplete
  rules; do not infer reaction chemistry from element proximity.
- Initially support pairwise connection rules. Multi-site concerted reactions,
  automatic mechanism discovery, and unrestricted reaction SMARTS execution
  are deferred. A generic connection rule can still create chains, branches,
  and crosslinks when its complete edits are supplied.
- Add `connect(state, rules, ...)` as a bounded operation. Require an explicit
  conversion target and candidate-attempt budget. Each committed event consumes
  one connector on each side; report conversion as consumed original connectors
  divided by the initial connector count, plus the number of new bonds.
- Select compatible candidates using periodic distances and a stable ordering
  by distance and atom IDs. For disconnected finite components, try bounded
  rigid alignment through the connector anchors and twelve seeded torsional
  orientations. Accept only candidates inside the rule's geometry limits that
  pass clash and local valence checks. Do not move an entire periodically
  percolating component as though it were a finite molecule.
- Use construction-specific geometric acceptance on the retained atoms after
  proposed deletions and bond edits. Require an explicit construction nonbonded
  minimum distance, recorded separately from the packing separation. Exempt
  exact bonded image pairs and pairs joined by two bonds in the edited periodic
  topology from this nonbonded cutoff; still apply severe-overlap, bond-length,
  and rule-geometry checks. All other nonbonded contacts, including periodic
  images, must pass the cutoff. Merging components must not exempt all contacts
  within the merged component. Existing `ValidationReport.is_valid` alone is
  insufficient because some nonbonded contacts are only warnings.
- Canonicalize periodic bond identity using both endpoints and the image offset,
  accounting for endpoint reversal. Determine percolation from the rank of
  independent net lattice translations around graph cycles: zero for finite
  components, and one, two, or three for periodic connectivity in that many
  independent directions. A finite molecule crossing a cell boundary is not
  necessarily percolating. Reconstruct a consistent unwrapped geometry for
  finite components before rigid placement.
- Same-component closures are opt-in (`allow_cycles=False` by default). When
  enabled, require the endpoints already to satisfy the geometry limits, reject
  duplicate periodic bonds, and enforce the rule's explicitly supplied minimum
  finite-cycle size. Evaluate ring size in the periodic topology with explicit
  atom images, not solely in the simple graph: a winding cycle can connect to
  another periodic image without forming a finite ring. Do not distort an
  existing network to force a closure.
- Apply edits transactionally on a copy. Retain an event ledger of consumed
  sites, new/changed/deleted bonds, removed atoms, and rejection reasons. Preserve
  untouched connectors as supplied; any termination/capping edits must be
  explicitly defined. Removed atoms remain accounted for in the ledger.
- Retain the best valid state when candidates are exhausted or the budget is
  reached. Report `partial` with the achieved conversion and diagnostic reason;
  do not claim success below the requested construction target. A connection
  can fail geometrically even when its graph edit is chemically plausible.
- Let `PopBuilder.build(...)` compose packing and connection steps. It requires
  connection rules; callers wanting packing alone must use `pack()`. External
  workflow code can instead call the stages individually and interleave external
  relaxation. `build()` does not launch MatKit or LAMMPS.

### 4. Box artifacts and external geometry updates

- Extend the existing periodic topology representation to retain bond orders
  without changing legacy constructors. Scope strict topology authority to
  `ConstructionState`: its periodic records generate a synchronized `Crystal`
  representation, including the simple NetworkX projection for existing searches.
  Native geometry refreshes must not infer bonds from proximity or reconcile
  topology from independent graph edits. Preserve legacy `Crystal` graph-edit
  synchronization and MOF replacement behavior; graph edits on a derived
  `Crystal` do not update the native state. Preserve all metadata through copy,
  subset, combination, wrapping, and native export.
- Write a versioned native box JSON alongside a P1 CIF and an XYZ interchange
  file. Native JSON contains species, coordinates, cell, periodicity, stable atom
  and component IDs, formal charges when known, bonds with image offsets and
  orders, connector state, event history, declared material state, initial packing
  density, current mass/volume/density, and provenance. Store experimental density
  references with their units and measurement definitions when supplied. CIF/XYZ
  alone are not lossless topology handoffs.
- Give the box state a separate versioned hash covering geometry, topology,
  connector state, and relevant metadata. Preserve legacy `structure_hash`
  semantics for existing artifacts; do not silently reinterpret old manifests.
  Link required files and hashes from the box manifest, publishing completion
  only after the complete set verifies. Retain inspectable failed attempts.
- Provide explicit save/load and `update_geometry` operations independent of
  MatKit. Geometry updates must identify the parent state and map every atom ID
  exactly once with unchanged species/topology. Accept new coordinates and
  supported cell lengths, update bond geometry, recompute current mass/volume/
  density, and invalidate prior geometric validation and derived results. Keep
  the original packing density in the history. Revalidate the updated geometry
  before further construction. Atom-deletion events likewise recompute current
  mass and density.
- Require unwrapped positions or explicit per-atom wrapping offsets when needed
  to preserve periodic bonds. Reject ambiguous image mappings, changed atom
  counts, and unannounced chemical edits. Loading an external coordinate file
  must never silently rebuild polymer connectivity. Geometry import preserves
  declared topology; it does not establish whether reactive MD changed chemical
  connectivity. Importing reactive topology changes requires a separate future
  interface and is outside this PR.
- An external caller may export, relax elsewhere, import updated geometry, and
  call `connect` again. The mofforge stage performs no simulation and does not
  require any particular calculator or orchestration framework.

### 5. Results, interfaces, and documentation

- Extend `POPResult` with the operation (`pack` or `connect`), a construction
  status (`completed`, `partial`, or `failed`), and the native state reference.
  Keep existing artifact, error, metadata, and validation fields. `success`
  requires completion of the requested operation and its geometric checks.
- Preserve status, achieved conversion, native-state references, artifact paths,
  and diagnostics in Python, CLI, and MCP responses for partial results, even
  though `success=False`. Failed attempts also expose saved diagnostics and any
  retained valid state when available. CLI exit status is zero only for completed
  operations and nonzero for partial or failed operations; machine-readable
  responses retain the distinct construction status.
- Packing reports requested/actual counts, composition, mass, box volume,
  initial packing density, current density, declared material state, and
  separation checks. Connection additionally reports conversion, component sizes,
  new bonds, and periodic connectivity rank. A connected graph is not by itself
  evidence of a three-dimensional percolating network. Construction completion
  does not label the current density as an equilibrated material property.
- Add a packing CLI command and matching MCP tool. Update the existing
  `polymerize` command/tool to require explicit connection definitions and expose
  the native implementation. Use one JSON configuration shape for component
  counts and connection rules across Python, CLI, and MCP.
- Gate chemistry inputs on RDKit availability and packing execution on Packmol
  availability. Reading saved boxes and inspecting results must not require
  external engines. Update `pop-doctor`, examples, capability descriptions, and
  `docs/polymerize.md` to describe initial box construction accurately.
- Since these interfaces originate in the unmerged PR, update its examples and
  tests together. Do not retain a nonfunctional pysimm compatibility backend or
  claim that old force-field options still work. Existing MOF APIs stay intact.

## Delivery order and acceptance

Deliver focused implementation commits in this order: native inputs/topology
and box persistence; direct Packmol packing; explicit geometric connections;
external geometry update support; CLI/MCP/docs and integration evidence.
Packing is the first independently usable milestone. PR #17 is complete only
after both the advertised packing and connection operations have real examples.

Required validation:

- Pack a single species, a counted mixture, and supplied oligomer chains using
  a real Packmol binary. Verify exact composition, unchanged intramolecular
  topology, periodic separation, and deterministic coordinates in the same
  pinned environment. Put this check in CI; absence of the binary in that job
  must fail rather than skip the only native integration test.
- Use generic bifunctional and trifunctional aromatic fragments with explicit
  replaceable hydrogens as connection fixtures. Check a finite chain, a branch,
  and an explicitly enabled periodic closure. Assert expected formulas, bond
  orders, connector consumption, image offsets, and geometric validity. These
  are construction fixtures, not a claimed polymerization mechanism.
- Add a seeded end-to-end example that performs real Packmol packing followed
  by multiple native connection events, verifies the saved bundle, reloads it,
  imports an explicitly mapped external geometry update, and resumes connection.
  Exercise `build()` as well as the individual stage interfaces. Independently
  positioned connection fixtures do not replace this combined acceptance test.
  The supplied geometry update requires no actual MatKit or MD execution.
- Reject coordinate-bearing inputs with unresolved implicit hydrogens and check
  stable labels for SMILES-generated hydrogens. Verify represented atom counts,
  molecular masses, and density-derived volumes agree. Test mutually exclusive
  sizing inputs, the removed `target_density` option, and density recomputation
  after atom deletions and cell changes while preserving initial packing history.
- Reject nonbonded construction contacts that the legacy validator only warns
  about, including a 1.5 Å unbonded C-C contact with a 2.0 Å construction cutoff.
  Test validation after leaving-atom deletion, topology-local exemptions, and
  contacts that become intracomponent after a merge.
- Distinguish wrapped finite molecules, finite rings, and networks with periodic
  connectivity ranks one and three. Test image-aware duplicate detection and
  ring-size checks independently of the simple graph projection.
- Run legacy MOF replacement regressions around graph/periodic-record
  synchronization, and check native states retain authoritative bond orders and
  images through geometry refreshes and derived `Crystal` representations.
- Exercise incomplete rules, invalid valence, attempted site reuse, bad counts,
  ambiguous XYZ topology, impossible packing, malformed output, timeout, and
  insufficient attainable conversion. Preserve diagnostics and partial states.
- Reproduce and prevent source-file overwrites and duplicate-name collisions.
  Verify that topology-only changes alter the native state identity, and that
  save/load and external geometry updates preserve stable IDs and periodic bonds.
- Test wrapping, cell changes, missing mappings, and interrupted artifact writes.
  A changed topology or incomplete bundle must not qualify for verified reload.
- Verify partial CLI/MCP responses retain achieved conversion, diagnostics, and
  a usable native-state reference; check CLI status codes for completed, partial,
  and failed operations.
- Run dependency-present and dependency-absent interface tests with explicit
  mocks; do not assert that pysimm or MCP happens to be absent on the developer's
  machine. Core installation and lock checks must pass without simulation extras.
- Run the affected polymer, periodic geometry, provenance/artifact, and CLI/MCP
  suites, then the existing full suite and Ruff. Report native execution evidence
  separately from unit-test counts. No MatKit execution is a prerequisite.

## Explicit defaults and exclusions

The scope includes explicit covalent connection construction as well as packing
supplied chains. Packing remains an independently usable milestone for callers
that already have finite polymer chains. Material-specific recipes and
experimental fitting are not required. Component counts, box lengths or initial
packing density, declared material state, and connection chemistry are user
inputs; no universal POP composition or density is invented. A generated root
seed is persisted when the caller omits one. Retain the 2.0 Å default packing
separation, 0.01 Å packing tolerance, 300-second Packmol timeout, twelve seeded
connection orientations, explicit conversion targets and attempt budgets, and
opt-in cycles. Construction nonbonded separation is an explicit separate input.

Packmol remains a packing dependency, RDKit supports molecular preparation and
local topology checks, and mofforge owns the construction state. Full Polymatic
integration, force-field typing, MD/MLIP/DFT execution, Nd adsorption, donor-site
placement, cluster extraction, pore-property prediction, broad artifact-system
rewrites, reactive-topology import, automatic solvent removal, workflow engines,
and unrelated review fixes are separate work.
