# Agent-Driven Linker Functionalization

mofforge lets an AI agent perform post-synthetic **linker functionalization** —
decorating a MOF's organic linker with a functional group (–NH₂, –NO₂, –F, …) —
without ever authoring 3-D geometry or SMILES.

The agent makes only **external decisions**:

1. **Which functional group** — chosen from a curated menu.
2. **Which site(s)** — chosen by integer index from a deterministic list.
3. **Coverage / concentration** — the fraction of linkers to modify.

All chemistry and geometry is generated deterministically by RDKit and then
validated by the existing find/replace/validate pipeline.

## Why not let the agent write SMILES?

Functionalization requires a correctly anchored replacement fragment with valid
3-D geometry (bond lengths, hybridisation, the exact attachment position).
LLM-authored SMILES/geometry is plausible-but-wrong often enough to silently
break frameworks. Instead, the "best position" is a first-class concept — the
aromatic C–H sites of the linker — and the agent selects among them by index.

## The workflow

```
MOFid ──► linker SMILES ──► find_sites ──► pick indices + group ──► functionalize ──► validate
```

### 1. Get the linker SMILES

Obtain the linker SMILES however you like (e.g. running MOFid on the CIF).
mofforge accepts the SMILES directly; it does not require MOFid to be installed.

### 2. Enumerate functionalizable sites

```python
from mofforge import find_functionalizable_sites

sites = find_functionalizable_sites("O=C(O)c1ccc(C(=O)O)cc1")  # BDC / IRMOF-1 linker
for s in sites:
    print(s.index, "symmetry class", s.symmetry_class, s.description)
```

Each site carries a **symmetry class**. This disambiguates an index selection:

* Indices in **different** classes (e.g. α vs β positions of a naphthalene
  linker) target genuinely different chemical environments.
* Multiple indices in the **same** class on one ring define a **substitution
  pattern** (mono- vs di-substitution — ortho/meta/para). Framework-wide loading
  is controlled separately by `coverage`.

Metal-binding carboxylate carbons carry no aromatic hydrogen and are never
returned, so a site can never be a connection point.

### 3. Choose a functional group

```python
from mofforge import available_groups
print(available_groups())
# ['Br', 'CH3', 'CN', 'COOH', 'Cl', 'F', 'H', 'NH2', 'NO2', 'OCH3', 'OH', 'acetamido']
```

### 4. Functionalize

```python
from mofforge import functionalize

result = functionalize(
    "IRMOF-1.cif",
    "O=C(O)c1ccc(C(=O)O)cc1",
    group="NO2",
    sites=0,          # index (or list of indices) from find_functionalizable_sites
    coverage=0.5,     # functionalize half of the matched linkers
    output_cif="IRMOF-1-NO2.cif",
)
print(result.n_matches, result.n_functionalized, result.is_valid, result.clashes)
```

`coverage` maps onto the number of linkers modified: `0.5` on a framework with
24 BDC linkers functionalizes 12 of them.

## Autonomous campaigns

`run_campaign` sweeps groups × coverages, validates each result, and returns
them ranked best-first (valid structures first, then fewest steric clashes):

```python
from mofforge import run_campaign

results = run_campaign(
    "IRMOF-1.cif",
    "O=C(O)c1ccc(C(=O)O)cc1",
    groups=["NH2", "F", "NO2"],
    coverages=[0.25, 0.5, 1.0],
    output_dir="campaign/",
)
best = results[0]
print(best.group, best.coverage, best.output_cif)
```

## MCP tools (for AI agents)

The stock server (`mofforge-mcp`) and the ChemGraph HPC server
(`mofforge-mcp-chemgraph`) expose:

| Tool | Purpose |
|------|---------|
| `mofforge_find_sites` | List functionalizable sites (index + symmetry class) for a linker SMILES |
| `mofforge_list_functional_groups` | The curated group menu |
| `mofforge_functionalize` | Functionalize a group at chosen site(s) with a coverage |
| `mofforge_functionalize_campaign` | Sweep groups × coverages, return ranked results |
| `mofforge_list_fragments` / `mofforge_get_fragment` | Discover the packaged moiety library (for manual `mofforge_replace`) |

On the ChemGraph server, `mofforge_functionalize` and
`mofforge_functionalize_campaign` run as backend tasks (RDKit embedding +
replacement is CPU-bound); site listing and the group menu run inline.

## Adding a new functional group

Groups are defined once, offline, in
`src/mofforge/functionalize/groups.py` as a SMILES fragment plus the index of
the attachment atom — no coordinates. Adding a group is a one-line edit; the
geometry is generated automatically at functionalization time.
