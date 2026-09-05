# Substructure search

These core-only lessons use the bundled IRMOF-1 structure. Searches compare
element-labeled connectivity; they do not optimize geometry or measure binding
energies. Begin with [structure preparation](../fundamentals/README.md).

## 1. Match an XYZ fragment

```bash
python examples/search/pattern_matching.py
python examples/search/pattern_matching.py --crystal examples/data/crystals/IRMOF-1.cif --query p-phenylene.xyz
```

[pattern_matching.py](pattern_matching.py) loads the periodic parent graph,
loads the fragment and its molecular graph, calls `find_pattern`, and extracts
the matched atoms to `matched_atoms.xyz`.

For the default parent/query, expect **24 locations**, **96 isomorphisms**, and
**four orientations per location**. A location identifies a set of parent atom
indices. An orientation maps query indices onto those atoms. The first mapping
in `summary.json` is `query index → parent index`; JSON represents its keys as
strings.

`--fragment-path` changes the fragment directory, while `--query` selects its
filename. `--crystal` accepts your own CIF. A zero-match search exits with an
explanation because this lesson promises a positive match. Missing explicit H
atoms or different bond inference can change a chemically plausible match.

**Try:** compare `p-phenylene.xyz` with a more specific tagged query and examine
orientation counts. The `!` marker is stripped for matching; its deletion
meaning is explained in the [modification guide](../modify/README.md).
The symmetry lesson remains at [symmetry_analysis.py](symmetry_analysis.py) and
is explained in the fundamentals guide.

Equivalent CLI:

```bash
mofforge search --parent examples/data/crystals/IRMOF-1.cif --query examples/data/moieties/p-phenylene.xyz
```

## 2. Search supported string patterns

```bash
python examples/search/string_pattern_search.py
```

[string_pattern_search.py](string_pattern_search.py) parses each pattern into
a graph and reports both locations and mappings:

| Pattern | Interpretation | Default locations |
| --- | --- | ---: |
| `[Zn]-[O]` | Zn–O contacts in the inferred graph | 128 |
| `[Zn]-[O]-C` | Metal–oxygen–carbon chains | 96 |
| `O-C-O` | Three-atom oxygen–carbon–oxygen paths | 48 |
| `[Zn]-[*]` | Any graph neighbor of Zn | 128 |
| `C1-C-C-C-C-C-1` | Six-membered carbon ring | 24 |
| `C-H` | Explicit carbon–hydrogen bonds | 96 |
| `[Xe]-[Xe]` | Intentional absent motif | 0 |

An unclosed `C1-C` ring is an expected parser error. The implementation is
**SMARTS-like**: it does not implement aromatic bond semantics, bond orders,
branching, stereochemistry, or general RDKit SMARTS predicates. Some unsupported
characters are ignored by tokenization rather than rejected; successful parsing
does not certify full SMARTS semantics.

**Try:** inspect the graph returned by `parse_smarts` before applying a new
pattern. Use an XYZ query when explicit graph geometry and anchor masks are
needed for replacement. String search results are for discovery; their dummy
query crystal is not a replacement scaffold.
