# MOFForge examples cookbook

These **36 executable lessons** teach the current Python API to MOF researchers
with basic Python experience. Each topic guide connects the code to a scientific
question, explains inputs and parameters, describes expected results, and offers
an extension exercise. The scripts execute real library calls and write
`summary.json` reports that the example tests inspect.

## Start here

From a source checkout, prepare the core environment:

```bash
uv sync --locked --extra dev
uv run --no-sync python examples/fundamentals/structure_io.py
uv run --no-sync python examples/search/pattern_matching.py
uv run --no-sync python examples/modify/linker_functionalization.py --nb-loc 6 --seed 42
```

Use the checkout rather than the installed wheel alone: the wheel includes
moiety assets, while this cookbook and its crystal inputs belong to the source
distribution. Commands below assume the repository root and an activated
environment, so `python` imports this checkout's installed mofforge. From another
directory, pass the absolute path to a script. Bundled input paths are anchored
to the script. User-supplied relative paths resolve from your working directory.

Follow **fundamentals → search → modification → validation → automation** first.
Then choose construction, databases, or adsorption and finish with a complete
workflow. The curriculum includes extra lessons on periodic geometry, custom
bond cutoffs, and conservative cleanup because these affect scientific results.

## Find a lesson

| Topic | Lessons | Guide |
| --- | ---: | --- |
| Fundamentals | 4 | [Walkthrough](fundamentals/README.md) |
| Search | 2 | [Walkthrough](search/README.md) |
| Fragment modification | 3 | [Walkthrough](modify/README.md) |
| SMILES functionalization | 2 | [Walkthrough](modify/functionalization.md) |
| Repair | 3 | [Walkthrough](repair/README.md) |
| Construction | 4 | [Walkthrough](build/README.md) |
| Adsorbates | 2 | [Walkthrough](adsorbate/README.md) |
| Validation | 1 | [Walkthrough](validation/README.md) |
| Visualization | 2 | [Walkthrough](vis/README.md) |
| Automation | 5 | [Walkthrough](pipeline/README.md) |
| Databases | 3 | [Walkthrough](databases/README.md) |
| MCP and ChemGraph | 3 | [Walkthrough](mcp/README.md) |
| Complete workflows | 2 | [Walkthrough](workflows/README.md) |

Every script, prerequisite, and guide is also recorded in
[catalog.json](catalog.json). The test suite checks that the catalog covers all
lesson scripts and that its links exist.

## Optional dependencies

Install only the extras for your chosen topic. Re-running `uv sync` replaces
the selected extras, so include every extra you want to keep.

| Capability | Checkout setup |
| --- | --- |
| Core, databases with demo data, HTML generation | `uv sync --locked --extra dev` |
| SMILES conversion and functionalization | `uv sync --locked --extra dev --extra chem` |
| Construction | `uv sync --locked --python 3.11 --extra dev --extra build` |
| MCP client/server | `uv sync --locked --extra dev --extra mcp` |
| PNG rendering | `uv sync --locked --extra dev --extra vis`, then `uv run --no-sync playwright install chromium` |
| Combined scientific examples | `uv sync --locked --python 3.11 --extra dev --extra all` |
| ChemGraph adapter lesson | MCP extra plus `pip install langchain-mcp-adapters` in the same environment |

Python 3.11 is the existing builder CI baseline. TOBACCO also needs its external
runtime data; see the [construction guide](build/README.md). PNG rendering and
viewing the generated HTML load 3Dmol.js from its CDN. Writing HTML itself is
offline. Full CSD examples require your licensed export; real CoRE examples
require separately obtained metadata and CIF files. The bundled database
demonstrations need neither a license nor downloads.

## Outputs and interpretation

Each script accepts `--help` and `--output-dir` (`-o`). Defaults go to the ignored
`examples/_outputs/<script stem>/` directory. Repeat runs may overwrite outputs
there; select a different directory to keep a comparison. Legacy `--output`
flags override the primary CIF path; summaries still go to `--output-dir`.
The MCP builder retains `--output-root` as an output-directory override.

Seeded lessons accept `--seed` (`-s`, default **42**). A builder lesson's linker
selection seed does not establish deterministic behavior inside an external
construction backend. Reports include installed package versions; reinstall
the checkout if editable-install version metadata is stale.

Read three independent kinds of result:

- **Execution:** the operation produced a result. Unexpected failures exit nonzero.
- **Validation:** geometry errors, warnings, and unperformed checks. A valid
  geometry report does not establish chemical stability or optimized energetics.
- **Verification:** recorded file integrity and, for batch outputs, complete input
  compatibility. A valid geometry can have incomplete artifact verification.

The validation and verification lessons deliberately include negative controls.
They exit successfully when those controls are correctly identified. Campaign
ranking uses validation and clashes; guest placement constructs starting
geometries and does not compute adsorption uptake or binding energies.

## Data and testing

Crystal and moiety inputs live in [data](data/). Database records are explicitly
fictional; [their README](data/databases/README.md) documents every geometry mapping.
Real input metadata is copied to the chosen output directory so SQLite caches
do not modify its source directory.

```bash
python -m pytest tests/test_examples.py
python -m ruff check examples tests/test_examples.py
MOFFORGE_RUN_BUILD_TESTS=1 python -m pytest tests/test_examples.py -k build
MOFFORGE_RUN_RENDER_TESTS=1 python -m pytest tests/test_examples.py -k png_rendering
```

Core lesson tests block network connections and optional-package imports, even
when those packages are installed. Missing optional dependencies produce
explained pytest skips. Installed optional
dependencies with runtime failures produce test failures. Consult
[verification notes](VERIFICATION.md) for the environment and optional runs
actually exercised during cookbook development.
