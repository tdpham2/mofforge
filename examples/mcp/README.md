# MCP clients and ChemGraph integration

Install the `mcp` extra. These examples launch a real stdio subprocess with
the current Python interpreter by default, initialize an MCP session, discover
tools, and call them. They do not require an LLM or an API key.

## 1. Discover and call a small tool catalog

```bash
python examples/mcp/local_client.py --seed 42
```

[local_client.py](local_client.py) starts the stock server with an explicit
three-tool allowlist:

1. Discover only `mofforge_search`, `mofforge_replace`, and `mofforge_validate`.
2. Search phenylene in the bundled IRMOF-1.
3. Replace two locations using the seed and save `modified.cif`.
4. Validate the returned structure and retain all JSON results.

`--server-python` selects another interpreter when the server dependencies live
in a different environment. Paths sent over MCP are absolute. Outputs go to the
configured lesson directory. Protocol errors and unsuccessful operations fail
the script; a successful validation call can still report geometric warnings
or unperformed checks.

The implementation accepts structured or text JSON results from the MCP SDK.
Inspect the advertised input schemas before adapting a call. Stock validation
can report skipped graph checks when its loaded input has no bond graph; read
`checks_performed` and `checks_skipped` rather than assuming API parity from a
success flag alone.

**Try:** run these startup configurations:

```bash
mofforge-mcp --tools mofforge_search,mofforge_validate
mofforge-mcp --available-tools-only
```

The first restricts registration; the second hides tools with unavailable
optional packages. Explicit unavailable or unknown tool names are startup
errors. Stop the server after examining it; the lesson manages its own process.

## 2. Construct through MCP

```bash
python examples/build/mcp_build_test.py
python examples/build/mcp_build_test.py --backend tobacco --tobacco-data /path/to/tobacco-data
```

[mcp_build_test.py](../build/mcp_build_test.py) requires both `mcp` and the selected
builder in the server environment. It defaults to Pormake's tested N108/E1/pcu
recipe. TOBACCO mode uses the dmc fixture from the construction lesson.

Only `mofforge_build` is exposed. The client inspects the operation's success,
returned CIF paths, atom count, and validation. Backend build accuracy is not
an argument of the current stock MCP build tool, so this lesson uses its
existing defaults rather than forwarding an unsupported setting.

The old `--server-python`, `--tobacco-data`, and `--output-root` flags remain;
there are no machine-specific interpreter or data defaults.

**Try:** give a compatible custom building-block path using `--node` and
`--edge`. Check errors in `summary.json` if the server environment cannot
import the requested backend.

## 3. Load tools for a ChemGraph agent

```bash
pip install langchain-mcp-adapters
python examples/mcp/chemgraph_integration.py
```

[chemgraph_integration.py](chemgraph_integration.py) uses the external-MCP path
described by mofforge's [ChemGraph guide](../../docs/chemgraph.md). It creates a
`MultiServerMCPClient`, retrieves LangChain-compatible tools, and directly invokes
the validation tool. Those same tool objects can be supplied to a ChemGraph
agent. No LLM is needed to verify this transport/adapter connection.

The adapter dependency is separate from mofforge's stock MCP extra. The checked
API is documented in the
[upstream adapter examples](https://github.com/langchain-ai/langchain-mcp-adapters#multiple-mcp-servers).
This local lesson tests the adapter path, not ChemGraph's remote scheduler.

**Try:** expose search and validation together and inspect their schemas before
binding them to an agent. Keep the scientific input/output paths explicit.

## Configure CGFastMCP for HPC

The optional `mofforge-mcp-chemgraph` entry point requires a ChemGraph build
providing `chemgraph.mcp.cg_fastmcp.CGFastMCP`. The repository integration guide
identifies the `dev-globus-hpc` source branch for that API. Install and configure
the execution backend according to that ChemGraph checkout and your site.

Before launch, configure paths visible to both the server and compute workers:

```bash
export MOFFORGE_COREMOF_DATA_PATH=/shared/coremof/metadata.csv
export MOFFORGE_COREMOF_STRUCTURES_PATH=/shared/coremof/structures
export MOFFORGE_LOG_DIR=/shared/results/mofforge
export MOFFORGE_MCP_JOBS_FILE=/shared/results/mofforge/jobs.json
mofforge-mcp-chemgraph
```

Create the output directory first and provide your site's backend configuration.
The entry point initializes the backend, registers job-management tools,
persists job tracking, and shuts the backend down when the server exits.
Use shared paths or ChemGraph's staging hooks for worker access.

The current fan-out operation accepts a **nested `params` object**. A client
call's arguments take this form:

```json
{
  "params": {
    "metal": "Zn",
    "pld_min": 3.8,
    "water_stability_min": 0.7,
    "adsorbate": "CO2",
    "n_adsorbates": 2,
    "limit": 2,
    "structures_dir": "/shared/coremof/structures",
    "output_dir": "/shared/results/mofforge"
  }
}
```

Pass it to `mofforge_screen_and_place`. One worker resolves and loads each
candidate structure. Inspect per-candidate missing-structure/failure results;
a submitted job is not a completed structure. This fan-out wrapper currently
does not expose a root placement seed; do not add unsupported screening keys
to `params` or claim seeded HPC repeatability.

Discover `check_job_status`, `get_job_results`, `list_jobs`, `cancel_job`, and
`check_endpoint_status` and use their advertised schemas for your installed
ChemGraph version. Preserve their names in an explicit allowlist if the client
needs to monitor submitted work. Poll status until completion, then retrieve and
inspect the scientific result. Actual scheduler submission is site-specific and
is not exercised by the local cookbook smoke tests.
