# ChemGraph Integration

[ChemGraph](https://github.com/argonne-lcf/ChemGraph) is an agentic framework
(LangGraph + ASE) from Argonne LCF that drives molecular-simulation workflows
with LLMs. mofforge exposes its MOF capabilities to ChemGraph as **MCP (Model
Context Protocol) tools**, so a ChemGraph agent can screen the CoRE MOF
database, fetch structures, place adsorbates, build MOFs, validate, and render —
then hand the resulting CIFs to ChemGraph's own ASE / gRASPA tools.

There are two ways to wire mofforge into ChemGraph.

---

## Path 1 — External stdio MCP server (no ChemGraph dependency)

mofforge ships a standalone [FastMCP](https://github.com/modelcontextprotocol)
server (`mofforge-mcp`) that any MCP client — including ChemGraph and OpenCode —
can launch as a child process over stdio. Tools run **inline** in the mofforge
process.

### Install

```bash
pip install "mofforge[mcp,chem]"
playwright install chromium     # only if you use mofforge_render
```

### Configure data paths

The database and structure tools resolve their data via environment variables
(or `mofforge.toml` / `set_paths(...)`):

| Variable | Purpose |
|----------|---------|
| `MOFFORGE_COREMOF_DATA_PATH` | CoRE MOF metadata CSV (properties) |
| `MOFFORGE_COREMOF_STRUCTURES_PATH` | Directory of CoRE MOF CIF files |
| `MOFFORGE_CSD_DATA_PATH` | CSD MOF subset export (TSV) |
| `MOFFORGE_LOG_DIR` | Base dir for relative output paths (CIFs, PNGs) |

> The CoRE MOF metadata CSV contains **properties only**. The CIF structure
> files are distributed separately (Zenodo: https://zenodo.org/records/14510695).
> Point `MOFFORGE_COREMOF_STRUCTURES_PATH` at the unpacked CIF directory so
> `mofforge_get_structure` can resolve files.

### Register with an MCP client

stdio config (OpenCode / Claude / `langchain-mcp-adapters` style):

```json
{
  "mcpServers": {
    "mofforge": {
      "command": "mofforge-mcp",
      "args": ["--transport", "stdio"],
      "env": {
        "MOFFORGE_COREMOF_DATA_PATH": "/data/coremof/CR_data_CSD_modified_20250227.csv",
        "MOFFORGE_COREMOF_STRUCTURES_PATH": "/data/coremof/structures",
        "MOFFORGE_CSD_DATA_PATH": "/data/csd/MOF_subset.tab",
        "MOFFORGE_LOG_DIR": "/scratch/mofforge_out"
      }
    }
  }
}
```

For an HTTP transport instead: `mofforge-mcp --transport streamable-http --port 9010`.

### Select exposed tools

The stock server exposes all 23 tools by default. Use a strict allowlist when an
agent should receive only a smaller catalog:

```bash
mofforge-mcp --tools mofforge_search_coremof,mofforge_get_structure,mofforge_validate
```

The same value can be supplied through `MOFFORGE_MCP_TOOLS`, which is convenient
in an MCP client's `env` configuration. Use `--available-tools-only` to omit
rendering, functionalization, or construction tools when their `vis`, `chem`, or
`build` dependencies are not installed. An explicitly requested unavailable
tool is treated as a startup configuration error.

In ChemGraph, load these as LangChain tools with `langchain-mcp-adapters`'
`MultiServerMCPClient` and bind them to your agent the same way ChemGraph loads
its own MCP servers (see ChemGraph `docs/mcp_servers.md`).

---

## Path 2 — ChemGraph CGFastMCP server (HPC backend execution)

ChemGraph's `dev-globus-hpc` branch defines `CGFastMCP` (a `FastMCP` subclass)
that submits tool calls to an execution backend (Parsl on Polaris/Aurora) with
async job tracking and ensemble fan-out. mofforge provides a `CGFastMCP` entry
point, `mofforge-mcp-chemgraph`, for this path.

### Install

```bash
pip install "mofforge[mcp,chem]"
pip install "chemgraph @ git+https://github.com/argonne-lcf/ChemGraph.git@dev-globus-hpc"
```

### Tool placement

| Tool | Execution |
|------|-----------|
| `mofforge_search_coremof`, `mofforge_screen_coremof`, `mofforge_search_csd`, `mofforge_lookup_mof`, `mofforge_get_structure`, `mofforge_list_adsorbates`, `mofforge_validate` | **inline** (fast, IO-bound) |
| `mofforge_build` (TOBACCO/pormake), `mofforge_render` (Playwright) | **backend task** (CPU/GPU-heavy, job-tracked) |
| `mofforge_screen_and_place` | **ensemble fan-out** — one backend task per screened MOF |

`CGFastMCP` also registers job-management tools (`check_job_status`,
`get_job_results`, `list_jobs`, `cancel_job`, `check_endpoint_status`) when the
backend is initialized. These names participate in the same strict
`MOFFORGE_MCP_TOOLS` allowlist as the mofforge-prefixed tools.

### Launch

```bash
mofforge-mcp-chemgraph                      # stdio (default), port 9011 for HTTP
MOFFORGE_MCP_TOOLS=mofforge_validate,mofforge_screen_coremof \
  mofforge-mcp-chemgraph                    # filtered catalog
```

`main()` wraps `init_backend()` / `run_mcp_server()` / `shutdown_backend()` just
like ChemGraph's own `*_mcp_hpc.py` servers. The job tracker persists to
`$MOFFORGE_MCP_JOBS_FILE` (default `~/.mofforge_mcp_jobs.json`).
Set `MOFFORGE_MCP_AVAILABLE_ONLY=1` to omit tools backed by unavailable optional
dependencies in this entry point.

> When tools run on compute nodes, the configured data paths
> (`MOFFORGE_COREMOF_*`) must be reachable from those nodes, or staged via a
> `set_pre_submit_hook`.

---

## Example agent workflow (adsorption screening)

The headline flow ChemGraph can drive end-to-end:

1. **Screen** — `mofforge_screen_coremof(pld_min=3.8, lcd_min=6.0, water_stability_min=0.7, metal="Cu", has_oms=True)`
   → shortlist of `coreid`s with pore size and stability metadata.
2. **Fetch** — `mofforge_get_structure(coreid)` → a local CIF path.
3. **Load** — hand the CIF to ChemGraph's `file_to_atomsdata` / `run_ase`.
4. **Place** — `mofforge_place_adsorbate(cif_path, adsorbate="CO2", n_adsorbates=10)`
   → a loaded structure for a gRASPA uptake run.
5. **Simulate & rank** — ChemGraph's `run_graspa` + `rank_mofs_performance`.

On HPC (Path 2), steps 1–4 collapse into a single fan-out call:

```text
mofforge_screen_and_place(
    pld_min=3.8, lcd_min=6.0, water_stability_min=0.7, metal="Cu",
    adsorbate="CO2", n_adsorbates=10, limit=50
)
```

which screens, then runs `get_structure` + `place_adsorbate` for every candidate
as independent backend tasks, returning per-MOF results ready for simulation.

---

## Tool reference

All tools return a JSON object with a `success` boolean; on failure the object
carries an `error` string (tools never raise out of the server).

| Tool | Summary |
|------|---------|
| `mofforge_search_coremof` | Search CoRE MOF (coreid/refcode/name/metal/topology) |
| `mofforge_screen_coremof` | Screen CoRE MOF by pore size, density, ASA, void fraction, stability, metal, topology, OMS |
| `mofforge_search_csd` | Search the CSD lookup table |
| `mofforge_lookup_mof` | Name → CSD → CoRE MOF bridge |
| `mofforge_get_structure` | Resolve coreid/refcode → local CIF path |
| `mofforge_place_adsorbate` | Place adsorbate molecule(s) into a MOF (stock server) |
| `mofforge_list_adsorbates` | List built-in adsorbates |
| `mofforge_validate` | Validate clashes / bonds / coordination |
| `mofforge_build` | Build a MOF from topology + building blocks |
| `mofforge_render` | Render a structure to PNG |
| `mofforge_search`, `mofforge_replace`, `mofforge_remove`, `mofforge_desolvate`, `mofforge_smarts_search`, `mofforge_list_topologies`, `mofforge_list_building_blocks` | Substructure search / fragment edit / construction helpers (stock server) |
| `mofforge_screen_and_place` | Ensemble fan-out: screen → get_structure → place (CGFastMCP only) |
