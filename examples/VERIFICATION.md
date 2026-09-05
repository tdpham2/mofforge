# Cookbook verification

Checked on **2026-09-05**, using the source checkout and an existing macOS ARM64
Python **3.12.2** environment. The source project declares 0.2.0; this environment's
older editable-install metadata reports 0.1.0. Reports retain that installed
metadata rather than silently relabeling it. CI creates fresh locked environments.

## Executed workflows

All **36 lesson recipes** were executed against their documented default inputs,
including both real construction backends, chemistry, real stdio MCP calls,
both complete workflows, and synchronous/asynchronous PNG rendering. TOBACCO
used an explicit path to prepared external runtime data.

Representative scientific checks included:

- IRMOF-1: 424 atoms, 512 inferred bonds, 24 phenylene locations.
- Missing-hydrogen repair: 328 → 424 atoms, adding 96 H.
- Disorder repair and acetylene removal: 200 → 136 → 104 atoms.
- Two CO2 guests: 424 → 430 atoms; desolvation restores 424.
- Pormake N108/E1/pcu: 123 atoms; one nitro substitution gives 125.
- TOBACCO dmc fixture: nine generated CIFs, each with a manifest.
- Batch descriptors verify completely; plain exports verify integrity only;
  an altered copy fails file integrity.
- Serial and two-worker batches produce matching run IDs and geometry.

These are reference observations, not promises of identical floating-point
coordinates across dependency versions. Some fixtures intentionally exercise
warnings, coordination heuristics, missing structures, or incomplete verification.
Do not convert those observations into claims of optimized chemistry.

Automated results:

- Full repository suite: **615 passed, 17 skipped**.
- Cookbook suite with network connections and optional imports blocked for core
  lessons: **74 passed, 8 skipped**.
- Additional Pormake, TOBACCO, MCP-builder, build-and-modify, adapter, and PNG
  lesson tests passed with their prerequisites enabled.
- Ruff lint and cookbook formatting checks passed.

## Optional environments

| Component | Version / setup exercised |
| --- | --- |
| pymatgen | 2026.5.4 |
| NumPy / SciPy / NetworkX | 2.4.4 / 1.17.1 / 3.6.1 |
| RDKit | 2026.3.3 |
| Pormake | 0.2.3 |
| tobacco3 | 3.1.0 with the existing verified external data cache |
| Stock MCP | 1.27.0 |
| Playwright | 1.58.0 with installed Chromium |
| LangChain adapter | langchain-mcp-adapters 0.3.2, langchain-core 1.6.2, MCP 1.29.1 in an isolated temporary dependency directory |

Chromium required permission to run outside the development sandbox. Both PNG
APIs then produced 640 × 480 images. The adapter lesson successfully discovered
and invoked the stock validation tool without an LLM.

Actual CGFastMCP scheduler submission and licensed CSD/full CoRE datasets were
**not exercised**. The HPC guide requires a site-configured backend; the optional
real-database test requires local dataset paths. Local adapter validation does
not imply that HPC scheduling was tested.

## Reproduce the checks

```bash
python -m ruff check src tests examples
python -m pytest tests/test_examples.py
MOFFORGE_RUN_BUILD_TESTS=1 python -m pytest tests/test_examples.py -m integration -k 'not external_database and not png_rendering'
MOFFORGE_RUN_RENDER_TESTS=1 python -m pytest tests/test_examples.py -k png_rendering
python examples/pipeline/batch_processing.py --compare-parallel
```

Normal CI executes core lessons and installed chemistry/MCP lessons. The builder
job additionally exercises real builder and MCP-builder lessons. Browser
execution is a manually dispatched CI job because it requires Chromium and the
3Dmol.js CDN. Missing packages are reported as skips; installed-but-broken
integrations fail when explicitly enabled.
