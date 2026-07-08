#!/usr/bin/env python3
"""End-to-end test: build several MOFs through the mofforge MCP server.

This is a *true* MCP client: it spawns the ``mofforge-mcp`` server as a stdio
subprocess and calls the ``mofforge_build`` tool over the MCP protocol, exactly
as an AI agent would.  For each spec it asserts the build succeeded and produced
a non-empty structure, then prints a summary table.

Environment prerequisites
-------------------------
The server process needs BOTH ``mofforge`` and ``tobacco3`` importable, plus the
TOBACCO data directories (the git-installed tobacco3 wheel does NOT bundle them).
A combined venv inside this repo satisfies that::

    cd /Users/tpham2/work/soft/mofforge
    python3 -m venv .venv
    .venv/bin/pip install -e ".[all]"          # mofforge + tobacco3 + mcp client SDK
    # tobacco3's data dirs live in the standalone checkout:
    export MOFFORGE_TOBACCO_DATA=/Users/tpham2/work/soft/tobacco_3.0
    .venv/bin/python -c "import mofforge, tobacco3; print('both OK')"

Run
---
    ./.venv/bin/python examples/build/mcp_build_test.py

Override the interpreter used to launch the server with ``--server-python`` and
the tobacco data dir with ``--tobacco-data`` (defaults below).
"""

from __future__ import annotations

import argparse
import asyncio
import json
import os
import sys
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

# --------------------------------------------------------------------------- #
# Defaults
# --------------------------------------------------------------------------- #
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SERVER_PYTHON = REPO_ROOT / ".venv" / "bin" / "python"
DEFAULT_TOBACCO_DATA = "/Users/tpham2/work/soft/tobacco_3.0"
DEFAULT_OUTPUT_ROOT = "/tmp/mofforge_out/mcp_build_test"

# Known-good (topology, node, edge) combos.  Node/edge are bare catalog names;
# the tobacco backend resolves them against nodes_database/ and edges_database/.
BUILD_SPECS = [
    {"name": "pcu_Zn", "topology": "pcu", "node": "6c_Zn_1_Ch", "edge": "1B_4H_Ch"},
    {"name": "pcu_Cu", "topology": "pcu", "node": "6c_Cu_1_Ch", "edge": "1B_4H_Ch"},
    {"name": "dia_Cd", "topology": "dia", "node": "4c_Cd_1_Ch", "edge": "1B_4H_Ch"},
]


def _tool_result_json(result) -> dict:
    """Extract the JSON payload from an MCP tool result's text content."""
    # Prefer structured content when the server provides it.
    structured = getattr(result, "structuredContent", None)
    if isinstance(structured, dict) and structured:
        # FastMCP wraps scalar returns as {"result": <value>}; unwrap if needed.
        inner = structured.get("result", structured)
        if isinstance(inner, str):
            return json.loads(inner)
        if isinstance(inner, dict):
            return inner
    for block in result.content:
        text = getattr(block, "text", None)
        if text:
            return json.loads(text)
    raise ValueError(f"No text content in tool result: {result!r}")


async def run(server_python: str, tobacco_data: str, output_root: str) -> int:
    env = dict(os.environ)
    env["MOFFORGE_TOBACCO_DATA"] = tobacco_data
    env.setdefault("MOFFORGE_LOG_DIR", output_root)

    # Launch the server via its module entry point so we don't depend on the
    # console-script being on PATH.
    params = StdioServerParameters(
        command=server_python,
        args=["-m", "mofforge.mcp.server"],
        env=env,
    )

    rows: list[dict] = []
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()

            tools = {t.name for t in (await session.list_tools()).tools}
            if "mofforge_build" not in tools:
                print("ERROR: mofforge_build tool not exposed by server", file=sys.stderr)
                print(f"  available: {sorted(tools)}", file=sys.stderr)
                return 1
            print(f"Connected. Server exposes {len(tools)} tools.\n")

            for spec in BUILD_SPECS:
                out_dir = str(Path(output_root) / spec["name"])
                print(f"Building {spec['name']}: topology={spec['topology']} "
                      f"node={spec['node']} edge={spec['edge']} ...")
                result = await session.call_tool(
                    "mofforge_build",
                    {
                        "topology": spec["topology"],
                        "backend": "tobacco",
                        "node_files": [spec["node"]],
                        "edge_files": [spec["edge"]],
                        "output_dir": out_dir,
                    },
                )
                payload = _tool_result_json(result)
                ok = bool(payload.get("success")) and (payload.get("atoms") or 0) > 0
                outputs = payload.get("output_paths") or []
                rows.append(
                    {
                        **spec,
                        "ok": ok,
                        "atoms": payload.get("atoms"),
                        "seconds": payload.get("elapsed_seconds"),
                        "output": outputs[0] if outputs else "",
                        "error": "" if ok else json.dumps(payload)[:200],
                    }
                )
                print("  -> {}  atoms={}  {}\n".format(
                    "OK" if ok else "FAIL", payload.get("atoms"),
                    rows[-1]["output"] or rows[-1]["error"],
                ))

    # ------------------------------------------------------------------ #
    # Summary
    # ------------------------------------------------------------------ #
    print("=" * 78)
    print(f"{'name':<10} {'topo':<6} {'node':<14} {'edge':<12} {'atoms':>6}  result")
    print("-" * 78)
    for r in rows:
        print(f"{r['name']:<10} {r['topology']:<6} {r['node']:<14} {r['edge']:<12} "
              f"{str(r['atoms'] or '-'):>6}  {'OK' if r['ok'] else 'FAIL'}")
    n_ok = sum(1 for r in rows if r["ok"])
    print("=" * 78)
    print(f"{n_ok}/{len(rows)} builds succeeded")
    return 0 if n_ok == len(rows) else 1


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--server-python", default=str(DEFAULT_SERVER_PYTHON),
                        help="Interpreter used to launch the MCP server "
                             f"(default: {DEFAULT_SERVER_PYTHON})")
    parser.add_argument("--tobacco-data", default=DEFAULT_TOBACCO_DATA,
                        help="TOBACCO data directory containing "
                             "template_database/ nodes_database/ edges_database/ "
                             f"(default: {DEFAULT_TOBACCO_DATA})")
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT,
                        help=f"Directory for output CIFs (default: {DEFAULT_OUTPUT_ROOT})")
    args = parser.parse_args()

    rc = asyncio.run(run(args.server_python, args.tobacco_data, args.output_root))
    sys.exit(rc)


if __name__ == "__main__":
    main()
