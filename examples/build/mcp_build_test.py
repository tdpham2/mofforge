#!/usr/bin/env python3
"""Build through a real stdio MCP client using a selected backend.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import example_parser, output_dir, require, tool_payload, write_report


async def run(args, out):
    import os

    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    env = {**os.environ, "MOFFORGE_LOG_DIR": str(out)}
    if args.tobacco_data:
        env["MOFFORGE_TOBACCO_DATA"] = str(args.tobacco_data.resolve())
    if args.backend == "tobacco":
        from mofforge.build.config import BuildConfig

        data = BuildConfig.load(
            **({"tobacco_data_dir": args.tobacco_data} if args.tobacco_data else {})
        ).resolve_tobacco_data_dir()
        fixture = data / "tests" / "fixtures"
        topology = args.topology or "dmc"
        nodes = args.node or [
            str(fixture / "nodes" / f"{n}.cif") for n in ["4c_1Zn_Ch", "triazole"]
        ]
        edges = args.edge or [
            str(fixture / "edges" / f"{n}.cif") for n in ["ntn_edge", "oxalic_edge", "squOxa_ch"]
        ]
    else:
        topology, nodes, edges = args.topology or "pcu", args.node or ["N108"], args.edge or ["E1"]
    params = StdioServerParameters(
        command=args.server_python,
        args=["-m", "mofforge.mcp.server", "--tools", "mofforge_build"],
        env=env,
        cwd=str(out),
    )
    async with stdio_client(params) as (read, write), ClientSession(read, write) as session:
        await session.initialize()
        result = tool_payload(
            await session.call_tool(
                "mofforge_build",
                {
                    "backend": args.backend,
                    "topology": topology,
                    "node_files": nodes,
                    "edge_files": edges,
                    "output_dir": str(out / "structures"),
                },
            )
        )
    write_report(out, result=result)
    if not result.get("success"):
        raise SystemExit("MCP construction failed; see summary.json.")


def main():
    parser = example_parser("mcp_build_test", __doc__)
    parser.add_argument("--server-python", default=sys.executable)
    parser.add_argument("--tobacco-data", type=Path)
    parser.add_argument("--output-root", type=Path, help="Legacy alias overriding --output-dir.")
    parser.add_argument("--backend", choices=["pormake", "tobacco"], default="pormake")
    parser.add_argument("--topology")
    parser.add_argument("--node", action="append")
    parser.add_argument("--edge", action="append")
    args = parser.parse_args()
    require("mcp", "mcp")
    if args.output_root:
        args.output_dir = args.output_root
    import asyncio

    asyncio.run(run(args, output_dir(args)))


if __name__ == "__main__":
    main()
