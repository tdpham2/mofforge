#!/usr/bin/env python3
"""Discover selected MCP tools and execute search, modification, and validation.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import (
    CRYSTALS,
    MOIETIES,
    example_parser,
    output_dir,
    require,
    tool_payload,
    write_report,
)


async def run(args, out):
    import os

    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    names = ["mofforge_search", "mofforge_replace", "mofforge_validate"]
    params = StdioServerParameters(
        command=args.server_python,
        args=["-m", "mofforge.mcp.server", "--tools", ",".join(names)],
        env={**os.environ, "MOFFORGE_LOG_DIR": str(out)},
        cwd=str(out),
    )
    async with stdio_client(params) as (read, write), ClientSession(read, write) as session:
        await session.initialize()
        catalog = await session.list_tools()
        advertised = sorted(tool.name for tool in catalog.tools)
        if advertised != sorted(names):
            raise RuntimeError("The server advertised tools outside the selected catalog.")
        search = tool_payload(
            await session.call_tool(
                "mofforge_search",
                {
                    "parent_cif": str(CRYSTALS / "IRMOF-1.cif"),
                    "query_xyz": str(MOIETIES / "p-phenylene.xyz"),
                },
            )
        )
        child_path = out / "modified.cif"
        replacement = tool_payload(
            await session.call_tool(
                "mofforge_replace",
                {
                    "parent_cif": str(CRYSTALS / "IRMOF-1.cif"),
                    "query_xyz": str(MOIETIES / "2-!-p-phenylene.xyz"),
                    "replacement_xyz": str(MOIETIES / "2-nitro-p-phenylene.xyz"),
                    "output_cif": str(child_path),
                    "nb_loc": 2,
                    "random_seed": args.seed,
                    "validate": True,
                },
            )
        )
        validation = tool_payload(
            await session.call_tool("mofforge_validate", {"cif_path": str(child_path)})
        )
    write_report(
        out, advertised=advertised, search=search, replacement=replacement, validation=validation
    )
    if any(not result.get("success") for result in [search, replacement, validation]):
        raise SystemExit("An MCP operation failed; see summary.json.")


def main():
    parser = example_parser("local_client", __doc__, seed=True)
    parser.add_argument("--server-python", default=sys.executable)
    args = parser.parse_args()
    require("mcp", "mcp")
    import asyncio

    asyncio.run(run(args, output_dir(args)))


if __name__ == "__main__":
    main()
