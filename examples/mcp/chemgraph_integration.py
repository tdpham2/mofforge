#!/usr/bin/env python3
"""Load MOF tools through the LangChain adapter used by ChemGraph.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, example_parser, output_dir, require, write_report


async def run(args, out):
    import json
    import os

    from langchain_mcp_adapters.client import MultiServerMCPClient

    # This external-server path needs no LLM or HPC backend.
    connection = {
        "command": args.server_python,
        "args": ["-m", "mofforge.mcp.server", "--tools", "mofforge_validate"],
        "transport": "stdio",
        "env": {**os.environ, "MOFFORGE_LOG_DIR": str(out)},
        "cwd": str(out),
    }
    client = MultiServerMCPClient({"mofforge": connection})
    tools = await client.get_tools()
    validate = next(tool for tool in tools if tool.name == "mofforge_validate")
    raw = await validate.ainvoke({"cif_path": str(CRYSTALS / "IRMOF-1.cif")})
    # Adapter calls return content blocks or text, rather than CallToolResult.
    payload = json.loads(
        raw
        if isinstance(raw, str)
        else next(block["text"] for block in raw if block.get("type") == "text")
    )
    write_report(
        out,
        tools=[tool.name for tool in tools],
        validation=payload,
        integration="External MCP tools ready to bind to a ChemGraph/LangGraph agent.",
    )
    if not payload.get("success"):
        raise SystemExit("The adapter's validation call failed.")


def main():
    parser = example_parser("chemgraph_integration", __doc__)
    parser.add_argument("--server-python", default=sys.executable)
    args = parser.parse_args()
    require("mcp", "mcp")
    import asyncio
    import importlib.util

    if importlib.util.find_spec("langchain_mcp_adapters") is None:
        raise SystemExit("Install this lesson's adapter: pip install langchain-mcp-adapters")
    asyncio.run(run(args, output_dir(args)))


if __name__ == "__main__":
    main()
