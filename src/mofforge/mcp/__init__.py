"""MCP (Model Context Protocol) server for mofforge.

Exposes mofforge's crystal structure manipulation capabilities as
tools for AI agents (e.g. ChemGraph, Claude Desktop, OpenCode).

Run the server::

    # stdio mode (default)
    mofforge-mcp

    # HTTP mode
    mofforge-mcp --transport streamable-http --port 9010

    # Or directly:
    python -m mofforge.mcp.server
"""
