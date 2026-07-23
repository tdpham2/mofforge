"""Tests for startup-time MCP tool selection."""

from __future__ import annotations

import asyncio

import pytest
from mcp.server.fastmcp.exceptions import ToolError

from mofforge.mcp import tool_selection
from mofforge.mcp.server import build_server
from mofforge.mcp.tool_selection import ToolSelectionError, parse_tool_list


def _tool_names(server) -> set[str]:
    return {tool.name for tool in asyncio.run(server.list_tools())}


def test_default_stock_server_keeps_full_catalog():
    assert len(_tool_names(build_server())) == 23


def test_stock_server_registers_only_allowlisted_tools():
    server = build_server({"mofforge_search", "mofforge_validate"})
    assert _tool_names(server) == {"mofforge_search", "mofforge_validate"}
    assert "functionalization" not in server.instructions


def test_hidden_stock_tool_cannot_be_called():
    server = build_server({"mofforge_validate"})
    with pytest.raises(ToolError, match="Unknown tool: mofforge_render"):
        asyncio.run(server.call_tool("mofforge_render", {}))


def test_unknown_explicit_tool_is_rejected():
    with pytest.raises(ToolSelectionError, match="does_not_exist"):
        build_server({"does_not_exist"})


def test_explicit_unavailable_tool_is_rejected(monkeypatch):
    monkeypatch.setattr(tool_selection, "capability_available", lambda capability: False)
    with pytest.raises(ToolSelectionError, match=r"mofforge_render \(vis\)"):
        build_server({"mofforge_render"})


def test_available_only_hides_unavailable_capabilities(monkeypatch):
    monkeypatch.setattr(
        tool_selection,
        "capability_available",
        lambda capability: capability == "vis",
    )
    names = _tool_names(build_server(available_only=True))
    assert "mofforge_render" in names
    assert "mofforge_build" not in names
    assert "mofforge_find_sites" not in names
    assert "mofforge_list_functional_groups" in names
    assert "mofforge_validate" in names


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        (None, None),
        ("", set()),
        ("mofforge_search, mofforge_validate", {"mofforge_search", "mofforge_validate"}),
    ],
)
def test_parse_tool_list(raw, expected):
    assert parse_tool_list(raw) == expected


def test_stock_cli_forwards_tool_selection(monkeypatch):
    from mofforge.mcp import server as server_module

    captured = {}

    class FakeServer:
        def run(self, *, transport):
            captured["transport"] = transport

    def fake_build_server(enabled_tools, available_only, *, port):
        captured.update(
            enabled_tools=enabled_tools,
            available_only=available_only,
            port=port,
        )
        return FakeServer()

    monkeypatch.setattr(server_module, "build_server", fake_build_server)
    monkeypatch.setattr(
        "sys.argv",
        [
            "mofforge-mcp",
            "--transport",
            "streamable-http",
            "--port",
            "9123",
            "--tools",
            "mofforge_search,mofforge_validate",
            "--available-tools-only",
        ],
    )

    server_module.main()

    assert captured == {
        "enabled_tools": {"mofforge_search", "mofforge_validate"},
        "available_only": True,
        "port": 9123,
        "transport": "streamable-http",
    }


def test_chemgraph_server_registers_only_allowlisted_tools():
    cs = pytest.importorskip("mofforge.mcp.chemgraph_server")
    try:
        server = cs.build_server({"mofforge_validate", "mofforge_screen_coremof"})
    except ImportError:
        pytest.skip("ChemGraph (CGFastMCP) not installed")
    assert _tool_names(server) == {"mofforge_validate", "mofforge_screen_coremof"}
    server.init_backend()
    assert _tool_names(server) == {"mofforge_validate", "mofforge_screen_coremof"}


def test_chemgraph_job_tools_can_be_allowlisted():
    cs = pytest.importorskip("mofforge.mcp.chemgraph_server")
    try:
        server = cs.build_server({"check_job_status"})
    except ImportError:
        pytest.skip("ChemGraph (CGFastMCP) not installed")
    assert _tool_names(server) == set()
    server.init_backend()
    assert _tool_names(server) == {"check_job_status"}
