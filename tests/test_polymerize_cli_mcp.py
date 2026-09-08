"""M2 tests: the polymerize MCP tools and CLI command.

Cover tool registration, capability gating, the ``_impl`` dict shape, and the
``polymerize`` / ``pop-doctor`` CLI commands.  None require the external
binaries; the polymerize path is exercised through its graceful missing-binary
failure so the plumbing (argument passing, JSON shape, exit codes) is verified.
"""

from __future__ import annotations

import json

import pytest
from click.testing import CliRunner

from mofforge.cli import main

# ---------------------------------------------------------------------------
# MCP tool registration + capability gating
# ---------------------------------------------------------------------------


def test_polymerize_tools_registered():
    from mofforge.mcp import server

    names = {r.name for r in server._TOOL_REGISTRY}
    assert "mofforge_polymerize" in names
    assert "mofforge_list_reactions" in names


def test_polymerize_tools_gated_on_pop_capability():
    from mofforge.mcp import server

    caps = {r.name: r.capability for r in server._TOOL_REGISTRY}
    assert caps["mofforge_polymerize"] == "pop"
    assert caps["mofforge_list_reactions"] == "pop"


def test_pop_capability_maps_to_pysimm():
    from mofforge.mcp.tool_selection import capability_available

    # pysimm is not installed in CI -> capability unavailable (no crash).
    assert capability_available("pop") is False


def test_pop_tools_excluded_when_unavailable():
    from mofforge.mcp.tool_selection import select_tool_names

    selected = select_tool_names(
        ["mofforge_polymerize", "mofforge_validate"],
        {"mofforge_polymerize": "pop", "mofforge_validate": None},
        available_only=True,
        availability=lambda cap: cap != "pop",
    )
    assert "mofforge_polymerize" not in selected
    assert "mofforge_validate" in selected


# ---------------------------------------------------------------------------
# _impl dict shape
# ---------------------------------------------------------------------------


def test_list_reactions_impl_shape():
    from mofforge.mcp._impl import list_reactions_impl

    result = list_reactions_impl()
    assert result["success"] is True
    assert any("imine" in r["reaction"] for r in result["reactions"])
    assert any(s["site_type"] == "amine" for s in result["site_types"])


def test_polymerize_impl_no_monomers():
    from mofforge.mcp._impl import polymerize_impl

    result = polymerize_impl([])
    assert result["success"] is False
    assert "monomer" in result["error"].lower()


def test_polymerize_impl_functionality_mismatch():
    from mofforge.mcp._impl import polymerize_impl

    result = polymerize_impl(["NCCN", "C=O"], functionality=[2])
    assert result["success"] is False
    assert "functionality" in result["error"].lower()


def test_polymerize_impl_missing_binary(tmp_path, monkeypatch):
    pytest.importorskip("rdkit")
    monkeypatch.setenv("MOFFORGE_PACKMOL_BIN", "/nonexistent/packmol")
    from mofforge.mcp._impl import polymerize_impl

    result = polymerize_impl(
        ["NCCN", "O=Cc1ccc(C=O)cc1"],
        output_dir=str(tmp_path),
        random_seed=42,
    )
    assert result["success"] is False
    # Reported as a clean errors list from the backend, not an exception.
    assert "errors" in result or "error" in result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def test_cli_polymerize_requires_monomer():
    runner = CliRunner()
    result = runner.invoke(main, ["polymerize"])
    assert result.exit_code != 0
    assert "monomer" in result.output.lower()


def test_cli_polymerize_functionality_mismatch():
    runner = CliRunner()
    result = runner.invoke(
        main, ["polymerize", "-m", "NCCN", "-m", "C=O", "--functionality", "2"]
    )
    assert result.exit_code == 1
    assert "once per" in result.output


def test_cli_polymerize_missing_binary(tmp_path, monkeypatch):
    pytest.importorskip("rdkit")
    monkeypatch.setenv("MOFFORGE_PACKMOL_BIN", "/nonexistent/packmol")
    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "polymerize",
            "-m",
            "NCCN",
            "-m",
            "O=Cc1ccc(C=O)cc1",
            "-o",
            str(tmp_path),
            "--random-seed",
            "42",
        ],
    )
    assert result.exit_code == 1
    assert "failed" in result.output.lower()


def test_cli_pop_doctor_json(monkeypatch):
    monkeypatch.setenv("MOFFORGE_PACKMOL_BIN", "/nonexistent/packmol")
    monkeypatch.setenv("MOFFORGE_LAMMPS_BIN", "/nonexistent/lmp")
    runner = CliRunner()
    result = runner.invoke(main, ["pop-doctor", "--as-json"])
    assert result.exit_code == 0
    report = json.loads(result.output)
    assert set(report) == {"pysimm", "packmol", "lammps"}
