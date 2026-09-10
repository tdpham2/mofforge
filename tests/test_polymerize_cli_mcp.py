"""Shared JSON interface contracts, including partial results and dependency gating."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from mofforge.cli import main
from mofforge.polymerize import POPResult


@pytest.mark.parametrize("capability_present", [False, True])
def test_pop_capability_checks_rdkit_explicitly(monkeypatch, capability_present):
    from mofforge.mcp import tool_selection

    monkeypatch.setattr(
        tool_selection.importlib.util,
        "find_spec",
        lambda name: object() if name == "rdkit" and capability_present else None,
    )
    assert tool_selection.capability_available("pop") is capability_present


def test_tools_registered_with_operation_specific_gates():
    pytest.importorskip("mcp")
    from mofforge.mcp import server

    caps = {tool.name: tool.capability for tool in server._TOOL_REGISTRY}
    assert caps["mofforge_pack"] == "pop"
    assert caps["mofforge_polymerize"] == "pop"
    assert caps["mofforge_list_reactions"] is None


def test_site_catalog_advertises_no_recipes():
    from mofforge.mcp._impl import list_reactions_impl

    result = list_reactions_impl()
    assert result["success"] and result["reactions"] == []
    assert result["site_types"]


@pytest.mark.parametrize("command", ["pack", "polymerize"])
def test_cli_requires_shared_config(command):
    result = CliRunner().invoke(main, [command])
    assert result.exit_code != 0 and "--config" in result.output


@pytest.mark.parametrize(
    "option,value,hint",
    [
        ("--target-density", "0.8", "initial_packing_density"),
        ("--forcefield", "gaff2", "MatKit"),
        ("--md-settings", "{}", "MatKit"),
        ("--n-monomers", "20", "count"),
    ],
)
def test_removed_cli_options_explain_migration(option, value, hint):
    result = CliRunner().invoke(main, ["polymerize", option, value])
    assert result.exit_code != 0
    assert hint in result.output


@pytest.mark.parametrize("status", ["completed", "partial", "failed"])
def test_cli_preserves_status_outputs_and_counts(tmp_path, monkeypatch, status):
    import mofforge.polymerize

    config = {
        "components": [{"source": "CC", "count": 3}],
        "packing": {"box_lengths": [10, 10, 10]},
    }
    path = tmp_path / "request.json"
    path.write_text(json.dumps(config))
    seen = []

    def run(request, **kwargs):
        seen.append((request, kwargs))
        return POPResult(
            status == "completed",
            status=status,
            operation="connect",
            output_paths=[Path("/saved/state.json")],
            state_path=Path("/saved/state.json"),
            metadata={"conversion": 0.5},
            errors=[] if status == "completed" else ["budget"],
        )

    monkeypatch.setattr(mofforge.polymerize, "run_config", run)
    result = CliRunner().invoke(main, ["polymerize", "--config", str(path), "--as-json"])
    assert result.exit_code == (0 if status == "completed" else 1)
    data = json.loads(result.output)
    assert data["status"] == status and data["metadata"]["conversion"] == 0.5
    assert data["state_path"] == "/saved/state.json"
    assert data["output_paths"] == ["/saved/state.json"]
    assert seen[0][0] == config and seen[0][1]["operation"] == "connect"


def test_mcp_partial_result_has_state_and_conversion(monkeypatch, tmp_path):
    import mofforge.polymerize
    from mofforge.mcp._impl import polymerize_impl

    monkeypatch.setattr(
        mofforge.polymerize,
        "run_config",
        lambda *a, **kw: POPResult(
            False,
            status="partial",
            operation="connect",
            state_path=tmp_path / "state.json",
            metadata={"conversion": 0.25},
            errors=["budget"],
            output_paths=[tmp_path / "state.json"],
        ),
    )
    result = polymerize_impl({}, output_dir=str(tmp_path))
    assert result["status"] == "partial" and not result["success"]
    assert result["metadata"]["conversion"] == 0.25
    assert result["state_path"] and result["output_paths"]


def test_mcp_rejects_removed_options():
    from mofforge.mcp._impl import polymerize_impl

    result = polymerize_impl({}, forcefield="gaff2")
    assert not result["success"] and "MatKit" in result["errors"][0]


def test_doctor_without_packmol(monkeypatch):
    monkeypatch.setenv("MOFFORGE_PACKMOL_BIN", "/missing/packmol")
    monkeypatch.delenv("MOFFORGE_LAMMPS_BIN", raising=False)
    result = CliRunner().invoke(main, ["pop-doctor", "--as-json"])
    assert result.exit_code == 0
    report = json.loads(result.output)
    assert set(report) == {"rdkit", "packmol"}
    assert not report["packmol"]["available"]
