"""Shared JSON interface contracts, including partial results and dependency gating."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from mofforge.cli import main
from mofforge.polymerize import POPResult


def test_site_catalog_advertises_no_recipes():
    from mofforge.mcp._impl import list_reactions_impl

    result = list_reactions_impl()
    assert result["success"] and result["reactions"] == []
    assert result["site_types"]


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


def test_doctor_without_packmol(monkeypatch):
    monkeypatch.setenv("MOFFORGE_PACKMOL_BIN", "/missing/packmol")
    result = CliRunner().invoke(main, ["pop-doctor", "--as-json"])
    assert result.exit_code == 0
    report = json.loads(result.output)
    assert set(report) == {"rdkit", "packmol"}
    assert not report["packmol"]["available"]


@pytest.mark.parametrize(
    "config,match",
    [
        ({"connection": {"output_dir": "x"}}, "Unknown connection options"),
        ({"components": [], "packing": {"rules": []}, "connection": {}}, "out of packing"),
        (
            {"components": [], "packing": {}, "connection": {"random_seed": 1}},
            "random_seed under packing",
        ),
    ],
)
def test_run_config_keeps_packing_and_connection_options_apart(config, match):
    from mofforge.polymerize import run_config

    with pytest.raises(ValueError, match=match):
        run_config(config, operation="connect")
