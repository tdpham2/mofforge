"""ChemGraph screening contracts without the optional HPC dependency."""

import inspect
import sys
from types import ModuleType
from unittest.mock import Mock

import pytest

from mofforge.mcp import _impl, chemgraph_server
from tests.test_coremof import _HEADER, _ROWS


@pytest.fixture
def screening_tools(monkeypatch):
    """Capture the real wrappers, substituting only ChemGraph registration."""
    class ToolRecorder:
        def __init__(self, **kwargs):
            self.functions = {}

        def add_tool(self, fn, name=None, **kwargs):
            self.functions[name or fn.__name__] = fn

        def tool(self, **kwargs):
            def register(fn):
                self.add_tool(fn, **kwargs)
                return fn
            return register

        schema_fanout_tool = tool

    for name in ("chemgraph", "chemgraph.mcp", "chemgraph.mcp.cg_fastmcp"):
        monkeypatch.setitem(sys.modules, name, ModuleType(name))
    monkeypatch.setattr(sys.modules["chemgraph.mcp.cg_fastmcp"], "CGFastMCP", ToolRecorder,
                        raising=False)
    server = chemgraph_server.build_server({
        "mofforge_screen_coremof", "mofforge_screen_and_place",
    })
    return server.functions


def _parameters(function):
    return {
        name: (parameter.annotation, parameter.default)
        for name, parameter in inspect.signature(function).parameters.items()
    }


def test_screening_signature_matches_shared_implementation(screening_tools):
    assert _parameters(screening_tools["mofforge_screen_coremof"]) == _parameters(
        _impl.screen_coremof_impl
    )


def test_screening_signature_matches_stock_server(screening_tools):
    pytest.importorskip("mcp")
    from mofforge.mcp.server import mofforge_screen_coremof

    assert _parameters(screening_tools["mofforge_screen_coremof"]) == _parameters(
        mofforge_screen_coremof
    )


@pytest.mark.parametrize("filters", [
    {},
    {"asa_min": 100.0, "asa_max": 2000.0, "void_fraction_max": 0.75, "extension": "ASR"},
    {"asa_min": 0.0, "asa_max": 0.0, "void_fraction_max": 0.0, "extension": None},
])
def test_screening_forwards_filters_and_defaults(screening_tools, monkeypatch, filters):
    defaults = {
        name: parameter.default
        for name, parameter in inspect.signature(_impl.screen_coremof_impl).parameters.items()
    }
    response = {"success": True, "records": []}
    screen = Mock(return_value=response)
    monkeypatch.setattr(_impl, "screen_coremof_impl", screen)

    result = screening_tools["mofforge_screen_coremof"](**filters)

    screen.assert_called_once_with(**{**defaults, **filters})
    assert result is response


@pytest.mark.parametrize(("filters", "expected_count"), [
    ({"asa_min": 1000.0}, 1),
    ({"asa_max": 100.0}, 5),
    ({"void_fraction_max": 0.5}, 5),
    ({"extension": "with ion"}, 1),
])
def test_screening_applies_filters_to_database(screening_tools, tmp_path, filters, expected_count):
    csv_file = tmp_path / "coremof.csv"
    csv_file.write_text(_HEADER + "\n" + "\n".join(_ROWS) + "\n", encoding="utf-8")

    result = screening_tools["mofforge_screen_coremof"](data_path=str(csv_file), **filters)

    assert result["success"]
    assert result["n_matches"] == expected_count


def test_screen_and_place_forwards_screening_filters(screening_tools, monkeypatch):
    filters = {"asa_min": 100.0, "asa_max": 2000.0, "void_fraction_max": 0.75, "extension": "ASR"}
    screen = Mock(return_value={"success": True, "records": [{"coreid": "candidate"}]})
    monkeypatch.setattr(_impl, "screen_coremof_impl", screen)
    params = {**filters, "adsorbate": "N2", "limit": 3, "data_path": "metadata.csv"}

    jobs = screening_tools["mofforge_screen_and_place"](params)

    screen.assert_called_once_with(limit=3, data_path="metadata.csv", **filters)
    assert len(jobs) == 1
    assert jobs[0]["coreid"] == "candidate"
    assert jobs[0]["adsorbate"] == "N2"
    assert params == {**filters, "adsorbate": "N2", "limit": 3, "data_path": "metadata.csv"}
