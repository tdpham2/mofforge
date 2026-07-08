"""Tests for the mofforge MCP layer (_impl functions and the stock server).

These exercise the backend-agnostic tool logic in ``mofforge.mcp._impl`` and the
stock FastMCP server registration. The CGFastMCP adapter is import-guarded and
only tested when ChemGraph is installed.
"""

from __future__ import annotations

import json

import pytest

from mofforge.mcp import _impl
from tests.test_coremof import _HEADER, _ROWS

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def coremof_csv(tmp_path):
    """A synthetic CoRE MOF CSV (reuses the test_coremof fixture rows)."""
    csv_file = tmp_path / "test_coremof.csv"
    csv_file.write_text(_HEADER + "\n" + "\n".join(_ROWS) + "\n", encoding="utf-8")
    return csv_file


@pytest.fixture
def host_cif(crystal_dir):
    """A real MOF CIF to use as an adsorbate host."""
    return str(crystal_dir / "IRMOF-1.cif")


# ---------------------------------------------------------------------------
# _impl: database tools
# ---------------------------------------------------------------------------


def test_search_coremof_impl(coremof_csv):
    out = _impl.search_coremof_impl("Cu", field="metal", data_path=str(coremof_csv))
    assert out["success"] is True
    assert out["n_matches"] >= 1
    assert all("coreid" in r for r in out["records"])


def test_screen_coremof_impl(coremof_csv):
    out = _impl.screen_coremof_impl(
        metal="Cu", lcd_min=8.0, has_oms=True, data_path=str(coremof_csv)
    )
    assert out["success"] is True
    assert out["n_matches"] >= 1
    for rec in out["records"]:
        assert rec["metal_types"].find("Cu") != -1
        assert rec["lcd"] >= 8.0


def test_screen_coremof_impl_no_match(coremof_csv):
    out = _impl.screen_coremof_impl(lcd_min=1000.0, data_path=str(coremof_csv))
    assert out["success"] is True
    assert out["n_matches"] == 0


def test_search_coremof_impl_bad_path():
    out = _impl.search_coremof_impl("Cu", data_path="/nonexistent/does_not_exist.csv")
    assert out["success"] is False
    assert "error" in out


# ---------------------------------------------------------------------------
# _impl: structure resolution
# ---------------------------------------------------------------------------


def test_get_structure_impl_unconfigured(monkeypatch):
    # Ensure no structures dir is configured anywhere.
    monkeypatch.delenv("MOFFORGE_COREMOF_STRUCTURES_PATH", raising=False)
    from mofforge.utils.config import config

    monkeypatch.setattr(config, "coremof_structures_path", None)
    out = _impl.get_structure_impl("ABACUF")
    assert out["success"] is False
    assert "error" in out


def test_get_structure_impl_found(tmp_path):
    # Lay down a fake CIF and resolve it.
    (tmp_path / "ABACUF_clean.cif").write_text("dummy", encoding="utf-8")
    out = _impl.get_structure_impl("ABACUF", structures_dir=str(tmp_path))
    assert out["success"] is True
    assert out["cif_path"].endswith("ABACUF_clean.cif")


def test_get_structure_impl_missing_dir(tmp_path):
    out = _impl.get_structure_impl("ABACUF", structures_dir=str(tmp_path / "nope"))
    assert out["success"] is False


# ---------------------------------------------------------------------------
# _impl: adsorbate tools
# ---------------------------------------------------------------------------


def test_list_adsorbates_impl():
    out = _impl.list_adsorbates_impl()
    assert out["success"] is True
    assert "CO2" in out["adsorbates"]


def test_place_adsorbate_impl(host_cif, tmp_path):
    out = _impl.place_adsorbate_impl(
        host_cif,
        adsorbate="CO2",
        n_adsorbates=1,
        output_cif=str(tmp_path / "out.cif"),
        random_seed=42,
        grid_spacing=1.5,  # coarse grid for test speed
    )
    assert out["success"] is True
    assert out["n_adsorbates_placed"] >= 1
    assert out["atoms_after"] > out["atoms_before"]
    assert (tmp_path / "out.cif").exists()


def test_place_adsorbate_impl_bad_input():
    out = _impl.place_adsorbate_impl("/nonexistent.cif")
    assert out["success"] is False
    assert "error" in out


# ---------------------------------------------------------------------------
# _impl: validation
# ---------------------------------------------------------------------------


def test_validate_impl(host_cif):
    out = _impl.validate_impl(host_cif)
    assert out["success"] is True
    assert "is_valid" in out


# ---------------------------------------------------------------------------
# Stock server registration + JSON wrapping
# ---------------------------------------------------------------------------


def test_stock_server_registers_new_tools():
    import asyncio

    from mofforge.mcp import server

    tools = asyncio.run(server.mcp.list_tools())
    names = {t.name for t in tools}
    for expected in (
        "mofforge_search_coremof",
        "mofforge_screen_coremof",
        "mofforge_search_csd",
        "mofforge_lookup_mof",
        "mofforge_get_structure",
        "mofforge_place_adsorbate",
        "mofforge_list_adsorbates",
    ):
        assert expected in names


def test_stock_server_tool_returns_json(coremof_csv):
    from mofforge.mcp import server

    raw = server.mofforge_screen_coremof(metal="Cu", data_path=str(coremof_csv))
    parsed = json.loads(raw)
    assert parsed["success"] is True


# ---------------------------------------------------------------------------
# Functionalization + fragment-discovery tools
# ---------------------------------------------------------------------------


def test_list_functional_groups_impl():
    out = _impl.list_functional_groups_impl()
    assert out["success"] is True
    assert out["count"] >= 1
    names = {g["name"] for g in out["groups"]}
    assert {"NH2", "NO2", "F"} <= names


def test_find_sites_impl():
    out = _impl.find_sites_impl("O=C(O)c1ccc(C(=O)O)cc1")
    assert out["success"] is True
    assert out["n_sites"] == 4
    assert out["n_symmetry_classes"] == 1


def test_functionalize_impl(host_cif, tmp_path):
    pytest.importorskip("rdkit")
    out = _impl.functionalize_impl(
        host_cif,
        "O=C(O)c1ccc(C(=O)O)cc1",
        "F",
        sites=[0],
        coverage=0.5,
        output_cif=str(tmp_path / "f.cif"),
        random_seed=1,
    )
    assert out["success"] is True
    assert out["n_matches"] == 24
    assert out["n_functionalized"] == 12


def test_functionalize_campaign_impl(host_cif, tmp_path):
    pytest.importorskip("rdkit")
    out = _impl.functionalize_campaign_impl(
        host_cif,
        "O=C(O)c1ccc(C(=O)O)cc1",
        ["F", "NH2"],
        coverages=[1.0],
        output_dir=str(tmp_path),
        random_seed=1,
    )
    assert out["success"] is True
    assert out["n_results"] == 2


def test_list_and_get_fragment_impl():
    listed = _impl.list_fragments_impl()
    assert listed["success"] is True
    assert listed["count"] >= 1
    name = listed["fragments"][0]
    got = _impl.get_fragment_impl(name)
    assert got["success"] is True
    assert got["path"].endswith(name)


def test_get_fragment_impl_unknown():
    out = _impl.get_fragment_impl("does_not_exist.xyz")
    assert out["success"] is False


def test_stock_server_registers_functionalization_tools():
    import asyncio

    from mofforge.mcp import server

    names = {t.name for t in asyncio.run(server.mcp.list_tools())}
    for expected in (
        "mofforge_list_functional_groups",
        "mofforge_find_sites",
        "mofforge_functionalize",
        "mofforge_functionalize_campaign",
        "mofforge_list_fragments",
        "mofforge_get_fragment",
    ):
        assert expected in names


# ---------------------------------------------------------------------------
# CGFastMCP adapter (only when ChemGraph is installed)
# ---------------------------------------------------------------------------


def test_chemgraph_server_importable():
    """The module must import even without ChemGraph installed."""
    import mofforge.mcp.chemgraph_server as cs

    assert callable(cs._screen_place_worker)
    assert callable(cs.build_server)


def test_chemgraph_build_server():
    cs = pytest.importorskip("mofforge.mcp.chemgraph_server")
    try:
        mcp = cs.build_server()
    except ImportError:
        pytest.skip("ChemGraph (CGFastMCP) not installed")
    import asyncio

    names = {t.name for t in asyncio.run(mcp.list_tools())}
    assert "mofforge_screen_and_place" in names
    assert "mofforge_screen_coremof" in names
