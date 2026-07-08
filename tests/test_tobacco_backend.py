"""Tests for the TOBACCO backend (importable ``tobacco3`` API).

The backend is a thin adapter over ``tobacco3.generate_mof``.  Most tests use a
temporary data directory of stub CIFs to exercise catalog discovery, name
resolution, and configuration.  One end-to-end test drives a real build against
the tobacco3 test fixtures when they are available.
"""

from pathlib import Path

import pytest

from mofforge.build.base import BuildingBlock, Topology
from mofforge.build.config import ConfigError

pytest.importorskip("tobacco3")


def _make_data_dir(tmp_path):
    """Create a minimal valid TOBACCO data directory."""
    root = tmp_path / "tobacco_data"
    root.mkdir()
    (root / "template_database").mkdir()
    (root / "nodes_database").mkdir()
    (root / "edges_database").mkdir()
    return root


def _tobacco3_fixtures():
    """Return the tobacco3 test fixture dir, or None if unavailable."""
    import tobacco3

    repo = Path(tobacco3.__file__).resolve().parent.parent
    fx = repo / "tests" / "fixtures"
    return fx if fx.is_dir() else None


class TestTobaccoBackendInit:
    """Tests for TobaccoBackend construction."""

    def test_valid_data_dir(self, tmp_path):
        """Backend initializes with a valid data directory."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        backend = TobaccoBackend(root)
        assert backend.name == "tobacco"

    def test_invalid_data_dir_raises(self, tmp_path):
        """Backend raises ConfigError for a directory with no data folders."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        with pytest.raises(ConfigError, match="Invalid TOBACCO data directory"):
            TobaccoBackend(tmp_path)


class TestTobaccoBackendTopologies:
    """Tests for topology listing."""

    def test_list_topologies_empty(self, tmp_path):
        """Empty database returns an empty list."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        backend = TobaccoBackend(root)
        assert backend.list_topologies() == []

    def test_list_topologies_with_cifs(self, tmp_path):
        """Database lists CIF files (union of database + active dirs)."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        (root / "template_database" / "pcu.cif").write_text("data_pcu")
        (root / "template_database" / "dia.cif").write_text("data_dia")
        (root / "template_database" / "readme.txt").write_text("ignored")

        backend = TobaccoBackend(root)
        assert backend.list_topologies() == ["dia.cif", "pcu.cif"]

    def test_describe_topology(self, tmp_path):
        """Describe returns info for an existing template."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        (root / "template_database" / "pcu.cif").write_text("data_pcu\n" * 10)

        backend = TobaccoBackend(root)
        assert "pcu.cif" in backend.describe_topology("pcu")

    def test_describe_topology_not_found(self, tmp_path):
        """Describe returns a message for a missing template."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        backend = TobaccoBackend(root)
        assert "not found" in backend.describe_topology("nonexistent")


class TestTobaccoBackendBuildingBlocks:
    """Tests for building-block discovery and validation."""

    def test_list_building_blocks(self, tmp_path):
        """Lists CIF files from the node/edge database directories."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        (root / "nodes_database" / "Zn_paddle.cif").write_text("data_zn")
        (root / "edges_database" / "BDC.cif").write_text("data_bdc")

        backend = TobaccoBackend(root)
        assert backend.list_building_blocks("node") == ["Zn_paddle.cif"]
        assert backend.list_building_blocks("edge") == ["BDC.cif"]

    def test_add_building_block_resolves_path(self, tmp_path):
        """add_building_block validates a real CIF path without copying."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        src = tmp_path / "my_node.cif"
        src.write_text("data_my_node")

        backend = TobaccoBackend(root)
        block = BuildingBlock(name="my_node", role="node", source=src)
        result = backend.add_building_block(block)
        assert result["success"]
        assert result["file"] == str(src.resolve())
        # Nothing copied into the data dir.
        assert backend.list_building_blocks("node") == []

    def test_add_building_block_resolves_catalog_name(self, tmp_path):
        """A bare catalog name resolves against the database directory."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        (root / "nodes_database" / "Zn_paddle.cif").write_text("data")

        backend = TobaccoBackend(root)
        block = BuildingBlock(name="Zn_paddle", role="node", source="Zn_paddle.cif")
        result = backend.add_building_block(block)
        assert result["success"]

    def test_add_building_block_not_cif(self, tmp_path):
        """Rejects non-CIF sources."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        src = tmp_path / "node.xyz"
        src.write_text("3\n\nC 0 0 0\n")

        backend = TobaccoBackend(root)
        block = BuildingBlock(name="node", role="node", source=src)
        result = backend.add_building_block(block)
        assert not result["success"]
        assert "CIF" in result["error"]

    def test_add_building_block_unresolvable(self, tmp_path):
        """Reports failure for a name that cannot be resolved."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        backend = TobaccoBackend(root)
        block = BuildingBlock(name="ghost", role="node", source="ghost.cif")
        result = backend.add_building_block(block)
        assert not result["success"]

    def test_copy_from_database_list(self, tmp_path):
        """List available files in the database directories."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        (root / "nodes_database" / "db_node.cif").write_text("data")

        backend = TobaccoBackend(root)
        result = backend.copy_from_database("node", names=None)
        assert result["success"]
        assert "db_node.cif" in result["available_in_database"]

    def test_copy_from_database_resolve(self, tmp_path):
        """Resolve specific names to their CIF paths."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        (root / "nodes_database" / "db_node.cif").write_text("data")

        backend = TobaccoBackend(root)
        result = backend.copy_from_database("node", names=["db_node.cif"], dry_run=False)
        assert result["success"]
        assert "db_node.cif" in result["resolved"]


class TestTobaccoBackendConfiguration:
    """Tests for configuration backed by TobaccoConfig."""

    def test_get_configuration(self, tmp_path):
        """Reads TobaccoConfig defaults."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        backend = TobaccoBackend(root)
        cfg = backend.get_configuration()
        assert cfg["success"]
        assert cfg["configuration"]["CHARGES"] is True

    def test_set_configuration(self, tmp_path):
        """Sets a TobaccoConfig field."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        backend = TobaccoBackend(root)
        result = backend.set_configuration("CHARGES", False)
        assert result["success"]
        assert backend.get_configuration()["configuration"]["CHARGES"] is False

    def test_set_configuration_unknown_key(self, tmp_path):
        """Setting an unknown key fails."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        backend = TobaccoBackend(root)
        result = backend.set_configuration("NONEXISTENT", True)
        assert not result["success"]


class TestTobaccoBackendStatus:
    """Tests for status reporting."""

    def test_status(self, tmp_path):
        """Status returns counts and readiness."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        (root / "template_database" / "pcu.cif").write_text("data")
        (root / "nodes_database" / "node.cif").write_text("data")
        (root / "edges_database" / "edge1.cif").write_text("data")
        (root / "edges_database" / "edge2.cif").write_text("data")

        backend = TobaccoBackend(root)
        s = backend.status()
        assert s["success"]
        assert s["templates_available"] == 1
        assert s["nodes_available"] == 1
        assert s["edges_available"] == 2
        assert s["ready"] is True


class TestTobaccoBackendBuild:
    """End-to-end build via generate_mof."""

    def test_build_requires_nodes_and_edges(self, tmp_path):
        """Build fails cleanly when no building blocks are provided."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        (root / "template_database" / "pcu.cif").write_text("data")
        backend = TobaccoBackend(root)
        result = backend.build(Topology(name="pcu"), [], [], tmp_path / "out")
        assert not result.success
        assert any("node" in e for e in result.errors)

    def test_build_unknown_topology(self, tmp_path):
        """Build fails when the topology cannot be resolved."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_data_dir(tmp_path)
        (root / "nodes_database" / "n.cif").write_text("data")
        (root / "edges_database" / "e.cif").write_text("data")
        backend = TobaccoBackend(root)
        nodes = [BuildingBlock(name="n", role="node", source="n.cif")]
        edges = [BuildingBlock(name="e", role="edge", source="e.cif")]
        result = backend.build(Topology(name="ghost"), nodes, edges, tmp_path / "out")
        assert not result.success
        assert any("not found" in e for e in result.errors)

    def test_build_end_to_end(self, tmp_path):
        """A real build produces structures and a loaded Crystal."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        fx = _tobacco3_fixtures()
        if fx is None:
            pytest.skip("tobacco3 test fixtures not available")

        # Use the installed tobacco3 repo as the data dir so name lookups work,
        # but pass fixture CIFs by explicit path.
        import tobacco3

        data_dir = Path(tobacco3.__file__).resolve().parent.parent
        backend = TobaccoBackend(data_dir)

        nodes = [
            BuildingBlock(name=n, role="node", source=str(fx / "nodes" / f"{n}.cif"))
            for n in ("4c_1Zn_Ch", "triazole")
        ]
        edges = [
            BuildingBlock(name=e, role="edge", source=str(fx / "edges" / f"{e}.cif"))
            for e in ("ntn_edge", "oxalic_edge", "squOxa_ch")
        ]
        topo = Topology(name="dmc", source=str(fx / "templates" / "dmc.cif"))

        result = backend.build(topo, nodes, edges, tmp_path / "out")
        assert result.success
        assert result.output_paths
        assert result.crystal is not None
        assert result.crystal.n_atoms > 0
        assert result.metadata["n_structures"] == len(result.output_paths)
