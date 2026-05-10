"""Tests for the MOFBuilder facade."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from mofforge.build.base import BuildingBlock, BuildResult, Topology
from mofforge.build.config import ConfigError


def _make_tobacco_dir(tmp_path):
    """Create a minimal valid TOBACCO directory structure."""
    root = tmp_path / "tobacco_3.0"
    root.mkdir()
    (root / "tobacco.py").write_text("def run_template(t): pass\n")
    (root / "configuration.py").write_text("CHARGES = True\n")
    (root / "templates").mkdir()
    (root / "nodes").mkdir()
    (root / "edges").mkdir()
    (root / "output_cifs").mkdir()
    return root


class TestMOFBuilderInit:
    """Tests for MOFBuilder construction."""

    def test_unknown_backend_raises(self):
        """Unknown backend name raises ValueError."""
        from mofforge.build.builder import MOFBuilder

        with pytest.raises(ValueError, match="Unknown backend"):
            MOFBuilder(backend="nonexistent")

    def test_tobacco_backend_without_path_raises(self):
        """TOBACCO backend without a configured path raises ConfigError."""
        from mofforge.build.builder import MOFBuilder

        with (
            patch("mofforge.build.config._find_toml", return_value=None),
            patch.dict("os.environ", {}, clear=True),
            pytest.raises(ConfigError),
        ):
            MOFBuilder(backend="tobacco")

    def test_tobacco_backend_with_path(self, tmp_path):
        """TOBACCO backend initializes with explicit path."""
        from mofforge.build.builder import MOFBuilder

        root = _make_tobacco_dir(tmp_path)
        builder = MOFBuilder(backend="tobacco", tobacco_path=str(root))
        assert builder.backend_name == "tobacco"

    def test_pormake_backend(self, tmp_path):
        """Pormake backend initializes with output_dir."""
        from mofforge.build.builder import MOFBuilder

        with patch("mofforge.build.config._find_toml", return_value=None):
            builder = MOFBuilder(backend="pormake", output_dir=str(tmp_path))
        assert builder.backend_name == "pormake"


class TestMOFBuilderBuildingBlocks:
    """Tests for building-block registration via the facade."""

    def test_add_node_cif(self, tmp_path):
        """Adding a CIF node registers it."""
        from mofforge.build.builder import MOFBuilder

        root = _make_tobacco_dir(tmp_path)
        cif = tmp_path / "my_node.cif"
        cif.write_text("data_node")

        builder = MOFBuilder(backend="tobacco", tobacco_path=str(root))
        builder.add_node(cif)

        assert "my_node.cif" in builder.list_nodes()

    def test_add_edge_cif(self, tmp_path):
        """Adding a CIF edge registers it."""
        from mofforge.build.builder import MOFBuilder

        root = _make_tobacco_dir(tmp_path)
        cif = tmp_path / "my_edge.cif"
        cif.write_text("data_edge")

        builder = MOFBuilder(backend="tobacco", tobacco_path=str(root))
        builder.add_edge(cif)

        assert "my_edge.cif" in builder.list_edges()

    def test_add_with_custom_name(self, tmp_path):
        """Building blocks can have custom names."""
        from mofforge.build.builder import MOFBuilder

        root = _make_tobacco_dir(tmp_path)
        cif = tmp_path / "block.cif"
        cif.write_text("data_block")

        builder = MOFBuilder(backend="tobacco", tobacco_path=str(root))
        builder.add_node(cif, name="custom_name")

        # The internal block list uses the custom name
        assert any(n.name == "custom_name" for n in builder._nodes)

    def test_remove_nodes(self, tmp_path):
        """Removing nodes works via the facade."""
        from mofforge.build.builder import MOFBuilder

        root = _make_tobacco_dir(tmp_path)
        cif = tmp_path / "removable.cif"
        cif.write_text("data_removable")

        builder = MOFBuilder(backend="tobacco", tobacco_path=str(root))
        builder.add_node(cif)
        assert "removable.cif" in builder.list_nodes()

        result = builder.remove_nodes(["removable.cif"], dry_run=False)
        assert result["success"]

    def test_clear_edges(self, tmp_path):
        """Clearing edges removes all edge blocks."""
        from mofforge.build.builder import MOFBuilder

        root = _make_tobacco_dir(tmp_path)
        for i in range(3):
            cif = tmp_path / f"edge_{i}.cif"
            cif.write_text(f"data_{i}")
            (root / "edges" / f"edge_{i}.cif").write_text(f"data_{i}")

        builder = MOFBuilder(backend="tobacco", tobacco_path=str(root))
        result = builder.clear_edges(dry_run=False)
        assert result["success"]
        assert builder.list_edges() == []


class TestMOFBuilderTopology:
    """Tests for topology listing via the facade."""

    def test_list_topologies(self, tmp_path):
        """Lists topologies from the backend."""
        from mofforge.build.builder import MOFBuilder

        root = _make_tobacco_dir(tmp_path)
        (root / "templates" / "pcu.cif").write_text("data_pcu")

        builder = MOFBuilder(backend="tobacco", tobacco_path=str(root))
        topos = builder.list_topologies()
        assert "pcu.cif" in topos

    def test_describe_topology(self, tmp_path):
        """Describes a topology from the backend."""
        from mofforge.build.builder import MOFBuilder

        root = _make_tobacco_dir(tmp_path)
        (root / "templates" / "dia.cif").write_text("data_dia\n" * 5)

        builder = MOFBuilder(backend="tobacco", tobacco_path=str(root))
        desc = builder.describe_topology("dia")
        assert "dia.cif" in desc


class TestMOFBuilderStatus:
    """Tests for status and configuration via the facade."""

    def test_status(self, tmp_path):
        """Status returns backend information."""
        from mofforge.build.builder import MOFBuilder

        root = _make_tobacco_dir(tmp_path)
        builder = MOFBuilder(backend="tobacco", tobacco_path=str(root))
        s = builder.status()
        assert s["success"]
        assert "project_root" in s

    def test_get_configuration(self, tmp_path):
        """Configuration is accessible."""
        from mofforge.build.builder import MOFBuilder

        root = _make_tobacco_dir(tmp_path)
        builder = MOFBuilder(backend="tobacco", tobacco_path=str(root))
        cfg = builder.get_configuration()
        assert cfg["success"]

    def test_copy_from_database(self, tmp_path):
        """copy_from_database delegates to backend."""
        from mofforge.build.builder import MOFBuilder

        root = _make_tobacco_dir(tmp_path)
        (root / "nodes_database").mkdir(exist_ok=True)
        (root / "nodes_database" / "db_node.cif").write_text("data")

        builder = MOFBuilder(backend="tobacco", tobacco_path=str(root))
        result = builder.copy_from_database("node", names=None)
        assert result["success"]
        assert "db_node.cif" in result["available_in_database"]


class TestMOFBuilderMakeBlock:
    """Tests for the _make_block helper."""

    def test_cif_source_infers_name(self):
        """Name is inferred from CIF filename."""
        from mofforge.build.builder import MOFBuilder

        block = MOFBuilder._make_block("node", "/path/to/Zn_paddle.cif", None, None)
        assert block.name == "Zn_paddle"
        assert block.role == "node"

    def test_smiles_source_generates_unique_name(self):
        """Name derived from SMILES includes a hash to avoid collisions."""
        from mofforge.build.builder import MOFBuilder

        smiles = "c1ccc(cc1)c1ccc(cc1)C(=O)O"
        block = MOFBuilder._make_block("edge", smiles, None, [0, 6])
        # Name is prefix + "_" + 8-char hash
        assert "_" in block.name
        assert block.connection_points == [0, 6]

        # Two different SMILES with the same 12-char prefix get different names
        smiles2 = "c1ccc(cc1)c1ccc(cc1)N"
        block2 = MOFBuilder._make_block("edge", smiles2, None, [0, 6])
        assert block.name != block2.name

    def test_explicit_name_used(self):
        """Explicit name overrides inference."""
        from mofforge.build.builder import MOFBuilder

        block = MOFBuilder._make_block("node", "node.cif", "my_custom_name", None)
        assert block.name == "my_custom_name"
