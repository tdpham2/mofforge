"""Tests for the MOFBuilder facade."""

from unittest.mock import patch

import pytest

from mofforge.build.config import ConfigError

pytest.importorskip("tobacco3")


def _make_data_dir(tmp_path):
    """Create a minimal valid TOBACCO data directory with a small catalog."""
    root = tmp_path / "tobacco_data"
    root.mkdir()
    (root / "template_database").mkdir()
    (root / "nodes_database").mkdir()
    (root / "edges_database").mkdir()
    (root / "template_database" / "pcu.cif").write_text("data_pcu\n" * 5)
    (root / "template_database" / "dia.cif").write_text("data_dia\n" * 5)
    (root / "nodes_database" / "Zn_paddle.cif").write_text("data_zn")
    (root / "edges_database" / "BDC.cif").write_text("data_bdc")
    return root


class TestMOFBuilderInit:
    """Tests for MOFBuilder construction."""

    def test_unknown_backend_raises(self):
        """Unknown backend name raises ValueError."""
        from mofforge.build.builder import MOFBuilder

        with pytest.raises(ValueError, match="Unknown backend"):
            MOFBuilder(backend="nonexistent")

    def test_tobacco_backend_without_data_dir_raises(self):
        """TOBACCO backend without a data dir (and no auto-detect) raises ConfigError."""
        from mofforge.build.builder import MOFBuilder

        with (
            patch("mofforge.build.config._find_toml", return_value=None),
            patch("mofforge.build.config._autodetect_tobacco_data_dir", return_value=None),
            patch(
                "mofforge.build.config._fetch_tobacco_data_from_github",
                return_value=None,
            ),
            patch.dict("os.environ", {}, clear=True),
            pytest.raises(ConfigError),
        ):
            MOFBuilder(backend="tobacco")

    def test_tobacco_backend_with_data_dir(self, tmp_path):
        """TOBACCO backend initializes with an explicit data dir."""
        from mofforge.build.builder import MOFBuilder

        root = _make_data_dir(tmp_path)
        builder = MOFBuilder(backend="tobacco", tobacco_data_dir=str(root))
        assert builder.backend_name == "tobacco"

    def test_pormake_backend(self, tmp_path):
        """Pormake backend initializes with output_dir."""
        from mofforge.build.builder import MOFBuilder

        pytest.importorskip("pormake")
        with patch("mofforge.build.config._find_toml", return_value=None):
            builder = MOFBuilder(backend="pormake", output_dir=str(tmp_path))
        assert builder.backend_name == "pormake"


class TestMOFBuilderBuildingBlocks:
    """Tests for building-block registration via the facade."""

    def test_add_node_cif_registers(self, tmp_path):
        """Adding a CIF node records it in the in-memory registry."""
        from mofforge.build.builder import MOFBuilder

        root = _make_data_dir(tmp_path)
        cif = tmp_path / "my_node.cif"
        cif.write_text("data_node")

        builder = MOFBuilder(backend="tobacco", tobacco_data_dir=str(root))
        builder.add_node(cif)

        assert any(n.name == "my_node" for n in builder._nodes)

    def test_add_edge_cif_registers(self, tmp_path):
        """Adding a CIF edge records it in the in-memory registry."""
        from mofforge.build.builder import MOFBuilder

        root = _make_data_dir(tmp_path)
        cif = tmp_path / "my_edge.cif"
        cif.write_text("data_edge")

        builder = MOFBuilder(backend="tobacco", tobacco_data_dir=str(root))
        builder.add_edge(cif)

        assert any(e.name == "my_edge" for e in builder._edges)

    def test_add_with_custom_name(self, tmp_path):
        """Building blocks can have custom names."""
        from mofforge.build.builder import MOFBuilder

        root = _make_data_dir(tmp_path)
        cif = tmp_path / "block.cif"
        cif.write_text("data_block")

        builder = MOFBuilder(backend="tobacco", tobacco_data_dir=str(root))
        builder.add_node(cif, name="custom_name")

        assert any(n.name == "custom_name" for n in builder._nodes)

    def test_add_node_by_catalog_name(self, tmp_path):
        """A bare catalog name resolves without an on-disk path."""
        from mofforge.build.builder import MOFBuilder

        root = _make_data_dir(tmp_path)
        builder = MOFBuilder(backend="tobacco", tobacco_data_dir=str(root))
        builder.add_node("Zn_paddle.cif")
        assert any(n.name == "Zn_paddle" for n in builder._nodes)

    def test_remove_nodes_clears_registry(self, tmp_path):
        """Removing nodes drops them from the in-memory registry."""
        from mofforge.build.builder import MOFBuilder

        root = _make_data_dir(tmp_path)
        cif = tmp_path / "removable.cif"
        cif.write_text("data_removable")

        builder = MOFBuilder(backend="tobacco", tobacco_data_dir=str(root))
        builder.add_node(cif)
        result = builder.remove_nodes(["removable"], dry_run=False)
        assert result["success"]
        assert not any(n.name == "removable" for n in builder._nodes)

    def test_clear_edges(self, tmp_path):
        """Clearing edges empties the in-memory registry."""
        from mofforge.build.builder import MOFBuilder

        root = _make_data_dir(tmp_path)
        for i in range(3):
            cif = tmp_path / f"edge_{i}.cif"
            cif.write_text(f"data_{i}")

        builder = MOFBuilder(backend="tobacco", tobacco_data_dir=str(root))
        for i in range(3):
            builder.add_edge(tmp_path / f"edge_{i}.cif")
        result = builder.clear_edges(dry_run=False)
        assert result["success"]
        assert builder._edges == []


class TestMOFBuilderTopology:
    """Tests for topology listing via the facade."""

    def test_list_topologies(self, tmp_path):
        """Lists topologies from the catalog."""
        from mofforge.build.builder import MOFBuilder

        root = _make_data_dir(tmp_path)
        builder = MOFBuilder(backend="tobacco", tobacco_data_dir=str(root))
        topos = builder.list_topologies()
        assert "pcu.cif" in topos
        assert "dia.cif" in topos

    def test_describe_topology(self, tmp_path):
        """Describes a topology from the catalog."""
        from mofforge.build.builder import MOFBuilder

        root = _make_data_dir(tmp_path)
        builder = MOFBuilder(backend="tobacco", tobacco_data_dir=str(root))
        desc = builder.describe_topology("dia")
        assert "dia.cif" in desc

    def test_list_nodes_edges_from_catalog(self, tmp_path):
        """Node/edge listings come from the catalog directories."""
        from mofforge.build.builder import MOFBuilder

        root = _make_data_dir(tmp_path)
        builder = MOFBuilder(backend="tobacco", tobacco_data_dir=str(root))
        assert "Zn_paddle.cif" in builder.list_nodes()
        assert "BDC.cif" in builder.list_edges()


class TestMOFBuilderStatus:
    """Tests for status and configuration via the facade."""

    def test_status(self, tmp_path):
        """Status returns backend information."""
        from mofforge.build.builder import MOFBuilder

        root = _make_data_dir(tmp_path)
        builder = MOFBuilder(backend="tobacco", tobacco_data_dir=str(root))
        s = builder.status()
        assert s["success"]
        assert "data_dir" in s
        assert s["templates_available"] == 2

    def test_get_configuration(self, tmp_path):
        """Configuration is accessible and reflects TobaccoConfig fields."""
        from mofforge.build.builder import MOFBuilder

        root = _make_data_dir(tmp_path)
        builder = MOFBuilder(backend="tobacco", tobacco_data_dir=str(root))
        cfg = builder.get_configuration()
        assert cfg["success"]
        assert "CHARGES" in cfg["configuration"]

    def test_copy_from_database_lists(self, tmp_path):
        """copy_from_database delegates to backend and lists the catalog."""
        from mofforge.build.builder import MOFBuilder

        root = _make_data_dir(tmp_path)
        builder = MOFBuilder(backend="tobacco", tobacco_data_dir=str(root))
        result = builder.copy_from_database("node", names=None)
        assert result["success"]
        assert "Zn_paddle.cif" in result["available_in_database"]


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
        assert "_" in block.name
        assert block.connection_points == [0, 6]

        smiles2 = "c1ccc(cc1)c1ccc(cc1)N"
        block2 = MOFBuilder._make_block("edge", smiles2, None, [0, 6])
        assert block.name != block2.name

    def test_explicit_name_used(self):
        """Explicit name overrides inference."""
        from mofforge.build.builder import MOFBuilder

        block = MOFBuilder._make_block("node", "node.cif", "my_custom_name", None)
        assert block.name == "my_custom_name"
