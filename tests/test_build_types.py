"""Tests for build data types (BuildResult, BuildingBlock, Topology)."""

from pathlib import Path

from mofforge.build.base import BuildingBlock, BuildResult, Topology


class TestTopology:
    """Tests for the Topology dataclass."""

    def test_basic_creation(self):
        """Create a topology with just a name."""
        t = Topology(name="pcu")
        assert t.name == "pcu"
        assert t.source is None

    def test_with_source(self):
        """Create a topology with a source path."""
        t = Topology(name="pcu", source="/path/to/pcu.cif")
        assert t.name == "pcu"
        assert t.source == "/path/to/pcu.cif"


class TestBuildingBlock:
    """Tests for the BuildingBlock dataclass."""

    def test_cif_block(self, tmp_path):
        """Create a CIF-based building block."""
        cif = tmp_path / "node.cif"
        cif.write_text("data_test")
        bb = BuildingBlock(name="test_node", role="node", source=cif)
        assert bb.name == "test_node"
        assert bb.role == "node"
        assert bb.connection_points is None
        assert not bb.is_smiles

    def test_smiles_block(self):
        """SMILES-based building blocks are detected."""
        bb = BuildingBlock(
            name="bdc",
            role="edge",
            source="c1ccc(cc1)C(=O)O",
            connection_points=[0, 6],
        )
        assert bb.is_smiles
        assert bb.connection_points == [0, 6]

    def test_xyz_block_not_smiles(self, tmp_path):
        """XYZ file is not detected as SMILES."""
        xyz = tmp_path / "block.xyz"
        xyz.write_text("3\n\nC 0 0 0\nO 1 0 0\nH 0 1 0\n")
        bb = BuildingBlock(name="block", role="edge", source=xyz)
        assert not bb.is_smiles

    def test_role_values(self):
        """Role must be 'node' or 'edge'."""
        bb_node = BuildingBlock(name="n", role="node", source="n.cif")
        bb_edge = BuildingBlock(name="e", role="edge", source="e.cif")
        assert bb_node.role == "node"
        assert bb_edge.role == "edge"


class TestBuildResult:
    """Tests for the BuildResult dataclass."""

    def test_success_result(self):
        """Create a successful build result."""
        r = BuildResult(
            success=True,
            output_paths=[Path("out.cif")],
            backend="tobacco",
            elapsed_seconds=1.5,
        )
        assert r.success
        assert r.crystal is None
        assert r.errors == []
        assert r.backend == "tobacco"

    def test_failure_result(self):
        """Create a failed build result."""
        r = BuildResult(
            success=False,
            errors=["Template not found"],
            backend="pormake",
        )
        assert not r.success
        assert "Template not found" in r.errors

    def test_defaults(self):
        """Default values are sensible."""
        r = BuildResult(success=True)
        assert r.output_paths == []
        assert r.crystal is None
        assert r.errors == []
        assert r.elapsed_seconds == 0.0
        assert r.backend == ""
        assert r.metadata == {}
