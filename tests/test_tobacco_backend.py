"""Tests for the TOBACCO backend.

Since TOBACCO is an external tool that may not be available in CI,
these tests mock the tobacco imports and focus on file management,
configuration, and integration logic.
"""

import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from mofforge.build.base import BuildingBlock, Topology
from mofforge.build.config import ConfigError


def _make_tobacco_dir(tmp_path):
    """Create a minimal valid TOBACCO directory structure."""
    root = tmp_path / "tobacco_3.0"
    root.mkdir()
    (root / "tobacco.py").write_text("def run_template(t): pass\n")
    (root / "configuration.py").write_text("CHARGES = True\nRUN_PARALLEL = False\n")
    (root / "templates").mkdir()
    (root / "nodes").mkdir()
    (root / "edges").mkdir()
    (root / "output_cifs").mkdir()
    # Databases
    (root / "nodes_database").mkdir()
    (root / "edges_database").mkdir()
    (root / "template_database").mkdir()
    return root


class TestTobaccoBackendInit:
    """Tests for TobaccoBackend construction."""

    def test_valid_path(self, tmp_path):
        """Backend initializes with a valid TOBACCO directory."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        backend = TobaccoBackend(root)
        assert backend.name == "tobacco"

    def test_invalid_path_raises(self, tmp_path):
        """Backend raises ConfigError for invalid directory."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        with pytest.raises(ConfigError, match="Invalid TOBACCO"):
            TobaccoBackend(tmp_path)


class TestTobaccoBackendTopologies:
    """Tests for topology listing."""

    def test_list_topologies_empty(self, tmp_path):
        """Empty templates directory returns empty list."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        backend = TobaccoBackend(root)
        assert backend.list_topologies() == []

    def test_list_topologies_with_cifs(self, tmp_path):
        """Templates directory lists CIF files."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        (root / "templates" / "pcu.cif").write_text("data_pcu")
        (root / "templates" / "dia.cif").write_text("data_dia")
        (root / "templates" / "readme.txt").write_text("ignored")

        backend = TobaccoBackend(root)
        topos = backend.list_topologies()
        assert topos == ["dia.cif", "pcu.cif"]

    def test_describe_topology(self, tmp_path):
        """Describe returns info for existing template."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        (root / "templates" / "pcu.cif").write_text("data_pcu\n" * 10)

        backend = TobaccoBackend(root)
        desc = backend.describe_topology("pcu")
        assert "pcu.cif" in desc

    def test_describe_topology_not_found(self, tmp_path):
        """Describe returns message for missing template."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        backend = TobaccoBackend(root)
        desc = backend.describe_topology("nonexistent")
        assert "not found" in desc


class TestTobaccoBackendBuildingBlocks:
    """Tests for building-block management."""

    def test_list_building_blocks(self, tmp_path):
        """Lists CIF files in nodes/ or edges/."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        (root / "nodes" / "Zn_paddle.cif").write_text("data_zn")
        (root / "edges" / "BDC.cif").write_text("data_bdc")

        backend = TobaccoBackend(root)
        assert backend.list_building_blocks("node") == ["Zn_paddle.cif"]
        assert backend.list_building_blocks("edge") == ["BDC.cif"]

    def test_add_building_block(self, tmp_path):
        """Copies a CIF into the right directory."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        src = tmp_path / "my_node.cif"
        src.write_text("data_my_node")

        backend = TobaccoBackend(root)
        block = BuildingBlock(name="my_node", role="node", source=src)
        result = backend.add_building_block(block)

        assert result["success"]
        assert (root / "nodes" / "my_node.cif").exists()

    def test_add_building_block_not_cif(self, tmp_path):
        """Rejects non-CIF sources."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        src = tmp_path / "node.xyz"
        src.write_text("3\n\nC 0 0 0\n")

        backend = TobaccoBackend(root)
        block = BuildingBlock(name="node", role="node", source=src)
        result = backend.add_building_block(block)
        assert not result["success"]
        assert "CIF" in result["error"]

    def test_remove_building_blocks_dry_run(self, tmp_path):
        """Dry-run removal reports what would be removed."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        (root / "nodes" / "a.cif").write_text("data_a")
        (root / "nodes" / "b.cif").write_text("data_b")

        backend = TobaccoBackend(root)
        result = backend.remove_building_blocks("node", ["a.cif"], dry_run=True)
        assert result["dry_run"]
        assert "a.cif" in result["would_remove"]
        assert (root / "nodes" / "a.cif").exists()  # Not actually removed

    def test_remove_building_blocks_force(self, tmp_path):
        """Force removal actually deletes files."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        (root / "nodes" / "a.cif").write_text("data_a")

        backend = TobaccoBackend(root)
        result = backend.remove_building_blocks("node", ["a.cif"], dry_run=False)
        assert result["success"]
        assert not (root / "nodes" / "a.cif").exists()

    def test_clear_building_blocks(self, tmp_path):
        """Clear removes all CIFs in a directory."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        (root / "edges" / "a.cif").write_text("data")
        (root / "edges" / "b.cif").write_text("data")

        backend = TobaccoBackend(root)
        result = backend.clear_building_blocks("edge", dry_run=False)
        assert result["success"]
        assert result["count"] == 2
        assert backend.list_building_blocks("edge") == []

    def test_copy_from_database_list(self, tmp_path):
        """List available files in the database directory."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        (root / "nodes_database" / "db_node.cif").write_text("data")

        backend = TobaccoBackend(root)
        result = backend.copy_from_database("node", names=None)
        assert result["success"]
        assert "db_node.cif" in result["available_in_database"]

    def test_copy_from_database_copy(self, tmp_path):
        """Copy files from database to active directory."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        (root / "nodes_database" / "db_node.cif").write_text("data_db_node")

        backend = TobaccoBackend(root)
        result = backend.copy_from_database("node", names=["db_node.cif"], dry_run=False)
        assert result["success"]
        assert (root / "nodes" / "db_node.cif").exists()


class TestTobaccoBackendConfiguration:
    """Tests for configuration reading/writing."""

    def test_get_configuration(self, tmp_path):
        """Reads configuration.py values."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        backend = TobaccoBackend(root)
        cfg = backend.get_configuration()
        assert cfg["success"]
        assert cfg["configuration"]["CHARGES"] is True
        assert cfg["configuration"]["RUN_PARALLEL"] is False

    def test_set_configuration(self, tmp_path):
        """Modifies a key in configuration.py."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        backend = TobaccoBackend(root)
        result = backend.set_configuration("CHARGES", "False")
        assert result["success"]

        # Verify it was written
        cfg = backend.get_configuration()
        assert cfg["configuration"]["CHARGES"] is False

    def test_set_configuration_unknown_key(self, tmp_path):
        """Setting an unknown key fails."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        backend = TobaccoBackend(root)
        result = backend.set_configuration("NONEXISTENT", "True")
        assert not result["success"]


class TestTobaccoBackendStatus:
    """Tests for status reporting."""

    def test_status(self, tmp_path):
        """Status returns counts and readiness."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        (root / "templates" / "pcu.cif").write_text("data")
        (root / "nodes" / "node.cif").write_text("data")
        (root / "edges" / "edge1.cif").write_text("data")
        (root / "edges" / "edge2.cif").write_text("data")

        backend = TobaccoBackend(root)
        s = backend.status()
        assert s["success"]
        assert s["templates_available"] == 1
        assert s["nodes_available"] == 1
        assert s["edges_available"] == 2
        assert s["ready"] is True
