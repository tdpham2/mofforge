"""Tests for TobaccoBackend building-block staging in build().

Validates that build() clears the nodes/ and edges/ directories and
copies in the explicitly passed building blocks before running TOBACCO,
so that TOBACCO sees exactly what the caller specified.
"""

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
    (root / "configuration.py").write_text("CHARGES = True\nRUN_PARALLEL = False\n")
    (root / "templates").mkdir()
    (root / "nodes").mkdir()
    (root / "edges").mkdir()
    (root / "output_cifs").mkdir()
    (root / "nodes_database").mkdir()
    (root / "edges_database").mkdir()
    (root / "template_database").mkdir()
    return root


def _make_cif(directory, name, content="data_test\n"):
    """Create a minimal CIF file and return its path."""
    path = directory / name
    path.write_text(content)
    return path


class TestStageBuildingBlocks:
    """Tests for _stage_building_blocks() directly."""

    def test_clears_existing_and_copies_new_nodes(self, tmp_path):
        """Staging replaces old nodes with the passed ones."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        # Pre-populate with a stale node
        _make_cif(root / "nodes", "old_node.cif")
        assert (root / "nodes" / "old_node.cif").exists()

        backend = TobaccoBackend(root)

        # Create a new node source file
        new_src = _make_cif(tmp_path, "new_node.cif", "data_new_node\n")
        nodes = [BuildingBlock(name="new_node", role="node", source=new_src)]

        errors = backend._stage_building_blocks(nodes, edges=[])

        assert errors == []
        # Old node is gone
        assert not (root / "nodes" / "old_node.cif").exists()
        # New node is present
        assert (root / "nodes" / "new_node.cif").exists()
        assert (root / "nodes" / "new_node.cif").read_text() == "data_new_node\n"

    def test_clears_existing_and_copies_new_edges(self, tmp_path):
        """Staging replaces old edges with the passed ones."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        _make_cif(root / "edges", "old_edge.cif")

        backend = TobaccoBackend(root)

        new_src = _make_cif(tmp_path, "new_edge.cif", "data_new_edge\n")
        edges = [BuildingBlock(name="new_edge", role="edge", source=new_src)]

        errors = backend._stage_building_blocks(nodes=[], edges=edges)

        assert errors == []
        assert not (root / "edges" / "old_edge.cif").exists()
        assert (root / "edges" / "new_edge.cif").exists()

    def test_stages_multiple_nodes_and_edges(self, tmp_path):
        """Multiple nodes and edges are all staged correctly."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        backend = TobaccoBackend(root)

        n1 = _make_cif(tmp_path, "node_a.cif", "data_a\n")
        n2 = _make_cif(tmp_path, "node_b.cif", "data_b\n")
        e1 = _make_cif(tmp_path, "edge_x.cif", "data_x\n")

        nodes = [
            BuildingBlock(name="node_a", role="node", source=n1),
            BuildingBlock(name="node_b", role="node", source=n2),
        ]
        edges = [BuildingBlock(name="edge_x", role="edge", source=e1)]

        errors = backend._stage_building_blocks(nodes, edges)

        assert errors == []
        assert sorted(f.name for f in (root / "nodes").iterdir()) == ["node_a.cif", "node_b.cif"]
        assert [f.name for f in (root / "edges").iterdir()] == ["edge_x.cif"]

    def test_skips_staging_when_both_empty(self, tmp_path):
        """When no nodes or edges are passed, existing files are preserved."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        _make_cif(root / "nodes", "existing.cif")

        backend = TobaccoBackend(root)
        errors = backend._stage_building_blocks(nodes=[], edges=[])

        assert errors == []
        # Existing file should still be there (staging was skipped)
        assert (root / "nodes" / "existing.cif").exists()

    def test_clears_both_dirs_even_if_only_nodes_passed(self, tmp_path):
        """Both dirs are cleared even when only nodes are provided."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        _make_cif(root / "edges", "stale_edge.cif")

        backend = TobaccoBackend(root)

        n1 = _make_cif(tmp_path, "my_node.cif")
        nodes = [BuildingBlock(name="my_node", role="node", source=n1)]

        errors = backend._stage_building_blocks(nodes, edges=[])

        assert errors == []
        # Stale edge should be cleared even though we only passed nodes
        assert not (root / "edges" / "stale_edge.cif").exists()
        assert (root / "nodes" / "my_node.cif").exists()

    def test_returns_error_for_missing_source(self, tmp_path):
        """Returns an error when a building block source file doesn't exist."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        backend = TobaccoBackend(root)

        ghost = tmp_path / "nonexistent.cif"
        nodes = [BuildingBlock(name="ghost", role="node", source=ghost)]

        errors = backend._stage_building_blocks(nodes, edges=[])

        assert len(errors) == 1
        assert "Failed to stage node 'ghost'" in errors[0]

    def test_returns_error_for_non_cif_source(self, tmp_path):
        """Returns an error when a building block source is not a CIF."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        backend = TobaccoBackend(root)

        xyz_src = _make_cif(tmp_path, "node.xyz", "3\n\nC 0 0 0\n")
        nodes = [BuildingBlock(name="node", role="node", source=xyz_src)]

        errors = backend._stage_building_blocks(nodes, edges=[])

        assert len(errors) == 1
        assert "CIF" in errors[0]

    def test_partial_failure_reports_all_errors(self, tmp_path):
        """When some blocks fail and others succeed, all errors are reported."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        backend = TobaccoBackend(root)

        good_src = _make_cif(tmp_path, "good.cif", "data_good\n")
        bad_src = tmp_path / "missing.cif"  # does not exist

        nodes = [
            BuildingBlock(name="good", role="node", source=good_src),
            BuildingBlock(name="bad", role="node", source=bad_src),
        ]

        errors = backend._stage_building_blocks(nodes, edges=[])

        # One error for the missing file, but the good one was still copied
        assert len(errors) == 1
        assert "bad" in errors[0]
        assert (root / "nodes" / "good.cif").exists()


def _write_capturing_tobacco_py(root, capture_file):
    """Write a tobacco.py that records nodes/ and edges/ contents when called.

    The capture file is written as ``node1.cif,node2.cif|edge1.cif,edge2.cif``.
    """
    (root / "tobacco.py").write_text(
        "import os, pathlib, json\n"
        f"_CAPTURE = {str(capture_file)!r}\n"
        "def run_template(t):\n"
        "    root = pathlib.Path(os.getcwd())\n"
        "    nodes = sorted(f.name for f in (root / 'nodes').iterdir() if f.suffix == '.cif')\n"
        "    edges = sorted(f.name for f in (root / 'edges').iterdir() if f.suffix == '.cif')\n"
        "    with open(_CAPTURE, 'a') as fh:\n"
        "        json.dump({'nodes': nodes, 'edges': edges}, fh)\n"
        "        fh.write('\\n')\n"
    )


def _read_captures(capture_file):
    """Read captured node/edge snapshots from the capture file."""
    import json

    results = []
    if capture_file.exists():
        for line in capture_file.read_text().strip().splitlines():
            results.append(json.loads(line))
    return results


class TestBuildStaging:
    """Tests that build() integrates staging correctly end-to-end."""

    def test_build_stages_blocks_before_running(self, tmp_path):
        """build() copies passed nodes/edges into TOBACCO dirs before execution."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        _make_cif(root / "templates", "pcu.cif", "data_pcu\n")
        # Pre-populate with stale files
        _make_cif(root / "nodes", "stale_node.cif")
        _make_cif(root / "edges", "stale_edge.cif")

        capture_file = tmp_path / "capture.jsonl"
        _write_capturing_tobacco_py(root, capture_file)

        backend = TobaccoBackend(root)

        node_src = _make_cif(tmp_path, "fresh_node.cif", "data_fresh\n")
        edge_src = _make_cif(tmp_path, "fresh_edge.cif", "data_fresh_edge\n")

        nodes = [BuildingBlock(name="fresh_node", role="node", source=node_src)]
        edges = [BuildingBlock(name="fresh_edge", role="edge", source=edge_src)]
        topo = Topology(name="pcu")

        result = backend.build(topo, nodes, edges, output_dir=tmp_path / "output")

        captures = _read_captures(capture_file)
        assert len(captures) == 1
        # TOBACCO saw only the fresh blocks, not the stale ones
        assert captures[0]["nodes"] == ["fresh_node.cif"]
        assert captures[0]["edges"] == ["fresh_edge.cif"]

    def test_build_fails_early_on_staging_error(self, tmp_path):
        """build() returns failure immediately if staging fails."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        _make_cif(root / "templates", "pcu.cif")

        backend = TobaccoBackend(root)

        ghost = tmp_path / "nonexistent.cif"
        nodes = [BuildingBlock(name="ghost", role="node", source=ghost)]
        topo = Topology(name="pcu")

        result = backend.build(topo, nodes, edges=[], output_dir=tmp_path / "output")

        assert not result.success
        assert any("Failed to stage" in e for e in result.errors)

    def test_build_with_empty_blocks_preserves_existing(self, tmp_path):
        """build() with empty nodes/edges skips staging, keeping existing files."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        _make_cif(root / "templates", "pcu.cif")
        _make_cif(root / "nodes", "preloaded.cif")

        capture_file = tmp_path / "capture.jsonl"
        _write_capturing_tobacco_py(root, capture_file)

        backend = TobaccoBackend(root)
        topo = Topology(name="pcu")

        result = backend.build(topo, nodes=[], edges=[], output_dir=tmp_path / "output")

        captures = _read_captures(capture_file)
        assert len(captures) == 1
        # Preloaded node should still be there since staging was skipped
        assert captures[0]["nodes"] == ["preloaded.cif"]

    def test_build_staging_clears_stale_from_previous_run(self, tmp_path):
        """Two successive builds don't leak building blocks between runs."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        _make_cif(root / "templates", "pcu.cif")

        capture_file = tmp_path / "capture.jsonl"
        _write_capturing_tobacco_py(root, capture_file)

        backend = TobaccoBackend(root)
        topo = Topology(name="pcu")

        # First build with node_a
        src_a = _make_cif(tmp_path, "node_a.cif", "data_a\n")
        nodes_a = [BuildingBlock(name="node_a", role="node", source=src_a)]
        backend.build(topo, nodes_a, edges=[], output_dir=tmp_path / "out1")

        # Second build with node_b
        src_b = _make_cif(tmp_path, "node_b.cif", "data_b\n")
        nodes_b = [BuildingBlock(name="node_b", role="node", source=src_b)]
        backend.build(topo, nodes_b, edges=[], output_dir=tmp_path / "out2")

        captures = _read_captures(capture_file)
        assert len(captures) == 2
        assert captures[0]["nodes"] == ["node_a.cif"]
        # node_a must NOT be present in the second run
        assert captures[1]["nodes"] == ["node_b.cif"]

    def test_build_result_metadata_still_populated(self, tmp_path):
        """build() still populates metadata correctly after staging."""
        from mofforge.build.tobacco_backend import TobaccoBackend

        root = _make_tobacco_dir(tmp_path)
        _make_cif(root / "templates", "pcu.cif")

        backend = TobaccoBackend(root)

        src = _make_cif(tmp_path, "node.cif")
        nodes = [BuildingBlock(name="node", role="node", source=src)]
        topo = Topology(name="pcu")

        result = backend.build(topo, nodes, edges=[], output_dir=tmp_path / "output")

        assert result.success
        assert result.backend == "tobacco"
        assert result.metadata["templates_requested"] == ["pcu.cif"]
        assert result.metadata["templates_processed"] == ["pcu.cif"]
