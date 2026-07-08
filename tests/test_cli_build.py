"""Tests for the CLI build subcommands."""

from unittest.mock import patch

import pytest
from click.testing import CliRunner

from mofforge.cli import main

pytest.importorskip("tobacco3")


def _make_data_dir(tmp_path):
    """Create a minimal valid TOBACCO data directory with a small catalog."""
    root = tmp_path / "tobacco_data"
    root.mkdir()
    (root / "template_database").mkdir()
    (root / "nodes_database").mkdir()
    (root / "edges_database").mkdir()
    return root


class TestBuildCommand:
    """Tests for 'mofforge build'."""

    def test_build_requires_topology(self):
        """Build command requires --topology."""
        runner = CliRunner()
        result = runner.invoke(main, ["build", "--backend", "tobacco"])
        assert result.exit_code != 0
        assert "Missing option" in result.output or "required" in result.output.lower()

    def test_build_invalid_data_dir(self, tmp_path):
        """Build command fails gracefully with an empty (invalid) data dir."""
        runner = CliRunner()
        empty = tmp_path / "empty"
        empty.mkdir()
        with patch("mofforge.build.config._find_toml", return_value=None):
            result = runner.invoke(
                main,
                [
                    "build",
                    "--backend",
                    "tobacco",
                    "--topology",
                    "pcu",
                    "--tobacco-data-dir",
                    str(empty),
                ],
            )
        assert result.exit_code != 0

    def test_build_missing_template(self, tmp_path):
        """Build reports an error when the topology cannot be resolved."""
        root = _make_data_dir(tmp_path)
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "build",
                "--backend",
                "tobacco",
                "--topology",
                "no_such_topology",
                "--tobacco-data-dir",
                str(root),
            ],
        )
        assert result.exit_code != 0
        assert "not found" in result.output.lower() or "failed" in result.output.lower()

    def test_build_legacy_tobacco_path_flag(self, tmp_path):
        """The deprecated --tobacco-path flag is still accepted."""
        root = _make_data_dir(tmp_path)
        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "build",
                "--backend",
                "tobacco",
                "--topology",
                "no_such_topology",
                "--tobacco-path",
                str(root),
            ],
        )
        # Resolves the data dir (legacy alias), then fails on the missing topology.
        assert "not found" in result.output.lower() or "failed" in result.output.lower()


class TestBuildStatusCommand:
    """Tests for 'mofforge build-status'."""

    def test_build_status_tobacco(self, tmp_path):
        """build-status shows tobacco status."""
        root = _make_data_dir(tmp_path)
        (root / "template_database" / "pcu.cif").write_text("data")

        runner = CliRunner()
        result = runner.invoke(
            main,
            ["build-status", "--backend", "tobacco", "--tobacco-data-dir", str(root)],
        )
        assert result.exit_code == 0
        assert "templates_available" in result.output


class TestBuildListCommand:
    """Tests for 'mofforge build-list'."""

    def test_build_list_topologies(self, tmp_path):
        """build-list shows available topologies."""
        root = _make_data_dir(tmp_path)
        (root / "template_database" / "pcu.cif").write_text("data_pcu")
        (root / "template_database" / "dia.cif").write_text("data_dia")

        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "build-list",
                "--backend",
                "tobacco",
                "--type",
                "topologies",
                "--tobacco-data-dir",
                str(root),
            ],
        )
        assert result.exit_code == 0
        assert "pcu.cif" in result.output
        assert "dia.cif" in result.output

    def test_build_list_nodes(self, tmp_path):
        """build-list shows available nodes."""
        root = _make_data_dir(tmp_path)
        (root / "nodes_database" / "Zn_paddle.cif").write_text("data")

        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "build-list",
                "--backend",
                "tobacco",
                "--type",
                "nodes",
                "--tobacco-data-dir",
                str(root),
            ],
        )
        assert result.exit_code == 0
        assert "Zn_paddle.cif" in result.output

    def test_build_list_edges(self, tmp_path):
        """build-list shows available edges."""
        root = _make_data_dir(tmp_path)
        (root / "edges_database" / "BDC.cif").write_text("data")

        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "build-list",
                "--backend",
                "tobacco",
                "--type",
                "edges",
                "--tobacco-data-dir",
                str(root),
            ],
        )
        assert result.exit_code == 0
        assert "BDC.cif" in result.output

    def test_build_list_requires_type(self, tmp_path):
        """build-list requires --type."""
        root = _make_data_dir(tmp_path)
        runner = CliRunner()
        result = runner.invoke(
            main,
            ["build-list", "--backend", "tobacco", "--tobacco-data-dir", str(root)],
        )
        assert result.exit_code != 0
