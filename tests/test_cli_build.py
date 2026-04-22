"""Tests for the CLI build subcommands."""

from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from mofforge.cli import main


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
    (root / "nodes_database").mkdir()
    (root / "edges_database").mkdir()
    (root / "template_database").mkdir()
    return root


class TestBuildCommand:
    """Tests for 'mofforge build'."""

    def test_build_requires_topology(self):
        """Build command requires --topology."""
        runner = CliRunner()
        result = runner.invoke(main, ["build", "--backend", "tobacco"])
        assert result.exit_code != 0
        assert "Missing option" in result.output or "required" in result.output.lower()

    def test_build_invalid_tobacco_path(self, tmp_path):
        """Build command fails gracefully with bad tobacco path."""
        runner = CliRunner()
        with patch("mofforge.build.config._find_toml", return_value=None):
            result = runner.invoke(
                main,
                [
                    "build",
                    "--backend",
                    "tobacco",
                    "--topology",
                    "pcu",
                    "--tobacco-path",
                    str(tmp_path),
                ],
            )
        assert result.exit_code != 0

    def test_build_missing_template(self, tmp_path):
        """Build reports error when template is missing."""
        root = _make_tobacco_dir(tmp_path)
        runner = CliRunner()
        result = runner.invoke(
            main,
            ["build", "--backend", "tobacco", "--topology", "pcu", "--tobacco-path", str(root)],
        )
        # Should fail because pcu.cif doesn't exist in templates/
        assert (
            result.exit_code != 0
            or "failed" in result.output.lower()
            or "not found" in result.output.lower()
        )


class TestBuildStatusCommand:
    """Tests for 'mofforge build-status'."""

    def test_build_status_tobacco(self, tmp_path):
        """build-status shows tobacco status."""
        root = _make_tobacco_dir(tmp_path)
        (root / "templates" / "pcu.cif").write_text("data")

        runner = CliRunner()
        result = runner.invoke(
            main,
            ["build-status", "--backend", "tobacco", "--tobacco-path", str(root)],
        )
        assert result.exit_code == 0
        assert "templates_available" in result.output


class TestBuildListCommand:
    """Tests for 'mofforge build-list'."""

    def test_build_list_topologies(self, tmp_path):
        """build-list shows available topologies."""
        root = _make_tobacco_dir(tmp_path)
        (root / "templates" / "pcu.cif").write_text("data_pcu")
        (root / "templates" / "dia.cif").write_text("data_dia")

        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "build-list",
                "--backend",
                "tobacco",
                "--type",
                "topologies",
                "--tobacco-path",
                str(root),
            ],
        )
        assert result.exit_code == 0
        assert "pcu.cif" in result.output
        assert "dia.cif" in result.output

    def test_build_list_nodes(self, tmp_path):
        """build-list shows available nodes."""
        root = _make_tobacco_dir(tmp_path)
        (root / "nodes" / "Zn_paddle.cif").write_text("data")

        runner = CliRunner()
        result = runner.invoke(
            main,
            ["build-list", "--backend", "tobacco", "--type", "nodes", "--tobacco-path", str(root)],
        )
        assert result.exit_code == 0
        assert "Zn_paddle.cif" in result.output

    def test_build_list_edges(self, tmp_path):
        """build-list shows available edges."""
        root = _make_tobacco_dir(tmp_path)
        (root / "edges" / "BDC.cif").write_text("data")

        runner = CliRunner()
        result = runner.invoke(
            main,
            ["build-list", "--backend", "tobacco", "--type", "edges", "--tobacco-path", str(root)],
        )
        assert result.exit_code == 0
        assert "BDC.cif" in result.output

    def test_build_list_requires_type(self, tmp_path):
        """build-list requires --type."""
        root = _make_tobacco_dir(tmp_path)
        runner = CliRunner()
        result = runner.invoke(
            main,
            ["build-list", "--backend", "tobacco", "--tobacco-path", str(root)],
        )
        assert result.exit_code != 0
