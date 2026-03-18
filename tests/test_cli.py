"""Tests for the CLI interface."""

from click.testing import CliRunner

from tests.conftest import CRYSTAL_DIR, MOIETY_DIR


class TestCLI:
    """Tests for CLI commands."""

    def test_version(self):
        """CLI should show version."""
        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        assert "0.1.0" in result.output

    def test_help(self):
        """CLI should show help."""
        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "mofforge" in result.output.lower()

    def test_search_help(self):
        """Search command should show help."""
        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["search", "--help"])
        assert result.exit_code == 0
        assert "parent" in result.output.lower()

    def test_replace_help(self):
        """Replace command should show help."""
        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["replace", "--help"])
        assert result.exit_code == 0
        assert "query" in result.output.lower()

    def test_validate_help(self):
        """Validate command should show help."""
        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["validate", "--help"])
        assert result.exit_code == 0


class TestCLIIntegration:
    """Integration tests for CLI commands with real data."""

    def test_search_command(self):
        """Search command should find matches in IRMOF-1."""
        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "search",
                "-p",
                str(CRYSTAL_DIR / "IRMOF-1.cif"),
                "-q",
                str(MOIETY_DIR / "p-phenylene.xyz"),
            ],
        )
        assert result.exit_code == 0
        assert "Locations:" in result.output

    def test_validate_command(self):
        """Validate command should run without error."""
        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(
            main,
            ["validate", str(CRYSTAL_DIR / "IRMOF-1.cif")],
        )
        assert result.exit_code == 0
        assert "Validation Report" in result.output

    def test_replace_command(self):
        """Replace command should produce output file."""
        import tempfile

        from mofforge.cli import main

        runner = CliRunner()
        with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as f:
            output_path = f.name

        try:
            result = runner.invoke(
                main,
                [
                    "replace",
                    "-p",
                    str(CRYSTAL_DIR / "IRMOF-1.cif"),
                    "-q",
                    str(MOIETY_DIR / "2-!-p-phenylene.xyz"),
                    "-r",
                    str(MOIETY_DIR / "2-acetylamido-p-phenylene.xyz"),
                    "-o",
                    output_path,
                    "--nb-loc",
                    "1",
                ],
            )
            assert result.exit_code == 0
            assert "Output written to" in result.output
        finally:
            import os

            os.unlink(output_path)

    def test_search_missing_parent(self):
        """Search with non-existent parent should fail gracefully."""
        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "search",
                "-p",
                "/nonexistent/file.cif",
                "-q",
                str(MOIETY_DIR / "p-phenylene.xyz"),
            ],
        )
        assert result.exit_code != 0
