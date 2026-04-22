"""Tests for build configuration (TOML loading, path validation)."""

import os
from pathlib import Path
from unittest.mock import patch

import pytest

from mofforge.build.config import BuildConfig, ConfigError, validate_tobacco_path


class TestValidateTobaccoPath:
    """Tests for validate_tobacco_path()."""

    def test_nonexistent_path(self, tmp_path):
        """Returns error for a path that doesn't exist."""
        errors = validate_tobacco_path(tmp_path / "nonexistent")
        assert len(errors) == 1
        assert "does not exist" in errors[0]

    def test_file_not_directory(self, tmp_path):
        """Returns error when path is a file, not a directory."""
        f = tmp_path / "afile.txt"
        f.write_text("hi")
        errors = validate_tobacco_path(f)
        assert len(errors) == 1
        assert "not a directory" in errors[0]

    def test_missing_files(self, tmp_path):
        """Returns errors for missing required files and dirs."""
        errors = validate_tobacco_path(tmp_path)
        # Should report missing tobacco.py, configuration.py, templates/, nodes/, edges/
        assert len(errors) == 5

    def test_valid_directory(self, tmp_path):
        """Returns no errors for a properly structured directory."""
        (tmp_path / "tobacco.py").write_text("")
        (tmp_path / "configuration.py").write_text("")
        (tmp_path / "templates").mkdir()
        (tmp_path / "nodes").mkdir()
        (tmp_path / "edges").mkdir()
        errors = validate_tobacco_path(tmp_path)
        assert errors == []

    def test_partial_structure(self, tmp_path):
        """Returns errors only for the missing pieces."""
        (tmp_path / "tobacco.py").write_text("")
        (tmp_path / "templates").mkdir()
        errors = validate_tobacco_path(tmp_path)
        # Missing: configuration.py, nodes/, edges/
        assert len(errors) == 3


class TestBuildConfig:
    """Tests for BuildConfig loading and resolution."""

    def test_default_config(self):
        """Default config has no tobacco_path and current dir for pormake."""
        with patch("mofforge.build.config._find_toml", return_value=None):
            cfg = BuildConfig.load()
        assert cfg.tobacco_path is None
        assert cfg.pormake_output_dir == Path(".")

    def test_explicit_kwargs_override(self, tmp_path):
        """Explicit kwargs take highest priority."""
        tobacco_dir = tmp_path / "tobacco"
        tobacco_dir.mkdir()
        with patch("mofforge.build.config._find_toml", return_value=None):
            cfg = BuildConfig.load(tobacco_path=str(tobacco_dir))
        assert cfg.tobacco_path == tobacco_dir.resolve()

    def test_env_var_override(self, tmp_path):
        """Environment variables override TOML values."""
        tobacco_dir = tmp_path / "tobacco_env"
        tobacco_dir.mkdir()
        with (
            patch("mofforge.build.config._find_toml", return_value=None),
            patch.dict(os.environ, {"MOFFORGE_TOBACCO_PATH": str(tobacco_dir)}),
        ):
            cfg = BuildConfig.load()
        assert cfg.tobacco_path == tobacco_dir.resolve()

    def test_toml_loading(self, tmp_path):
        """Config loads from a TOML file."""
        tobacco_dir = tmp_path / "tobacco_toml"
        tobacco_dir.mkdir()
        toml_file = tmp_path / "mofforge.toml"
        toml_file.write_text(
            f'[backends.tobacco]\npath = "{tobacco_dir}"\n'
            f'[backends.pormake]\noutput_dir = "/tmp/pm_out"\n'
        )
        with patch("mofforge.build.config._find_toml", return_value=toml_file):
            cfg = BuildConfig.load()
        assert cfg.tobacco_path == tobacco_dir.resolve()
        assert cfg.pormake_output_dir == Path("/tmp/pm_out")

    def test_kwargs_override_toml(self, tmp_path):
        """Explicit kwargs beat TOML values."""
        toml_dir = tmp_path / "from_toml"
        toml_dir.mkdir()
        kwarg_dir = tmp_path / "from_kwarg"
        kwarg_dir.mkdir()

        toml_file = tmp_path / "mofforge.toml"
        toml_file.write_text(f'[backends.tobacco]\npath = "{toml_dir}"\n')

        with patch("mofforge.build.config._find_toml", return_value=toml_file):
            cfg = BuildConfig.load(tobacco_path=str(kwarg_dir))
        assert cfg.tobacco_path == kwarg_dir.resolve()

    def test_resolve_tobacco_path_raises_when_none(self):
        """resolve_tobacco_path raises ConfigError when no path configured."""
        with patch("mofforge.build.config._find_toml", return_value=None):
            cfg = BuildConfig.load()
        with pytest.raises(ConfigError, match="not configured"):
            cfg.resolve_tobacco_path()

    def test_resolve_tobacco_path_raises_for_invalid(self, tmp_path):
        """resolve_tobacco_path raises ConfigError for invalid directory."""
        # Directory exists but is empty
        with patch("mofforge.build.config._find_toml", return_value=None):
            cfg = BuildConfig.load(tobacco_path=str(tmp_path))
        with pytest.raises(ConfigError, match="Invalid TOBACCO"):
            cfg.resolve_tobacco_path()

    def test_resolve_tobacco_path_succeeds(self, tmp_path):
        """resolve_tobacco_path returns path for valid directory."""
        (tmp_path / "tobacco.py").write_text("")
        (tmp_path / "configuration.py").write_text("")
        (tmp_path / "templates").mkdir()
        (tmp_path / "nodes").mkdir()
        (tmp_path / "edges").mkdir()

        with patch("mofforge.build.config._find_toml", return_value=None):
            cfg = BuildConfig.load(tobacco_path=str(tmp_path))
        result = cfg.resolve_tobacco_path()
        assert result == tmp_path.resolve()

    def test_validate_tobacco_convenience(self, tmp_path):
        """validate_tobacco() is a convenience wrapper for resolve_tobacco_path."""
        with patch("mofforge.build.config._find_toml", return_value=None):
            cfg = BuildConfig.load(tobacco_path=str(tmp_path))
        with pytest.raises(ConfigError):
            cfg.validate_tobacco()

    def test_env_var_pormake_output(self, tmp_path):
        """MOFFORGE_PORMAKE_OUTPUT_DIR env var sets pormake output dir."""
        with (
            patch("mofforge.build.config._find_toml", return_value=None),
            patch.dict(os.environ, {"MOFFORGE_PORMAKE_OUTPUT_DIR": str(tmp_path)}),
        ):
            cfg = BuildConfig.load()
        assert cfg.pormake_output_dir == tmp_path
