"""Tests for batch processing."""

import tempfile
from pathlib import Path

import pytest
import yaml

from mofforge.batch import BatchConfig, BatchResult, _resolve_parent_paths
from tests.conftest import CRYSTAL_DIR, MOIETY_DIR


class TestBatchConfig:
    """Tests for BatchConfig YAML parsing and validation."""

    def test_from_yaml_basic(self):
        """Load a basic YAML config."""
        config_data = {
            "parents": [{"path": str(CRYSTAL_DIR / "IRMOF-1.cif")}],
            "operations": [{"type": "validate"}],
            "output": {
                "directory": "test_output",
                "format": "cif",
                "naming": "{parent_name}_test",
            },
            "parallel": 0,
            "moiety_path": str(MOIETY_DIR),
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(config_data, f)
            filepath = f.name

        try:
            config = BatchConfig.from_yaml(filepath)
            assert len(config.parent_paths) == 1
            assert config.output_format == "cif"
            assert config.naming == "{parent_name}_test"
            assert config.moiety_path == str(MOIETY_DIR)
        finally:
            Path(filepath).unlink(missing_ok=True)

    def test_from_yaml_empty_file(self):
        """Empty YAML file should not crash (returns empty config)."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            f.write("")  # empty file
            filepath = f.name

        try:
            config = BatchConfig.from_yaml(filepath)
            assert config.parent_paths == []
            assert config.operations == []
        finally:
            Path(filepath).unlink(missing_ok=True)

    def test_from_yaml_string_parents(self):
        """Parents can be plain strings (not dicts)."""
        config_data = {
            "parents": ["path/to/file1.cif", "path/to/file2.cif"],
            "operations": [],
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(config_data, f)
            filepath = f.name

        try:
            config = BatchConfig.from_yaml(filepath)
            assert len(config.parent_paths) == 2
            assert config.parent_paths[0] == "path/to/file1.cif"
        finally:
            Path(filepath).unlink(missing_ok=True)

    def test_invalid_output_format_raises(self):
        """Invalid output format should raise ValueError."""
        with pytest.raises(ValueError, match="Unsupported output format"):
            BatchConfig(output_format="pdb")

    def test_valid_output_formats(self):
        """Valid output formats should not raise."""
        config_cif = BatchConfig(output_format="cif")
        assert config_cif.output_format == "cif"

        config_xyz = BatchConfig(output_format="xyz")
        assert config_xyz.output_format == "xyz"

    def test_fragment_path_deprecation(self):
        """The 'fragment_path' YAML key should work with a deprecation warning."""
        config_data = {
            "parents": [],
            "operations": [],
            "fragment_path": "/some/path",
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(config_data, f)
            filepath = f.name

        try:
            with pytest.warns(DeprecationWarning, match="fragment_path"):
                config = BatchConfig.from_yaml(filepath)
            assert config.moiety_path == "/some/path"
        finally:
            Path(filepath).unlink(missing_ok=True)

    def test_moiety_path_preferred_over_fragment_path(self):
        """If both keys are present, moiety_path takes precedence."""
        config_data = {
            "parents": [],
            "operations": [],
            "moiety_path": "/preferred",
            "fragment_path": "/deprecated",
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(config_data, f)
            filepath = f.name

        try:
            config = BatchConfig.from_yaml(filepath)
            assert config.moiety_path == "/preferred"
        finally:
            Path(filepath).unlink(missing_ok=True)


class TestResolveParentPaths:
    """Tests for _resolve_parent_paths."""

    def test_exact_path(self):
        """Exact path should resolve correctly."""
        cif_path = str(CRYSTAL_DIR / "IRMOF-1.cif")
        paths = _resolve_parent_paths([cif_path])
        assert len(paths) == 1
        assert paths[0].name == "IRMOF-1.cif"

    def test_glob_pattern(self):
        """Glob pattern should resolve to multiple files."""
        pattern = str(CRYSTAL_DIR / "*.cif")
        paths = _resolve_parent_paths([pattern])
        assert len(paths) > 0
        assert all(p.suffix == ".cif" for p in paths)

    def test_nonexistent_pattern_warns(self):
        """Non-existent pattern should log warning and return empty."""
        paths = _resolve_parent_paths(["/nonexistent/path/*.cif"])
        assert len(paths) == 0


class TestBatchResult:
    """Tests for BatchResult dataclass."""

    def test_default_success(self):
        """BatchResult defaults to success."""
        result = BatchResult(parent_name="test")
        assert result.success is True
        assert result.error is None
        assert result.output_path is None
        assert result.validation is None

    def test_failure_result(self):
        """BatchResult can represent failure."""
        result = BatchResult(parent_name="test", success=False, error="something went wrong")
        assert result.success is False
        assert result.error == "something went wrong"
