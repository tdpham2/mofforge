"""Tests for build configuration (TOML loading, data-dir resolution)."""

import os
from pathlib import Path
from unittest.mock import patch

import pytest

from mofforge.build.config import (
    BuildConfig,
    ConfigError,
    validate_tobacco_data_dir,
)


def _make_data_dir(root: Path) -> Path:
    """Create a minimal valid TOBACCO data directory."""
    root.mkdir(parents=True, exist_ok=True)
    (root / "template_database").mkdir()
    (root / "nodes_database").mkdir()
    (root / "edges_database").mkdir()
    return root


class TestValidateTobaccoDataDir:
    """Tests for validate_tobacco_data_dir()."""

    def test_nonexistent_path(self, tmp_path):
        """Returns error for a path that doesn't exist."""
        errors = validate_tobacco_data_dir(tmp_path / "nonexistent")
        assert len(errors) == 1
        assert "does not exist" in errors[0]

    def test_file_not_directory(self, tmp_path):
        """Returns error when path is a file, not a directory."""
        f = tmp_path / "afile.txt"
        f.write_text("hi")
        errors = validate_tobacco_data_dir(f)
        assert len(errors) == 1
        assert "not a directory" in errors[0]

    def test_missing_data_dirs(self, tmp_path):
        """Returns an error when no data subdirectories are present."""
        errors = validate_tobacco_data_dir(tmp_path)
        assert len(errors) == 1
        assert "No TOBACCO data directories" in errors[0]

    def test_valid_directory(self, tmp_path):
        """Returns no errors when at least one data subdirectory exists."""
        _make_data_dir(tmp_path)
        assert validate_tobacco_data_dir(tmp_path) == []

    def test_active_templates_only_is_valid(self, tmp_path):
        """A bare ``templates/`` folder is sufficient."""
        (tmp_path / "templates").mkdir()
        assert validate_tobacco_data_dir(tmp_path) == []


class TestBuildConfig:
    """Tests for BuildConfig loading and resolution."""

    def test_default_config(self):
        """Default config has no data dir and current dir for pormake."""
        with patch("mofforge.build.config._find_toml", return_value=None):
            cfg = BuildConfig.load()
        assert cfg.tobacco_data_dir is None
        assert cfg.pormake_output_dir == Path(".")

    def test_explicit_kwargs_override(self, tmp_path):
        """Explicit kwargs take highest priority."""
        data_dir = _make_data_dir(tmp_path / "tobacco")
        with patch("mofforge.build.config._find_toml", return_value=None):
            cfg = BuildConfig.load(tobacco_data_dir=str(data_dir))
        assert cfg.tobacco_data_dir == data_dir.resolve()

    def test_legacy_tobacco_path_kwarg(self, tmp_path):
        """The legacy ``tobacco_path`` kwarg still maps to the data dir."""
        data_dir = _make_data_dir(tmp_path / "tobacco_legacy")
        with patch("mofforge.build.config._find_toml", return_value=None):
            cfg = BuildConfig.load(tobacco_path=str(data_dir))
        assert cfg.tobacco_data_dir == data_dir.resolve()

    def test_env_var_override(self, tmp_path):
        """MOFFORGE_TOBACCO_DATA overrides TOML values."""
        data_dir = tmp_path / "tobacco_env"
        data_dir.mkdir()
        with (
            patch("mofforge.build.config._find_toml", return_value=None),
            patch.dict(os.environ, {"MOFFORGE_TOBACCO_DATA": str(data_dir)}),
        ):
            cfg = BuildConfig.load()
        assert cfg.tobacco_data_dir == data_dir.resolve()

    def test_legacy_env_var(self, tmp_path):
        """The legacy MOFFORGE_TOBACCO_PATH env var still works."""
        data_dir = tmp_path / "tobacco_env_legacy"
        data_dir.mkdir()
        with (
            patch("mofforge.build.config._find_toml", return_value=None),
            patch.dict(os.environ, {"MOFFORGE_TOBACCO_PATH": str(data_dir)}, clear=False),
        ):
            os.environ.pop("MOFFORGE_TOBACCO_DATA", None)
            cfg = BuildConfig.load()
        assert cfg.tobacco_data_dir == data_dir.resolve()

    def test_toml_loading(self, tmp_path):
        """Config loads from a TOML file (new data_dir key)."""
        data_dir = tmp_path / "tobacco_toml"
        data_dir.mkdir()
        toml_file = tmp_path / "mofforge.toml"
        toml_file.write_text(
            f'[backends.tobacco]\ndata_dir = "{data_dir}"\n'
            f'[backends.pormake]\noutput_dir = "/tmp/pm_out"\n'
        )
        with patch("mofforge.build.config._find_toml", return_value=toml_file):
            cfg = BuildConfig.load()
        assert cfg.tobacco_data_dir == data_dir.resolve()
        assert cfg.pormake_output_dir == Path("/tmp/pm_out")

    def test_toml_legacy_path_key(self, tmp_path):
        """The legacy ``path`` TOML key still maps to the data dir."""
        data_dir = tmp_path / "tobacco_toml_legacy"
        data_dir.mkdir()
        toml_file = tmp_path / "mofforge.toml"
        toml_file.write_text(f'[backends.tobacco]\npath = "{data_dir}"\n')
        with patch("mofforge.build.config._find_toml", return_value=toml_file):
            cfg = BuildConfig.load()
        assert cfg.tobacco_data_dir == data_dir.resolve()

    def test_kwargs_override_toml(self, tmp_path):
        """Explicit kwargs beat TOML values."""
        toml_dir = tmp_path / "from_toml"
        toml_dir.mkdir()
        kwarg_dir = _make_data_dir(tmp_path / "from_kwarg")

        toml_file = tmp_path / "mofforge.toml"
        toml_file.write_text(f'[backends.tobacco]\ndata_dir = "{toml_dir}"\n')

        with patch("mofforge.build.config._find_toml", return_value=toml_file):
            cfg = BuildConfig.load(tobacco_data_dir=str(kwarg_dir))
        assert cfg.tobacco_data_dir == kwarg_dir.resolve()

    def test_resolve_raises_when_unconfigured_and_no_autodetect(self):
        """resolve raises ConfigError when nothing configured, no auto-detect, no fetch."""
        with (
            patch("mofforge.build.config._find_toml", return_value=None),
            patch("mofforge.build.config._autodetect_tobacco_data_dir", return_value=None),
            patch(
                "mofforge.build.config._fetch_tobacco_data_from_github",
                return_value=None,
            ),
        ):
            cfg = BuildConfig.load()
            with pytest.raises(ConfigError, match="not configured"):
                cfg.resolve_tobacco_data_dir()

    def test_resolve_fetches_from_github_as_last_resort(self, tmp_path):
        """resolve falls back to the GitHub fetch when unconfigured + no autodetect."""
        data_dir = _make_data_dir(tmp_path / "fetched")
        with (
            patch("mofforge.build.config._find_toml", return_value=None),
            patch("mofforge.build.config._autodetect_tobacco_data_dir", return_value=None),
            patch(
                "mofforge.build.config._fetch_tobacco_data_from_github",
                return_value=data_dir,
            ) as fetch,
        ):
            cfg = BuildConfig.load()
            assert cfg.resolve_tobacco_data_dir() == data_dir.resolve()
        fetch.assert_called_once()

    def test_configured_dir_skips_github_fetch(self, tmp_path):
        """An explicit data dir short-circuits before any network fetch."""
        data_dir = _make_data_dir(tmp_path / "local")
        with (
            patch("mofforge.build.config._find_toml", return_value=None),
            patch(
                "mofforge.build.config._fetch_tobacco_data_from_github",
            ) as fetch,
        ):
            cfg = BuildConfig.load(tobacco_data_dir=str(data_dir))
            assert cfg.resolve_tobacco_data_dir() == data_dir.resolve()
        fetch.assert_not_called()

    def test_resolve_raises_for_invalid(self, tmp_path):
        """resolve raises ConfigError for a directory with no data folders."""
        with patch("mofforge.build.config._find_toml", return_value=None):
            cfg = BuildConfig.load(tobacco_data_dir=str(tmp_path))
        with pytest.raises(ConfigError, match="Invalid TOBACCO data directory"):
            cfg.resolve_tobacco_data_dir()

    def test_resolve_succeeds(self, tmp_path):
        """resolve returns the path for a valid data directory."""
        data_dir = _make_data_dir(tmp_path)
        with patch("mofforge.build.config._find_toml", return_value=None):
            cfg = BuildConfig.load(tobacco_data_dir=str(data_dir))
        assert cfg.resolve_tobacco_data_dir() == data_dir.resolve()

    def test_resolve_autodetects(self, tmp_path):
        """resolve falls back to auto-detection when unconfigured."""
        data_dir = _make_data_dir(tmp_path)
        with (
            patch("mofforge.build.config._find_toml", return_value=None),
            patch(
                "mofforge.build.config._autodetect_tobacco_data_dir",
                return_value=data_dir,
            ),
        ):
            cfg = BuildConfig.load()
            assert cfg.resolve_tobacco_data_dir() == data_dir.resolve()

    def test_env_var_pormake_output(self, tmp_path):
        """MOFFORGE_PORMAKE_OUTPUT_DIR env var sets pormake output dir."""
        with (
            patch("mofforge.build.config._find_toml", return_value=None),
            patch.dict(os.environ, {"MOFFORGE_PORMAKE_OUTPUT_DIR": str(tmp_path)}),
        ):
            cfg = BuildConfig.load()
        assert cfg.pormake_output_dir == tmp_path


def _make_tobacco_tarball(dest: Path, tag: str = "data-v1") -> Path:
    """Build a GitHub-style tarball whose single root holds TOBACCO data dirs."""
    import tarfile

    root_name = f"tobacco_3.0-{tag}"
    src = dest / "src" / root_name
    _make_data_dir(src)
    (src / "template_database" / "pcu.cif").write_text("pcu\n")

    tarball = dest / f"{tag}.tar.gz"
    with tarfile.open(tarball, "w:gz") as tf:
        tf.add(src, arcname=root_name)
    return tarball


class TestFetchTobaccoDataFromGithub:
    """Tests for the GitHub download-and-cache fallback."""

    def test_downloads_extracts_and_caches(self, tmp_path):
        """First call downloads + extracts; second call is a cache hit."""
        from mofforge.build import config as cfgmod

        cache_root = tmp_path / "cache"
        tarball = _make_tobacco_tarball(tmp_path)

        calls = {"n": 0}

        def fake_urlretrieve(url, filename):
            calls["n"] += 1
            import shutil

            shutil.copy(tarball, filename)
            return filename, None

        with (
            patch.dict(os.environ, {"XDG_CACHE_HOME": str(cache_root)}),
            patch("urllib.request.urlretrieve", side_effect=fake_urlretrieve),
        ):
            result = cfgmod._fetch_tobacco_data_from_github()
            assert result is not None
            assert validate_tobacco_data_dir(result) == []
            assert (result / "template_database" / "pcu.cif").is_file()
            assert calls["n"] == 1

            # Second call: cached, no new download.
            result2 = cfgmod._fetch_tobacco_data_from_github()
            assert result2 == result
            assert calls["n"] == 1

    def test_returns_none_on_download_failure(self, tmp_path):
        """A failed download yields None (caller raises a helpful error)."""
        from mofforge.build import config as cfgmod

        cache_root = tmp_path / "cache"
        with (
            patch.dict(os.environ, {"XDG_CACHE_HOME": str(cache_root)}),
            patch(
                "urllib.request.urlretrieve",
                side_effect=OSError("network down"),
            ),
        ):
            assert cfgmod._fetch_tobacco_data_from_github() is None

    def test_respects_env_repo_and_tag(self, tmp_path):
        """Repo/tag env overrides shape the download URL and cache subdir."""
        from mofforge.build import config as cfgmod

        cache_root = tmp_path / "cache"
        tarball = _make_tobacco_tarball(tmp_path, tag="custom-tag")
        seen = {}

        def fake_urlretrieve(url, filename):
            seen["url"] = url
            import shutil

            shutil.copy(tarball, filename)
            return filename, None

        with (
            patch.dict(
                os.environ,
                {
                    "XDG_CACHE_HOME": str(cache_root),
                    "MOFFORGE_TOBACCO_DATA_REPO": "me/myfork",
                    "MOFFORGE_TOBACCO_DATA_TAG": "custom-tag",
                },
            ),
            patch("urllib.request.urlretrieve", side_effect=fake_urlretrieve),
        ):
            result = cfgmod._fetch_tobacco_data_from_github()
            assert result is not None
            assert result.name == "custom-tag"
            assert "me/myfork" in seen["url"]
            assert "custom-tag.tar.gz" in seen["url"]
