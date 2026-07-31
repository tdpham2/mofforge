"""Build configuration: loads settings from toml, env vars, and kwargs (highest priority wins)."""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger("mofforge")

# A usable TOBACCO data directory holds at least one of these subdirectories
# of topology / building-block CIFs.  ``tobacco3`` ships as an importable
# package but does *not* bundle these data folders, so mofforge locates them
# separately (auto-detected next to the installed package, or configured).
_TOBACCO_DATA_DIRS = (
    "template_database",
    "templates",
    "nodes_database",
    "edges_database",
)

# The ``tobacco3`` package ships code only; its topology / building-block data
# folders live in the source repo.  When the data cannot be found locally,
# mofforge fetches a pinned tarball of that repo from GitHub and caches it.  The
# repo is owned by the mofforge author, so availability and versioning are under
# our control.  The data is GPLv3 (ToBaCCo, (c) 2019 Ryther Anderson); mofforge
# only fetches and uses it at runtime, which is aggregation, so mofforge itself
# stays MIT.  Both values are overridable via environment variables.
_TOBACCO_DATA_REPO = "tdpham2/tobacco_3.0"
_TOBACCO_DATA_TAG = "data-v1"


class ConfigError(Exception):
    """Raised when build configuration is invalid or incomplete."""


def _load_toml(path: Path) -> dict[str, Any]:
    """Load a TOML file and return it as a nested dict."""
    try:
        import tomllib  # Python 3.11+
    except ModuleNotFoundError:
        try:
            import tomli as tomllib  # type: ignore[no-redef]
        except ModuleNotFoundError:
            logger.debug("Neither tomllib nor tomli available; cannot read %s", path)
            return {}

    try:
        with open(path, "rb") as fh:
            return tomllib.load(fh)
    except Exception as exc:
        logger.warning("Failed to read %s: %s", path, exc)
        return {}


def _find_toml() -> Path | None:
    """Search for ``mofforge.toml`` in the standard locations."""
    candidates = [
        Path.cwd() / "mofforge.toml",
        Path.home() / ".mofforge.toml",
    ]
    for p in candidates:
        if p.is_file():
            logger.debug("Found config file: %s", p)
            return p
    return None


def validate_tobacco() -> list[str]:
    """Validate that the importable ``tobacco3`` package is available."""
    errors: list[str] = []
    try:
        import tobacco3
    except ImportError:
        errors.append(
            "tobacco3 is not installed.  Install it with:  "
            "pip install 'tobacco3 @ git+https://github.com/tdpham2/tobacco_3.0.git'"
        )
        return errors

    if not hasattr(tobacco3, "generate_mof"):
        errors.append(
            "tobacco3 is installed but does not expose generate_mof; "
            "an incompatible (pre-3.1) version may be present."
        )
    return errors


def validate_tobacco_data_dir(path: Path) -> list[str]:
    """Validate that *path* holds TOBACCO topology / building-block data."""
    errors: list[str] = []
    path = Path(path).resolve()

    if not path.exists():
        errors.append(f"Path does not exist: {path}")
        return errors

    if not path.is_dir():
        errors.append(f"Path is not a directory: {path}")
        return errors

    if not any((path / d).is_dir() for d in _TOBACCO_DATA_DIRS):
        errors.append(
            f"No TOBACCO data directories found under {path}; "
            f"expected at least one of: {', '.join(_TOBACCO_DATA_DIRS)}"
        )

    return errors


def _autodetect_tobacco_data_dir() -> Path | None:
    """Locate the TOBACCO data directory next to the installed ``tobacco3`` package.

    ``tobacco3`` lives at ``<repo>/tobacco3/``; the data folders
    (``template_database``, ``nodes_database``, ...) sit alongside it at
    ``<repo>/``.  Return that repo root if it looks valid, else *None*.
    """
    try:
        import tobacco3
    except ImportError:
        return None

    pkg_file = getattr(tobacco3, "__file__", None)
    if not pkg_file:
        return None

    repo_root = Path(pkg_file).resolve().parent.parent
    if not validate_tobacco_data_dir(repo_root):
        return repo_root
    return None


def _tobacco_cache_root() -> Path:
    """Return the base cache directory for downloaded TOBACCO data."""
    base = os.environ.get("XDG_CACHE_HOME")
    cache_home = Path(base) if base else Path.home() / ".cache"
    return cache_home / "mofforge" / "tobacco-data"


def _fetch_tobacco_data_from_github() -> Path | None:
    """Download and cache the TOBACCO data repo, returning the local data dir.

    Fetches a pinned tarball of ``_TOBACCO_DATA_REPO`` at ``_TOBACCO_DATA_TAG``
    (both overridable via ``MOFFORGE_TOBACCO_DATA_REPO`` /
    ``MOFFORGE_TOBACCO_DATA_TAG``) into ``~/.cache/mofforge/tobacco-data/<tag>``.
    A cached copy that already validates is reused without re-downloading.
    Returns the cached directory, or *None* if the download/extract fails (the
    caller then raises a helpful :class:`ConfigError`).
    """
    import shutil
    import tarfile
    import tempfile
    import urllib.request

    repo = os.environ.get("MOFFORGE_TOBACCO_DATA_REPO", _TOBACCO_DATA_REPO)
    tag = os.environ.get("MOFFORGE_TOBACCO_DATA_TAG", _TOBACCO_DATA_TAG)

    dest = _tobacco_cache_root() / tag
    # Cache hit: a previously extracted, still-valid copy.
    if dest.is_dir() and not validate_tobacco_data_dir(dest):
        logger.debug("Using cached TOBACCO data at %s", dest)
        return dest

    url = f"https://github.com/{repo}/archive/refs/tags/{tag}.tar.gz"
    logger.warning(
        "TOBACCO data not found locally; downloading %s@%s (~27 MB, one-time) "
        "from %s into %s",
        repo,
        tag,
        url,
        dest,
    )

    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        with tempfile.TemporaryDirectory(dir=dest.parent) as tmp:
            tmp_path = Path(tmp)
            tarball = tmp_path / "data.tar.gz"
            urllib.request.urlretrieve(url, tarball)

            extract_dir = tmp_path / "extract"
            extract_dir.mkdir()
            with tarfile.open(tarball) as tf:
                # filter="data" rejects unsafe members (absolute paths, links
                # escaping the target) and silences the 3.14 deprecation.
                tf.extractall(extract_dir, filter="data")

            # A GitHub tag tarball extracts to a single ``<repo>-<tag>/`` root.
            roots = [p for p in extract_dir.iterdir() if p.is_dir()]
            data_root = roots[0] if len(roots) == 1 else extract_dir
            if validate_tobacco_data_dir(data_root):
                logger.warning(
                    "Downloaded archive from %s is missing expected TOBACCO "
                    "data directories",
                    url,
                )
                return None

            # Atomic-ish move into place (temp dir is on the same filesystem).
            if dest.exists():
                shutil.rmtree(dest)
            shutil.move(str(data_root), str(dest))
    except Exception as exc:
        logger.warning("Failed to fetch TOBACCO data from %s: %s", url, exc)
        return None

    return dest if dest.is_dir() else None


@dataclass
class BuildConfig:
    """Configuration for MOF builder backends."""

    tobacco_data_dir: Path | None = None
    pormake_output_dir: Path = field(default_factory=lambda: Path("."))

    @classmethod
    def load(cls, **overrides: Any) -> BuildConfig:
        """Build a :class:`BuildConfig` by merging all configuration sources."""
        # --- 1. Start from TOML file (lowest priority) ---
        toml_path = _find_toml()
        toml_data: dict[str, Any] = {}
        if toml_path is not None:
            toml_data = _load_toml(toml_path)
            logger.debug("Loaded build config from %s", toml_path)

        backends = toml_data.get("backends", {})
        tobacco_cfg = backends.get("tobacco", {})
        pormake_cfg = backends.get("pormake", {})

        # ``data_dir`` is the new key; ``path`` is the legacy alias.
        tobacco_data_dir: str | Path | None = tobacco_cfg.get("data_dir") or tobacco_cfg.get(
            "path"
        )
        pormake_output_dir: str | Path | None = pormake_cfg.get("output_dir")

        # --- 2. Environment variables (medium priority) ---
        # MOFFORGE_TOBACCO_DATA is preferred; MOFFORGE_TOBACCO_PATH is the
        # legacy alias kept for back-compat with the old file-based backend.
        env_tobacco = os.environ.get("MOFFORGE_TOBACCO_DATA") or os.environ.get(
            "MOFFORGE_TOBACCO_PATH"
        )
        if env_tobacco:
            tobacco_data_dir = env_tobacco

        env_pormake_out = os.environ.get("MOFFORGE_PORMAKE_OUTPUT_DIR")
        if env_pormake_out:
            pormake_output_dir = env_pormake_out

        # --- 3. Explicit kwargs (highest priority) ---
        # Accept both ``tobacco_data_dir`` and the legacy ``tobacco_path``.
        override_data_dir = overrides.get("tobacco_data_dir") or overrides.get("tobacco_path")
        if override_data_dir is not None:
            tobacco_data_dir = override_data_dir
        if "pormake_output_dir" in overrides and overrides["pormake_output_dir"] is not None:
            pormake_output_dir = overrides["pormake_output_dir"]

        return cls(
            tobacco_data_dir=Path(tobacco_data_dir).resolve() if tobacco_data_dir else None,
            pormake_output_dir=Path(pormake_output_dir) if pormake_output_dir else Path("."),
        )

    def resolve_tobacco_data_dir(self) -> Path:
        """Return a validated TOBACCO data directory, or raise :class:`ConfigError`.

        Resolution order: an explicitly configured ``tobacco_data_dir``, then
        auto-detection next to the installed ``tobacco3`` package, then a
        one-time download of the pinned data tarball from GitHub (cached under
        ``~/.cache/mofforge/tobacco-data``).
        """
        candidate = self.tobacco_data_dir
        if candidate is None:
            candidate = _autodetect_tobacco_data_dir()
        if candidate is None:
            candidate = _fetch_tobacco_data_from_github()

        if candidate is None:
            raise ConfigError(
                "TOBACCO data directory is not configured, could not be "
                "auto-detected next to the installed tobacco3 package, and could "
                "not be downloaded from GitHub.  Set MOFFORGE_TOBACCO_DATA to a "
                "directory containing template_database/, nodes_database/, and "
                "edges_database/, or check your network connection."
            )

        errors = validate_tobacco_data_dir(candidate)
        if errors:
            details = "; ".join(errors)
            raise ConfigError(f"Invalid TOBACCO data directory at {candidate}: {details}")

        return Path(candidate).resolve()

    def validate_tobacco(self) -> None:
        """Raise :class:`ConfigError` if ``tobacco3`` or its data dir is unusable."""
        errors = validate_tobacco()
        if errors:
            raise ConfigError("; ".join(errors))
        self.resolve_tobacco_data_dir()


def validate_pormake() -> list[str]:
    """Validate that pormake is importable and return any errors."""
    errors: list[str] = []
    try:
        import pormake
    except ImportError:
        errors.append(
            "pormake is not installed.  Install it with:  pip install pormake"
        )
        return errors

    try:
        db = pormake.Database()
        _ = db.topology_list
    except Exception as exc:
        errors.append(f"pormake is installed but database init failed: {exc}")

    return errors
