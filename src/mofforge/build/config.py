"""Build configuration: loads settings from toml, env vars, and kwargs (highest priority wins)."""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger("mofforge")

# Required entries inside a valid TOBACCO installation directory.
_TOBACCO_REQUIRED_FILES = ("tobacco.py", "configuration.py")
_TOBACCO_REQUIRED_DIRS = ("templates", "nodes", "edges")


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


def validate_tobacco_path(path: Path) -> list[str]:
    """Validate that *path* points to a usable TOBACCO installation."""
    errors: list[str] = []
    path = Path(path).resolve()

    if not path.exists():
        errors.append(f"Path does not exist: {path}")
        return errors

    if not path.is_dir():
        errors.append(f"Path is not a directory: {path}")
        return errors

    for fname in _TOBACCO_REQUIRED_FILES:
        if not (path / fname).is_file():
            errors.append(f"Missing required file: {path / fname}")

    for dname in _TOBACCO_REQUIRED_DIRS:
        if not (path / dname).is_dir():
            errors.append(f"Missing required directory: {path / dname}")

    return errors


@dataclass
class BuildConfig:
    """Configuration for MOF builder backends."""

    tobacco_path: Path | None = None
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

        tobacco_path: str | Path | None = tobacco_cfg.get("path")
        pormake_output_dir: str | Path | None = pormake_cfg.get("output_dir")

        # --- 2. Environment variables (medium priority) ---
        env_tobacco = os.environ.get("MOFFORGE_TOBACCO_PATH")
        if env_tobacco:
            tobacco_path = env_tobacco

        env_pormake_out = os.environ.get("MOFFORGE_PORMAKE_OUTPUT_DIR")
        if env_pormake_out:
            pormake_output_dir = env_pormake_out

        # --- 3. Explicit kwargs (highest priority) ---
        if "tobacco_path" in overrides and overrides["tobacco_path"] is not None:
            tobacco_path = overrides["tobacco_path"]
        if "pormake_output_dir" in overrides and overrides["pormake_output_dir"] is not None:
            pormake_output_dir = overrides["pormake_output_dir"]

        return cls(
            tobacco_path=Path(tobacco_path).resolve() if tobacco_path else None,
            pormake_output_dir=Path(pormake_output_dir) if pormake_output_dir else Path("."),
        )

    def resolve_tobacco_path(self) -> Path:
        """Return a validated TOBACCO path, or raise :class:`ConfigError`."""
        if self.tobacco_path is None:
            raise ConfigError("TOBACCO path is not configured")

        errors = validate_tobacco_path(self.tobacco_path)
        if errors:
            details = "; ".join(errors)
            raise ConfigError(
                f"Invalid TOBACCO installation at {self.tobacco_path}: {details}"
            )

        return self.tobacco_path

    def validate_tobacco(self) -> None:
        """Raise :class:`ConfigError` if ``tobacco_path`` is invalid."""
        self.resolve_tobacco_path()


def validate_pormake() -> list[str]:
    """Validate that pormake is importable and return any errors."""
    errors: list[str] = []
    try:
        import pormake  # noqa: F401
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
