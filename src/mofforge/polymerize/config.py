"""Configuration for the POP (porous organic polymer) subsystem.

Loads settings from ``mofforge.toml``, environment variables, and kwargs
(highest priority wins), and resolves the external binaries the pysimm engine
drives — **Packmol** and **LAMMPS** — the same way :mod:`mofforge.build.config`
resolves the TOBACCO data directory: explicit config, then environment, then a
PATH search.  The binaries are never installed by mofforge; a missing one raises
an actionable :class:`ConfigError` only when a step actually needs it.
"""

from __future__ import annotations

import logging
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

logger = logging.getLogger("mofforge")

# Candidate executable names, in preference order.
_PACKMOL_NAMES = ("packmol",)
_LAMMPS_NAMES = ("lmp", "lmp_serial", "lmp_mpi", "lammps")


class ConfigError(Exception):
    """Raised when POP configuration is invalid or a required tool is missing."""


def _load_toml(path: Path) -> dict[str, Any]:
    """Load a TOML file and return it as a nested dict (mirrors build.config)."""
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
    for p in (Path.cwd() / "mofforge.toml", Path.home() / ".mofforge.toml"):
        if p.is_file():
            return p
    return None


def validate_pysimm() -> list[str]:
    """Validate that the importable ``pysimm`` package is available."""
    import importlib.util

    if importlib.util.find_spec("pysimm") is None:
        return ["pysimm is not installed.  Install it with:  pip install 'mofforge[pop]'"]
    return []


def _resolve_binary(
    names: tuple[str, ...],
    configured: str | Path | None,
    env_var: str,
    label: str,
) -> Path:
    """Resolve an external executable: config -> env -> PATH search.

    Raises :class:`ConfigError` with an actionable message if none is found.
    """
    # 1. Explicit configuration (kwargs / toml), highest priority.
    candidates: list[str] = []
    if configured:
        candidates.append(str(configured))
    # 2. Environment variable.
    env_value = os.environ.get(env_var)
    if env_value:
        candidates.append(env_value)

    for candidate in candidates:
        # A configured value may be an absolute path or a bare name on PATH.
        resolved = shutil.which(candidate) or (
            candidate if Path(candidate).is_file() else None
        )
        if resolved:
            return Path(resolved).resolve()

    # 3. PATH search over known executable names.
    for name in names:
        found = shutil.which(name)
        if found:
            return Path(found).resolve()

    raise ConfigError(
        f"{label} executable not found.  Install {label} and ensure it is on "
        f"PATH, set {env_var} to its path, or configure it in mofforge.toml "
        f"([backends.pop] section).  Searched names: {', '.join(names)}."
    )


@dataclass
class PopBuildConfig:
    """Resolved locations of the external tools the pysimm engine drives."""

    packmol_bin: str | Path | None = None
    lammps_bin: str | Path | None = None

    @classmethod
    def load(cls, **overrides: Any) -> PopBuildConfig:
        """Merge TOML file, environment, and kwargs (kwargs win)."""
        toml_path = _find_toml()
        pop_cfg: dict[str, Any] = {}
        if toml_path is not None:
            pop_cfg = _load_toml(toml_path).get("backends", {}).get("pop", {})

        packmol = overrides.get("packmol_bin") or pop_cfg.get("packmol_bin")
        lammps = overrides.get("lammps_bin") or pop_cfg.get("lammps_bin")
        return cls(packmol_bin=packmol, lammps_bin=lammps)

    def resolve_packmol_binary(self) -> Path:
        """Return the Packmol executable path, or raise :class:`ConfigError`."""
        return _resolve_binary(
            _PACKMOL_NAMES, self.packmol_bin, "MOFFORGE_PACKMOL_BIN", "Packmol"
        )

    def resolve_lammps_binary(self) -> Path:
        """Return the LAMMPS executable path, or raise :class:`ConfigError`."""
        return _resolve_binary(
            _LAMMPS_NAMES, self.lammps_bin, "MOFFORGE_LAMMPS_BIN", "LAMMPS"
        )


def doctor() -> dict[str, Any]:
    """Report which POP dependencies and binaries are available.

    Useful for a CLI ``doctor``/status command: never raises, returns a plain
    dict describing pysimm, Packmol, and LAMMPS availability.
    """
    cfg = PopBuildConfig.load()
    report: dict[str, Any] = {}

    report["pysimm"] = {"available": not validate_pysimm(), "errors": validate_pysimm()}

    for label, resolver in (
        ("packmol", cfg.resolve_packmol_binary),
        ("lammps", cfg.resolve_lammps_binary),
    ):
        try:
            report[label] = {"available": True, "path": str(resolver())}
        except ConfigError as exc:
            report[label] = {"available": False, "error": str(exc)}
    return report
