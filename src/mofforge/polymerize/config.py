"""Packmol executable resolution and feature detection; no simulation backend."""

from __future__ import annotations

import importlib.util
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from mofforge.polymerize.base import reject_removed


class ConfigError(ValueError):
    """Invalid POP configuration or unavailable Packmol executable."""


def _find_toml():
    for path in (Path.cwd() / "mofforge.toml", Path.home() / ".mofforge.toml"):
        if path.is_file():
            return path
    return None


def _load_toml(path):
    try:
        import tomllib
    except ImportError:
        import tomli as tomllib
    with Path(path).open("rb") as handle:
        return tomllib.load(handle)


@dataclass
class PopBuildConfig:
    packmol_bin: str | Path | None = None

    @classmethod
    def load(cls, **overrides):
        path = _find_toml()
        config = _load_toml(path).get("backends", {}).get("pop", {}) if path else {}
        reject_removed(config)
        reject_removed(overrides)
        if os.environ.get("MOFFORGE_LAMMPS_BIN"):
            raise ConfigError("MOFFORGE_LAMMPS_BIN was removed; configure LAMMPS in MatKit.")
        unknown = (set(config) | set(overrides)) - {"packmol_bin"}
        if unknown:
            raise ConfigError(f"Unknown POP backend settings: {sorted(unknown)}")
        binary = overrides.get("packmol_bin")
        return cls(config.get("packmol_bin") if binary is None else binary)

    def resolve_packmol_binary(self):
        configured = self.packmol_bin or os.environ.get("MOFFORGE_PACKMOL_BIN")
        if configured:
            found = shutil.which(str(configured))
            if not found:
                raise ConfigError(
                    f"Configured Packmol executable not found or not executable: {configured}"
                )
        else:
            found = shutil.which("packmol")
            if not found:
                found = shutil.which(str(Path(sys.prefix) / "bin" / "packmol"))
        if not found:
            raise ConfigError(
                "Packmol executable not found. Install Packmol >=20.15.0; "
                "set MOFFORGE_PACKMOL_BIN or [backends.pop].packmol_bin."
            )
        return Path(found).resolve()


def probe_packmol(binary):
    # Packmol has no portable --version flag; even empty stdin prints its banner.
    result = subprocess.run([str(binary)], input="", capture_output=True, text=True, timeout=10)
    transcript = result.stdout + result.stderr
    match = re.search(r"Version\s+(\d+)\.(\d+)\.(\d+)", transcript, re.IGNORECASE)
    if not match or tuple(map(int, match.groups())) < (20, 15, 0):
        raise ConfigError(
            "Packmol >=20.15.0 with orthorhombic pbc support is required. "
            f"Version probe output: {transcript[:1000]}"
        )
    return ".".join(match.groups()), transcript


def doctor():
    report = {"rdkit": {"available": importlib.util.find_spec("rdkit") is not None}}
    try:
        binary = PopBuildConfig.load().resolve_packmol_binary()
        version, _ = probe_packmol(binary)
        report["packmol"] = {"available": True, "path": str(binary), "version": version}
    except (ValueError, TypeError, OSError, subprocess.SubprocessError) as exc:
        report["packmol"] = {"available": False, "error": str(exc)}
    return report
