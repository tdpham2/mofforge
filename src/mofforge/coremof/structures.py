"""Resolve CoRE MOF structure (CIF) files from a local structures directory.

The CoRE MOF metadata CSV bundled with mofforge contains *properties only*; the
actual crystal structure files (CIFs) are distributed separately (see the Zenodo
record referenced by :data:`mofforge.coremof.database.ZENODO_URL`). This module
maps a ``coreid`` / refcode to a CIF file inside a user-configured structures
directory so downstream tools (adsorbate placement, ASE, gRASPA) can load it.

The structures directory is resolved in order of priority:

1. Explicit *structures_dir* argument.
2. ``config.coremof_structures_path`` (set via :func:`mofforge.utils.config.set_paths`).
3. ``MOFFORGE_COREMOF_STRUCTURES_PATH`` environment variable.
4. ``[coremof] structures_path`` in ``mofforge.toml``.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

logger = logging.getLogger("mofforge")

__all__ = ["resolve_structure_path", "resolve_structures_dir"]


def _load_structures_path_from_toml() -> Path | None:
    """Try to read ``[coremof] structures_path`` from ``mofforge.toml``."""
    from mofforge.build.config import _find_toml, _load_toml

    toml_path = _find_toml()
    if toml_path is None:
        return None
    data = _load_toml(toml_path)
    section = data.get("coremof", {})
    raw = section.get("structures_path")
    return Path(raw) if raw else None


def resolve_structures_dir(structures_dir: str | Path | None = None) -> Path | None:
    """Resolve the CoRE MOF structures directory (without requiring it to exist).

    Returns ``None`` if no directory is configured by any mechanism.
    """
    if structures_dir is not None:
        return Path(structures_dir)

    from mofforge.utils.config import config

    if config.coremof_structures_path is not None:
        return config.coremof_structures_path

    env = os.environ.get("MOFFORGE_COREMOF_STRUCTURES_PATH")
    if env:
        return Path(env)

    return _load_structures_path_from_toml()


def _candidate_names(identifier: str) -> list[str]:
    """Generate candidate CIF filenames for a coreid / refcode.

    CoRE MOF CIFs are commonly named ``<coreid>.cif`` (e.g.
    ``ABAVIJ_clean.cif``). We also try the bare identifier and a couple of
    common processing suffixes seen across CoRE MOF releases.
    """
    ident = identifier.strip()
    stem = ident[:-4] if ident.endswith(".cif") else ident
    names = [f"{stem}.cif"]
    # Common processed variants when only a base refcode is known.
    for suffix in ("_clean", "_ASR", "_FSR", "_ION"):
        names.append(f"{stem}{suffix}.cif")
    # De-duplicate while preserving order.
    seen: set[str] = set()
    unique: list[str] = []
    for n in names:
        if n not in seen:
            seen.add(n)
            unique.append(n)
    return unique


def resolve_structure_path(
    identifier: str,
    structures_dir: str | Path | None = None,
) -> Path | None:
    """Return the local CIF path for a CoRE MOF ``coreid`` or refcode.

    Parameters
    ----------
    identifier : str
        A CoRE MOF ``coreid`` (preferred) or refcode.
    structures_dir : str or Path, optional
        Directory containing CoRE MOF CIF files. Resolved via
        :func:`resolve_structures_dir` when omitted.

    Returns
    -------
    Path or None
        The resolved CIF path if a matching file is found, else ``None``.
    """
    base = resolve_structures_dir(structures_dir)
    if base is None or not base.exists():
        return None

    # 1. Direct candidate filenames at the top level.
    for name in _candidate_names(identifier):
        candidate = base / name
        if candidate.is_file():
            return candidate

    # 2. Recursive search (datasets are sometimes split into subdirectories).
    ident = identifier.strip()
    stem = ident[:-4] if ident.endswith(".cif") else ident
    matches = sorted(base.rglob(f"{stem}*.cif"))
    if matches:
        # Prefer an exact stem match, otherwise the first lexical match.
        for m in matches:
            if m.stem == stem:
                return m
        return matches[0]

    return None
