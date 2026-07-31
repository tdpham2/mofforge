"""Packaged data assets shipped with mofforge.

Currently this holds the curated moiety (fragment) library under
``moieties/``: a set of ``!``-anchor-tagged XYZ fragments used as query and
replacement patterns for find-and-replace functionalization.  These are the
same fragments as in ``examples/data/moieties`` but installed with the package
so they are discoverable at runtime (including by the MCP server).
"""

from __future__ import annotations

from pathlib import Path

# Directory containing the packaged moiety XYZ fragments.
MOIETIES_DIR = Path(__file__).parent / "moieties"


def moieties_dir() -> Path:
    """Return the path to the packaged moiety library directory."""
    return MOIETIES_DIR


def list_moieties() -> list[str]:
    """Return the sorted names of packaged moiety XYZ files."""
    if not MOIETIES_DIR.is_dir():
        return []
    return sorted(p.name for p in MOIETIES_DIR.glob("*.xyz"))


def moiety_path(name: str) -> Path:
    """Return the full path to a packaged moiety by file name.

    Raises
    ------
    FileNotFoundError
        If no packaged moiety with that name exists.
    """
    path = MOIETIES_DIR / name
    if not path.is_file():
        raise FileNotFoundError(
            f"No packaged moiety named {name!r}. "
            f"Available: {', '.join(list_moieties())}"
        )
    return path
