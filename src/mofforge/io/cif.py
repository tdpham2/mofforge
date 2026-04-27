"""CIF file I/O using pymatgen."""

from __future__ import annotations

import warnings
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.io.cif import CifParser, CifWriter


def read_cif(filepath: str | Path, primitive: bool = False) -> Structure:
    """Read a CIF file and return a pymatgen Structure."""
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"CIF file not found: {filepath}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parser = CifParser(str(filepath))
        structures = parser.parse_structures(primitive=primitive)

    if not structures:
        raise ValueError(f"No structures found in CIF file: {filepath}")

    return structures[0]


def write_cif(structure: Structure, filepath: str | Path) -> None:
    """Write a pymatgen Structure to a CIF file."""
    filepath = Path(filepath)
    writer = CifWriter(structure)
    writer.write_file(str(filepath))
