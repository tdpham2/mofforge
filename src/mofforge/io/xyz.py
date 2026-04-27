"""XYZ file I/O with R-group tag support."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def read_xyz(filepath: str | Path) -> tuple[list[str], np.ndarray]:
    """Read an XYZ file, returning species labels and Cartesian coordinates."""
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"XYZ file not found: {filepath}")

    with open(filepath, encoding="utf-8") as f:
        lines = f.readlines()

    try:
        n_atoms = int(lines[0].strip())
    except (ValueError, IndexError) as e:
        raise ValueError(f"Invalid XYZ header in {filepath}: {e}") from e

    # Line 1 is the comment line (skipped)
    if len(lines) < 2 + n_atoms:
        raise ValueError(
            f"XYZ file {filepath} declares {n_atoms} atoms but has only "
            f"{len(lines)} lines (expected at least {2 + n_atoms})"
        )

    if n_atoms == 0:
        return [], np.empty((0, 3), dtype=np.float64)

    species = []
    coords = []
    for i in range(2, 2 + n_atoms):
        parts = lines[i].split()
        if len(parts) < 4:
            raise ValueError(
                f"Malformed XYZ data on line {i + 1} of {filepath}: "
                f"expected at least 4 fields, got {len(parts)}"
            )
        species.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    return species, np.array(coords, dtype=np.float64)


def write_xyz(
    species: list[str],
    coords: np.ndarray,
    filepath: str | Path,
    comment: str = "",
) -> None:
    """Write an XYZ file."""
    filepath = Path(filepath)
    n_atoms = len(species)
    if coords.shape != (n_atoms, 3):
        raise ValueError(f"coords shape {coords.shape} does not match {n_atoms} atoms")

    with open(filepath, "w", encoding="utf-8") as f:
        f.write(f"{n_atoms}\n")
        f.write(f"{comment}\n")
        for i in range(n_atoms):
            f.write(
                f"{species[i]:8s} {coords[i, 0]:14.5f} {coords[i, 1]:14.5f} {coords[i, 2]:14.5f}\n"
            )
