"""Periodic boundary condition utilities for fractional coordinate handling."""

from __future__ import annotations

import numpy as np
from pymatgen.core import Lattice


def nearest_image(dx: np.ndarray) -> np.ndarray:
    """Shift a fractional displacement vector to its nearest image.

    Applies the minimum image convention: subtracts round(dx) so that
    each component is in [-0.5, 0.5).

    Args:
        dx: Fractional displacement vector(s), shape (3,) or (N, 3).

    Returns:
        Shifted displacement vector(s) with the same shape.
    """
    return dx - np.round(dx)


def wrap_coords(xf: np.ndarray) -> np.ndarray:
    """Wrap fractional coordinates into [0, 1).

    Args:
        xf: Fractional coordinates, shape (N, 3) or (3,).

    Returns:
        Wrapped coordinates with the same shape.
    """
    return xf % 1.0


def min_image_distance(
    xf1: np.ndarray,
    xf2: np.ndarray,
    lattice: Lattice,
) -> float:
    """Compute the minimum image distance between two fractional coordinate points.

    Args:
        xf1: Fractional coordinates of point 1, shape (3,).
        xf2: Fractional coordinates of point 2, shape (3,).
        lattice: pymatgen Lattice for converting to Cartesian.

    Returns:
        Minimum image distance in Angstroms.
    """
    dx_frac = nearest_image(xf1 - xf2)
    dx_cart = lattice.get_cartesian_coords(dx_frac)
    return float(np.linalg.norm(dx_cart))


def is_cross_boundary(
    xf1: np.ndarray,
    xf2: np.ndarray,
    lattice: Lattice,
    bond_distance: float,
    atol: float = 0.1,
) -> bool:
    """Determine if a bond between two atoms crosses a periodic boundary.

    Compares the direct Cartesian distance to the known bond distance.
    If they differ significantly, the bond crosses the periodic boundary.

    Args:
        xf1: Fractional coordinates of atom 1.
        xf2: Fractional coordinates of atom 2.
        lattice: pymatgen Lattice.
        bond_distance: The actual bond distance (minimum image).
        atol: Tolerance in Angstroms.

    Returns:
        True if the bond crosses a periodic boundary.
    """
    dx_frac = xf1 - xf2
    dx_cart = lattice.get_cartesian_coords(dx_frac)
    direct_distance = float(np.linalg.norm(dx_cart))
    return not np.isclose(bond_distance, direct_distance, atol=atol)
