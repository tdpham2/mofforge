"""Periodic boundary condition utilities for fractional coordinate handling."""

from __future__ import annotations

import numpy as np
from pymatgen.core import Lattice


def nearest_image(dx: np.ndarray) -> np.ndarray:
    """Shift a fractional displacement vector to its nearest image in [-0.5, 0.5)."""
    return dx - np.round(dx)


def wrap_coords(xf: np.ndarray) -> np.ndarray:
    """Wrap fractional coordinates into [0, 1)."""
    return xf % 1.0


def min_image_distance(
    xf1: np.ndarray,
    xf2: np.ndarray,
    lattice: Lattice,
) -> float:
    """Compute the minimum image distance between two fractional coordinate points."""
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
    """Determine if a bond between two atoms crosses a periodic boundary."""
    # Wrap to [0, 1) so the direct distance is well-defined within one unit cell,
    # even when inputs have not been wrapped yet (e.g. after reassemble()).
    dx_frac = wrap_coords(xf1) - wrap_coords(xf2)
    dx_cart = lattice.get_cartesian_coords(dx_frac)
    direct_distance = float(np.linalg.norm(dx_cart))
    return not np.isclose(bond_distance, direct_distance, atol=atol)
