#!/usr/bin/env python3
"""Measure periodic distances in a skewed cell and inspect wrapping.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import example_parser, output_dir, write_report


def main():
    args = example_parser("periodic_geometry", __doc__).parse_args()
    out = output_dir(args)
    import numpy as np
    from pymatgen.core import Lattice

    from mofforge import Crystal, min_image_distance

    lattice = Lattice([[8, 0, 0], [6, 5, 0], [1, 2, 9]])
    fractional = np.array([[0.98, 0.5, 0.5], [0.02, 0.5, 0.5]])
    crystal = Crystal.from_xyz(
        ["He", "He"],
        lattice.get_cartesian_coords(fractional),
        lattice=lattice,
        name="periodic_distance_demo",
    )
    direct = np.linalg.norm(crystal.cart_coords[1] - crystal.cart_coords[0])
    # Use the full lattice metric; Cartesian separation alone misses the neighbor.
    periodic = lattice.get_distance_and_image(fractional[0], fractional[1])[0]
    shifted = crystal.copy()
    shifted.set_frac_coords(fractional + np.array([1, -1, 0]))
    wrapped = shifted.wrap()
    wrapped.write_cif(out / "wrapped.cif")
    write_report(
        out,
        direct_distance_angstrom=direct,
        periodic_distance_angstrom=periodic,
        library_distance_angstrom=min_image_distance(fractional[0], fractional[1], lattice),
        original_fractional=fractional,
        wrapped_fractional=wrapped.frac_coords,
    )


if __name__ == "__main__":
    main()
