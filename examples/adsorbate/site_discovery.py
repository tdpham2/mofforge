#!/usr/bin/env python3
"""Inspect void sites and the coordination heuristic for open metal sites.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, example_parser, output_dir, write_report


def main():
    args = example_parser("site_discovery", __doc__).parse_args()
    out = output_dir(args)
    from dataclasses import asdict

    import numpy as np
    from pymatgen.core import Lattice

    from mofforge import Crystal, find_adsorption_sites, infer_bonds

    parent = infer_bonds(Crystal.from_cif(CRYSTALS / "IRMOF-1.cif"))
    voids = find_adsorption_sites(parent, strategy="void", grid_spacing=1.0, max_sites=5)
    closed = find_adsorption_sites(parent, strategy="open_metal")
    # Explicitly illustrative CuO3 motif, not an experimental MOF or optimized complex.
    motif = Crystal.from_xyz(
        ["Cu", "O", "O", "O"],
        np.array([[10, 10, 10], [12, 10, 10], [10, 12, 10], [8, 10, 10]], dtype=float),
        lattice=Lattice.cubic(20),
        name="illustrative_CuO3",
    )
    motif = infer_bonds(motif)
    open_sites = find_adsorption_sites(motif, strategy="open_metal")
    motif.write_cif(out / "illustrative_CuO3.cif")
    # An empty site list is a legitimate result, not proof of nonporosity.
    unavailable = find_adsorption_sites(parent, min_distance=100, grid_spacing=1.0)
    write_report(
        out,
        void_sites=[asdict(site) for site in voids],
        framework_open_metal_count=len(closed),
        illustrative_open_sites=[asdict(site) for site in open_sites],
        impossible_clearance_count=len(unavailable),
    )


if __name__ == "__main__":
    main()
