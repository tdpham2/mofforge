#!/usr/bin/env python3
"""Place built-in and custom guests with seeded or fixed orientations.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, example_parser, output_dir, write_report


def main():
    parser = example_parser("placement", __doc__, seed=True)
    parser.add_argument("--count", type=int, default=2)
    args = parser.parse_args()
    out = output_dir(args)
    from mofforge import (
        Crystal,
        available_molecules,
        find_adsorption_sites,
        get_molecule,
        infer_bonds,
        place_adsorbate,
    )

    parent = infer_bonds(Crystal.from_cif(CRYSTALS / "IRMOF-1.cif"))
    # The actual number may be smaller after filtering distances between site centers.
    result = place_adsorbate(
        parent,
        "CO2",
        n_adsorbates=args.count,
        grid_spacing=1.0,
        min_intermolecular_dist=3.0,
        random_seed=args.seed,
    )
    result.crystal.write_cif(out / "CO2_loaded.cif")
    # Custom molecules use a Crystal and an explicit site; each placement starts from parent.
    species, coordinates = get_molecule("H2O")
    custom = Crystal.from_xyz(species, coordinates, name="custom_water")
    site = find_adsorption_sites(parent, grid_spacing=1.0, max_sites=1)[0]
    water = place_adsorbate(parent, custom, site=site, orient="fixed", random_seed=args.seed)
    water.crystal.write_cif(out / "water_loaded.cif")
    write_report(
        out,
        seed=args.seed,
        available=available_molecules(),
        requested=args.count,
        placed=result.n_adsorbates,
        parent_atoms=parent.n_atoms,
        loaded_atoms=result.crystal.n_atoms,
        adsorbate_indices=result.adsorbate_indices,
        validation=result.validation.to_dict(),
        custom_validation=water.validation.to_dict(),
    )


if __name__ == "__main__":
    main()
