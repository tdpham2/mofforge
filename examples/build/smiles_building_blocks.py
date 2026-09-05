#!/usr/bin/env python3
"""Convert a SMILES linker into both supported building-block formats.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import BDC_SMILES, example_parser, output_dir, require, write_report


def main():
    parser = example_parser("smiles_building_blocks", __doc__, seed=True)
    parser.add_argument("--smiles", default=BDC_SMILES)
    parser.add_argument(
        "--mode", choices=["auto", "carboxylate", "direct", "carboxylic"], default="auto"
    )
    args = parser.parse_args()
    require("rdkit", "chem")
    out = output_dir(args)
    from dataclasses import asdict

    from mofforge.build.smiles_to_bb import (
        detect_connection_points,
        smiles_to_pormake_edge_xyz,
        smiles_to_tobacco_edge_cif,
    )

    detected = detect_connection_points(args.smiles)
    tobacco = out / "linker.cif"
    pormake = out / "linker.xyz"
    smiles_to_tobacco_edge_cif(args.smiles, tobacco, mode=args.mode, random_seed=args.seed)
    smiles_to_pormake_edge_xyz(args.smiles, pormake, mode=args.mode, random_seed=args.seed)
    write_report(
        out,
        seed=args.seed,
        auto_detection=asdict(detected),
        requested_mode=args.mode,
        tobacco_cif=tobacco,
        pormake_xyz=pormake,
    )


if __name__ == "__main__":
    main()
