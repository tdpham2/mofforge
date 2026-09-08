#!/usr/bin/env python3
"""Prepare an amorphous porous organic polymer (POP) build.

Shows the parts of the POP workflow that run with only rdkit: inspect the
curated reactions, detect reactive sites on monomer SMILES, generate 3-D monomer
geometry, and size the packing box for a target density.  The final packing +
bonding + MD step (pysimm driving Packmol and LAMMPS) is described here and run
via ``mofforge polymerize`` / ``mofforge_polymerize`` once those binaries are
installed.  Run this file directly; see the topic README for the walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import example_parser, output_dir, require, write_report


def main():
    parser = example_parser("amorphous_pop", __doc__, seed=True)
    parser.add_argument(
        "--monomer",
        dest="monomers",
        action="append",
        help="Monomer SMILES (repeatable). Defaults to an imine CMP pair.",
    )
    parser.add_argument("--target-density", type=float, default=0.8)
    args = parser.parse_args()
    require("rdkit", "pop")
    out = output_dir(args)

    from mofforge.polymerize import Monomer
    from mofforge.polymerize.config import doctor
    from mofforge.polymerize.monomer import box_length_for_density, prepare_monomer
    from mofforge.polymerize.reactions import (
        available_reactions,
        detect_reactive_sites,
        reaction_name,
    )

    # A diamine + a dialdehyde condense into an imine (Schiff base) network.
    monomers = args.monomers or ["NCCN", "O=Cc1ccc(C=O)cc1"]

    # 1. Detect reactive sites on each monomer.
    site_report = {}
    for smiles in monomers:
        sites = detect_reactive_sites(smiles)
        site_report[smiles] = [{"atom_idx": s.atom_idx, "site_type": s.site_type} for s in sites]

    # 2. Confirm the two monomers carry a compatible pair of reactive groups.
    types = sorted({s["site_type"] for sites in site_report.values() for s in sites})
    linkage = reaction_name(types[0], types[-1]) if len(types) >= 2 else None

    # 3. Generate 3-D geometry and size a cubic box for the target density.
    counts = [20] * len(monomers)
    prepared = [
        prepare_monomer(Monomer(name=f"m{i}", source=s), out, random_seed=args.seed)
        for i, s in enumerate(monomers)
    ]
    box_length = box_length_for_density(prepared, counts, args.target_density)

    write_report(
        out,
        seed=args.seed,
        monomers=monomers,
        reactive_sites=site_report,
        detected_linkage=linkage,
        available_reactions=available_reactions(),
        monomer_masses={p.name: round(p.molar_mass, 3) for p in prepared},
        monomer_counts=dict(zip([p.name for p in prepared], counts, strict=True)),
        target_density=args.target_density,
        box_length_angstrom=round(box_length, 3),
        tool_availability=doctor(),
        next_step=(
            "Install Packmol and LAMMPS, then run: "
            "mofforge polymerize "
            + " ".join(f"-m '{s}'" for s in monomers)
            + f" --target-density {args.target_density} -o {out}"
        ),
    )


if __name__ == "__main__":
    main()
