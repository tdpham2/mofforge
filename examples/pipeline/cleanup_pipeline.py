#!/usr/bin/env python3
"""Chain disorder repair, explicit guest removal, desolvation, and validation.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, MOIETIES, example_parser, output_dir, write_report


def main():
    args = example_parser("cleanup_pipeline", __doc__, seed=True).parse_args()
    out = output_dir(args)
    from mofforge import Pipeline

    pipeline = Pipeline(
        CRYSTALS / "SIFSIX-2-Cu-i.cif", fragment_path=MOIETIES, random_seed=args.seed
    )
    pipeline.replace("disordered_ligand!.xyz", "4-pyridyl.xyz")
    pipeline.remove("acetylene.xyz")  # Complete disconnected molecules by default.
    # Preserve seven-atom SiF6 units even if inferred metal coordination misses them.
    pipeline.desolvate(min_atoms=7, keep_metal_containing=True)
    pipeline.validate()
    child = pipeline.build(name="cleaned_SIFSIX")
    child.write_cif(out / "cleaned.cif")
    write_report(
        out,
        seed=args.seed,
        child_atoms=child.n_atoms,
        validation=[r.to_dict() for r in pipeline.validation_reports],
        provenance=child.provenance.to_dict(),
    )


if __name__ == "__main__":
    main()
