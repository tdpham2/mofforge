#!/usr/bin/env python3
"""Inspect provenance and reproduce a replacement using its recorded seed.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, MOIETIES, example_parser, output_dir, write_report


def main():
    args = example_parser("provenance_replay", __doc__, seed=True).parse_args()
    out = output_dir(args)
    import numpy as np

    from mofforge import Crystal, Provenance, find_pattern, fragment, infer_bonds, replace_pattern

    parent = infer_bonds(Crystal.from_cif(CRYSTALS / "IRMOF-1.cif"))
    query = fragment("2-!-p-phenylene.xyz", fragment_path=MOIETIES)
    replacement = fragment("2-nitro-p-phenylene.xyz", fragment_path=MOIETIES)
    matches = find_pattern(query, parent)
    child = replace_pattern(matches, replacement, nb_loc=3, random_seed=args.seed)
    child.provenance.to_json(out / "provenance.json")
    record = Provenance.from_json(out / "provenance.json")
    # Replay starts from the original parent and fragments, not the exported child.
    replay = replace_pattern(
        matches, replacement, nb_loc=3, random_seed=record.parameters["random_seed"]
    )
    np.testing.assert_allclose(child.frac_coords, replay.frac_coords, atol=1e-10)
    if child.species != replay.species:
        raise RuntimeError("Replay changed species.")
    child.write_cif(out / "original.cif")
    replay.write_cif(out / "replayed.cif")
    write_report(
        out,
        same_geometry=True,
        seed=record.parameters["random_seed"],
        locations=record.parameters["locations"],
        provenance=record.to_dict(),
    )


if __name__ == "__main__":
    main()
