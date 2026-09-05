#!/usr/bin/env python3
"""Sweep functional groups and coverage, retaining validation and errors.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import (
    BDC_SMILES,
    CRYSTALS,
    example_parser,
    output_dir,
    require,
    write_report,
)


def main():
    parser = example_parser("functionalization_campaign", __doc__, seed=True)
    parser.add_argument("--groups", nargs="+", default=["NH2", "F"])
    parser.add_argument("--coverages", nargs="+", type=float, default=[0.0, 0.5])
    args = parser.parse_args()
    require("rdkit", "chem")
    out = output_dir(args)
    from mofforge import run_campaign

    # Zero coverage is a control. Ranking uses validity and clashes, not adsorption.
    results = run_campaign(
        str(CRYSTALS / "IRMOF-1.cif"),
        BDC_SMILES,
        groups=args.groups,
        coverages=args.coverages,
        output_dir=str(out / "structures"),
        random_seed=args.seed,
    )
    rows = [{k: v for k, v in vars(result).items() if k != "crystal"} for result in results]
    write_report(out, seed=args.seed, results=rows)
    if any(result.error for result in results):
        raise SystemExit("One or more campaign members failed; see summary.json.")


if __name__ == "__main__":
    main()
