#!/usr/bin/env python3
"""Enumerate aromatic C-H sites and functionalize selected linker fractions.

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
    output_file,
    require,
    write_report,
)


def main():
    parser = example_parser("agent_functionalization", __doc__, seed=True)
    parser.add_argument("--group", default="NO2")
    parser.add_argument("--site", type=int, default=0)
    parser.add_argument("--coverage", type=float, default=0.5)
    parser.add_argument("--output")
    parser.add_argument("--campaign", action="store_true")
    args = parser.parse_args()
    require("rdkit", "chem")
    out = output_dir(args)
    from dataclasses import asdict

    from mofforge import available_groups, find_functionalizable_sites, functionalize, run_campaign

    sites = find_functionalizable_sites(BDC_SMILES)
    source = str(CRYSTALS / "IRMOF-1.cif")
    if args.campaign:
        results = run_campaign(
            source,
            BDC_SMILES,
            groups=["NH2", "F"],
            coverages=[0.25, 0.5],
            output_dir=str(out / "campaign"),
            random_seed=args.seed,
        )
    else:
        results = [
            functionalize(
                source,
                BDC_SMILES,
                args.group,
                sites=args.site,
                coverage=args.coverage,
                random_seed=args.seed,
                output_cif=str(output_file(args, "functionalized.cif")),
            )
        ]
    rows = [{k: v for k, v in vars(result).items() if k != "crystal"} for result in results]
    write_report(
        out,
        seed=args.seed,
        linker_smiles=BDC_SMILES,
        sites=[asdict(site) for site in sites],
        groups=available_groups(),
        results=rows,
    )
    if any(result.error for result in results):
        raise SystemExit("Functionalization failed; see results[].error in summary.json.")


if __name__ == "__main__":
    main()
