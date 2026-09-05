#!/usr/bin/env python3
"""Inspect default and explicit distance cutoffs without changing global settings.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, example_parser, output_dir, write_report


def main():
    args = example_parser("bonding_rules", __doc__).parse_args()
    out = output_dir(args)
    from mofforge import BondingRule, Crystal, infer_bonds, remove_bonds
    from mofforge.core.bonding import default_bonding_rules

    raw = Crystal.from_cif(CRYSTALS / "IRMOF-1.cif")
    baseline = infer_bonds(raw, periodic=True)
    # A deliberately short Zn-O cutoff illustrates how inferred connectivity changes.
    rules = [
        rule for rule in default_bonding_rules() if {rule.species_i, rule.species_j} != {"Zn", "O"}
    ]
    restricted = infer_bonds(raw, bonding_rules=[*rules, BondingRule("Zn", "O", 1.0)])
    nonperiodic = infer_bonds(raw, periodic=False)
    unbonded = remove_bonds(baseline)
    write_report(
        out,
        default_bonds=baseline.n_bonds,
        restricted_bonds=restricted.n_bonds,
        nonperiodic_bonds=nonperiodic.n_bonds,
        removed_bonds=unbonded.n_bonds,
        default_zinc_coordination=baseline.coordination_number(raw.species.index("Zn")),
        restricted_zinc_coordination=restricted.coordination_number(raw.species.index("Zn")),
    )


if __name__ == "__main__":
    main()
