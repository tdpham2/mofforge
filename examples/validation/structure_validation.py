#!/usr/bin/env python3
"""Separate fatal geometry errors, warnings, and unperformed checks.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, example_parser, output_dir, write_report


def main():
    parser = example_parser("structure_validation", __doc__)
    parser.add_argument(
        "structure",
        type=Path,
        nargs="?",
        help="CIF to diagnose; omit to inspect all bundled structures.",
    )
    args = parser.parse_args()
    out = output_dir(args)
    import numpy as np

    from mofforge import Crystal, infer_bonds, validate_structure

    sources = [args.structure] if args.structure else sorted(CRYSTALS.glob("*.cif"))
    reports = {}
    for source in sources:
        parent = Crystal.from_cif(source)
        # Diagnose partial occupancy without forcing an ordered bonding graph.
        if parent.structure.is_ordered:
            parent = infer_bonds(parent)
        report = validate_structure(parent, check_charges=True)
        reports[source.stem] = report.to_dict()
        print(
            f"{source.stem}: geometry valid={report.is_valid}; "
            f"warnings={len(report.warnings)}; skipped={list(report.checks_skipped)}"
        )
    # A positive parser result alone is not a valid geometry.
    overlap = Crystal.from_xyz(["C", "C"], np.zeros((2, 3)), name="overlap_control")
    invalid = validate_structure(overlap)
    unchecked = validate_structure(Crystal.empty(), check_clashes=False)
    write_report(
        out,
        reports=reports,
        expected_invalid=invalid.to_dict(),
        expected_unchecked=unchecked.to_dict(),
    )
    if invalid.is_valid or unchecked.is_valid:
        raise SystemExit("A negative validation control unexpectedly passed.")


if __name__ == "__main__":
    main()
