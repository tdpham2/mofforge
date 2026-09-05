#!/usr/bin/env python3
"""Search chains, wildcards, and rings using the supported string syntax.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, example_parser, output_dir, write_report


def main():
    args = example_parser("string_pattern_search", __doc__).parse_args()
    out = output_dir(args)
    from mofforge import Crystal, infer_bonds, parse_smarts, smarts_search

    parent = infer_bonds(Crystal.from_cif(CRYSTALS / "IRMOF-1.cif"))
    patterns = ["[Zn]-[O]", "[Zn]-[O]-C", "O-C-O", "[Zn]-[*]", "C1-C-C-C-C-C-1", "C-H", "[Xe]-[Xe]"]
    rows = {}
    for pattern in patterns:
        graph = parse_smarts(pattern)
        matches = smarts_search(pattern, parent)
        rows[pattern] = {
            "query_atoms": graph.number_of_nodes(),
            "locations": matches.nb_locations(),
            "isomorphisms": matches.nb_isomorphisms(),
        }
    # This parser implements a subset, not the full RDKit SMARTS language.
    try:
        parse_smarts("C1-C")
    except ValueError as exc:
        invalid_pattern = str(exc)
    else:
        raise RuntimeError("The deliberately malformed pattern unexpectedly parsed.")
    write_report(out, patterns=rows, expected_parse_error=invalid_pattern)


if __name__ == "__main__":
    main()
