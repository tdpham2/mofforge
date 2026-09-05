#!/usr/bin/env python3
"""Construct a framework, modify a phenylene site, validate, and optionally render.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import MOIETIES, example_parser, output_dir, require, write_report


def main():
    parser = example_parser("build_and_modify", __doc__, seed=True)
    parser.add_argument("--render", action="store_true", help="Also render a PNG (vis extra).")
    args = parser.parse_args()
    require("pormake", "build")
    if args.render:
        require("playwright", "vis")
    out = output_dir(args)
    from mofforge import (
        MOFBuilder,
        find_pattern,
        fragment,
        infer_bonds,
        replace_pattern,
        validate_structure,
    )

    builder = MOFBuilder("pormake", output_dir=out / "build")
    builder.add_node("N108")
    builder.add_edge("E1")
    built = builder.build("pcu", output_dir=out / "build", accuracy=3)
    if not built.success:
        raise SystemExit(f"Construction failed: {built.errors}")
    parent = infer_bonds(built.crystal)
    query = fragment("2-!-p-phenylene.xyz", fragment_path=MOIETIES)
    replacement = fragment("2-nitro-p-phenylene.xyz", fragment_path=MOIETIES)
    match = find_pattern(query, parent)
    if not match.nb_locations():
        raise SystemExit("Constructed framework has no matching phenylene C-H motif.")
    child = replace_pattern(match, replacement, nb_loc=1, random_seed=args.seed)
    child.write_cif(out / "modified.cif")
    validation = validate_structure(infer_bonds(child))
    image = None
    if args.render:
        from mofforge import render_to_png

        image = render_to_png(
            child, output_file=str(out / "modified.png"), label_mode="none", show_unit_cell=True
        )
    write_report(
        out,
        seed=args.seed,
        build_outputs=built.output_paths,
        parent_atoms=parent.n_atoms,
        child_atoms=child.n_atoms,
        matched_locations=match.nb_locations(),
        image=image,
        validation=validation.to_dict(),
        provenance=child.provenance.to_dict(),
    )


if __name__ == "__main__":
    main()
