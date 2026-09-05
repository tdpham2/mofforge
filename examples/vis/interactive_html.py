#!/usr/bin/env python3
"""Export browser-viewable structures with labels and representations.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, example_parser, output_dir, write_report


def main():
    args = example_parser("interactive_html", __doc__).parse_args()
    out = output_dir(args)
    from mofforge import Crystal, build_html

    parent = Crystal.from_cif(CRYSTALS / "IRMOF-1.cif")
    paths = []
    for representation in ["ball_stick", "stick", "sphere"]:
        html = build_html(
            parent,
            representation=representation,
            label_mode="none",
            show_unit_cell=True,
            rotate=(15, 30, 0),
        )
        destination = out / f"{representation}.html"
        destination.write_text(html, encoding="utf-8")
        paths.append(destination)
    write_report(
        out,
        html_files=paths,
        viewing_requirement="A browser with network access to the 3Dmol.js CDN.",
    )


if __name__ == "__main__":
    main()
