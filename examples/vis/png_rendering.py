#!/usr/bin/env python3
"""Render PNG files with synchronous and asynchronous APIs.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, example_parser, output_dir, require, write_report


def main():
    args = example_parser("png_rendering", __doc__).parse_args()
    require("playwright", "vis")
    out = output_dir(args)
    import asyncio

    from mofforge import Crystal, async_render_to_png, render_to_png

    parent = Crystal.from_cif(CRYSTALS / "IRMOF-1.cif")
    # Run the synchronous API outside an event loop.
    sync_path = render_to_png(
        parent,
        output_file=str(out / "sync.png"),
        show_unit_cell=True,
        label_mode="none",
        width=640,
        height=480,
    )
    async_path = asyncio.run(
        async_render_to_png(
            parent,
            output_file=str(out / "async.png"),
            representation="stick",
            show_unit_cell=True,
            label_mode="none",
            width=640,
            height=480,
        )
    )
    write_report(out, images=[sync_path, async_path], width=640, height=480)


if __name__ == "__main__":
    main()
