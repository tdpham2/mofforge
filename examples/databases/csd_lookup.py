#!/usr/bin/env python3
"""Query a CSD-format table by name, refcode, DOI, formula, and CCDC number.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import example_parser, output_dir, stage_demo_data, write_report


def main():
    parser = example_parser("csd_lookup", __doc__)
    parser.add_argument("--data-path", type=Path, help="Your licensed ConQuest TSV export.")
    parser.add_argument("--query", default="DEMOZN")
    parser.add_argument("--field", default="auto")
    args = parser.parse_args()
    out = output_dir(args)
    import shutil
    from dataclasses import asdict

    from mofforge.csd import CSDDatabase

    if args.data_path:
        source = out / "csd_input.tab"
        shutil.copyfile(args.data_path, source)
    else:
        source, _, _ = stage_demo_data(out)
    db = CSDDatabase(source)
    try:
        result = db.search(args.query, field=args.field, limit=5)
        rows = [asdict(rec) for rec in result.records]
        # These extra queries deliberately use fictional identifiers in demo mode.
        examples = {}
        if not args.data_path:
            for query, field in [
                ("Fictional zinc", "name"),
                ("10.0000/fictional-zinc", "doi"),
                ("999991", "ccdc"),
                ("Zn", "formula"),
            ]:
                examples[field] = [r.refcode for r in db.search(query, field=field).records]
        write_report(
            out,
            synthetic_data=args.data_path is None,
            records=rows,
            detected_field=result.field,
            searches=examples,
            cached_metadata=source,
        )
    finally:
        db.close()


if __name__ == "__main__":
    main()
