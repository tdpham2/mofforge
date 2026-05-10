#!/usr/bin/env python3
"""Extract descriptive metadata for all topologies in the pormake database.

Produces a structured output file (JSONL, JSON, or CSV) with enriched
information for each RCSR topology: spacegroup, node/edge counts,
coordination numbers, node types, and edge types.

Usage:
    python scripts/describe_topologies.py
    python scripts/describe_topologies.py --format csv -o topologies.csv
    python scripts/describe_topologies.py --limit 50        # first 50 only
    python scripts/describe_topologies.py --names pcu dia tbo  # specific ones

Requirements:
    This script requires a separate Python environment with compatible
    package versions.  pormake 0.2.x depends on pymatgen<2024 which
    requires numpy<2.  See scripts/requirements.txt for details.

    Install into a dedicated environment:
        pip install -r scripts/requirements.txt

    Or manually:
        pip install pormake 'pymatgen>=2023.8.10,<2024.0.0' 'numpy<2'
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import time
from pathlib import Path

_COMPAT_MSG = (
    "Error: Failed to import pormake.\n\n"
    "These scripts require a Python environment with compatible packages.\n"
    "Install with:\n"
    "  pip install -r scripts/requirements.txt\n\n"
    "The key constraint is numpy<2 (pormake's pymatgen<2024 dependency\n"
    "uses APIs removed in numpy 2.0).  The mofforge .venv may have\n"
    "incompatible versions -- use a separate environment for these scripts."
)


def get_database():
    """Return a pormake Database instance.

    Catches import failures caused by version conflicts (e.g. numpy>=2
    with pymatgen<2024) and prints an actionable error message.
    """
    try:
        import pormake as pm
    except (ImportError, AttributeError) as exc:
        print(f"{_COMPAT_MSG}\n\nUnderlying error: {exc}", file=sys.stderr)
        sys.exit(1)

    return pm.Database()


def describe_one_topology(db, name: str) -> dict | None:
    """Extract metadata for a single topology.

    Returns a dict on success or None on failure (with a message
    printed to stderr).
    """
    try:
        topo = db.get_topo(name)
    except Exception as exc:
        print(f"  SKIP {name}: {exc}", file=sys.stderr)
        return None

    # Capture the full describe() text for reference.
    buf = io.StringIO()
    topo.describe(file=buf)
    describe_text = buf.getvalue()

    # Build per-node-type detail list.
    import numpy as np

    node_type_details = []
    for t, cn in zip(topo.unique_node_types, topo.unique_cn):
        indices = np.argwhere(topo.node_types == t).reshape(-1).tolist()
        node_type_details.append({
            "type": int(t),
            "coordination_number": int(cn),
            "count": len(indices),
        })

    # Build per-edge-type detail list.
    edge_type_details = []
    for t in topo.unique_edge_types:
        condition = np.all(topo.edge_types == t, axis=1)
        indices = np.argwhere(condition).reshape(-1).tolist()
        edge_type_details.append({
            "type": (int(t[0]), int(t[1])),
            "count": len(indices),
        })

    return {
        "name": topo.name,
        "spacegroup": topo.spacegroup,
        "n_nodes": int(topo.n_nodes),
        "n_edges": int(topo.n_edges),
        "n_slots": int(topo.n_slots),
        "n_node_types": int(topo.n_node_types),
        "n_edge_types": int(topo.n_edge_types),
        "coordination_numbers": [int(c) for c in topo.unique_cn],
        "node_type_details": node_type_details,
        "edge_type_details": edge_type_details,
        "describe_text": describe_text,
    }


# ------------------------------------------------------------------ #
# Writers
# ------------------------------------------------------------------ #


def write_jsonl(records: list[dict], path: Path) -> None:
    with path.open("w") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")


def write_json(records: list[dict], path: Path) -> None:
    with path.open("w") as f:
        json.dump(records, f, indent=2, ensure_ascii=False)


def write_csv(records: list[dict], path: Path) -> None:
    """Write a flat CSV (nested fields are JSON-encoded)."""
    if not records:
        return

    flat_fields = [
        "name",
        "spacegroup",
        "n_nodes",
        "n_edges",
        "n_slots",
        "n_node_types",
        "n_edge_types",
        "coordination_numbers",
        "node_type_details",
        "edge_type_details",
    ]

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=flat_fields)
        writer.writeheader()
        for rec in records:
            row = {}
            for k in flat_fields:
                v = rec[k]
                if isinstance(v, (list, dict)):
                    row[k] = json.dumps(v)
                else:
                    row[k] = v
            writer.writerow(row)


# ------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------ #


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract descriptive metadata for pormake topologies.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output file path. Default: scripts/output/topologies.<format>",
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=["jsonl", "json", "csv"],
        default="jsonl",
        help="Output format (default: jsonl).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N topologies (for testing).",
    )
    parser.add_argument(
        "--names",
        nargs="+",
        default=None,
        help="Process only these specific topology names.",
    )
    args = parser.parse_args()

    # Resolve output path.
    if args.output is None:
        out_dir = Path(__file__).resolve().parent / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        args.output = out_dir / f"topologies.{args.format}"

    print(f"Loading pormake database...", file=sys.stderr)
    db = get_database()

    # Determine which topologies to process.
    if args.names:
        names = args.names
    else:
        names = sorted(db.topology_list)
        if args.limit:
            names = names[: args.limit]

    total = len(names)
    print(f"Processing {total} topologies...", file=sys.stderr)

    records = []
    skipped = 0
    t_start = time.time()

    for i, name in enumerate(names, start=1):
        rec = describe_one_topology(db, name)
        if rec is not None:
            records.append(rec)
        else:
            skipped += 1

        # Progress every 50 items or at the end.
        if i % 50 == 0 or i == total:
            elapsed = time.time() - t_start
            rate = i / elapsed if elapsed > 0 else 0
            eta = (total - i) / rate if rate > 0 else 0
            print(
                f"\r  [{i}/{total}] {len(records)} ok, {skipped} skipped "
                f"({elapsed:.1f}s elapsed, ~{eta:.0f}s remaining)",
                end="",
                file=sys.stderr,
            )

    print(file=sys.stderr)  # newline after progress

    # Write output.
    writers = {"jsonl": write_jsonl, "json": write_json, "csv": write_csv}
    writers[args.format](records, args.output)

    elapsed_total = time.time() - t_start
    print(
        f"Done. Wrote {len(records)} topologies to {args.output} "
        f"({skipped} skipped, {elapsed_total:.1f}s total).",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
