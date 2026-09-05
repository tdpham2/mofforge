#!/usr/bin/env python3
"""Build a topology with Pormake and inspect execution and validation separately.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import example_parser, output_dir, require, write_report


def main():
    parser = example_parser("pormake_build", __doc__)
    parser.add_argument("--topology", "-t", default="pcu")
    parser.add_argument("--node", "-n", action="append")
    parser.add_argument("--edge", "-e", action="append")
    parser.add_argument("--bb-dir", type=Path)
    parser.add_argument("--list-topologies", action="store_true")
    parser.add_argument("--list-bbs", action="store_true")
    parser.add_argument("--describe-topology")
    parser.add_argument("--status", action="store_true")
    parser.add_argument("--accuracy", type=float, default=3)
    args = parser.parse_args()
    require("pormake", "build")
    out = output_dir(args)
    from mofforge import MOFBuilder

    builder = MOFBuilder("pormake", output_dir=out, bb_dir=args.bb_dir)
    if args.status or args.list_topologies or args.list_bbs or args.describe_topology:
        write_report(
            out,
            status=builder.status(),
            topologies=builder.list_topologies() if args.list_topologies else [],
            nodes=builder.copy_from_database("node") if args.list_bbs else {},
            edges=builder.copy_from_database("edge") if args.list_bbs else {},
            topology=builder.describe_topology(args.describe_topology)
            if args.describe_topology
            else None,
        )
        return
    # The default N108/E1/pcu combination is also exercised by backend integration tests.
    for node in args.node or ["N108"]:
        builder.add_node(node)
    for edge in args.edge if args.edge is not None else ([] if args.node else ["E1"]):
        builder.add_edge(edge)
    result = builder.build(args.topology, output_dir=out, accuracy=args.accuracy)
    write_report(
        out,
        success=result.success,
        errors=result.errors,
        outputs=result.output_paths,
        atoms=result.crystal.n_atoms if result.crystal else None,
        validation=result.validation.to_dict() if result.validation else None,
    )
    if not result.success:
        raise SystemExit("Pormake did not produce usable outputs; see summary.json.")


if __name__ == "__main__":
    main()
