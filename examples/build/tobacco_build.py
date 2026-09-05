#!/usr/bin/env python3
"""Construct with TOBACCO using a tested dmc fixture or explicit catalog blocks.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import example_parser, output_dir, require, write_report


def main():
    parser = example_parser("tobacco_build", __doc__)
    parser.add_argument("--tobacco-path", "--tobacco-data", type=Path)
    parser.add_argument("--topology", "-t", default="dmc")
    parser.add_argument("--node", "-n", action="append")
    parser.add_argument("--edge", "-e", action="append")
    parser.add_argument("--parallel", action="store_true")
    for flag in ["list-nodes", "list-edges", "list-topologies", "status"]:
        parser.add_argument("--" + flag, action="store_true")
    args = parser.parse_args()
    require("tobacco3", "build")
    out = output_dir(args)
    from mofforge import MOFBuilder

    options = {"tobacco_data_dir": args.tobacco_path} if args.tobacco_path else {}
    builder = MOFBuilder("tobacco", **options)
    if args.status or args.list_nodes or args.list_edges or args.list_topologies:
        write_report(
            out,
            status=builder.status(),
            nodes=builder.copy_from_database("node") if args.list_nodes else {},
            edges=builder.copy_from_database("edge") if args.list_edges else {},
            topologies=builder.list_topologies() if args.list_topologies else [],
        )
        return
    # Keep runtime data external. This fixture matches test_build_integration.py.
    fixture = Path(builder.status()["data_dir"]) / "tests" / "fixtures"
    nodes = args.node
    edges = args.edge
    if args.topology == "dmc" and nodes is None and edges is None:
        nodes = [fixture / "nodes" / f"{name}.cif" for name in ["4c_1Zn_Ch", "triazole"]]
        edges = [
            fixture / "edges" / f"{name}.cif" for name in ["ntn_edge", "oxalic_edge", "squOxa_ch"]
        ]
    else:
        nodes = nodes or ["6c_Zn_1_Ch.cif"]
        edges = edges or ["1B_4H_Ch.cif"]
    for node in nodes:
        builder.add_node(node)
    for edge in edges:
        builder.add_edge(edge)
    result = builder.build(args.topology, output_dir=out, parallel=args.parallel)
    write_report(
        out,
        success=result.success,
        errors=result.errors,
        outputs=result.output_paths,
        atoms=result.crystal.n_atoms if result.crystal else None,
        validation=result.validation.to_dict() if result.validation else None,
    )
    if not result.success:
        raise SystemExit("TOBACCO did not produce usable outputs; see summary.json.")


if __name__ == "__main__":
    main()
