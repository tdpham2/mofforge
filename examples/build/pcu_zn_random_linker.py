#!/usr/bin/env python3
"""Select a topology-compatible Zn node and a reproducible random linker.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import example_parser, output_dir, require, write_report


def main():
    parser = example_parser("pcu_zn_random_linker", __doc__, seed=True)
    parser.add_argument("--topology", "-t", default="pcu")
    parser.add_argument("--node", "-n")
    parser.add_argument("--edge", "-e")
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--accuracy", type=float, default=3)
    parser.add_argument("--list-nodes", action="store_true")
    parser.add_argument("--list-edges", action="store_true")
    parser.add_argument("--describe-topology", action="store_true")
    args = parser.parse_args()
    require("pormake", "build")
    if args.retries < 1:
        parser.error("--retries must be positive.")
    out = output_dir(args)
    import random

    import pormake as pm

    from mofforge import MOFBuilder

    db = pm.Database()
    topology = db.get_topo(args.topology)
    required = {int(cn) for cn in topology.unique_cn}
    if len(required) != 1:
        parser.error("This single-node lesson requires a topology with one coordination number.")
    connectivity = next(iter(required))
    nodes, edges, unreadable = [], [], []
    for name in sorted(db.bb_list):
        try:
            block = db.get_bb(name)
        except (OSError, ValueError, KeyError) as exc:
            unreadable.append({"name": name, "error": str(exc)})
            continue
        if block.is_node:
            if (
                "Zn" in block.atoms.get_chemical_symbols()
                and block.n_connection_points == connectivity
            ):
                nodes.append(name)
        elif block.n_connection_points == 2:
            edges.append(name)
    if args.list_nodes or args.list_edges or args.describe_topology:
        write_report(
            out,
            required_connectivity=connectivity,
            compatible_zinc_nodes=nodes,
            two_connected_edges=edges,
            unreadable=unreadable,
        )
        return
    if not nodes:
        raise SystemExit("No compatible Zn node found for this topology.")
    node = args.node or ("N108" if "N108" in nodes else nodes[0])
    if node not in nodes:
        parser.error(f"{node} is not a Zn node with {connectivity} connection points.")
    if args.edge and args.edge not in edges:
        parser.error("The selected edge is not a two-connected catalog linker.")
    rng = random.Random(args.seed)
    candidates = [args.edge] if args.edge else rng.sample(edges, min(args.retries, len(edges)))
    attempts = []
    for index, edge in enumerate(candidates):
        builder = MOFBuilder("pormake", output_dir=out / f"attempt_{index + 1}")
        builder.add_node(node)
        builder.add_edge(edge)
        result = builder.build(
            args.topology, output_dir=out / f"attempt_{index + 1}", accuracy=args.accuracy
        )
        attempts.append(
            {
                "edge": edge,
                "success": result.success,
                "errors": result.errors,
                "outputs": result.output_paths,
                "validation": result.validation.to_dict() if result.validation else None,
            }
        )
        if result.success:
            break
    write_report(
        out,
        seed=args.seed,
        node=node,
        required_connectivity=connectivity,
        attempts=attempts,
        unreadable=unreadable,
    )
    if not attempts or not attempts[-1]["success"]:
        raise SystemExit("Every selected linker failed to build; see recorded attempts.")


if __name__ == "__main__":
    main()
