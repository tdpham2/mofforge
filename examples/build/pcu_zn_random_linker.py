#!/usr/bin/env python3
"""Build a pcu MOF with a Zn paddlewheel node and a random linker.

Demonstrates how to:
  1. Query the pormake database for available building blocks
  2. Auto-detect a Zn paddlewheel node (or accept one via CLI)
  3. Pick a random edge (linker) from the database
  4. Build a MOF with the pcu topology via MOFBuilder

Usage:
    python build/pcu_zn_random_linker.py
    python build/pcu_zn_random_linker.py --seed 42
    python build/pcu_zn_random_linker.py --node N109 --seed 7
    python build/pcu_zn_random_linker.py --list-edges
    python build/pcu_zn_random_linker.py --list-nodes
    python build/pcu_zn_random_linker.py --describe-topology
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path


def _get_pormake_database():
    """Return a pormake Database instance."""
    try:
        import pormake as pm
    except ImportError:
        print(
            "Error: pormake is required.\n"
            "Install with:  pip install pormake",
            file=sys.stderr,
        )
        sys.exit(1)
    return pm.Database()


def _classify_building_blocks(db) -> tuple[list[str], list[str]]:
    """Split pormake building blocks into nodes and edges.

    Pormake distinguishes nodes from edges via the ``is_node`` attribute
    on each BuildingBlock object.  As a fast heuristic the database also
    uses the naming convention N* = node, E* = edge.
    """
    nodes: list[str] = []
    edges: list[str] = []
    for name in sorted(db.bb_list):
        try:
            bb = db.get_bb(name)
            if bb.is_node:
                nodes.append(name)
            else:
                edges.append(name)
        except Exception:
            # Some entries may fail to load -- skip them.
            continue
    return nodes, edges


def _find_zn_paddlewheel(db, node_names: list[str]) -> str | None:
    """Find a Zn-containing paddlewheel node in the database.

    Strategy:
      1. Look for nodes whose chemical formula contains Zn.
      2. Among those, prefer 4-connected nodes (paddlewheel coordination).
      3. Return the first match, or None.
    """
    zn_candidates: list[tuple[str, int]] = []
    for name in node_names:
        try:
            bb = db.get_bb(name)
            symbols = list(bb.atoms.symbols)
            if "Zn" in symbols:
                n_cp = int(bb.n_connection_points)
                zn_candidates.append((name, n_cp))
        except Exception:
            continue

    if not zn_candidates:
        return None

    # Prefer 4-connected (classic paddlewheel), then 6-connected.
    for target_cp in (4, 6):
        for name, n_cp in zn_candidates:
            if n_cp == target_cp:
                return name

    # Fall back to the first Zn node found.
    return zn_candidates[0][0]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a pcu MOF with a Zn paddlewheel node and a random linker",
    )
    parser.add_argument(
        "--topology", "-t",
        default="pcu",
        help="RCSR topology code (default: pcu)",
    )
    parser.add_argument(
        "--node", "-n",
        default=None,
        help="Node building block: pormake DB name (e.g. N109). "
             "If omitted, auto-detects a Zn paddlewheel from the database.",
    )
    parser.add_argument(
        "--edge", "-e",
        default=None,
        help="Edge building block: pormake DB name. "
             "If omitted, one is chosen at random.",
    )
    parser.add_argument(
        "-o", "--output-dir",
        default="./pormake_output",
        help="Directory for output CIF files (default: ./pormake_output)",
    )
    parser.add_argument(
        "--seed", "-s",
        type=int,
        default=None,
        help="Random seed for reproducible linker selection",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Number of random edges to try if the build fails (default: 3)",
    )
    parser.add_argument(
        "--list-nodes",
        action="store_true",
        help="List available node building blocks and exit",
    )
    parser.add_argument(
        "--list-edges",
        action="store_true",
        help="List available edge building blocks and exit",
    )
    parser.add_argument(
        "--describe-topology",
        action="store_true",
        help="Describe the selected topology and exit",
    )
    parser.add_argument(
        "--accuracy",
        type=float,
        default=None,
        help="Pormake build accuracy parameter",
    )
    args = parser.parse_args()

    # ------------------------------------------------------------------ #
    # Load the pormake database and classify building blocks
    # ------------------------------------------------------------------ #
    db = _get_pormake_database()
    print("Loading pormake database...")
    node_names, edge_names = _classify_building_blocks(db)
    print(f"  Found {len(node_names)} nodes and {len(edge_names)} edges")

    # ------------------------------------------------------------------ #
    # Info-only modes
    # ------------------------------------------------------------------ #
    if args.list_nodes:
        print(f"\nAvailable node building blocks ({len(node_names)}):")
        for name in node_names:
            try:
                bb = db.get_bb(name)
                formula = bb.atoms.get_chemical_formula()
                n_cp = bb.n_connection_points
                print(f"  {name:<10s}  {formula:<20s}  {n_cp} connection points")
            except Exception:
                print(f"  {name:<10s}  (failed to load)")
        return

    if args.list_edges:
        print(f"\nAvailable edge building blocks ({len(edge_names)}):")
        for name in edge_names:
            try:
                bb = db.get_bb(name)
                formula = bb.atoms.get_chemical_formula()
                n_cp = bb.n_connection_points
                print(f"  {name:<10s}  {formula:<20s}  {n_cp} connection points")
            except Exception:
                print(f"  {name:<10s}  (failed to load)")
        return

    if args.describe_topology:
        from mofforge.build import MOFBuilder

        builder = MOFBuilder(backend="pormake", output_dir=args.output_dir)
        desc = builder.describe_topology(args.topology)
        print(f"\nTopology: {args.topology}")
        print(desc)
        return

    # ------------------------------------------------------------------ #
    # Select the node building block
    # ------------------------------------------------------------------ #
    if args.node is not None:
        node_name = args.node
        print(f"\nUsing user-specified node: {node_name}")
    else:
        print("\nSearching for a Zn paddlewheel node in the database...")
        node_name = _find_zn_paddlewheel(db, node_names)
        if node_name is None:
            print(
                "Error: Could not find a Zn-containing node in the pormake database.\n"
                "Specify one manually with --node <NAME>.\n"
                "Use --list-nodes to see available options.",
                file=sys.stderr,
            )
            sys.exit(1)
        try:
            bb = db.get_bb(node_name)
            formula = bb.atoms.get_chemical_formula()
            n_cp = bb.n_connection_points
            print(f"  Found: {node_name}  ({formula}, {n_cp} connection points)")
        except Exception:
            print(f"  Found: {node_name}")

    # ------------------------------------------------------------------ #
    # Select the edge building block (random or user-specified)
    # ------------------------------------------------------------------ #
    if args.seed is not None:
        random.seed(args.seed)
        print(f"\nRandom seed: {args.seed}")

    if args.edge is not None:
        selected_edges = [args.edge]
        print(f"\nUsing user-specified edge: {args.edge}")
    else:
        if not edge_names:
            print("Error: No edge building blocks found in the database.", file=sys.stderr)
            sys.exit(1)

        # Shuffle and pick up to `retries` edges to try.
        candidates = list(edge_names)
        random.shuffle(candidates)
        selected_edges = candidates[: max(1, args.retries)]
        print(f"\nRandomly selected edge candidate(s): {', '.join(selected_edges)}")

    # ------------------------------------------------------------------ #
    # Build with MOFBuilder, retrying with different edges on failure
    # ------------------------------------------------------------------ #
    from mofforge.build import MOFBuilder

    for attempt, edge_name in enumerate(selected_edges, start=1):
        print(f"\n{'=' * 60}")
        print(f"  Attempt {attempt}/{len(selected_edges)}")
        print(f"  Topology: {args.topology}")
        print(f"  Node:     {node_name}")
        print(f"  Edge:     {edge_name}")
        print(f"  Output:   {args.output_dir}")
        print(f"{'=' * 60}")

        try:
            builder = MOFBuilder(backend="pormake", output_dir=args.output_dir)

            # Register the node (using database name as source).
            builder.add_node(node_name, name=node_name)

            # Register the edge.
            builder.add_edge(edge_name, name=edge_name)

            # Build.
            build_kwargs: dict[str, object] = {}
            if args.accuracy is not None:
                build_kwargs["accuracy"] = args.accuracy

            print(f"\nBuilding MOF with topology '{args.topology}' ...")
            result = builder.build(
                topology=args.topology,
                output_dir=args.output_dir,
                **build_kwargs,
            )

            if result.success:
                print(f"\nBuild SUCCEEDED in {result.elapsed_seconds}s")
                print(f"  Output files: {len(result.output_paths)}")
                for p in result.output_paths:
                    print(f"    {p}")
                if result.crystal:
                    print(f"  Crystal loaded: {result.crystal.n_atoms} atoms")
                if result.errors:
                    print(f"  Warnings ({len(result.errors)}):")
                    for err in result.errors:
                        print(f"    {err}")

                print("\nMetadata:")
                print(json.dumps(result.metadata, indent=2, default=str))
                return  # success -- done

            else:
                print(f"\nBuild FAILED")
                for err in result.errors:
                    print(f"  {err}")
                if attempt < len(selected_edges):
                    print("  Retrying with a different edge...")

        except Exception as exc:
            print(f"\nError during build: {exc}", file=sys.stderr)
            if attempt < len(selected_edges):
                print("  Retrying with a different edge...")

    # All attempts exhausted.
    print("\nAll attempts failed. Try a different node or topology.", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
