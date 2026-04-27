#!/usr/bin/env python3
"""Pormake MOF Construction

Demonstrates the Pormake backend for MOF construction via the MOFBuilder
facade.

Uses the MOFBuilder facade to:
  1. Register node (and optionally edge) building blocks
  2. Build a MOF with a specified RCSR topology
  3. Print the result summary

Building blocks can be specified as file paths to XYZ files (with ``X``
atoms marking connection points) or as names from the pormake database.

Usage:
    python build/pormake_build.py --node node.xyz --topology pcu
    python build/pormake_build.py --node N109 --edge E41 --topology pcu
    python build/pormake_build.py --node node.xyz --edge edge.xyz -t pcu --accuracy 10
    python build/pormake_build.py --list-topologies
    python build/pormake_build.py --list-bbs
    python build/pormake_build.py --describe-topology pcu
    python build/pormake_build.py --status
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
EXAMPLES_DIR = SCRIPT_DIR.parent


def _resolve_source(name: str) -> str | Path:
    """Return a Path if *name* points to an existing file, else the raw name.

    This lets users pass either a file path (``./bbs/node.xyz``) or a
    pormake database name (``N109``).
    """
    p = Path(name)
    if p.is_file():
        return p.resolve()
    return name


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the Pormake backend via mofforge",
    )
    parser.add_argument(
        "--topology", "-t",
        default="pcu",
        help="RCSR topology code, e.g. 'pcu', 'dia' (default: pcu)",
    )
    parser.add_argument(
        "--node", "-n",
        action="append",
        default=None,
        help="Node building block: path to XYZ file or pormake DB name "
             "(repeatable, at least one required for building)",
    )
    parser.add_argument(
        "--edge", "-e",
        action="append",
        default=None,
        help="Edge building block: path to XYZ file or pormake DB name "
             "(repeatable, optional)",
    )
    parser.add_argument(
        "-o", "--output-dir",
        default="./pormake_output",
        help="Directory for output CIF files (default: ./pormake_output)",
    )
    parser.add_argument(
        "--bb-dir",
        default=None,
        help="Optional directory for on-disk building-block file storage",
    )
    parser.add_argument(
        "--list-topologies",
        action="store_true",
        help="List available RCSR topologies and exit",
    )
    parser.add_argument(
        "--list-bbs",
        action="store_true",
        help="List available building blocks in the pormake database and exit",
    )
    parser.add_argument(
        "--describe-topology",
        metavar="NAME",
        default=None,
        help="Describe a topology (node/edge types, coordination) and exit",
    )
    parser.add_argument(
        "--status",
        action="store_true",
        help="Show backend status and exit",
    )
    parser.add_argument(
        "--accuracy",
        type=float,
        default=None,
        help="Pormake build accuracy parameter (passed to build_by_type)",
    )
    args = parser.parse_args()

    # ------------------------------------------------------------------ #
    # Initialize the builder
    # ------------------------------------------------------------------ #
    from mofforge.build import MOFBuilder

    kwargs: dict[str, object] = {"output_dir": args.output_dir}
    if args.bb_dir:
        kwargs["bb_dir"] = args.bb_dir

    try:
        builder = MOFBuilder(backend="pormake", **kwargs)
    except Exception as exc:
        print(f"Error initializing Pormake backend: {exc}", file=sys.stderr)
        sys.exit(1)

    # ------------------------------------------------------------------ #
    # Info-only modes
    # ------------------------------------------------------------------ #
    if args.status:
        status = builder.status()
        print(json.dumps(status, indent=2, default=str))
        return

    if args.list_topologies:
        topologies = builder.list_topologies()
        print(f"Available RCSR topologies ({len(topologies)}):")
        for t in topologies:
            print(f"  {t}")
        return

    if args.list_bbs:
        # Pormake's database doesn't strongly distinguish node vs edge;
        # role is determined by connection-point count at build time.
        backend = builder.backend
        bbs = backend.list_building_blocks("node")
        print(f"Available building blocks in pormake database ({len(bbs)}):")
        for b in bbs:
            print(f"  {b}")
        return

    if args.describe_topology is not None:
        description = builder.describe_topology(args.describe_topology)
        print(description)
        return

    # ------------------------------------------------------------------ #
    # Validate building-block arguments
    # ------------------------------------------------------------------ #
    if not args.node:
        print(
            "Error: At least one --node / -n is required for building.\n"
            "Provide a path to an XYZ file or a pormake database name.\n"
            "Use --list-bbs to see available database entries.",
            file=sys.stderr,
        )
        sys.exit(1)

    node_names = args.node
    edge_names = args.edge or []
    topology = args.topology

    print(f"Topology:  {topology}")
    print(f"Nodes:     {', '.join(node_names)}")
    if edge_names:
        print(f"Edges:     {', '.join(edge_names)}")
    else:
        print("Edges:     (none)")
    print(f"Output:    {args.output_dir}")
    if args.accuracy is not None:
        print(f"Accuracy:  {args.accuracy}")
    print()

    # ------------------------------------------------------------------ #
    # Register building blocks
    # ------------------------------------------------------------------ #
    for name in node_names:
        source = _resolve_source(name)
        bb_name = Path(name).stem if Path(name).is_file() else name
        try:
            builder.add_node(source, name=bb_name)
        except Exception as exc:
            print(f"Error registering node '{name}': {exc}", file=sys.stderr)
            sys.exit(1)
    print(f"Registered {len(node_names)} node(s)")

    for name in edge_names:
        source = _resolve_source(name)
        bb_name = Path(name).stem if Path(name).is_file() else name
        try:
            builder.add_edge(source, name=bb_name)
        except Exception as exc:
            print(f"Error registering edge '{name}': {exc}", file=sys.stderr)
            sys.exit(1)
    if edge_names:
        print(f"Registered {len(edge_names)} edge(s)")

    # ------------------------------------------------------------------ #
    # Build
    # ------------------------------------------------------------------ #
    print()
    print(f"Building MOF with topology '{topology}' ...")

    build_kwargs: dict[str, object] = {}
    if args.accuracy is not None:
        build_kwargs["accuracy"] = args.accuracy

    result = builder.build(
        topology=topology, output_dir=args.output_dir, **build_kwargs
    )

    # ------------------------------------------------------------------ #
    # Report results
    # ------------------------------------------------------------------ #
    print()
    if result.success:
        print(f"Build succeeded in {result.elapsed_seconds}s")
        print(f"  Output files: {len(result.output_paths)}")
        for p in result.output_paths:
            print(f"    {p}")
        if result.crystal:
            print(f"  Crystal loaded: {result.crystal.n_atoms} atoms")
        if result.errors:
            print(f"  Warnings ({len(result.errors)}):")
            for err in result.errors:
                print(f"    {err}")
    else:
        print("Build FAILED")
        for err in result.errors:
            print(f"  {err}")
        sys.exit(1)

    print()
    print("Metadata:")
    print(json.dumps(result.metadata, indent=2, default=str))


if __name__ == "__main__":
    main()
