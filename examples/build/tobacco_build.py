#!/usr/bin/env python3
"""TOBACCO MOF Construction

Demonstrates the TOBACCO backend with building-block staging via the
MOFBuilder facade.

Uses the MOFBuilder facade to:
  1. Copy a node and edge from the TOBACCO database into the active dirs
  2. Copy a template from the template database
  3. Build a MOF with the pcu topology
  4. Print the result summary

Usage:
    python build/tobacco_build.py
    python build/tobacco_build.py --topology pcu --node 6c_Zn_1_Ch.cif --edge 1B_4H_Ch.cif
    python build/tobacco_build.py --list-nodes
    python build/tobacco_build.py --list-edges
    python build/tobacco_build.py --list-topologies
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
EXAMPLES_DIR = SCRIPT_DIR.parent


def _find_in_database(db_dir: Path, name: str) -> Path | None:
    """Find a CIF file by name in a database directory (recursive)."""
    matches = list(db_dir.rglob(name))
    return matches[0] if matches else None


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the TOBACCO backend via mofforge")
    parser.add_argument(
        "--tobacco-path",
        default=None,
        help="Path to the TOBACCO 3.0 directory (reads mofforge.toml if omitted)",
    )
    parser.add_argument(
        "--topology", "-t",
        default="pcu",
        help="Topology name, e.g. 'pcu' (default: pcu)",
    )
    parser.add_argument(
        "--node", "-n",
        action="append",
        default=None,
        help="Node CIF filename from the database (repeatable). "
             "Default: 6c_Zn_1_Ch.cif",
    )
    parser.add_argument(
        "--edge", "-e",
        action="append",
        default=None,
        help="Edge CIF filename from the database (repeatable). "
             "Default: 1B_4H_Ch.cif",
    )
    parser.add_argument(
        "-o", "--output-dir",
        default="./tobacco_output",
        help="Directory for output CIF files (default: ./tobacco_output)",
    )
    parser.add_argument(
        "--list-nodes",
        action="store_true",
        help="List available nodes in the database and exit",
    )
    parser.add_argument(
        "--list-edges",
        action="store_true",
        help="List available edges in the database and exit",
    )
    parser.add_argument(
        "--list-topologies",
        action="store_true",
        help="List available templates in the database and exit",
    )
    parser.add_argument(
        "--status",
        action="store_true",
        help="Show backend status and exit",
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Run TOBACCO in parallel mode",
    )
    args = parser.parse_args()

    # ------------------------------------------------------------------ #
    # Initialize the builder
    # ------------------------------------------------------------------ #
    from mofforge.build import MOFBuilder

    kwargs = {}
    if args.tobacco_path:
        kwargs["tobacco_path"] = args.tobacco_path

    try:
        builder = MOFBuilder(backend="tobacco", **kwargs)
    except Exception as exc:
        print(f"Error initializing TOBACCO backend: {exc}", file=sys.stderr)
        sys.exit(1)

    # ------------------------------------------------------------------ #
    # Info-only modes
    # ------------------------------------------------------------------ #
    if args.status:
        status = builder.status()
        print(json.dumps(status, indent=2, default=str))
        return

    if args.list_nodes:
        result = builder.copy_from_database("node", names=None)
        nodes = result.get("available_in_database", [])
        print(f"Available nodes in database ({len(nodes)}):")
        for n in nodes:
            print(f"  {n}")
        return

    if args.list_edges:
        result = builder.copy_from_database("edge", names=None)
        edges = result.get("available_in_database", [])
        print(f"Available edges in database ({len(edges)}):")
        for e in edges:
            print(f"  {e}")
        return

    if args.list_topologies:
        # List from the template_database (not the active templates/ dir)
        backend = builder.backend
        db_dir = backend._db_dir_for_role("template")
        if db_dir.is_dir():
            templates = sorted(f.name for f in db_dir.iterdir() if f.suffix == ".cif")
        else:
            templates = builder.list_topologies()
        print(f"Available topologies ({len(templates)}):")
        for t in templates:
            print(f"  {t}")
        return

    # ------------------------------------------------------------------ #
    # Set up building blocks
    # ------------------------------------------------------------------ #
    node_names = args.node or ["6c_Zn_1_Ch.cif"]
    edge_names = args.edge or ["1B_4H_Ch.cif"]
    topology = args.topology

    print(f"Topology:  {topology}")
    print(f"Nodes:     {', '.join(node_names)}")
    print(f"Edges:     {', '.join(edge_names)}")
    print(f"Output:    {args.output_dir}")
    print()

    # Copy template from database into templates/
    backend = builder.backend
    template_name = topology if topology.endswith(".cif") else topology + ".cif"
    template_result = backend.copy_from_database(
        "template", names=[template_name], dry_run=False  # type: ignore[arg-type]
    )
    if not template_result.get("success") and not template_result.get("copied"):
        # Template might already be in templates/ -- check
        if template_name not in builder.list_topologies():
            print(f"Error: Could not find template '{template_name}' in database or templates/")
            print(f"  Result: {json.dumps(template_result, indent=2, default=str)}")
            sys.exit(1)
    print(f"Template '{template_name}' ready")

    # Resolve database source paths for each building block.
    # We point add_node/add_edge at the *database* copies so that
    # _stage_building_blocks() can clear nodes/ and edges/ and then
    # copy from these stable source paths.
    node_db = backend._db_dir_for_role("node")
    edge_db = backend._db_dir_for_role("edge")

    for name in node_names:
        src = _find_in_database(node_db, name)
        if src is None:
            print(f"Error: Node '{name}' not found in database ({node_db})")
            sys.exit(1)
        builder.add_node(src, name=Path(name).stem)
    print(f"Registered {len(node_names)} node(s)")

    for name in edge_names:
        src = _find_in_database(edge_db, name)
        if src is None:
            print(f"Error: Edge '{name}' not found in database ({edge_db})")
            sys.exit(1)
        builder.add_edge(src, name=Path(name).stem)
    print(f"Registered {len(edge_names)} edge(s)")

    # ------------------------------------------------------------------ #
    # Build
    # ------------------------------------------------------------------ #
    print()
    print(f"Building MOF with topology '{topology}' ...")

    build_kwargs = {}
    if args.parallel:
        build_kwargs["parallel"] = True

    result = builder.build(topology=topology, output_dir=args.output_dir, **build_kwargs)

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
