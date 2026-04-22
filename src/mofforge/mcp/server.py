"""mofforge MCP Server.

Exposes mofforge's crystal structure manipulation capabilities as
Model Context Protocol (MCP) tools for AI agents.

Available tools
---------------
- ``mofforge_search`` -- Find substructure patterns via VF2 isomorphism.
- ``mofforge_replace`` -- Replace molecular fragments in a crystal.
- ``mofforge_remove`` -- Remove guest molecules from a crystal.
- ``mofforge_validate`` -- Validate a crystal structure.
- ``mofforge_render`` -- Render a structure file to a labelled PNG image.
- ``mofforge_smarts_search`` -- Search using SMARTS-like pattern strings.
- ``mofforge_build`` -- Build a MOF from topology + building blocks.
- ``mofforge_list_topologies`` -- List available topologies.
- ``mofforge_list_building_blocks`` -- List available node/edge building blocks.

Usage
-----
::

    # stdio mode (for LangGraph / OpenCode / Claude Desktop)
    mofforge-mcp

    # Or run directly:
    python -m mofforge.mcp.server
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path

from mcp.server.fastmcp import FastMCP

logger = logging.getLogger("mofforge.mcp")

# ---------------------------------------------------------------------------
# FastMCP server instance
# ---------------------------------------------------------------------------

mcp = FastMCP(
    name="mofforge",
    instructions=(
        "mofforge tools for building and modifying atomistic crystal structure "
        "models, especially Metal-Organic Frameworks (MOFs).  Capabilities "
        "include substructure search (VF2 graph isomorphism), find-and-replace "
        "of molecular fragments, MOF construction from topology + building "
        "blocks, structure validation, SMARTS-like pattern matching, and "
        "structure rendering to PNG images for visual inspection.\n\n"
        "General guidance:\n"
        "- File paths should be absolute.\n"
        "- Crystal structures are read from CIF files.\n"
        "- Fragments are read from XYZ files with optional '!' anchor tags.\n"
        "- Distances are in Angstroms.\n"
        "- Call mofforge_validate after modifications to check for issues."
    ),
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _resolve_output(path: str) -> str:
    """Ensure output path is absolute and parent directory exists."""
    p = Path(path)
    if not p.is_absolute():
        log_dir = os.environ.get("MOFFORGE_LOG_DIR", ".")
        p = Path(log_dir) / p
    p.parent.mkdir(parents=True, exist_ok=True)
    return str(p.resolve())


def _load_crystal(cif_path: str, with_bonds: bool = True):
    """Load a Crystal from CIF, optionally inferring bonds."""
    from mofforge.core.bonding import infer_bonds
    from mofforge.core.crystal import Crystal

    xtal = Crystal.from_cif(cif_path)
    if with_bonds:
        xtal = infer_bonds(xtal, periodic=True)
    return xtal


def _load_fragment(xyz_path: str):
    """Load a fragment from an XYZ file."""
    from mofforge.core.moiety import fragment as load_fragment

    p = Path(xyz_path)
    return load_fragment(p.name, fragment_path=str(p.parent))


# ---------------------------------------------------------------------------
# MCP tool definitions
# ---------------------------------------------------------------------------


@mcp.tool(
    name="mofforge_search",
    description=(
        "Search for a molecular substructure pattern within a crystal "
        "structure using VF2 graph isomorphism.  Returns the number of "
        "matches and their locations."
    ),
)
def mofforge_search(
    parent_cif: str,
    query_xyz: str,
    disconnected: bool = False,
) -> str:
    """Find pattern matches in a crystal.

    Parameters
    ----------
    parent_cif : str
        Absolute path to the parent crystal CIF file.
    query_xyz : str
        Absolute path to the query fragment XYZ file.
    disconnected : bool
        If True, search for isolated (guest) molecules only.

    Returns
    -------
    str
        JSON summary of search results.
    """
    from mofforge.search.search import find_pattern

    xtal = _load_crystal(parent_cif)
    query = _load_fragment(query_xyz)

    result = find_pattern(query, xtal, disconnected_component=disconnected)

    return json.dumps(
        {
            "success": True,
            "parent": parent_cif,
            "query": query_xyz,
            "total_isomorphisms": result.nb_isomorphisms(),
            "locations": result.nb_locations(),
            "orientations_per_location": result.nb_ori_at_loc(),
        },
        indent=2,
    )


@mcp.tool(
    name="mofforge_replace",
    description=(
        "Find a molecular pattern in a crystal and replace it with a "
        "different fragment.  Writes the modified structure to a CIF file."
    ),
)
def mofforge_replace(
    parent_cif: str,
    query_xyz: str,
    replacement_xyz: str,
    output_cif: str = "modified.cif",
    nb_loc: int = 0,
    use_random: bool = False,
    validate: bool = False,
) -> str:
    """Find and replace a fragment in a crystal.

    Parameters
    ----------
    parent_cif : str
        Absolute path to the parent crystal CIF file.
    query_xyz : str
        Absolute path to the query fragment XYZ file.
    replacement_xyz : str
        Absolute path to the replacement fragment XYZ file.
    output_cif : str
        Output CIF file path.
    nb_loc : int
        Number of random locations to replace (0 = all).
    use_random : bool
        If True, use random orientations at each location.
    validate : bool
        If True, validate the output structure.

    Returns
    -------
    str
        JSON summary of the replacement result.
    """
    from mofforge.replace.replace import replace_pattern
    from mofforge.search.search import find_pattern

    xtal = _load_crystal(parent_cif)
    query = _load_fragment(query_xyz)
    replacement = _load_fragment(replacement_xyz)

    result = find_pattern(query, xtal)

    kwargs = {"verbose": False}
    if nb_loc > 0:
        kwargs["nb_loc"] = nb_loc
    if use_random:
        kwargs["random"] = True

    child = replace_pattern(result, replacement, **kwargs)

    output_path = _resolve_output(output_cif)
    child.write_cif(output_path)

    response = {
        "success": True,
        "output_cif": output_path,
        "parent_atoms": xtal.n_atoms,
        "child_atoms": child.n_atoms,
        "matches_found": result.nb_isomorphisms(),
        "locations_found": result.nb_locations(),
    }

    if validate:
        from mofforge.core.bonding import infer_bonds
        from mofforge.validation import validate_structure

        child = infer_bonds(child, periodic=True)
        report = validate_structure(child)
        response["validation"] = {
            "is_valid": report.is_valid,
            "steric_clashes": len(report.steric_clashes),
            "unusual_bonds": len(report.unusual_bonds),
            "coordination_issues": len(report.coordination_issues),
            "summary": report.summary(),
        }

    return json.dumps(response, indent=2)


@mcp.tool(
    name="mofforge_remove",
    description=(
        "Remove guest molecules from a crystal structure.  Uses "
        "disconnected-component search to find isolated molecules."
    ),
)
def mofforge_remove(
    parent_cif: str,
    guest_xyz: str,
    output_cif: str = "clean.cif",
) -> str:
    """Remove guest molecules from a crystal.

    Parameters
    ----------
    parent_cif : str
        Absolute path to the parent crystal CIF file.
    guest_xyz : str
        Absolute path to the guest molecule XYZ file.
    output_cif : str
        Output CIF file path.

    Returns
    -------
    str
        JSON summary of the removal result.
    """
    from mofforge.replace.replace import replace_pattern
    from mofforge.search.search import find_pattern

    xtal = _load_crystal(parent_cif)
    guest = _load_fragment(guest_xyz)

    result = find_pattern(guest, xtal, disconnected_component=True)
    child = replace_pattern(result, None)

    output_path = _resolve_output(output_cif)
    child.write_cif(output_path)

    return json.dumps(
        {
            "success": True,
            "output_cif": output_path,
            "guests_removed": result.nb_locations(),
            "atoms_before": xtal.n_atoms,
            "atoms_after": child.n_atoms,
            "atoms_removed": xtal.n_atoms - child.n_atoms,
        },
        indent=2,
    )


@mcp.tool(
    name="mofforge_validate",
    description=(
        "Validate a crystal structure for steric clashes, unusual bond "
        "lengths, and metal coordination number issues."
    ),
)
def mofforge_validate(
    cif_path: str,
    check_clashes: bool = True,
    check_bonds: bool = True,
    check_coordination: bool = True,
) -> str:
    """Validate a crystal structure.

    Parameters
    ----------
    cif_path : str
        Absolute path to the CIF file.
    check_clashes : bool
        Check for steric clashes.
    check_bonds : bool
        Check for unusual bond lengths.
    check_coordination : bool
        Check metal coordination numbers.

    Returns
    -------
    str
        JSON validation report.
    """
    from mofforge.validation import validate_structure

    xtal = _load_crystal(cif_path, with_bonds=True)
    report = validate_structure(
        xtal,
        check_clashes=check_clashes,
        check_bonds=check_bonds,
        check_coordination=check_coordination,
    )

    return json.dumps(
        {
            "success": True,
            "is_valid": report.is_valid,
            "steric_clashes": len(report.steric_clashes),
            "unusual_bonds": len(report.unusual_bonds),
            "coordination_issues": len(report.coordination_issues),
            "warnings": report.warnings,
            "summary": report.summary(),
        },
        indent=2,
    )


@mcp.tool(
    name="mofforge_render",
    description=(
        "Render a crystal structure file (CIF or XYZ) to a labelled PNG "
        "image.  Useful for visual inspection and vision-LLM analysis.  "
        "Supports unit cell rendering for periodic structures."
    ),
)
def mofforge_render(
    input_file: str,
    output_file: str = "structure.png",
    label_mode: str = "sequential",
    representation: str = "ball_stick",
    show_unit_cell: bool = False,
    width: int = 800,
    height: int = 600,
) -> str:
    """Render a structure to PNG.

    Parameters
    ----------
    input_file : str
        Path to the input CIF or XYZ file.
    output_file : str
        Path for the output PNG file.
    label_mode : str
        Atom label mode: "sequential", "per_element", or "none".
    representation : str
        Visual style: "ball_stick", "stick", or "sphere".
    show_unit_cell : bool
        Draw unit cell edges (auto-enabled for CIF files).
    width, height : int
        Image dimensions in pixels.

    Returns
    -------
    str
        JSON with the path to the rendered PNG.
    """
    from mofforge.vis.render import render_file_to_png

    output_path = _resolve_output(output_file)

    rendered = render_file_to_png(
        input_file=input_file,
        output_file=output_path,
        label_mode=label_mode,
        representation=representation,
        show_unit_cell=show_unit_cell,
        width=width,
        height=height,
    )

    return json.dumps(
        {
            "success": True,
            "input_file": input_file,
            "output_png": rendered,
        },
        indent=2,
    )


@mcp.tool(
    name="mofforge_smarts_search",
    description=(
        "Search for a SMARTS-like pattern string in a crystal structure.  "
        "Supports element symbols, wildcards [*], bonds (-), and ring "
        "closures (e.g. 'C1-C-C-C-C-C-1' for a 6-ring)."
    ),
)
def mofforge_smarts_search(
    cif_path: str,
    pattern: str,
) -> str:
    """Search with a SMARTS-like pattern.

    Parameters
    ----------
    cif_path : str
        Absolute path to the CIF file.
    pattern : str
        SMARTS-like pattern (e.g. "Zn-O-C", "[*]-N-[*]").

    Returns
    -------
    str
        JSON summary of matches found.
    """
    from mofforge.smarts import smarts_search

    xtal = _load_crystal(cif_path, with_bonds=True)
    result = smarts_search(pattern, xtal)

    return json.dumps(
        {
            "success": True,
            "pattern": pattern,
            "cif_path": cif_path,
            "total_isomorphisms": result.nb_isomorphisms(),
            "locations": result.nb_locations(),
            "orientations_per_location": result.nb_ori_at_loc(),
        },
        indent=2,
    )


@mcp.tool(
    name="mofforge_build",
    description=(
        "Build a MOF structure from a topology and building blocks using "
        "TOBACCO or Pormake backends.  Node and edge building blocks are "
        "specified as file paths."
    ),
)
def mofforge_build(
    topology: str,
    backend: str = "tobacco",
    node_files: list[str] | None = None,
    edge_files: list[str] | None = None,
    output_dir: str = ".",
) -> str:
    """Build a MOF from topology + building blocks.

    Parameters
    ----------
    topology : str
        Topology name (e.g. "pcu", "pcu.cif").
    backend : str
        Builder backend: "tobacco" or "pormake".
    node_files : list[str] or None
        Paths to node building block files.
    edge_files : list[str] or None
        Paths to edge building block files.
    output_dir : str
        Output directory for CIF files.

    Returns
    -------
    str
        JSON summary of the build result.
    """
    from mofforge.build import MOFBuilder

    try:
        builder = MOFBuilder(backend=backend)
    except Exception as exc:
        return json.dumps(
            {
                "success": False,
                "error": f"Failed to initialize {backend} backend: {exc}",
            },
            indent=2,
        )

    for n in node_files or []:
        builder.add_node(n)
    for e in edge_files or []:
        builder.add_edge(e)

    out = _resolve_output(os.path.join(output_dir, "placeholder.cif"))
    out_dir = str(Path(out).parent)

    result = builder.build(topology=topology, output_dir=out_dir)

    if result.success:
        response = {
            "success": True,
            "topology": topology,
            "backend": backend,
            "elapsed_seconds": result.elapsed_seconds,
            "output_paths": [str(p) for p in result.output_paths],
            "atoms": result.crystal.n_atoms if result.crystal else None,
        }
    else:
        response = {
            "success": False,
            "topology": topology,
            "backend": backend,
            "errors": result.errors,
        }

    return json.dumps(response, indent=2)


@mcp.tool(
    name="mofforge_list_topologies",
    description="List available topology templates for MOF construction.",
)
def mofforge_list_topologies(
    backend: str = "tobacco",
) -> str:
    """List available topologies.

    Parameters
    ----------
    backend : str
        Builder backend: "tobacco" or "pormake".

    Returns
    -------
    str
        JSON list of topology names.
    """
    from mofforge.build import MOFBuilder

    try:
        builder = MOFBuilder(backend=backend)
        topologies = builder.list_topologies()
        return json.dumps(
            {
                "success": True,
                "backend": backend,
                "count": len(topologies),
                "topologies": [str(t) for t in topologies],
            },
            indent=2,
        )
    except Exception as exc:
        return json.dumps(
            {
                "success": False,
                "error": str(exc),
            },
            indent=2,
        )


@mcp.tool(
    name="mofforge_list_building_blocks",
    description=("List available building blocks (nodes or edges) for MOF construction."),
)
def mofforge_list_building_blocks(
    block_type: str = "nodes",
    backend: str = "tobacco",
) -> str:
    """List available building blocks.

    Parameters
    ----------
    block_type : str
        Type of building blocks: "nodes" or "edges".
    backend : str
        Builder backend: "tobacco" or "pormake".

    Returns
    -------
    str
        JSON list of building block names.
    """
    from mofforge.build import MOFBuilder

    try:
        builder = MOFBuilder(backend=backend)

        if block_type == "nodes":
            blocks = builder.list_nodes()
        elif block_type == "edges":
            blocks = builder.list_edges()
        else:
            return json.dumps(
                {
                    "success": False,
                    "error": f"Unknown block_type: {block_type!r}. Use 'nodes' or 'edges'.",
                },
                indent=2,
            )

        return json.dumps(
            {
                "success": True,
                "backend": backend,
                "block_type": block_type,
                "count": len(blocks),
                "blocks": [str(b) for b in blocks],
            },
            indent=2,
        )
    except Exception as exc:
        return json.dumps(
            {
                "success": False,
                "error": str(exc),
            },
            indent=2,
        )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main():
    """Run the mofforge MCP server (stdio transport)."""
    import argparse

    parser = argparse.ArgumentParser(description="mofforge MCP server")
    parser.add_argument(
        "--transport",
        choices=["stdio", "streamable-http"],
        default="stdio",
        help="MCP transport mode (default: stdio).",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=9010,
        help="HTTP port (only for streamable-http transport).",
    )
    args = parser.parse_args()

    if args.transport == "stdio":
        mcp.run(transport="stdio")
    else:
        mcp.run(transport="streamable-http", port=args.port)


if __name__ == "__main__":
    main()
