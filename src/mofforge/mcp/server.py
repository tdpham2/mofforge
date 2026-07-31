"""mofforge MCP Server exposing crystal manipulation tools for AI agents."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path

from mcp.server.fastmcp import FastMCP

logger = logging.getLogger("mofforge.mcp")

mcp = FastMCP(
    name="mofforge",
    instructions=(
        "mofforge tools for building and modifying atomistic crystal structure "
        "models, especially Metal-Organic Frameworks (MOFs).  Capabilities "
        "include substructure search (VF2 graph isomorphism), find-and-replace "
        "of molecular fragments, MOF construction from topology + building "
        "blocks, structure validation, SMARTS-like pattern matching, "
        "structure rendering to PNG images for visual inspection, CoRE MOF / "
        "CSD database search and screening, structure-file retrieval, and "
        "adsorbate placement.\n\n"
        "General guidance:\n"
        "- File paths should be absolute.\n"
        "- Crystal structures are read from CIF files.\n"
        "- Fragments are read from XYZ files with optional '!' anchor tags.\n"
        "- Distances are in Angstroms.\n"
        "- Typical screening workflow: mofforge_screen_coremof to find "
        "candidate MOFs by pore size / stability, mofforge_get_structure to "
        "fetch a CIF, then mofforge_place_adsorbate before a simulation.\n"
        "- Call mofforge_validate after modifications to check for issues."
    ),
)



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
    name="mofforge_desolvate",
    description=(
        "Automatically identify and remove all uncoordinated solvent/guest "
        "molecules from a MOF crystal structure. Unlike mofforge_remove, "
        "this does not require specifying which solvent to remove."
    ),
)
def mofforge_desolvate(
    parent_cif: str,
    output_cif: str = "desolvated.cif",
    min_atoms: int = 1,
    keep_metal_containing: bool = False,
    n_framework_components: int | None = None,
) -> str:
    """Automatically remove solvent from a crystal structure.

    Parameters
    ----------
    parent_cif : str
        Absolute path to the parent crystal CIF file.
    output_cif : str
        Output CIF file path.
    min_atoms : int
        Minimum atoms for a component to be considered framework.
    keep_metal_containing : bool
        If True, keep components containing metal atoms.
    n_framework_components : int | None
        Number of framework components to keep (None = auto-detect).
    """
    from collections import Counter

    from mofforge.solvent.removal import remove_solvent

    xtal = _load_crystal(parent_cif)
    result = remove_solvent(
        xtal,
        min_atoms=min_atoms,
        keep_metal_containing=keep_metal_containing,
        n_framework_components=n_framework_components,
    )

    output_path = _resolve_output(output_cif)
    result.crystal.write_cif(output_path)

    formula_counts = Counter(m.formula for m in result.removed_molecules)

    return json.dumps(
        {
            "success": True,
            "output_cif": output_path,
            "atoms_before": result.n_atoms_original,
            "atoms_after": result.crystal.n_atoms,
            "atoms_removed": result.n_atoms_removed,
            "molecules_removed": result.n_components_removed,
            "framework_components": result.n_framework_components,
            "removed_species": dict(formula_counts),
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



@mcp.tool(
    name="mofforge_search_coremof",
    description=(
        "Search the CoRE MOF database of simulation-ready MOF structures.  "
        "Auto-detects query type (coreid, refcode, name, DOI, metal, "
        "topology).  Returns matching entries with pore-size, density, "
        "stability, and open-metal-site metadata."
    ),
)
def mofforge_search_coremof(
    query: str,
    field: str = "auto",
    limit: int = 25,
    data_path: str | None = None,
) -> str:
    """Search the CoRE MOF database.

    Parameters
    ----------
    query : str
        Search term (coreid, refcode, name, DOI, metal symbol, or topology).
    field : str
        One of "auto", "coreid", "refcode", "name", "doi", "metal", "topology".
    limit : int
        Maximum number of results.
    data_path : str or None
        Optional path to the CoRE MOF metadata CSV (else resolved from config).
    """
    from mofforge.mcp._impl import search_coremof_impl

    return json.dumps(
        search_coremof_impl(query, field=field, limit=limit, data_path=data_path),
        indent=2,
    )


@mcp.tool(
    name="mofforge_screen_coremof",
    description=(
        "Screen the CoRE MOF database by property ranges and filters, e.g. "
        "pore size (lcd/pld), density, accessible surface area, void "
        "fraction, water/thermal stability, metal, topology, and open metal "
        "sites.  Use this to shortlist MOFs before simulation."
    ),
)
def mofforge_screen_coremof(
    lcd_min: float | None = None,
    lcd_max: float | None = None,
    pld_min: float | None = None,
    pld_max: float | None = None,
    density_min: float | None = None,
    density_max: float | None = None,
    asa_min: float | None = None,
    asa_max: float | None = None,
    void_fraction_min: float | None = None,
    void_fraction_max: float | None = None,
    water_stability_min: float | None = None,
    thermal_stability_min: float | None = None,
    metal: str | None = None,
    topology: str | None = None,
    has_oms: bool | None = None,
    extension: str | None = None,
    limit: int = 50,
    data_path: str | None = None,
) -> str:
    """Screen CoRE MOFs by property ranges and categorical filters.

    Parameters
    ----------
    lcd_min, lcd_max : float or None
        Largest cavity diameter range (Angstroms).
    pld_min, pld_max : float or None
        Pore limiting diameter range (Angstroms).
    density_min, density_max : float or None
        Crystal density range (g/cm^3).
    asa_min, asa_max : float or None
        Accessible surface area range (m^2/g).
    void_fraction_min, void_fraction_max : float or None
        Void fraction range (0-1).
    water_stability_min, thermal_stability_min : float or None
        Minimum stability scores.
    metal : str or None
        Restrict to MOFs containing this metal element.
    topology : str or None
        Restrict to a topology (single-node), e.g. "pcu".
    has_oms : bool or None
        Filter on presence of open metal sites.
    extension : str or None
        CoRE MOF processing variant (e.g. "ASR", "FSR", "ION").
    limit : int
        Maximum number of results.
    data_path : str or None
        Optional path to the CoRE MOF metadata CSV.
    """
    from mofforge.mcp._impl import screen_coremof_impl

    return json.dumps(
        screen_coremof_impl(
            lcd_min=lcd_min,
            lcd_max=lcd_max,
            pld_min=pld_min,
            pld_max=pld_max,
            density_min=density_min,
            density_max=density_max,
            asa_min=asa_min,
            asa_max=asa_max,
            void_fraction_min=void_fraction_min,
            void_fraction_max=void_fraction_max,
            water_stability_min=water_stability_min,
            thermal_stability_min=thermal_stability_min,
            metal=metal,
            topology=topology,
            has_oms=has_oms,
            extension=extension,
            limit=limit,
            data_path=data_path,
        ),
        indent=2,
    )


@mcp.tool(
    name="mofforge_search_csd",
    description=(
        "Search the Cambridge Structural Database (CSD) lookup table.  "
        "Auto-detects query type (refcode, name, DOI, formula, CCDC number)."
    ),
)
def mofforge_search_csd(
    query: str,
    field: str = "auto",
    limit: int = 25,
    data_path: str | None = None,
) -> str:
    """Search the CSD database.

    Parameters
    ----------
    query : str
        Search term (refcode, name, DOI, formula, or CCDC number).
    field : str
        One of "auto", "refcode", "name", "doi", "formula", "ccdc".
    limit : int
        Maximum number of results.
    data_path : str or None
        Optional path to the CSD data file (else resolved from config).
    """
    from mofforge.mcp._impl import search_csd_impl

    return json.dumps(
        search_csd_impl(query, field=field, limit=limit, data_path=data_path),
        indent=2,
    )


@mcp.tool(
    name="mofforge_lookup_mof",
    description=(
        "Look up a MOF by name in the CSD and return the corresponding "
        "simulation-ready CoRE MOF entries in one step (CSD->CoRE bridge)."
    ),
)
def mofforge_lookup_mof(
    name: str,
    limit: int = 25,
    csd_data_path: str | None = None,
    coremof_data_path: str | None = None,
) -> str:
    """Bridge a MOF name from CSD to CoRE MOF entries.

    Parameters
    ----------
    name : str
        MOF name to search in the CSD (substring match).
    limit : int
        Maximum number of CSD records to consider.
    csd_data_path : str or None
        Optional path to the CSD data file.
    coremof_data_path : str or None
        Optional path to the CoRE MOF metadata CSV.
    """
    from mofforge.mcp._impl import lookup_mof_impl

    return json.dumps(
        lookup_mof_impl(
            name,
            limit=limit,
            csd_data_path=csd_data_path,
            coremof_data_path=coremof_data_path,
        ),
        indent=2,
    )


@mcp.tool(
    name="mofforge_get_structure",
    description=(
        "Resolve a CoRE MOF coreid or refcode to a local CIF structure file "
        "path, suitable for loading into ASE / gRASPA / adsorbate placement.  "
        "Requires a configured CoRE MOF structures directory (CIFs are "
        "distributed separately from the metadata)."
    ),
)
def mofforge_get_structure(
    identifier: str,
    structures_dir: str | None = None,
) -> str:
    """Resolve a CoRE MOF identifier to a local CIF path.

    Parameters
    ----------
    identifier : str
        A CoRE MOF coreid (preferred) or refcode.
    structures_dir : str or None
        Directory containing CoRE MOF CIF files (else resolved from config /
        MOFFORGE_COREMOF_STRUCTURES_PATH).
    """
    from mofforge.mcp._impl import get_structure_impl

    return json.dumps(
        get_structure_impl(identifier, structures_dir=structures_dir),
        indent=2,
    )


@mcp.tool(
    name="mofforge_place_adsorbate",
    description=(
        "Place one or more adsorbate molecules (e.g. CO2, H2, CH4, N2, H2O) "
        "into a MOF crystal structure and write the result to a CIF file.  "
        "Auto-detects void or open-metal adsorption sites."
    ),
)
def mofforge_place_adsorbate(
    cif_path: str,
    adsorbate: str = "CO2",
    n_adsorbates: int = 1,
    strategy: str = "void",
    output_cif: str = "with_adsorbate.cif",
    validate: bool = True,
    random_seed: int | None = None,
) -> str:
    """Place adsorbate(s) into a MOF.

    Parameters
    ----------
    cif_path : str
        Absolute path to the host MOF CIF file.
    adsorbate : str
        Built-in adsorbate name (see mofforge_list_adsorbates).
    n_adsorbates : int
        Number of adsorbate molecules to place.
    strategy : str
        Site-detection strategy: "void" or "oms" (open metal sites).
    output_cif : str
        Output CIF file path.
    validate : bool
        Check for steric clashes after placement.
    random_seed : int or None
        Seed for reproducible orientations/placement.
    """
    from mofforge.mcp._impl import place_adsorbate_impl

    return json.dumps(
        place_adsorbate_impl(
            cif_path,
            adsorbate=adsorbate,
            n_adsorbates=n_adsorbates,
            strategy=strategy,
            output_cif=output_cif,
            validate=validate,
            random_seed=random_seed,
        ),
        indent=2,
    )


@mcp.tool(
    name="mofforge_list_adsorbates",
    description="List the built-in adsorbate molecules available for placement.",
)
def mofforge_list_adsorbates() -> str:
    """List built-in adsorbate molecule names."""
    from mofforge.mcp._impl import list_adsorbates_impl

    return json.dumps(list_adsorbates_impl(), indent=2)


def main():
    """Run the mofforge MCP server (stdio transport by default).

    ``--transport`` accepts both ``streamable-http`` (mofforge legacy) and
    ``streamable_http`` (ChemGraph/swarm supervisor convention) so this
    binary works under either launcher.
    """
    import argparse

    parser = argparse.ArgumentParser(description="mofforge MCP server")
    parser.add_argument(
        "--transport",
        choices=["stdio", "streamable-http", "streamable_http"],
        default="stdio",
        help="MCP transport mode (default: stdio).",
    )
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="Host for HTTP transport (default: 127.0.0.1).",
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
        # FastMCP's transport enum expects the dash form; normalize.
        mcp.settings.host = args.host
        mcp.settings.port = args.port
        mcp.run(transport="streamable-http")


if __name__ == "__main__":
    main()
