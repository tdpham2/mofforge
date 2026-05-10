"""Command-line interface for mofforge."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import click

logger = logging.getLogger("mofforge")


def _setup_logging(verbose: bool) -> None:
    """Configure logging for CLI usage."""
    level = logging.DEBUG if verbose else logging.INFO
    root_logger = logging.getLogger("mofforge")
    root_logger.setLevel(level)
    # Avoid adding duplicate handlers on repeated calls (e.g. in tests)
    if not root_logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
        root_logger.addHandler(handler)


def _load_parent(parent_path: str) -> "Crystal":  # noqa: F821
    """Load a parent crystal from CIF and infer bonds."""
    from mofforge.core.bonding import infer_bonds
    from mofforge.core.crystal import Crystal

    xtal = Crystal.from_cif(parent_path)
    xtal = infer_bonds(xtal, periodic=True)
    return xtal


def _load_frag(filepath: str, fragment_path: str | None = None) -> "Crystal":  # noqa: F821
    """Load a fragment from an XYZ file, resolving the path."""
    from mofforge.core.moiety import fragment as load_fragment

    p = Path(filepath)
    mp = fragment_path if fragment_path else str(p.parent)
    return load_fragment(p.name, fragment_path=mp)


@click.group()
@click.version_option(package_name="mofforge", prog_name="mofforge")
def main():
    """mofforge: Find-and-replace tool for crystal structures."""
    pass


@main.command()
@click.option("-p", "--parent", required=True, help="Path to parent CIF file.")
@click.option("-q", "--query", required=True, help="Path to query fragment XYZ file.")
@click.option("--disconnected", is_flag=True, help="Search for isolated components only.")
@click.option("--fragment-path", default=None, help="Directory containing fragment XYZ files.")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
def search(parent, query, disconnected, fragment_path, verbose):
    """Search for a pattern in a crystal."""
    _setup_logging(verbose)

    from mofforge.search.search import find_pattern

    click.echo(f"Loading parent: {parent}")
    xtal = _load_parent(parent)

    q = _load_frag(query, fragment_path=fragment_path)

    result = find_pattern(q, xtal, disconnected_component=disconnected)

    click.echo(f"Results: {result}")
    click.echo(f"  Isomorphisms: {result.nb_isomorphisms()}")
    click.echo(f"  Locations: {result.nb_locations()}")
    click.echo(f"  Orientations per location: {result.nb_ori_at_loc()}")


@main.command("replace")
@click.option("-p", "--parent", required=True, help="Path to parent CIF file.")
@click.option("-q", "--query", required=True, help="Path to query fragment XYZ file.")
@click.option("-r", "--replacement", required=True, help="Path to replacement fragment XYZ file.")
@click.option("-o", "--output", default="new_xtal.cif", help="Output CIF file path.")
@click.option("--nb-loc", default=0, type=int, help="Number of random locations.")
@click.option("--random", "use_random", is_flag=True, help="Use random orientations.")
@click.option("--validate", "do_validate", is_flag=True, help="Validate output structure.")
@click.option("--fragment-path", default=None, help="Directory containing fragment XYZ files.")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
def replace_cmd(
    parent, query, replacement, output, nb_loc, use_random, do_validate, fragment_path, verbose
):
    """Find and replace a pattern in a crystal."""
    _setup_logging(verbose)

    from mofforge.replace.replace import replace_pattern
    from mofforge.search.search import find_pattern

    click.echo(f"Loading parent: {parent}")
    xtal = _load_parent(parent)

    q = _load_frag(query, fragment_path=fragment_path)
    r = _load_frag(replacement, fragment_path=fragment_path)

    result = find_pattern(q, xtal)
    click.echo(f"Found {result.nb_isomorphisms()} matches at {result.nb_locations()} locations")

    kwargs = {"verbose": True}
    if nb_loc > 0:
        kwargs["nb_loc"] = nb_loc
    if use_random:
        kwargs["random"] = True

    child = replace_pattern(result, r, **kwargs)
    child.write_cif(output)
    click.echo(f"Output written to: {output}")
    click.echo(f"  Atoms: {child.n_atoms}, Bonds: {child.n_bonds}")

    # Validate if requested
    if do_validate:
        from mofforge.core.bonding import infer_bonds
        from mofforge.validation import validate_structure

        child = infer_bonds(child, periodic=True)
        report = validate_structure(child)
        click.echo(report.summary())


@main.command("remove")
@click.option("-p", "--parent", required=True, help="Path to parent CIF file.")
@click.option("-g", "--guest", required=True, help="Path to guest fragment XYZ file.")
@click.option("-o", "--output", default="clean.cif", help="Output CIF file path.")
@click.option("--fragment-path", default=None, help="Directory containing fragment XYZ files.")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
def remove_cmd(parent, guest, output, fragment_path, verbose):
    """Remove guest molecules from a crystal structure."""
    _setup_logging(verbose)

    from mofforge.replace.replace import replace_pattern
    from mofforge.search.search import find_pattern

    click.echo(f"Loading parent: {parent}")
    xtal = _load_parent(parent)

    g = _load_frag(guest, fragment_path=fragment_path)

    result = find_pattern(g, xtal, disconnected_component=True)
    click.echo(f"Found {result.nb_locations()} guest molecule(s)")

    child = replace_pattern(result, None)
    child.write_cif(output)
    click.echo(f"Output written to: {output}")
    click.echo(f"  Atoms: {child.n_atoms} (removed {xtal.n_atoms - child.n_atoms})")


@main.command("desolvate")
@click.option("-p", "--parent", required=True, help="Path to parent CIF file.")
@click.option("-o", "--output", default="desolvated.cif", help="Output CIF file path.")
@click.option("--min-atoms", default=1, type=int, help="Min atoms for a component to be framework.")
@click.option("--keep-metals", is_flag=True, help="Keep metal-containing components.")
@click.option(
    "--n-frameworks", default=None, type=int, help="Number of framework components to keep."
)
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
def desolvate_cmd(parent, output, min_atoms, keep_metals, n_frameworks, verbose):
    """Remove all uncoordinated solvent molecules from a crystal structure."""
    _setup_logging(verbose)

    from collections import Counter

    from mofforge.solvent.removal import remove_solvent

    click.echo(f"Loading parent: {parent}")
    xtal = _load_parent(parent)

    result = remove_solvent(
        xtal,
        min_atoms=min_atoms,
        keep_metal_containing=keep_metals,
        n_framework_components=n_frameworks,
    )

    result.crystal.write_cif(output)
    click.echo(f"Output written to: {output}")
    click.echo(f"  Atoms: {result.crystal.n_atoms} (removed {result.n_atoms_removed})")
    click.echo(f"  Molecules removed: {result.n_components_removed}")
    click.echo(f"  Framework components: {result.n_framework_components}")

    if verbose and result.removed_molecules:
        click.echo("  Removed species:")
        formula_counts = Counter(m.formula for m in result.removed_molecules)
        for formula, count in formula_counts.most_common():
            click.echo(f"    {count}x {formula}")


@main.command("validate")
@click.argument("structure")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
def validate_cmd(structure, verbose):
    """Validate a crystal structure."""
    _setup_logging(verbose)

    from mofforge.validation import validate_structure

    click.echo(f"Loading structure: {structure}")
    xtal = _load_parent(structure)

    report = validate_structure(xtal)
    click.echo(report.summary())


@main.command("batch")
@click.option("-c", "--config", "config_path", required=True, help="Path to YAML config file.")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
def batch_cmd(config_path, verbose):
    """Run batch processing from a YAML config file."""
    _setup_logging(verbose)

    from mofforge.batch import run_batch

    results = run_batch(config_path)

    click.echo(f"\nBatch Results ({len(results)} structures):")
    for r in results:
        status = "OK" if r.success else f"FAILED: {r.error}"
        click.echo(f"  {r.parent_name}: {status}")
        if r.output_path:
            click.echo(f"    -> {r.output_path}")


@main.command("render")
@click.option(
    "-i", "--input", "input_file", required=True, help="Path to input structure file (CIF or XYZ)."
)
@click.option(
    "-o", "--output", "output_file", default="structure.png", help="Output PNG file path."
)
@click.option(
    "--label-mode",
    type=click.Choice(["sequential", "per_element", "none"]),
    default="sequential",
    help="Atom label mode.",
)
@click.option(
    "--representation",
    type=click.Choice(["ball_stick", "stick", "sphere"]),
    default="ball_stick",
    help="Visual representation style.",
)
@click.option("--width", default=800, type=int, help="Image width in pixels.")
@click.option("--height", default=600, type=int, help="Image height in pixels.")
@click.option(
    "--show-unit-cell/--no-unit-cell",
    default=None,
    help="Draw unit cell edges (auto for CIF files).",
)
@click.option("--show-formula/--no-formula", default=True, help="Show chemical formula overlay.")
@click.option("--bg-color", default="white", help="Background color.")
@click.option("--label-size", default=14, type=int, help="Font size for atom labels.")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
def render_cmd(
    input_file,
    output_file,
    label_mode,
    representation,
    width,
    height,
    show_unit_cell,
    show_formula,
    bg_color,
    label_size,
    verbose,
):
    """Render a structure file to a PNG image."""
    _setup_logging(verbose)

    from mofforge.vis.render import render_file_to_png

    click.echo(f"Rendering: {input_file}")

    # If show_unit_cell is None (not explicitly set), let render_file_to_png
    # auto-detect based on file extension.
    kwargs = {}
    if show_unit_cell is not None:
        kwargs["show_unit_cell"] = show_unit_cell

    output_path = render_file_to_png(
        input_file=input_file,
        output_file=output_file,
        label_mode=label_mode,
        width=width,
        height=height,
        representation=representation,
        label_size=label_size,
        show_formula=show_formula,
        bg_color=bg_color,
        **kwargs,
    )
    click.echo(f"Output written to: {output_path}")


@main.command("build")
@click.option(
    "--backend",
    "-b",
    type=click.Choice(["tobacco", "pormake"]),
    default="tobacco",
    help="Builder backend to use.",
)
@click.option("--topology", "-t", required=True, help="Topology name (e.g. 'pcu' or 'pcu.cif').")
@click.option("--node", "-n", multiple=True, help="Node building-block file(s).")
@click.option("--edge", "-e", multiple=True, help="Edge building-block file(s).")
@click.option("-o", "--output", "output_dir", default=".", help="Output directory for CIF files.")
@click.option("--tobacco-path", default=None, help="Path to TOBACCO 3.0 directory.")
@click.option("--parallel", is_flag=True, help="Run TOBACCO in parallel mode.")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
def build_cmd(backend, topology, node, edge, output_dir, tobacco_path, parallel, verbose):
    """Build a MOF structure from topology + building blocks."""
    _setup_logging(verbose)

    from mofforge.build import MOFBuilder

    kwargs = {}
    if tobacco_path:
        kwargs["tobacco_path"] = tobacco_path

    try:
        builder = MOFBuilder(backend=backend, **kwargs)
    except Exception as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    for n in node:
        builder.add_node(n)
    for e in edge:
        builder.add_edge(e)

    click.echo(f"Building MOF with {backend} backend (topology={topology})")
    build_kwargs = {}
    if parallel:
        build_kwargs["parallel"] = True

    result = builder.build(topology=topology, output_dir=output_dir, **build_kwargs)

    if result.success:
        click.echo(f"Build succeeded in {result.elapsed_seconds}s")
        for p in result.output_paths:
            click.echo(f"  Output: {p}")
        if result.crystal:
            click.echo(f"  Atoms: {result.crystal.n_atoms}")
    else:
        click.echo("Build failed:")
        for err in result.errors:
            click.echo(f"  {err}")
        sys.exit(1)


@main.command("build-status")
@click.option(
    "--backend",
    "-b",
    type=click.Choice(["tobacco", "pormake"]),
    default="tobacco",
    help="Builder backend.",
)
@click.option("--tobacco-path", default=None, help="Path to TOBACCO 3.0 directory.")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
def build_status_cmd(backend, tobacco_path, verbose):
    """Show status of a builder backend."""
    _setup_logging(verbose)

    import json

    from mofforge.build import MOFBuilder

    kwargs = {}
    if tobacco_path:
        kwargs["tobacco_path"] = tobacco_path

    try:
        builder = MOFBuilder(backend=backend, **kwargs)
    except Exception as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    status = builder.status()
    click.echo(json.dumps(status, indent=2, default=str))


@main.command("build-list")
@click.option(
    "--backend",
    "-b",
    type=click.Choice(["tobacco", "pormake"]),
    default="tobacco",
    help="Builder backend.",
)
@click.option(
    "--type",
    "list_type",
    type=click.Choice(["topologies", "nodes", "edges"]),
    required=True,
    help="What to list.",
)
@click.option("--tobacco-path", default=None, help="Path to TOBACCO 3.0 directory.")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
def build_list_cmd(backend, list_type, tobacco_path, verbose):
    """List available topologies or building blocks."""
    _setup_logging(verbose)

    from mofforge.build import MOFBuilder

    kwargs = {}
    if tobacco_path:
        kwargs["tobacco_path"] = tobacco_path

    try:
        builder = MOFBuilder(backend=backend, **kwargs)
    except Exception as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    if list_type == "topologies":
        items = builder.list_topologies()
        click.echo(f"Available topologies ({len(items)}):")
    elif list_type == "nodes":
        items = builder.list_nodes()
        click.echo(f"Available nodes ({len(items)}):")
    elif list_type == "edges":
        items = builder.list_edges()
        click.echo(f"Available edges ({len(items)}):")

    for item in items:
        click.echo(f"  {item}")


@main.command("csd")
@click.argument("query")
@click.option(
    "--field",
    "-f",
    type=click.Choice(["auto", "refcode", "name", "doi", "formula", "ccdc"]),
    default="auto",
    help="Field to search (default: auto-detect).",
)
@click.option("--limit", "-n", default=50, type=int, help="Max results to display (default: 50).")
@click.option("--data-path", default=None, help="Path to CSD TSV file (overrides config).")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output.")
def csd_cmd(query, field, limit, data_path, verbose):
    """Look up MOF entries in the CSD database by name, DOI, or REFcode."""
    _setup_logging(verbose)

    from mofforge.csd import get_database
    from mofforge.utils.config import set_paths

    if data_path:
        set_paths(csd_data=data_path)

    try:
        db = get_database()
    except FileNotFoundError as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    # Get total count, then limited results
    all_result = db.search(query, field=field)
    total = all_result.n_matches
    result = db.search(query, field=field, limit=limit)

    click.echo(f"CSD lookup: {total} match(es) for '{query}' (field: {result.field})")
    if total > limit:
        click.echo(f"  (showing first {limit} of {total}; use -n {total} to see all)")
    for rec in result.records:
        click.echo(f"  {rec.refcode}: {rec.chemical_name_common or rec.chemical_name_systematic}")
        if rec.doi:
            click.echo(f"    DOI: {rec.doi}")
        if rec.ccdc_number:
            click.echo(f"    CCDC: {rec.ccdc_number}")
        if verbose:
            click.echo(f"    Formula: {rec.chemical_formula_moiety}")
            click.echo(f"    Journal: {rec.journal} ({rec.year})")
            click.echo(f"    Space group: {rec.space_group}")


@main.command("coremof")
@click.argument("query")
@click.option(
    "--field",
    "-f",
    type=click.Choice(["auto", "coreid", "refcode", "name", "doi", "metal", "topology"]),
    default="auto",
    help="Field to search (default: auto-detect).",
)
@click.option("--limit", "-n", default=50, type=int, help="Max results to display (default: 50).")
@click.option("--data-path", default=None, help="Path to CoRE MOF CSV file (overrides config).")
@click.option("--bridge", is_flag=True, help="Treat QUERY as CSD refcode, return CoreMOF entries.")
@click.option("-v", "--verbose", is_flag=True, help="Show extended properties.")
def coremof_cmd(query, field, limit, data_path, bridge, verbose):
    """Search the CoRE MOF database for simulation-ready structures."""
    _setup_logging(verbose)

    from mofforge.coremof import get_database, csd_to_coremof
    from mofforge.utils.config import set_paths

    if data_path:
        set_paths(coremof_data=data_path)

    try:
        db = get_database()
    except FileNotFoundError as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    if bridge:
        records = csd_to_coremof(query, db=db)
        click.echo(f"CoreMOF bridge: {len(records)} match(es) for CSD refcode '{query}'")
        for rec in records:
            click.echo(f"  {rec.coreid}")
            click.echo(f"    Refcode:   {rec.refcode}")
            click.echo(f"    Extension: {rec.extension}")
            click.echo(f"    Metals:    {rec.metal_types}")
            if verbose:
                click.echo(rec.properties_summary())
        return

    # Get total count, then limited results
    all_result = db.search(query, field=field)
    total = all_result.n_matches
    result = db.search(query, field=field, limit=limit)

    click.echo(
        f"CoreMOF lookup: {total} match(es) "
        f"for '{query}' (field: {result.field})"
    )
    if total > limit:
        click.echo(f"  (showing first {limit} of {total}; use -n {total} to see all)")
    for rec in result.records:
        click.echo(f"  {rec.summary()}")
        if verbose:
            click.echo(rec.properties_summary())
            click.echo()


@main.command("lookup")
@click.argument("name")
@click.option("--limit", "-n", default=50, type=int, help="Max CSD results to display (default: 50).")
@click.option("--csd-data-path", default=None, help="Path to CSD TSV file.")
@click.option("--coremof-data-path", default=None, help="Path to CoRE MOF CSV file.")
@click.option("-v", "--verbose", is_flag=True, help="Show extended CoreMOF properties.")
def lookup_cmd(name, limit, csd_data_path, coremof_data_path, verbose):
    """Search CSD by MOF name and return CoreMOF simulation IDs.

    Chains a CSD name search with CoreMOF refcode lookup: for each CSD
    record found, displays all corresponding CoreMOF entries with their
    coreid identifiers for simulation.
    """
    _setup_logging(verbose)

    from mofforge.coremof import search_csd_name
    from mofforge.coremof import get_database as get_coremof_db
    from mofforge.csd import get_database as get_csd_db
    from mofforge.utils.config import set_paths

    if csd_data_path:
        set_paths(csd_data=csd_data_path)
    if coremof_data_path:
        set_paths(coremof_data=coremof_data_path)

    try:
        csd_db = get_csd_db(data_path=csd_data_path) if csd_data_path else get_csd_db()
        coremof_db = get_coremof_db(data_path=coremof_data_path) if coremof_data_path else get_coremof_db()
        # Run without limit to get total count
        all_results = search_csd_name(name, coremof_db=coremof_db, csd_db=csd_db)
    except FileNotFoundError as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    total_csd = len(all_results)
    results = all_results[:limit]
    total_coremof = sum(len(br.coremof_records) for br in results)
    click.echo(
        f"Lookup '{name}': {total_csd} total CSD match(es), "
        f"showing {len(results)}, "
        f"{total_coremof} CoreMOF entry(ies)"
    )
    if total_csd > limit:
        click.echo(
            f"  (showing first {limit} of {total_csd} CSD matches; "
            f"use -n {total_csd} to see all)"
        )
    click.echo()

    for br in results:
        click.echo(br.summary())
        if verbose and br.coremof_records:
            for rec in br.coremof_records:
                click.echo(rec.properties_summary())
        click.echo()


if __name__ == "__main__":
    main()
