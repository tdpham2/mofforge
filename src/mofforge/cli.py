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

    # Load parent
    click.echo(f"Loading parent: {parent}")
    xtal = _load_parent(parent)

    # Load query
    q = _load_frag(query, fragment_path=fragment_path)

    # Search
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

    # Load parent
    click.echo(f"Loading parent: {parent}")
    xtal = _load_parent(parent)

    # Load fragments
    q = _load_frag(query, fragment_path=fragment_path)
    r = _load_frag(replacement, fragment_path=fragment_path)

    # Search
    result = find_pattern(q, xtal)
    click.echo(f"Found {result.nb_isomorphisms()} matches at {result.nb_locations()} locations")

    # Replace
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

    # Load parent
    click.echo(f"Loading parent: {parent}")
    xtal = _load_parent(parent)

    # Load guest
    g = _load_frag(guest, fragment_path=fragment_path)

    # Search for guest molecules
    result = find_pattern(g, xtal, disconnected_component=True)
    click.echo(f"Found {result.nb_locations()} guest molecule(s)")

    # Remove
    child = replace_pattern(result, None)
    child.write_cif(output)
    click.echo(f"Output written to: {output}")
    click.echo(f"  Atoms: {child.n_atoms} (removed {xtal.n_atoms - child.n_atoms})")


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

    # Print summary
    click.echo(f"\nBatch Results ({len(results)} structures):")
    for r in results:
        status = "OK" if r.success else f"FAILED: {r.error}"
        click.echo(f"  {r.parent_name}: {status}")
        if r.output_path:
            click.echo(f"    -> {r.output_path}")


if __name__ == "__main__":
    main()
