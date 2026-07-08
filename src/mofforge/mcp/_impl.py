"""Backend-agnostic implementations of mofforge MCP tools.

Each function here returns a plain ``dict`` (never a JSON string) and never
raises out to the caller: failures are reported as ``{"success": False,
"error": ...}``. This keeps the logic reusable by both:

* the stock FastMCP server (``mofforge.mcp.server``), which JSON-encodes the
  return value, and
* the ChemGraph ``CGFastMCP`` adapter (``mofforge.mcp.chemgraph_server``),
  which may submit these functions to an HPC execution backend.

Because the CGFastMCP backend pickles tool functions *by reference*, every
function here is a module-level callable with no closures over server state.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any

logger = logging.getLogger("mofforge.mcp")

__all__ = [
    "build_impl",
    "find_sites_impl",
    "functionalize_campaign_impl",
    "functionalize_impl",
    "get_fragment_impl",
    "get_structure_impl",
    "list_adsorbates_impl",
    "list_fragments_impl",
    "list_functional_groups_impl",
    "load_crystal",
    "load_fragment",
    "lookup_mof_impl",
    "place_adsorbate_impl",
    "render_impl",
    "resolve_output",
    "screen_coremof_impl",
    "search_coremof_impl",
    "search_csd_impl",
    "validate_impl",
]


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def resolve_output(path: str) -> str:
    """Ensure an output path is absolute and its parent directory exists."""
    p = Path(path)
    if not p.is_absolute():
        log_dir = os.environ.get("MOFFORGE_LOG_DIR", ".")
        p = Path(log_dir) / p
    p.parent.mkdir(parents=True, exist_ok=True)
    return str(p.resolve())


def load_crystal(cif_path: str, with_bonds: bool = True):
    """Load a Crystal from CIF, optionally inferring bonds."""
    from mofforge.core.bonding import infer_bonds
    from mofforge.core.crystal import Crystal

    xtal = Crystal.from_cif(cif_path)
    if with_bonds:
        xtal = infer_bonds(xtal, periodic=True)
    return xtal


def load_fragment(xyz_path: str):
    """Load a fragment from an XYZ file."""
    from mofforge.core.moiety import fragment as load_fragment_fn

    p = Path(xyz_path)
    return load_fragment_fn(p.name, fragment_path=str(p.parent))


def _coremof_record_dict(rec) -> dict[str, Any]:
    """Serialize the ChemGraph-relevant fields of a CoreMOFRecord."""
    return {
        "coreid": rec.coreid,
        "refcode": rec.refcode,
        "base_refcode": rec.base_refcode,
        "name": rec.name,
        "metal_types": rec.metal_types,
        "topology": rec.topology_single,
        "lcd": rec.lcd,
        "pld": rec.pld,
        "density": rec.density,
        "asa": rec.asa,
        "pore_volume": rec.pore_volume,
        "void_fraction": rec.void_fraction,
        "has_oms": rec.has_oms,
        "water_stability": rec.water_stability,
        "thermal_stability": rec.thermal_stability,
        "extension": rec.extension,
        "natoms": rec.natoms,
        "doi": rec.doi,
        "year": rec.year,
    }


# ---------------------------------------------------------------------------
# Database tools
# ---------------------------------------------------------------------------


def search_coremof_impl(
    query: str,
    field: str = "auto",
    limit: int = 25,
    data_path: str | None = None,
) -> dict[str, Any]:
    """Search the CoRE MOF database (auto-detects coreid/refcode/name/etc.)."""
    try:
        from mofforge.coremof import get_database

        db = get_database(data_path=data_path)
        result = db.search(query, field=field, limit=limit)
        return {
            "success": True,
            "query": query,
            "field": result.field,
            "n_matches": result.n_matches,
            "records": [_coremof_record_dict(r) for r in result.records],
        }
    except Exception as exc:
        logger.warning("search_coremof failed", exc_info=True)
        return {"success": False, "error": str(exc)}


def screen_coremof_impl(
    *,
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
) -> dict[str, Any]:
    """Screen CoRE MOFs by property ranges and categorical filters."""
    try:
        from mofforge.coremof import get_database

        db = get_database(data_path=data_path)
        records = db.screen(
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
        )
        return {
            "success": True,
            "n_matches": len(records),
            "records": [_coremof_record_dict(r) for r in records],
        }
    except Exception as exc:
        logger.warning("screen_coremof failed", exc_info=True)
        return {"success": False, "error": str(exc)}


def search_csd_impl(
    query: str,
    field: str = "auto",
    limit: int = 25,
    data_path: str | None = None,
) -> dict[str, Any]:
    """Search the CSD database (auto-detects refcode/name/doi/formula/etc.)."""
    try:
        from mofforge.csd import get_database

        db = get_database(data_path=data_path)
        result = db.search(query, field=field, limit=limit)
        records = [
            {
                "refcode": r.refcode,
                "name": r.chemical_name_common or r.chemical_name_systematic,
                "formula": r.chemical_formula_moiety,
                "doi": r.doi,
                "year": r.year,
                "space_group": r.space_group,
            }
            for r in result.records
        ]
        return {
            "success": True,
            "query": query,
            "field": result.field,
            "n_matches": len(records),
            "records": records,
        }
    except Exception as exc:
        logger.warning("search_csd failed", exc_info=True)
        return {"success": False, "error": str(exc)}


def lookup_mof_impl(
    name: str,
    limit: int = 25,
    csd_data_path: str | None = None,
    coremof_data_path: str | None = None,
) -> dict[str, Any]:
    """Search a MOF name in CSD and return its simulation-ready CoRE MOF entries."""
    try:
        from mofforge.coremof import get_database as get_coremof_db
        from mofforge.coremof import search_csd_name
        from mofforge.csd import get_database as get_csd_db

        csd_db = get_csd_db(data_path=csd_data_path)
        coremof_db = get_coremof_db(data_path=coremof_data_path)
        bridges = search_csd_name(
            name, coremof_db=coremof_db, csd_db=csd_db, limit=limit
        )
        results = []
        for b in bridges:
            results.append(
                {
                    "csd_refcode": b.csd_record.refcode,
                    "csd_name": (
                        b.csd_record.chemical_name_common
                        or b.csd_record.chemical_name_systematic
                    ),
                    "has_coremof": b.has_coremof,
                    "coremof": [_coremof_record_dict(r) for r in b.coremof_records],
                }
            )
        return {
            "success": True,
            "query": name,
            "n_csd_matches": len(results),
            "results": results,
        }
    except Exception as exc:
        logger.warning("lookup_mof failed", exc_info=True)
        return {"success": False, "error": str(exc)}


def get_structure_impl(
    identifier: str,
    structures_dir: str | None = None,
) -> dict[str, Any]:
    """Resolve a CoRE MOF coreid/refcode to a local CIF file path.

    The CoRE MOF metadata CSV is properties-only; CIF structure files are
    distributed separately. Configure their location via
    ``MOFFORGE_COREMOF_STRUCTURES_PATH`` (or ``set_paths(coremof_structures=...)``,
    or ``[coremof] structures_path`` in mofforge.toml).
    """
    try:
        from mofforge.coremof import resolve_structure_path, resolve_structures_dir
        from mofforge.coremof.database import ZENODO_URL

        path = resolve_structure_path(identifier, structures_dir=structures_dir)
        if path is not None:
            return {
                "success": True,
                "identifier": identifier,
                "cif_path": str(path),
            }

        base = resolve_structures_dir(structures_dir)
        if base is None:
            error = (
                "No CoRE MOF structures directory configured. Set "
                "MOFFORGE_COREMOF_STRUCTURES_PATH, call "
                "set_paths(coremof_structures=...), or add [coremof] "
                "structures_path to mofforge.toml. "
                f"Download CIF structures from: {ZENODO_URL}"
            )
        elif not base.exists():
            error = (
                f"Configured CoRE MOF structures directory does not exist: {base}. "
                f"Download CIF structures from: {ZENODO_URL}"
            )
        else:
            error = (
                f"No CIF file found for {identifier!r} under {base}. "
                f"Download CIF structures from: {ZENODO_URL}"
            )
        return {"success": False, "identifier": identifier, "error": error}
    except Exception as exc:
        logger.warning("get_structure failed", exc_info=True)
        return {"success": False, "error": str(exc)}


# ---------------------------------------------------------------------------
# Adsorbate tools
# ---------------------------------------------------------------------------


def list_adsorbates_impl() -> dict[str, Any]:
    """List built-in adsorbate molecule names."""
    try:
        from mofforge.adsorbate import available_molecules

        names = available_molecules()
        return {"success": True, "count": len(names), "adsorbates": names}
    except Exception as exc:
        logger.warning("list_adsorbates failed", exc_info=True)
        return {"success": False, "error": str(exc)}


def place_adsorbate_impl(
    cif_path: str,
    adsorbate: str = "CO2",
    n_adsorbates: int = 1,
    strategy: str = "void",
    output_cif: str = "with_adsorbate.cif",
    validate: bool = True,
    random_seed: int | None = None,
    **site_kwargs: Any,
) -> dict[str, Any]:
    """Place adsorbate molecule(s) into a MOF and write the result to CIF.

    Extra keyword arguments (e.g. ``grid_spacing``) are forwarded to the
    adsorption-site detector.
    """
    try:
        from mofforge.adsorbate import place_adsorbate

        xtal = load_crystal(cif_path, with_bonds=True)
        result = place_adsorbate(
            xtal,
            adsorbate,
            n_adsorbates=n_adsorbates,
            strategy=strategy,
            validate=validate,
            random_seed=random_seed,
            **site_kwargs,
        )

        out_path = resolve_output(output_cif)
        result.crystal.write_cif(out_path)

        return {
            "success": True,
            "input_cif": cif_path,
            "output_cif": out_path,
            "adsorbate": result.adsorbate_name,
            "n_adsorbates_placed": result.n_adsorbates,
            "sites": [
                {"site_type": s.site_type, "cart_coords": [float(c) for c in s.cart_coords]}
                for s in result.sites
            ],
            "atoms_before": xtal.n_atoms,
            "atoms_after": result.crystal.n_atoms,
            "clashes": result.clashes,
        }
    except Exception as exc:
        logger.warning("place_adsorbate failed", exc_info=True)
        return {"success": False, "error": str(exc)}


# ---------------------------------------------------------------------------
# Validation / construction / rendering
# ---------------------------------------------------------------------------


def validate_impl(
    cif_path: str,
    check_clashes: bool = True,
    check_bonds: bool = True,
    check_coordination: bool = True,
) -> dict[str, Any]:
    """Validate a crystal structure for clashes, bonds, and coordination."""
    try:
        from mofforge.validation import validate_structure

        xtal = load_crystal(cif_path, with_bonds=True)
        report = validate_structure(
            xtal,
            check_clashes=check_clashes,
            check_bonds=check_bonds,
            check_coordination=check_coordination,
        )
        return {
            "success": True,
            "is_valid": report.is_valid,
            "steric_clashes": len(report.steric_clashes),
            "unusual_bonds": len(report.unusual_bonds),
            "coordination_issues": len(report.coordination_issues),
            "warnings": report.warnings,
            "summary": report.summary(),
        }
    except Exception as exc:
        logger.warning("validate failed", exc_info=True)
        return {"success": False, "error": str(exc)}


def build_impl(
    topology: str,
    backend: str = "tobacco",
    node_files: list[str] | None = None,
    edge_files: list[str] | None = None,
    output_dir: str = ".",
) -> dict[str, Any]:
    """Build a MOF from a topology and building blocks."""
    try:
        from mofforge.build import MOFBuilder

        try:
            builder = MOFBuilder(backend=backend)
        except Exception as exc:
            return {
                "success": False,
                "error": f"Failed to initialize {backend} backend: {exc}",
            }

        for n in node_files or []:
            builder.add_node(n)
        for e in edge_files or []:
            builder.add_edge(e)

        out = resolve_output(os.path.join(output_dir, "placeholder.cif"))
        out_dir = str(Path(out).parent)

        result = builder.build(topology=topology, output_dir=out_dir)

        if result.success:
            return {
                "success": True,
                "topology": topology,
                "backend": backend,
                "elapsed_seconds": result.elapsed_seconds,
                "output_paths": [str(p) for p in result.output_paths],
                "atoms": result.crystal.n_atoms if result.crystal else None,
            }
        return {
            "success": False,
            "topology": topology,
            "backend": backend,
            "errors": result.errors,
        }
    except Exception as exc:
        logger.warning("build failed", exc_info=True)
        return {"success": False, "error": str(exc)}


# ---------------------------------------------------------------------------
# Packaged moiety library discovery
# ---------------------------------------------------------------------------


def list_fragments_impl() -> dict[str, Any]:
    """List the packaged moiety (fragment) XYZ files bundled with mofforge."""
    try:
        from mofforge.data import list_moieties, moieties_dir

        names = list_moieties()
        return {
            "success": True,
            "count": len(names),
            "directory": str(moieties_dir()),
            "fragments": names,
        }
    except Exception as exc:
        logger.warning("list_fragments failed", exc_info=True)
        return {"success": False, "error": str(exc)}


def get_fragment_impl(name: str) -> dict[str, Any]:
    """Resolve a packaged moiety name to its absolute XYZ file path."""
    try:
        from mofforge.data import moiety_path

        path = moiety_path(name)
        return {"success": True, "name": name, "path": str(path)}
    except Exception as exc:
        logger.warning("get_fragment failed", exc_info=True)
        return {"success": False, "error": str(exc)}


# ---------------------------------------------------------------------------
# Functionalization tools
# ---------------------------------------------------------------------------


def list_functional_groups_impl() -> dict[str, Any]:
    """List the curated functional groups an agent may graft onto a linker."""
    try:
        from mofforge.functionalize.groups import _GROUPS, available_groups

        return {
            "success": True,
            "count": len(_GROUPS),
            "groups": [
                {"name": g.name, "description": g.description}
                for g in (_GROUPS[name] for name in available_groups())
            ],
        }
    except Exception as exc:
        logger.warning("list_functional_groups failed", exc_info=True)
        return {"success": False, "error": str(exc)}


def find_sites_impl(linker_smiles: str) -> dict[str, Any]:
    """Enumerate functionalizable aromatic C-H sites on a linker SMILES.

    Each site has a stable ``index`` (what the agent selects) and a
    ``symmetry_class``: equal classes are chemically equivalent positions;
    distinct classes are different chemical environments.
    """
    try:
        from mofforge.functionalize.sites import find_functionalizable_sites

        sites = find_functionalizable_sites(linker_smiles)
        n_classes = len({s.symmetry_class for s in sites})
        return {
            "success": True,
            "linker_smiles": linker_smiles,
            "n_sites": len(sites),
            "n_symmetry_classes": n_classes,
            "sites": [
                {
                    "index": s.index,
                    "symmetry_class": s.symmetry_class,
                    "element": s.element,
                    "ring_id": s.ring_id,
                    "description": s.description,
                }
                for s in sites
            ],
        }
    except Exception as exc:
        logger.warning("find_sites failed", exc_info=True)
        return {"success": False, "error": str(exc)}


def functionalize_impl(
    parent_cif: str,
    linker_smiles: str,
    group: str,
    sites: int | list[int] = 0,
    coverage: float = 1.0,
    output_cif: str = "functionalized.cif",
    validate: bool = True,
    random_seed: int | None = None,
) -> dict[str, Any]:
    """Functionalize a MOF linker with a chosen group at chosen site(s)."""
    try:
        from mofforge.functionalize.campaign import functionalize

        out_path = resolve_output(output_cif)
        res = functionalize(
            parent_cif,
            linker_smiles,
            group,
            sites=sites,
            coverage=coverage,
            output_cif=out_path,
            validate=validate,
            random_seed=random_seed,
        )
        if res.error is not None:
            return {"success": False, "error": res.error}
        return {
            "success": True,
            "group": res.group,
            "sites": res.sites,
            "coverage": res.coverage,
            "n_matches": res.n_matches,
            "n_functionalized": res.n_functionalized,
            "output_cif": res.output_cif,
            "is_valid": res.is_valid,
            "clashes": res.clashes,
            "validation_summary": res.validation_summary,
        }
    except Exception as exc:
        logger.warning("functionalize failed", exc_info=True)
        return {"success": False, "error": str(exc)}


def functionalize_campaign_impl(
    parent_cif: str,
    linker_smiles: str,
    groups: list[str],
    coverages: list[float] | None = None,
    sites: int | list[int] = 0,
    output_dir: str = "functionalization_campaign",
    validate: bool = True,
    random_seed: int | None = None,
) -> dict[str, Any]:
    """Sweep groups x coverages, functionalize each, return ranked results."""
    try:
        from mofforge.functionalize.campaign import run_campaign

        out_dir = resolve_output(os.path.join(output_dir, "placeholder"))
        out_dir = str(Path(out_dir).parent)

        results = run_campaign(
            parent_cif,
            linker_smiles,
            groups=groups,
            coverages=coverages,
            sites=sites,
            output_dir=out_dir,
            validate=validate,
            random_seed=random_seed,
        )
        return {
            "success": True,
            "n_results": len(results),
            "results": [
                {
                    "group": r.group,
                    "sites": r.sites,
                    "coverage": r.coverage,
                    "n_matches": r.n_matches,
                    "n_functionalized": r.n_functionalized,
                    "output_cif": r.output_cif,
                    "is_valid": r.is_valid,
                    "clashes": r.clashes,
                    "error": r.error,
                }
                for r in results
            ],
        }
    except Exception as exc:
        logger.warning("functionalize_campaign failed", exc_info=True)
        return {"success": False, "error": str(exc)}


def render_impl(
    input_file: str,
    output_file: str = "structure.png",
    label_mode: str = "sequential",
    representation: str = "ball_stick",
    show_unit_cell: bool = False,
    width: int = 800,
    height: int = 600,
) -> dict[str, Any]:
    """Render a crystal structure file (CIF or XYZ) to a PNG image."""
    try:
        from mofforge.vis.render import render_file_to_png

        out_path = resolve_output(output_file)
        rendered = render_file_to_png(
            input_file=input_file,
            output_file=out_path,
            label_mode=label_mode,
            representation=representation,
            show_unit_cell=show_unit_cell,
            width=width,
            height=height,
        )
        return {
            "success": True,
            "input_file": input_file,
            "output_png": rendered,
        }
    except Exception as exc:
        logger.warning("render failed", exc_info=True)
        return {"success": False, "error": str(exc)}
