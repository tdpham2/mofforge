"""ChemGraph-native MCP server for mofforge, built on ``CGFastMCP``.

This is the *HPC* entry point. Unlike the stock server
(:mod:`mofforge.mcp.server`), which runs every tool inline, this server uses
ChemGraph's :class:`CGFastMCP` (from the ``dev-globus-hpc`` branch) so that
heavy and fan-out workloads are submitted to an execution backend (e.g. Parsl
on Polaris/Aurora) with async job tracking.

Tool placement:

* **Inline** (run in-process): database search/screen/lookup, adsorbate
  listing, and validation — fast and IO-bound.
* **Backend** (submitted as tasks): MOF construction (TOBACCO/pormake) and PNG
  rendering — CPU/GPU-heavy.
* **Ensemble fan-out**: ``mofforge_screen_and_place`` screens the CoRE MOF
  database and, for every hit, resolves a structure and places an adsorbate —
  one backend task per candidate.

``chemgraph`` is imported lazily inside :func:`build_server` so importing this
module does not hard-require ChemGraph; only running the server does.

Tool logic lives in :mod:`mofforge.mcp._impl` (plain, top-level functions) so it
is shared with the stock server and pickles by reference for remote workers.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Callable, Collection
from pathlib import Path
from typing import Any

from mofforge.mcp import _impl
from mofforge.mcp.tool_selection import parse_tool_list, select_tool_names

logger = logging.getLogger("mofforge.mcp.chemgraph")

# Resource hints for backend tasks. Conservative defaults; override per site.
_BUILD_TASK = {"num_nodes": 1, "processes_per_node": 1}
_RENDER_TASK = {"num_nodes": 1, "processes_per_node": 1}

_TOOL_REQUIREMENTS: dict[str, str | None] = {
    "mofforge_search_coremof": None,
    "mofforge_screen_coremof": None,
    "mofforge_search_csd": None,
    "mofforge_lookup_mof": None,
    "mofforge_get_structure": None,
    "mofforge_list_adsorbates": None,
    "mofforge_validate": None,
    "mofforge_list_fragments": None,
    "mofforge_get_fragment": None,
    "mofforge_list_functional_groups": None,
    "mofforge_find_sites": "chem",
    "mofforge_functionalize": "chem",
    "mofforge_functionalize_campaign": "chem",
    "mofforge_build": "build",
    "mofforge_render": "vis",
    "mofforge_screen_and_place": None,
    "check_job_status": None,
    "get_job_results": None,
    "list_jobs": None,
    "cancel_job": None,
    "check_endpoint_status": None,
}


# ---------------------------------------------------------------------------
# Fan-out: screen the CoRE MOF DB, then place an adsorbate in each hit.
# These are module-level so CGFastMCP can pickle them by reference.
# ---------------------------------------------------------------------------


def _screen_place_worker(item: dict[str, Any]) -> dict[str, Any]:
    """Backend worker: resolve a structure and place an adsorbate for one hit.

    Runs on a compute node, so the configured CoRE MOF structures directory
    must be reachable there (or staged by a pre-submit hook).
    """
    coreid = item.get("coreid", "")
    structures_dir = item.get("structures_dir")
    cif_path = item.get("remote_structure_file") or item.get("cif_path")

    if not cif_path:
        resolved = _impl.get_structure_impl(coreid, structures_dir=structures_dir)
        if not resolved.get("success"):
            return {"coreid": coreid, "status": "no_structure", **resolved}
        cif_path = resolved["cif_path"]

    out_name = item.get("output_cif") or f"{coreid}_{item.get('adsorbate', 'CO2')}.cif"
    placement = _impl.place_adsorbate_impl(
        cif_path,
        adsorbate=item.get("adsorbate", "CO2"),
        n_adsorbates=int(item.get("n_adsorbates", 1)),
        strategy=item.get("strategy", "void"),
        output_cif=out_name,
        validate=bool(item.get("validate", True)),
        random_seed=item.get("random_seed"),
    )
    status = "success" if placement.get("success") else "failed"
    return {"coreid": coreid, "status": status, **placement}


def _server_instructions(enabled_names: Collection[str]) -> str:
    """Build ChemGraph instructions that match the selected catalog."""
    enabled = set(enabled_names)
    descriptions: list[str] = []
    if enabled & {
        "mofforge_search_coremof",
        "mofforge_screen_coremof",
        "mofforge_search_csd",
        "mofforge_lookup_mof",
        "mofforge_get_structure",
    }:
        descriptions.append("database discovery")
    if "mofforge_validate" in enabled:
        descriptions.append("structure validation")
    if "mofforge_build" in enabled:
        descriptions.append("backend MOF construction")
    if "mofforge_render" in enabled:
        descriptions.append("backend rendering")
    if "mofforge_functionalize" in enabled:
        descriptions.append("backend linker functionalization")
    if "mofforge_screen_and_place" in enabled:
        descriptions.append("screen-and-place ensemble fan-out")
    summary = ", ".join(descriptions) if descriptions else "the configured tool set"
    return (
        "mofforge MOF tools with HPC backend execution. Exposed capabilities: "
        + summary
        + ". Lightweight database operations run inline; selected heavy and "
        "ensemble tools are submitted to the configured execution backend."
    )


def build_server(
    enabled_tools: Collection[str] | None = None,
    available_only: bool = False,
):
    """Construct and return the ``CGFastMCP`` server instance.

    Raises
    ------
    ImportError
        If ChemGraph (with ``CGFastMCP``) is not installed.
    """
    try:
        from chemgraph.mcp.cg_fastmcp import CGFastMCP
    except ImportError as exc:  # pragma: no cover - depends on optional dep
        raise ImportError(
            "mofforge.mcp.chemgraph_server requires ChemGraph with CGFastMCP "
            "(the 'dev-globus-hpc' branch). Install ChemGraph, or use the stock "
            "server 'mofforge-mcp' instead."
        ) from exc

    selected = select_tool_names(
        _TOOL_REQUIREMENTS,
        _TOOL_REQUIREMENTS,
        requested=enabled_tools,
        available_only=available_only,
    )

    class SelectedCGFastMCP(CGFastMCP):
        """CGFastMCP variant that registers only this process's selected names."""

        def add_tool(self, fn, name=None, **kwargs):
            if (name or fn.__name__) not in selected:
                return None
            return super().add_tool(fn, name=name, **kwargs)

    mcp = SelectedCGFastMCP(
        name="mofforge HPC Tools",
        instructions=_server_instructions(selected),
    )

    def selected_tool(**kwargs: Any) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
        """Register a normal CGFastMCP tool only when selected."""
        if kwargs["name"] in selected:
            return mcp.tool(**kwargs)
        return lambda function: function

    def selected_fanout_tool(
        **kwargs: Any,
    ) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
        """Register a CGFastMCP fan-out tool only when selected."""
        if kwargs["name"] in selected:
            return mcp.schema_fanout_tool(**kwargs)
        return lambda function: function

    # ---- Inline (lightweight) tools ------------------------------------- #
    @selected_tool(
        name="mofforge_search_coremof",
        description="Search the CoRE MOF database (coreid/refcode/name/metal/topology).",
    )
    def mofforge_search_coremof(
        query: str, field: str = "auto", limit: int = 25, data_path: str | None = None
    ) -> dict:
        return _impl.search_coremof_impl(
            query, field=field, limit=limit, data_path=data_path
        )

    @selected_tool(
        name="mofforge_screen_coremof",
        description="Screen CoRE MOFs by pore size, density, stability, metal, topology, OMS.",
    )
    def mofforge_screen_coremof(
        lcd_min: float | None = None,
        lcd_max: float | None = None,
        pld_min: float | None = None,
        pld_max: float | None = None,
        density_min: float | None = None,
        density_max: float | None = None,
        void_fraction_min: float | None = None,
        water_stability_min: float | None = None,
        thermal_stability_min: float | None = None,
        metal: str | None = None,
        topology: str | None = None,
        has_oms: bool | None = None,
        limit: int = 50,
        data_path: str | None = None,
    ) -> dict:
        return _impl.screen_coremof_impl(
            lcd_min=lcd_min,
            lcd_max=lcd_max,
            pld_min=pld_min,
            pld_max=pld_max,
            density_min=density_min,
            density_max=density_max,
            void_fraction_min=void_fraction_min,
            water_stability_min=water_stability_min,
            thermal_stability_min=thermal_stability_min,
            metal=metal,
            topology=topology,
            has_oms=has_oms,
            limit=limit,
            data_path=data_path,
        )

    @selected_tool(
        name="mofforge_search_csd",
        description="Search the CSD lookup table (refcode/name/doi/formula/ccdc).",
    )
    def mofforge_search_csd(
        query: str, field: str = "auto", limit: int = 25, data_path: str | None = None
    ) -> dict:
        return _impl.search_csd_impl(query, field=field, limit=limit, data_path=data_path)

    @selected_tool(
        name="mofforge_lookup_mof",
        description="Bridge a MOF name from CSD to simulation-ready CoRE MOF entries.",
    )
    def mofforge_lookup_mof(
        name: str,
        limit: int = 25,
        csd_data_path: str | None = None,
        coremof_data_path: str | None = None,
    ) -> dict:
        return _impl.lookup_mof_impl(
            name,
            limit=limit,
            csd_data_path=csd_data_path,
            coremof_data_path=coremof_data_path,
        )

    @selected_tool(
        name="mofforge_get_structure",
        description="Resolve a CoRE MOF coreid/refcode to a local CIF file path.",
    )
    def mofforge_get_structure(identifier: str, structures_dir: str | None = None) -> dict:
        return _impl.get_structure_impl(identifier, structures_dir=structures_dir)

    @selected_tool(
        name="mofforge_list_adsorbates",
        description="List the built-in adsorbate molecules available for placement.",
    )
    def mofforge_list_adsorbates() -> dict:
        return _impl.list_adsorbates_impl()

    @selected_tool(
        name="mofforge_validate",
        description="Validate a crystal structure for clashes, bonds, and coordination.",
    )
    def mofforge_validate(
        cif_path: str,
        check_clashes: bool = True,
        check_bonds: bool = True,
        check_coordination: bool = True,
    ) -> dict:
        return _impl.validate_impl(
            cif_path,
            check_clashes=check_clashes,
            check_bonds=check_bonds,
            check_coordination=check_coordination,
        )

    @selected_tool(
        name="mofforge_list_fragments",
        description="List packaged moiety fragment XYZ files bundled with mofforge.",
    )
    def mofforge_list_fragments() -> dict:
        return _impl.list_fragments_impl()

    @selected_tool(
        name="mofforge_get_fragment",
        description="Resolve a packaged moiety fragment name to an absolute XYZ path.",
    )
    def mofforge_get_fragment(name: str) -> dict:
        return _impl.get_fragment_impl(name)

    @selected_tool(
        name="mofforge_list_functional_groups",
        description="List curated functional groups for linker functionalization.",
    )
    def mofforge_list_functional_groups() -> dict:
        return _impl.list_functional_groups_impl()

    @selected_tool(
        name="mofforge_find_sites",
        description=(
            "Enumerate functionalizable aromatic C-H sites on a linker SMILES; "
            "each has an 'index' to select and a 'symmetry_class'."
        ),
    )
    def mofforge_find_sites(linker_smiles: str) -> dict:
        return _impl.find_sites_impl(linker_smiles)

    # ---- Backend (heavy) tools ----------------------------------------- #
    @selected_tool(
        name="mofforge_functionalize",
        description=(
            "Functionalize a MOF linker with a chosen group at chosen site(s); "
            "geometry generated automatically, coverage sets concentration."
        ),
        **_BUILD_TASK,
    )
    def mofforge_functionalize(
        parent_cif: str,
        linker_smiles: str,
        group: str,
        sites: list[int] | None = None,
        coverage: float = 1.0,
        output_cif: str = "functionalized.cif",
        validate: bool = True,
        random_seed: int | None = None,
    ) -> dict:
        return _impl.functionalize_impl(
            parent_cif,
            linker_smiles,
            group,
            sites=sites if sites is not None else 0,
            coverage=coverage,
            output_cif=output_cif,
            validate=validate,
            random_seed=random_seed,
        )

    @selected_tool(
        name="mofforge_functionalize_campaign",
        description=(
            "Sweep functional groups x coverages on a linker, validate each, "
            "and return results ranked best-first (valid, then fewest clashes)."
        ),
        **_BUILD_TASK,
    )
    def mofforge_functionalize_campaign(
        parent_cif: str,
        linker_smiles: str,
        groups: list[str],
        coverages: list[float] | None = None,
        sites: list[int] | None = None,
        output_dir: str = "functionalization_campaign",
        validate: bool = True,
        random_seed: int | None = None,
    ) -> dict:
        return _impl.functionalize_campaign_impl(
            parent_cif,
            linker_smiles,
            groups,
            coverages=coverages,
            sites=sites if sites is not None else 0,
            output_dir=output_dir,
            validate=validate,
            random_seed=random_seed,
        )

    # ---- Backend (heavy) tools ----------------------------------------- #
    @selected_tool(
        name="mofforge_build",
        description="Build a MOF from a topology and building blocks (TOBACCO/pormake).",
        **_BUILD_TASK,
    )
    def mofforge_build(
        topology: str,
        backend: str = "tobacco",
        node_files: list[str] | None = None,
        edge_files: list[str] | None = None,
        output_dir: str = ".",
    ) -> dict:
        return _impl.build_impl(
            topology,
            backend=backend,
            node_files=node_files,
            edge_files=edge_files,
            output_dir=output_dir,
        )

    @selected_tool(
        name="mofforge_render",
        description="Render a crystal structure file (CIF/XYZ) to a PNG image.",
        **_RENDER_TASK,
    )
    def mofforge_render(
        input_file: str,
        output_file: str = "structure.png",
        label_mode: str = "sequential",
        representation: str = "ball_stick",
        show_unit_cell: bool = False,
        width: int = 800,
        height: int = 600,
    ) -> dict:
        return _impl.render_impl(
            input_file,
            output_file=output_file,
            label_mode=label_mode,
            representation=representation,
            show_unit_cell=show_unit_cell,
            width=width,
            height=height,
        )

    # ---- Ensemble fan-out: screen -> (per hit) get_structure + place ---- #
    @selected_fanout_tool(
        name="mofforge_screen_and_place",
        description=(
            "Screen the CoRE MOF database by property ranges, then for each "
            "matching MOF resolve its structure and place an adsorbate. One "
            "backend task per candidate; returns per-MOF placement results "
            "ready for gRASPA/ASE."
        ),
        worker=_screen_place_worker,
        **_BUILD_TASK,
    )
    def mofforge_screen_and_place(params: dict) -> list[dict]:
        """Expander: turn screening criteria into one job per matching MOF.

        ``params`` is the ensemble schema (a dict). Recognized keys: the same
        screening filters as ``mofforge_screen_coremof`` plus ``adsorbate``,
        ``n_adsorbates``, ``strategy``, ``structures_dir``, ``output_dir``,
        ``data_path``, and ``limit``.
        """
        p = dict(params or {})
        adsorbate = p.pop("adsorbate", "CO2")
        n_adsorbates = p.pop("n_adsorbates", 1)
        strategy = p.pop("strategy", "void")
        structures_dir = p.pop("structures_dir", None)
        output_dir = p.pop("output_dir", ".")
        data_path = p.pop("data_path", None)
        limit = p.pop("limit", 50)

        screen = _impl.screen_coremof_impl(limit=limit, data_path=data_path, **p)
        if not screen.get("success"):
            return []

        jobs: list[dict] = []
        for rec in screen.get("records", []):
            coreid = rec.get("coreid", "")
            jobs.append(
                {
                    "coreid": coreid,
                    "structures_dir": structures_dir,
                    "adsorbate": adsorbate,
                    "n_adsorbates": n_adsorbates,
                    "strategy": strategy,
                    "output_cif": os.path.join(
                        output_dir, f"{coreid}_{adsorbate}.cif"
                    ),
                }
            )
        return jobs

    return mcp


def main() -> None:
    """Launch the mofforge CGFastMCP server with backend lifecycle management."""
    from chemgraph.mcp.server_utils import run_mcp_server

    mcp = build_server(
        enabled_tools=parse_tool_list(os.environ.get("MOFFORGE_MCP_TOOLS")),
        available_only=os.environ.get("MOFFORGE_MCP_AVAILABLE_ONLY", "").lower()
        in {"1", "true", "yes", "on"},
    )

    jobs_file = os.environ.get(
        "MOFFORGE_MCP_JOBS_FILE", str(Path.home() / ".mofforge_mcp_jobs.json")
    )
    mcp.init_backend(tracker_kwargs={"persist_file": jobs_file})
    try:
        run_mcp_server(mcp, default_port=9011)
    finally:
        mcp.shutdown_backend()


if __name__ == "__main__":
    main()
