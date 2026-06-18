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
from pathlib import Path
from typing import Any

from mofforge.mcp import _impl

logger = logging.getLogger("mofforge.mcp.chemgraph")

# Resource hints for backend tasks. Conservative defaults; override per site.
_BUILD_TASK = {"num_nodes": 1, "processes_per_node": 1}
_RENDER_TASK = {"num_nodes": 1, "processes_per_node": 1}


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


def build_server():
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

    mcp = CGFastMCP(
        name="mofforge HPC Tools",
        instructions=(
            "mofforge MOF tools with HPC backend execution. Lightweight "
            "database search/screening, structure lookup, and validation run "
            "inline; MOF construction, rendering, and the screen+place "
            "ensemble are submitted to the execution backend with job "
            "tracking. Typical flow: screen the CoRE MOF database by pore size "
            "/ stability, then fan out structure retrieval + adsorbate "
            "placement across all candidates for downstream gRASPA/ASE runs."
        ),
    )

    # ---- Inline (lightweight) tools ------------------------------------- #
    @mcp.tool(
        name="mofforge_search_coremof",
        description="Search the CoRE MOF database (coreid/refcode/name/metal/topology).",
    )
    def mofforge_search_coremof(
        query: str, field: str = "auto", limit: int = 25, data_path: str | None = None
    ) -> dict:
        return _impl.search_coremof_impl(
            query, field=field, limit=limit, data_path=data_path
        )

    @mcp.tool(
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

    @mcp.tool(
        name="mofforge_search_csd",
        description="Search the CSD lookup table (refcode/name/doi/formula/ccdc).",
    )
    def mofforge_search_csd(
        query: str, field: str = "auto", limit: int = 25, data_path: str | None = None
    ) -> dict:
        return _impl.search_csd_impl(query, field=field, limit=limit, data_path=data_path)

    @mcp.tool(
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

    @mcp.tool(
        name="mofforge_get_structure",
        description="Resolve a CoRE MOF coreid/refcode to a local CIF file path.",
    )
    def mofforge_get_structure(identifier: str, structures_dir: str | None = None) -> dict:
        return _impl.get_structure_impl(identifier, structures_dir=structures_dir)

    @mcp.tool(
        name="mofforge_list_adsorbates",
        description="List the built-in adsorbate molecules available for placement.",
    )
    def mofforge_list_adsorbates() -> dict:
        return _impl.list_adsorbates_impl()

    @mcp.tool(
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

    # ---- Backend (heavy) tools ----------------------------------------- #
    @mcp.tool(
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

    @mcp.tool(
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
    @mcp.schema_fanout_tool(
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

    mcp = build_server()

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
