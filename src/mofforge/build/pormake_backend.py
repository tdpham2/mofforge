"""Pormake builder backend."""

from __future__ import annotations

import copy
import io
import logging
import shutil
from pathlib import Path
from typing import Any, Literal

from mofforge.build.base import BuildResult, BuildingBlock, Timer, Topology
from mofforge.core.crystal import Crystal

logger = logging.getLogger("mofforge")

_pm = None  # cached pormake module


def _get_pormake():  # noqa: ANN202
    global _pm  # noqa: PLW0603
    if _pm is None:
        try:
            import pormake as pm

            _pm = pm
        except ImportError as exc:
            raise ImportError(
                "pormake is required for the pormake backend.  "
                "Install it with:  pip install pormake"
            ) from exc
    return _pm


class PormakeBackend:
    """Backend that delegates MOF construction to Pormake."""

    name: str = "pormake"

    def __init__(
        self,
        output_dir: Path | str = ".",
        bb_dir: Path | str | None = None,
    ) -> None:
        self._output_dir = Path(output_dir).resolve()
        self._output_dir.mkdir(parents=True, exist_ok=True)

        # Building-block directory (optional on-disk storage)
        self._bb_dir: Path | None = None
        if bb_dir is not None:
            self._bb_dir = Path(bb_dir).resolve()
            self._bb_dir.mkdir(parents=True, exist_ok=True)

        # In-memory registries (name -> BuildingBlock metadata)
        self._registered_nodes: dict[str, BuildingBlock] = {}
        self._registered_edges: dict[str, BuildingBlock] = {}

    def _get_database(self):  # noqa: ANN202
        """Return a pormake.Database instance."""
        pm = _get_pormake()
        return pm.Database()

    def _load_pormake_bb(self, block: BuildingBlock):  # noqa: ANN202
        """Convert a :class:`BuildingBlock` into a ``pormake.BuildingBlock``."""
        pm = _get_pormake()
        src = Path(str(block.source))

        # Try as a file path first
        if src.is_file():
            return pm.BuildingBlock(str(src))

        # Try as a database name
        try:
            db = self._get_database()
            return db.get_bb(block.name)
        except Exception:
            logger.debug("Failed to load '%s' from pormake database", block.name, exc_info=True)

        raise FileNotFoundError(
            f"Building block '{block.name}' not found as file ({src}) or in the pormake database."
        )

    def build(
        self,
        topology: Topology,
        nodes: list[BuildingBlock],
        edges: list[BuildingBlock],
        output_dir: Path,
        **options: Any,
    ) -> BuildResult:
        """Build a MOF using pormake."""
        pm = _get_pormake()
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

        errors: list[str] = []
        timer = Timer()

        with timer:
            # --- Load topology ---
            try:
                db = self._get_database()
                topo = db.get_topo(topology.name)
            except Exception as exc:
                return BuildResult(
                    success=False,
                    errors=[f"Failed to load topology '{topology.name}': {exc}"],
                    backend=self.name,
                )

            # --- Load building blocks ---
            pm_node_bbs: dict[int, Any] = {}
            pm_edge_bbs: dict[tuple[int, ...], Any] = {}

            # Map node building blocks to node types.
            # If a single node is provided, assign it to all node types.
            unique_node_types = list(topo.unique_node_types)
            if len(nodes) == 1:
                try:
                    bb = self._load_pormake_bb(nodes[0])
                    for nt in unique_node_types:
                        pm_node_bbs[nt] = copy.deepcopy(bb)
                except Exception as exc:
                    errors.append(f"Failed to load node '{nodes[0].name}': {exc}")
            else:
                for i, node_block in enumerate(nodes):
                    if i >= len(unique_node_types):
                        errors.append(
                            f"More node building blocks ({len(nodes)}) than "
                            f"node types ({len(unique_node_types)}) in topology "
                            f"'{topology.name}'"
                        )
                        break
                    try:
                        bb = self._load_pormake_bb(node_block)
                        pm_node_bbs[unique_node_types[i]] = bb
                    except Exception as exc:
                        errors.append(f"Failed to load node '{node_block.name}': {exc}")

            # Map edge building blocks to edge types.
            unique_edge_types = [tuple(et) for et in topo.unique_edge_types]
            if len(edges) == 1:
                try:
                    bb = self._load_pormake_bb(edges[0])
                    for et in unique_edge_types:
                        pm_edge_bbs[et] = copy.deepcopy(bb)
                except Exception as exc:
                    errors.append(f"Failed to load edge '{edges[0].name}': {exc}")
            elif len(edges) > 0:
                for i, edge_block in enumerate(edges):
                    if i >= len(unique_edge_types):
                        errors.append(
                            f"More edge building blocks ({len(edges)}) than "
                            f"edge types ({len(unique_edge_types)}) in topology "
                            f"'{topology.name}'"
                        )
                        break
                    try:
                        bb = self._load_pormake_bb(edge_block)
                        pm_edge_bbs[unique_edge_types[i]] = bb
                    except Exception as exc:
                        errors.append(f"Failed to load edge '{edge_block.name}': {exc}")

            if errors:
                return BuildResult(
                    success=False,
                    errors=errors,
                    backend=self.name,
                )

            # --- Build the framework ---
            try:
                builder = pm.Builder()
                framework = builder.build_by_type(
                    topo,
                    node_bbs=pm_node_bbs,
                    edge_bbs=pm_edge_bbs if pm_edge_bbs else None,
                    **options,
                )
            except Exception as exc:
                return BuildResult(
                    success=False,
                    errors=[f"Pormake build failed: {exc}"],
                    backend=self.name,
                )

            # --- Write output CIF ---
            node_names = "_".join(n.name for n in nodes)
            edge_names = "_".join(e.name for e in edges) if edges else "no_edge"
            cif_name = f"{topology.name}_{node_names}_{edge_names}.cif"
            cif_path = output_dir / cif_name

            try:
                framework.write_cif(str(cif_path))
            except Exception as exc:
                return BuildResult(
                    success=False,
                    errors=[f"Failed to write CIF: {exc}"],
                    backend=self.name,
                )

        # --- Load as Crystal ---
        crystal: Crystal | None = None
        output_paths: list[Path] = []
        if cif_path.is_file():
            output_paths.append(cif_path)
            try:
                from mofforge.core.bonding import infer_bonds

                crystal = Crystal.from_cif(str(cif_path))
                crystal = infer_bonds(crystal, periodic=True)
            except Exception as exc:
                logger.warning("Could not load output CIF as Crystal: %s", exc)

        return BuildResult(
            success=True,
            output_paths=output_paths,
            crystal=crystal,
            errors=errors,
            elapsed_seconds=round(timer.elapsed, 2),
            backend=self.name,
            metadata={
                "topology": topology.name,
                "node_bbs": [n.name for n in nodes],
                "edge_bbs": [e.name for e in edges],
                "framework_info": {
                    k: str(v)
                    for k, v in getattr(framework, "info", {}).items()
                    if k in ("relax_obj", "max_rmsd", "mean_rmsd")
                },
            },
        )

    def list_topologies(self) -> list[str]:
        """List available RCSR topology names from the pormake database."""
        try:
            db = self._get_database()
            return sorted(db.topology_list)
        except Exception as exc:
            logger.warning("Failed to list topologies: %s", exc)
            return []

    def describe_topology(self, name: str) -> str:
        """Return the ``topo.describe()`` output for a given topology."""
        try:
            db = self._get_database()
            topo = db.get_topo(name)
            buf = io.StringIO()
            topo.describe(file=buf)
            return buf.getvalue()
        except Exception as exc:
            return f"Failed to describe topology '{name}': {exc}"

    def list_building_blocks(self, role: Literal["node", "edge"]) -> list[str]:
        """List registered building blocks for the given role."""
        registry = self._registered_nodes if role == "node" else self._registered_edges
        names = list(registry.keys())

        # Also include pormake database entries
        try:
            db = self._get_database()
            for bb_name in db.bb_list:
                if bb_name not in names:
                    names.append(bb_name)
        except Exception:
            logger.debug("Failed to list pormake database entries", exc_info=True)

        return sorted(names)

    def add_building_block(self, block: BuildingBlock) -> dict[str, Any]:
        """Register a building block."""
        registry = self._registered_nodes if block.role == "node" else self._registered_edges
        registry[block.name] = block

        # Optionally copy file to bb_dir
        if self._bb_dir is not None:
            src = Path(str(block.source))
            if src.is_file():
                dst = self._bb_dir / src.name
                shutil.copy2(src, dst)
                return {"success": True, "file": str(dst), "registered": block.name}

        return {"success": True, "registered": block.name}

    def remove_building_blocks(
        self,
        role: Literal["node", "edge"],
        names: list[str],
        dry_run: bool = True,
    ) -> dict[str, Any]:
        """Remove building blocks from the in-memory registry."""
        registry = self._registered_nodes if role == "node" else self._registered_edges
        existing = [n for n in names if n in registry]
        not_found = [n for n in names if n not in registry]

        if dry_run:
            return {
                "success": True,
                "dry_run": True,
                "would_remove": existing,
                "not_found": not_found if not_found else None,
                "count": len(existing),
            }

        for n in existing:
            del registry[n]

        return {
            "success": True,
            "removed": existing,
            "count": len(existing),
            "dry_run": False,
        }

    def clear_building_blocks(
        self,
        role: Literal["node", "edge"],
        dry_run: bool = True,
    ) -> dict[str, Any]:
        """Clear all registered building blocks for the given role."""
        registry = self._registered_nodes if role == "node" else self._registered_edges
        names = list(registry.keys())

        if dry_run:
            return {
                "success": True,
                "dry_run": True,
                "would_remove": names,
                "count": len(names),
            }

        registry.clear()
        return {
            "success": True,
            "removed": names,
            "count": len(names),
            "dry_run": False,
        }

    def copy_from_database(
        self,
        role: Literal["node", "edge"],
        names: list[str] | None = None,
        source: Path | None = None,
        dry_run: bool = True,
    ) -> dict[str, Any]:
        """List or copy building blocks from the pormake database."""
        try:
            db = self._get_database()
            available = sorted(db.bb_list)
        except Exception as exc:
            return {"success": False, "error": str(exc)}

        if names is None:
            return {
                "success": True,
                "available_in_database": available,
                "count": len(available),
            }

        available_set = set(available)
        to_copy = [n for n in names if n in available_set]
        not_found = [n for n in names if n not in available_set]

        if dry_run:
            return {
                "success": True,
                "dry_run": True,
                "would_register": to_copy,
                "not_found": not_found if not_found else None,
                "count": len(to_copy),
            }

        registered: list[str] = []
        errors: list[dict[str, str]] = []
        registry = self._registered_nodes if role == "node" else self._registered_edges

        for bb_name in to_copy:
            try:
                # Load from database to verify it works
                db.get_bb(bb_name)
                registry[bb_name] = BuildingBlock(
                    name=bb_name,
                    role=role,
                    source=bb_name,  # Will be resolved via database at build time
                )
                registered.append(bb_name)
            except Exception as exc:
                errors.append({"name": bb_name, "error": str(exc)})

        return {
            "success": len(errors) == 0,
            "registered": registered,
            "count": len(registered),
            "errors": errors if errors else None,
            "dry_run": False,
        }

    def get_configuration(self) -> dict[str, Any]:
        """Return current pormake backend configuration."""
        return {
            "success": True,
            "configuration": {
                "output_dir": str(self._output_dir),
                "bb_dir": str(self._bb_dir) if self._bb_dir else None,
                "registered_nodes": len(self._registered_nodes),
                "registered_edges": len(self._registered_edges),
            },
        }

    def set_configuration(self, key: str, value: Any) -> dict[str, Any]:
        """Set a configuration value."""
        if key == "output_dir":
            self._output_dir = Path(str(value)).resolve()
            self._output_dir.mkdir(parents=True, exist_ok=True)
            return {"success": True, "message": f"Set output_dir = {self._output_dir}"}
        if key == "bb_dir":
            self._bb_dir = Path(str(value)).resolve()
            self._bb_dir.mkdir(parents=True, exist_ok=True)
            return {"success": True, "message": f"Set bb_dir = {self._bb_dir}"}
        return {"success": False, "error": f"Unknown configuration key: {key}"}

    def status(self) -> dict[str, Any]:
        """Return overall status of the pormake backend."""
        try:
            db = self._get_database()
            n_topos = len(db.topology_list)
            n_bbs = len(db.bb_list)
            db_available = True
        except Exception:
            logger.debug("Pormake database unavailable", exc_info=True)
            n_topos = 0
            n_bbs = 0
            db_available = False

        return {
            "success": True,
            "backend": self.name,
            "database_available": db_available,
            "topologies_available": n_topos,
            "database_building_blocks": n_bbs,
            "registered_nodes": len(self._registered_nodes),
            "registered_edges": len(self._registered_edges),
            "output_dir": str(self._output_dir),
            "ready": db_available,
        }

    def visualize_molecule(
        self,
        smiles: str,
        labelsize: int = 10,
        png_filename: str = "molecule.png",
        width: int = 800,
        height: int = 800,
        device_scale_factor: int = 2,
    ) -> str:
        """Visualize a molecule with py3Dmol and export a PNG."""
        import os
        import tempfile

        import py3Dmol
        from architector.io_molecule import convert_io_molecule
        from architector import io_ptable
        from playwright.sync_api import sync_playwright

        mol = convert_io_molecule(smiles)
        metal_ind = [i for i, x in enumerate(mol.ase_atoms) if x.symbol in io_ptable.all_metals]

        with tempfile.TemporaryDirectory() as tmpdir:
            mol2_path = os.path.join(tmpdir, "molecule.mol2")
            coords = mol.write_mol2(mol2_path, writestring=True)

            view_ats = py3Dmol.view(width=width, height=height)
            view_ats.addModel(coords.replace("un", "1"), "mol2")
            view_ats.addStyle({"sphere": {"colorscheme": "Jmol", "scale": 0.3}})

            msyms = [mol.ase_atoms.get_chemical_symbols()[x] for x in metal_ind]
            for ms in set(msyms):
                view_ats.setStyle(
                    {"elem": ms},
                    {"sphere": {"colorscheme": "Jmol", "scale": 0.75}},
                )

            view_ats.addStyle({"stick": {"colorscheme": "Jmol", "radius": 0.25}})

            inds = list(range(len(mol.ase_atoms)))
            for p, i in enumerate(inds):
                atom_posit = mol.ase_atoms.positions[p]
                view_ats.addLabel(
                    str(i),
                    {
                        "position": {
                            "x": float(atom_posit[0]),
                            "y": float(atom_posit[1]),
                            "z": float(atom_posit[2]),
                        },
                        "backgroundColor": "black",
                        "backgroundOpacity": 0.4,
                        "fontOpacity": 1.0,
                        "fontSize": int(labelsize),
                        "fontColor": "white",
                        "inFront": True,
                    },
                )

            view_ats.zoomTo()
            view_ats.render()
            html_str = view_ats._make_html()
            html_path = os.path.join(tmpdir, "view.html")
            with open(html_path, "w") as f:
                f.write(html_str)

            with sync_playwright() as pw:
                browser = pw.chromium.launch(headless=True)
                context = browser.new_context(
                    viewport={"width": width, "height": height},
                    device_scale_factor=device_scale_factor,
                )
                page = context.new_page()
                page.goto("file://" + html_path)
                page.wait_for_timeout(1000)
                canvas = page.locator("canvas")
                canvas.screenshot(path=png_filename)
                browser.close()

        return f"The file is saved under {os.path.abspath(png_filename)}"

    def image_to_connection_points(
        self,
        image_path: str,
        smiles: str,
        model_name: str = "argo:gpt-4o-mini",
        base_url: str | None = None,
    ) -> dict[str, Any]:
        """Use a multimodal LLM to infer connection points from an image."""
        import base64
        import json

        from chemgraph.tools.argoapi_loader import load_argoapi_model

        with open(image_path, "rb") as f:
            image_bytes = f.read()
        image_b64 = base64.b64encode(image_bytes).decode("utf-8")

        vision_llm = load_argoapi_model(
            model_name=model_name, temperature=0, base_url=base_url, api_key=None
        )

        messages = [
            {
                "role": "system",
                "content": (
                    "You are a chemistry vision model specialized in analyzing "
                    "ligand structures. Given a 2D molecule diagram, identify "
                    "atom indices that can coordinate to a metal. Return ONLY "
                    "a JSON dictionary with the fields: "
                    "{smiles: str, connection_points: [int], reasoning: str}. "
                    "If the molecule only contains 1 atom, the output is [0]"
                ),
            },
            {
                "role": "user",
                "content": [
                    {
                        "type": "text",
                        "text": (
                            "Here is a ligand image. Identify all atoms that "
                            "are likely to serve as metal coordination points. "
                            "Use the index labels visible in the diagram."
                        ),
                    },
                    {
                        "type": "image_url",
                        "image_url": {"url": f"data:image/png;base64,{image_b64}"},
                    },
                ],
            },
        ]

        response = vision_llm.invoke(messages)
        raw_output = response.content

        try:
            parsed = json.loads(raw_output)
        except Exception:
            try:
                start = raw_output.index("{")
                end = raw_output.rindex("}") + 1
                parsed = json.loads(raw_output[start:end])
            except Exception:
                parsed = {
                    "smiles": smiles,
                    "connection_points": [],
                    "reasoning": f"Failed to parse model output: {raw_output}",
                }

        return parsed
