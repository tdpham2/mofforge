"""TOBACCO 3.0 builder backend (importable ``tobacco3`` API)."""

from __future__ import annotations

import logging
from dataclasses import fields
from pathlib import Path
from typing import Any, Literal

from mofforge.build.base import BuildingBlock, BuildResult, Timer, Topology
from mofforge.build.config import ConfigError, validate_tobacco, validate_tobacco_data_dir
from mofforge.core.crystal import Crystal

logger = logging.getLogger("mofforge")


class TobaccoBackend:
    """Backend that delegates MOF construction to the ``tobacco3`` package.

    Uses the filesystem-independent ``tobacco3.generate_mof`` API: building
    blocks and templates are fed in as in-memory objects and produced CIFs are
    returned as :class:`tobacco3.MofResult` objects, so no working-directory
    juggling or output scraping is required.  Topology and building-block CIFs
    are discovered from the TOBACCO *data directory* (``template_database/``,
    ``nodes_database/``, ``edges_database/``, and the active ``templates/``,
    ``nodes/``, ``edges/`` folders).
    """

    name: str = "tobacco"

    # Map a logical role to (active dir, database dir).  The active dirs are the
    # small curated defaults shipped in the tobacco repo; the *_database dirs are
    # the full catalog.  Both are searched when resolving a block/topology name.
    _ROLE_DIRS: dict[str, tuple[str, str]] = {
        "node": ("nodes", "nodes_database"),
        "edge": ("edges", "edges_database"),
        "template": ("templates", "template_database"),
    }

    def __init__(self, data_dir: Path | str) -> None:
        errors = validate_tobacco()
        if errors:
            raise ConfigError("; ".join(errors))

        self._data_dir = Path(data_dir).resolve()
        data_errors = validate_tobacco_data_dir(self._data_dir)
        if data_errors:
            details = "\n  ".join(data_errors)
            raise ConfigError(
                f"Invalid TOBACCO data directory at {self._data_dir}:\n  {details}"
            )

        import tobacco3

        self._tobacco3 = tobacco3
        self._cfg = tobacco3.TobaccoConfig()

    # ------------------------------------------------------------------ #
    # Data-directory helpers
    # ------------------------------------------------------------------ #
    def _search_dirs(self, role: Literal["node", "edge", "template"]) -> list[Path]:
        """Return existing (database, active) directories to search for *role*."""
        active, database = self._ROLE_DIRS[role]
        dirs = [self._data_dir / database, self._data_dir / active]
        return [d for d in dirs if d.is_dir()]

    def _list_cifs(self, role: Literal["node", "edge", "template"]) -> list[str]:
        """Return the sorted union of ``.cif`` names available for *role*."""
        names: set[str] = set()
        for d in self._search_dirs(role):
            names.update(f.name for f in d.iterdir() if f.suffix == ".cif")
        return sorted(names)

    def _resolve_cif(
        self, role: Literal["node", "edge", "template"], name: str
    ) -> Path | None:
        """Resolve a bare block/topology *name* to a concrete CIF path.

        Accepts either an existing filesystem path or a name looked up in the
        role's database / active directories (``.cif`` suffix optional).
        """
        p = Path(name)
        if p.is_file():
            return p.resolve()

        cif_name = name if name.endswith(".cif") else f"{name}.cif"
        for d in self._search_dirs(role):
            candidate = d / cif_name
            if candidate.is_file():
                return candidate.resolve()
        return None

    def _to_tobacco_bb(self, block: BuildingBlock) -> Any:
        """Convert a mofforge :class:`BuildingBlock` to a ``tobacco3.BuildingBlock``.

        The block's *name* preserves the CIF basename because the TOBACCO
        pipeline keys some behavior off it (e.g. the ``ntn_edge.cif`` special
        case), so identity must survive the round-trip.
        """
        path = self._resolve_cif(block.role, str(block.source))
        if path is None:
            raise FileNotFoundError(
                f"Could not resolve {block.role} '{block.source}' to a CIF file "
                f"(searched {[str(d) for d in self._search_dirs(block.role)]})"
            )
        name = path.name if path.name.endswith(".cif") else f"{path.name}.cif"
        return self._tobacco3.BuildingBlock.from_cif(str(path), name=name)

    # ------------------------------------------------------------------ #
    # Build
    # ------------------------------------------------------------------ #
    def build(
        self,
        topology: Topology,
        nodes: list[BuildingBlock],
        edges: list[BuildingBlock],
        output_dir: Path,
        *,
        verbose: bool = False,
        **options: Any,
    ) -> BuildResult:
        """Run ``tobacco3.generate_mof`` to construct MOF structure(s)."""
        output_dir = Path(output_dir).resolve()

        if not nodes:
            return BuildResult(
                success=False,
                errors=["No node building blocks provided."],
                backend=self.name,
            )
        if not edges:
            return BuildResult(
                success=False,
                errors=["No edge building blocks provided."],
                backend=self.name,
            )

        # Resolve template.  Prefer an explicit source path, else look up the
        # topology name in template_database / templates.
        template_ref = topology.source or topology.name
        template_path = self._resolve_cif("template", str(template_ref))
        if template_path is None:
            return BuildResult(
                success=False,
                errors=[
                    f"Topology '{topology.name}' not found in "
                    f"{[str(d) for d in self._search_dirs('template')]}"
                ],
                backend=self.name,
            )

        try:
            template_obj = self._tobacco3.Template.from_cif(str(template_path))
            node_objs = [self._to_tobacco_bb(n) for n in nodes]
            edge_objs = [self._to_tobacco_bb(e) for e in edges]
        except Exception as exc:
            return BuildResult(
                success=False,
                errors=[f"Failed to load building blocks: {exc}"],
                backend=self.name,
            )

        results: list[Any] = []
        errors: list[str] = []
        timer = Timer()
        with timer:
            try:
                results = self._tobacco3.generate_mof(
                    template_obj,
                    node_objs,
                    edge_objs,
                    config=self._cfg,
                    quiet=not verbose,
                )
            except Exception as exc:
                logger.warning("tobacco3.generate_mof failed", exc_info=True)
                return BuildResult(
                    success=False,
                    errors=[f"generate_mof failed: {exc}"],
                    elapsed_seconds=round(timer.elapsed, 2),
                    backend=self.name,
                )

        if not results:
            return BuildResult(
                success=False,
                errors=[
                    "generate_mof produced no structures "
                    "(no valid vertex/edge assignment for this topology + building blocks)."
                ],
                elapsed_seconds=round(timer.elapsed, 2),
                backend=self.name,
            )

        # Persist every produced structure.
        output_dir.mkdir(parents=True, exist_ok=True)
        copied: list[Path] = []
        for r in results:
            try:
                copied.append(Path(r.write(str(output_dir))))
            except Exception as exc:
                errors.append(f"Failed to write {r.cifname}: {exc}")

        # Load the first structure as a Crystal for downstream inspection.
        crystal: Crystal | None = None
        if copied:
            try:
                from mofforge.core.bonding import infer_bonds

                crystal = Crystal.from_cif(str(copied[0]))
                crystal = infer_bonds(crystal, periodic=True)
            except Exception as exc:
                logger.warning("Could not load output CIF as Crystal: %s", exc)

        return BuildResult(
            success=bool(copied),
            output_paths=copied,
            crystal=crystal,
            errors=errors,
            elapsed_seconds=round(timer.elapsed, 2),
            backend=self.name,
            metadata={
                "topology": topology.name,
                "n_structures": len(results),
                "structures": [
                    {
                        "cifname": r.cifname,
                        "n_atoms": r.n_atoms,
                        "net_charge": r.net_charge,
                        "bond_check_passed": r.bond_check_passed,
                    }
                    for r in results
                ],
            },
        )

    # ------------------------------------------------------------------ #
    # Topology / building-block discovery
    # ------------------------------------------------------------------ #
    def list_topologies(self) -> list[str]:
        """List available topology template CIF names."""
        return self._list_cifs("template")

    def describe_topology(self, name: str) -> str:
        """Return basic information about a topology template."""
        path = self._resolve_cif("template", name)
        if path is None:
            return f"Template '{name}' not found."
        return f"Template: {path.name} ({path.stat().st_size} bytes) [{path.parent.name}]"

    def list_building_blocks(self, role: Literal["node", "edge"]) -> list[str]:
        """List available node or edge building-block CIF names."""
        return self._list_cifs(role)

    def add_building_block(self, block: BuildingBlock) -> dict[str, Any]:
        """Validate that *block* resolves to a usable CIF.

        Registration is tracked in-memory by :class:`~mofforge.build.builder.MOFBuilder`;
        the block is passed directly to ``generate_mof`` at build time, so nothing
        is copied onto disk here.
        """
        src = str(block.source)
        if src.endswith((".xyz", ".mol2")):
            return {"success": False, "error": "TOBACCO building blocks must be CIF files"}

        path = self._resolve_cif(block.role, src)
        if path is None:
            return {
                "success": False,
                "error": f"Could not resolve {block.role} '{src}' to a CIF file",
            }
        if path.suffix != ".cif":
            return {"success": False, "error": "TOBACCO building blocks must be CIF files"}
        return {"success": True, "file": str(path)}

    def remove_building_blocks(
        self,
        role: Literal["node", "edge"],
        names: list[str],
        dry_run: bool = True,
    ) -> dict[str, Any]:
        """Building blocks are registered in-memory; removal is a caller-side op.

        Kept for backend-protocol compatibility.  Reports which of *names* are
        known in the catalog but does not touch any files.
        """
        available = set(self._list_cifs(role))
        matched = [n for n in names if n in available or f"{n}.cif" in available]
        return {
            "success": True,
            "dry_run": dry_run,
            "removed": [] if dry_run else matched,
            "would_remove": matched if dry_run else [],
            "count": len(matched),
            "note": "TOBACCO building blocks are registered in-memory; no files removed.",
        }

    def clear_building_blocks(
        self,
        role: Literal["node", "edge"],
        dry_run: bool = True,
    ) -> dict[str, Any]:
        """No-op for the importable backend (registration is in-memory)."""
        return {
            "success": True,
            "dry_run": dry_run,
            "removed": [],
            "count": 0,
            "note": "TOBACCO building blocks are registered in-memory; nothing to clear.",
        }

    def copy_from_database(
        self,
        role: Literal["node", "edge"],
        names: list[str] | None = None,
        source: Path | None = None,
        dry_run: bool = True,
    ) -> dict[str, Any]:
        """List catalog entries or resolve names to CIF paths.

        With the in-memory API there is no active directory to copy into, so this
        reduces to catalog discovery: ``names=None`` lists everything available;
        otherwise it resolves each requested name to its CIF path.
        """
        role_key: Literal["node", "edge", "template"] = role  # type: ignore[assignment]

        if source is not None:
            source = Path(source).resolve()
            if not source.exists():
                return {"success": False, "error": f"Source path not found: {source}"}
            if source.is_file():
                return {
                    "success": True,
                    "resolved": {source.name: str(source)},
                    "count": 1,
                }
            available_map = {f.name: f for f in source.rglob("*.cif")}
        else:
            available_map = {}
            for d in self._search_dirs(role_key):
                for f in d.iterdir():
                    if f.suffix == ".cif":
                        available_map.setdefault(f.name, f)

        available = sorted(available_map)

        if names is None:
            return {
                "success": True,
                "available_in_database": available,
                "count": len(available),
            }

        resolved: dict[str, str] = {}
        not_found: list[str] = []
        for n in names:
            key = n if n.endswith(".cif") else f"{n}.cif"
            if key in available_map:
                resolved[key] = str(available_map[key])
            else:
                not_found.append(n)

        return {
            "success": not not_found,
            "resolved": resolved,
            "not_found": not_found or None,
            "count": len(resolved),
        }

    # ------------------------------------------------------------------ #
    # Configuration (backed by tobacco3.TobaccoConfig)
    # ------------------------------------------------------------------ #
    def get_configuration(self) -> dict[str, Any]:
        """Return the current :class:`tobacco3.TobaccoConfig` as a dict."""
        return {"success": True, "configuration": self._cfg.as_dict()}

    def set_configuration(self, key: str, value: Any) -> dict[str, Any]:
        """Set a single :class:`tobacco3.TobaccoConfig` field."""
        valid = {f.name for f in fields(self._cfg)}
        if key not in valid:
            return {
                "success": False,
                "error": f"Unknown configuration key '{key}'. Valid keys: {sorted(valid)}",
            }
        setattr(self._cfg, key, value)
        return {"success": True, "message": f"Set {key} = {value!r}"}

    # ------------------------------------------------------------------ #
    # Status
    # ------------------------------------------------------------------ #
    def status(self) -> dict[str, Any]:
        """Return overall status of the TOBACCO data directory and config."""
        return {
            "success": True,
            "backend": self.name,
            "tobacco3_version": getattr(self._tobacco3, "__version__", "unknown"),
            "data_dir": str(self._data_dir),
            "templates_available": len(self._list_cifs("template")),
            "nodes_available": len(self._list_cifs("node")),
            "edges_available": len(self._list_cifs("edge")),
            "configuration": self._cfg.as_dict(),
            "ready": True,
        }
