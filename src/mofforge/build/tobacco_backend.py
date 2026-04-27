"""TOBACCO 3.0 builder backend."""

from __future__ import annotations

import ast
import importlib
import logging
import os
import re
import shutil
import sys
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Literal

from mofforge.build.base import BuildResult, BuildingBlock, Timer, Topology
from mofforge.build.config import ConfigError, validate_tobacco_path
from mofforge.core.crystal import Crystal

logger = logging.getLogger("mofforge")

# Module names that TOBACCO uses internally.  Cached here to avoid
# duplicating the list in _tobacco_context()'s save and restore blocks.
_TOBACCO_MODULES = (
    "tobacco",
    "configuration",
    "ciftemplate2graph",
    "vertex_edge_assign",
    "cycle_cocyle",
    "bbcif_properties",
    "SBU_geometry",
    "scale",
    "scaled_embedding2coords",
    "place_bbs",
    "remove_net_charge",
    "remove_dummy_atoms",
    "adjust_edges",
    "write_cifs",
    "scale_animation",
)


class TobaccoBackend:
    """Backend that delegates MOF construction to TOBACCO 3.0."""

    name: str = "tobacco"

    def __init__(self, tobacco_path: Path | str) -> None:
        self._root = Path(tobacco_path).resolve()
        errors = validate_tobacco_path(self._root)
        if errors:
            details = "\n  ".join(errors)
            raise ConfigError(f"Invalid TOBACCO installation at {self._root}:\n  {details}")

    @contextmanager
    def _tobacco_context(self) -> Iterator[None]:
        """Temporarily switch to TOBACCO's directory for execution."""
        saved_cwd = os.getcwd()
        saved_path = sys.path[:]
        # Remove any previously cached tobacco modules so a fresh
        # import picks up the current configuration.py on disk.
        tobacco_modules = [k for k in sys.modules if k in _TOBACCO_MODULES]
        saved_modules = {k: sys.modules.pop(k) for k in tobacco_modules}
        try:
            os.chdir(self._root)
            sys.path.insert(0, str(self._root))
            yield
        finally:
            os.chdir(saved_cwd)
            sys.path[:] = saved_path
            # Restore previously cached modules (if any) so we don't
            # permanently pollute the module namespace.
            for k in list(sys.modules):
                if k in _TOBACCO_MODULES:
                    del sys.modules[k]
            sys.modules.update(saved_modules)

    def _dir_for_role(self, role: Literal["node", "edge", "template"]) -> Path:
        """Map a logical role to the corresponding TOBACCO subdirectory."""
        mapping = {
            "node": "nodes",
            "edge": "edges",
            "template": "templates",
        }
        return self._root / mapping[role]

    def _db_dir_for_role(self, role: Literal["node", "edge", "template"]) -> Path:
        """Map a logical role to the corresponding TOBACCO database directory."""
        mapping = {
            "node": "nodes_database",
            "edge": "edges_database",
            "template": "template_database",
        }
        return self._root / mapping[role]

    @staticmethod
    def _list_cifs(directory: Path) -> list[str]:
        """Return sorted CIF filenames in *directory*."""
        if not directory.is_dir():
            return []
        return sorted(f.name for f in directory.iterdir() if f.suffix == ".cif")

    def _collect_outputs(self) -> dict[str, list[str]]:
        """Return ``{subdirectory: [cif_names]}`` from ``output_cifs/``."""
        output_dir = self._root / "output_cifs"
        results: dict[str, list[str]] = {}
        if not output_dir.is_dir():
            return results
        top_level_cifs: list[str] = []
        for entry in sorted(output_dir.iterdir()):
            if entry.is_file() and entry.suffix == ".cif":
                top_level_cifs.append(entry.name)
            elif entry.is_dir():
                cifs = [f.name for f in entry.iterdir() if f.suffix == ".cif"]
                if cifs:
                    results[entry.name] = sorted(cifs)
        if top_level_cifs:
            results["."] = sorted(top_level_cifs)
        return results

    def build(
        self,
        topology: Topology,
        nodes: list[BuildingBlock],
        edges: list[BuildingBlock],
        output_dir: Path,
        *,
        parallel: bool = False,
        ignore_errors: bool = True,
        **options: Any,
    ) -> BuildResult:
        """Run TOBACCO to generate MOF structures."""
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

        template_name = topology.name
        if not template_name.endswith(".cif"):
            template_name += ".cif"

        # Stage building blocks: clear existing files and copy in
        # the explicitly passed nodes/edges so TOBACCO sees exactly
        # what the caller specified.
        stage_errors = self._stage_building_blocks(nodes, edges)
        if stage_errors:
            return BuildResult(
                success=False,
                errors=stage_errors,
                backend=self.name,
            )

        # Snapshot existing outputs so we can detect new ones.
        initial_outputs = self._snapshot_outputs()

        errors: list[str] = []
        processed: list[str] = []
        templates: list[str] = []

        timer = Timer()
        with timer, self._tobacco_context():
            try:
                tobacco_mod = importlib.import_module("tobacco")
                run_template = tobacco_mod.run_template
            except Exception as exc:
                return BuildResult(
                    success=False,
                    errors=[f"Failed to import tobacco: {exc}"],
                    backend=self.name,
                )

            # Clean stray .DS_Store files
            for d in ("templates", "nodes", "edges"):
                ds = Path(d) / ".DS_Store"
                if ds.exists():
                    ds.unlink()

            # Determine template list
            templates_dir = Path("templates")
            if template_name == "all.cif" or template_name == "all":
                templates = sorted(f.name for f in templates_dir.iterdir() if f.suffix == ".cif")
            else:
                if not (templates_dir / template_name).exists():
                    return BuildResult(
                        success=False,
                        errors=[
                            f"Template '{template_name}' not found in {self._root / 'templates'}"
                        ],
                        backend=self.name,
                    )
                templates = [template_name]

            if parallel:
                try:
                    run_parallel = tobacco_mod.run_tobacco_parallel
                    cfg = importlib.import_module("configuration")
                    run_parallel(templates, cfg.CHARGES)
                    processed = list(templates)
                except Exception as exc:
                    errors.append(f"Parallel execution failed: {exc}")
            else:
                for tmpl in templates:
                    try:
                        run_template(tmpl)
                        processed.append(tmpl)
                    except Exception as exc:
                        msg = f"Template '{tmpl}' failed: {exc}"
                        if not ignore_errors:
                            errors.append(msg)
                            break
                        errors.append(msg)
                        logger.warning(msg)

        # Detect new output CIFs
        current_outputs = self._snapshot_outputs()
        new_paths = sorted(current_outputs - initial_outputs)

        # Copy new outputs to the requested output_dir
        copied: list[Path] = []
        for p in new_paths:
            dst = output_dir / p.name
            try:
                shutil.copy2(p, dst)
                copied.append(dst)
            except Exception as exc:
                errors.append(f"Failed to copy {p} -> {dst}: {exc}")

        # Load the first output as a Crystal
        crystal: Crystal | None = None
        if copied:
            try:
                from mofforge.core.bonding import infer_bonds

                crystal = Crystal.from_cif(str(copied[0]))
                crystal = infer_bonds(crystal, periodic=True)
            except Exception as exc:
                logger.warning("Could not load output CIF as Crystal: %s", exc)

        success = len(errors) == 0 or (ignore_errors and len(processed) > 0)

        return BuildResult(
            success=success,
            output_paths=copied,
            crystal=crystal,
            errors=errors,
            elapsed_seconds=round(timer.elapsed, 2),
            backend=self.name,
            metadata={
                "templates_requested": templates,
                "templates_processed": processed,
                "parallel": parallel,
            },
        )

    def _snapshot_outputs(self) -> set[Path]:
        """Return the set of all CIF paths currently in ``output_cifs/``."""
        output_dir = self._root / "output_cifs"
        paths: set[Path] = set()
        if output_dir.is_dir():
            for entry in output_dir.iterdir():
                if entry.is_file() and entry.suffix == ".cif":
                    paths.add(entry.resolve())
                elif entry.is_dir():
                    for f in entry.iterdir():
                        if f.suffix == ".cif":
                            paths.add(f.resolve())
        return paths

    def _stage_building_blocks(
        self,
        nodes: list[BuildingBlock],
        edges: list[BuildingBlock],
    ) -> list[str]:
        """Clear TOBACCO's ``nodes/`` and ``edges/`` dirs and copy in the given blocks."""
        errors: list[str] = []

        if not nodes and not edges:
            logger.debug("build() called with no nodes or edges; skipping staging")
            return errors

        # Clear existing CIFs in nodes/ and edges/
        for role in ("node", "edge"):
            self.clear_building_blocks(role, dry_run=False)

        # Copy in the passed building blocks
        for block in nodes:
            result = self.add_building_block(block)
            if not result.get("success"):
                errors.append(
                    f"Failed to stage node '{block.name}': {result.get('error', 'unknown error')}"
                )

        for block in edges:
            result = self.add_building_block(block)
            if not result.get("success"):
                errors.append(
                    f"Failed to stage edge '{block.name}': {result.get('error', 'unknown error')}"
                )

        if nodes or edges:
            staged_nodes = self._list_cifs(self._dir_for_role("node"))
            staged_edges = self._list_cifs(self._dir_for_role("edge"))
            logger.debug(
                "Staged %d node(s) and %d edge(s) for TOBACCO build",
                len(staged_nodes),
                len(staged_edges),
            )

        return errors

    def list_topologies(self) -> list[str]:
        """List available template CIF files."""
        return self._list_cifs(self._dir_for_role("template"))

    def describe_topology(self, name: str) -> str:
        """Return basic information about a template."""
        if not name.endswith(".cif"):
            name += ".cif"
        path = self._dir_for_role("template") / name
        if not path.exists():
            return f"Template '{name}' not found."
        size = path.stat().st_size
        return f"Template: {name} ({size} bytes)"

    def list_building_blocks(self, role: Literal["node", "edge"]) -> list[str]:
        """List CIF files in the ``nodes/`` or ``edges/`` directory."""
        return self._list_cifs(self._dir_for_role(role))

    def add_building_block(self, block: BuildingBlock) -> dict[str, Any]:
        """Copy a CIF file into the appropriate TOBACCO directory."""
        src = Path(block.source).resolve()
        if not src.is_file():
            return {"success": False, "error": f"Source file not found: {src}"}
        if not src.suffix == ".cif":
            return {"success": False, "error": "TOBACCO building blocks must be CIF files"}

        dest_dir = self._dir_for_role(block.role)
        dest = dest_dir / (block.name + ".cif" if not block.name.endswith(".cif") else block.name)
        try:
            shutil.copy2(src, dest)
            return {"success": True, "file": str(dest)}
        except Exception as exc:
            return {"success": False, "error": str(exc)}

    def remove_building_blocks(
        self,
        role: Literal["node", "edge"],
        names: list[str],
        dry_run: bool = True,
    ) -> dict[str, Any]:
        """Remove specific CIF files from a TOBACCO directory."""
        target_dir = self._dir_for_role(role)
        available = set(self._list_cifs(target_dir))

        to_remove = [n for n in names if n in available]
        not_found = [n for n in names if n not in available]

        if dry_run:
            return {
                "success": True,
                "dry_run": True,
                "would_remove": to_remove,
                "not_found": not_found if not_found else None,
                "count": len(to_remove),
            }

        removed: list[str] = []
        errors: list[dict[str, str]] = []
        for f in to_remove:
            try:
                (target_dir / f).unlink()
                removed.append(f)
            except Exception as exc:
                errors.append({"file": f, "error": str(exc)})

        return {
            "success": len(errors) == 0,
            "removed": removed,
            "count": len(removed),
            "errors": errors if errors else None,
            "dry_run": False,
        }

    def clear_building_blocks(
        self,
        role: Literal["node", "edge"],
        dry_run: bool = True,
    ) -> dict[str, Any]:
        """Remove **all** CIF files for the given role."""
        target_dir = self._dir_for_role(role)
        cifs = self._list_cifs(target_dir)

        if dry_run:
            return {
                "success": True,
                "dry_run": True,
                "would_remove": cifs,
                "count": len(cifs),
            }

        removed: list[str] = []
        errors: list[dict[str, str]] = []
        for f in cifs:
            try:
                (target_dir / f).unlink()
                removed.append(f)
            except Exception as exc:
                errors.append({"file": f, "error": str(exc)})

        return {
            "success": len(errors) == 0,
            "removed": removed,
            "count": len(removed),
            "errors": errors if errors else None,
            "dry_run": False,
        }

    def copy_from_database(
        self,
        role: Literal["node", "edge"],
        names: list[str] | None = None,
        source: Path | None = None,
        dry_run: bool = True,
    ) -> dict[str, Any]:
        """Copy building blocks from a database directory into the active directory."""
        # Also support "template" role for internal use
        role_key: Literal["node", "edge", "template"] = role  # type: ignore[assignment]
        target_dir = self._dir_for_role(role_key)

        if source is not None:
            source = Path(source).resolve()
            if not source.exists():
                return {"success": False, "error": f"Source path not found: {source}"}
            if source.is_file():
                # Copy a single file
                if not source.suffix == ".cif":
                    return {"success": False, "error": "Source file must be a .cif file"}
                if dry_run:
                    return {
                        "success": True,
                        "dry_run": True,
                        "would_copy": [source.name],
                        "source": str(source),
                        "destination": str(target_dir),
                        "count": 1,
                    }
                dst = target_dir / source.name
                shutil.copy2(source, dst)
                return {
                    "success": True,
                    "copied": [source.name],
                    "source": str(source),
                    "destination": str(target_dir),
                    "count": 1,
                    "dry_run": False,
                }
            db_dir = source
        else:
            db_dir = self._db_dir_for_role(role_key)
            if not db_dir.is_dir():
                return {
                    "success": False,
                    "error": f"Database directory not found: {db_dir}",
                }

        # Collect available CIFs in the database (recursively)
        available: list[str] = sorted({f.name for f in db_dir.rglob("*.cif")})

        if names is None:
            result: dict[str, Any] = {
                "success": True,
                "available_in_database": available,
                "count": len(available),
            }
            if source is not None:
                result["source"] = str(source)
            return result

        # Determine which files to copy
        available_set = set(available)
        to_copy = [n for n in names if n in available_set]
        not_found = [n for n in names if n not in available_set]

        if dry_run:
            result = {
                "success": True,
                "dry_run": True,
                "would_copy": to_copy,
                "not_found": not_found if not_found else None,
                "count": len(to_copy),
            }
            if source is not None:
                result["source"] = str(source)
            return result

        # Actually copy
        copied: list[str] = []
        errors: list[dict[str, str]] = []
        for fname in to_copy:
            # Find the file in the database tree
            matches = list(db_dir.rglob(fname))
            if matches:
                dst = target_dir / fname
                try:
                    shutil.copy2(matches[0], dst)
                    copied.append(fname)
                except Exception as exc:
                    errors.append({"file": fname, "error": str(exc)})

        result = {
            "success": len(errors) == 0,
            "copied": copied,
            "count": len(copied),
            "errors": errors if errors else None,
            "dry_run": False,
        }
        if source is not None:
            result["source"] = str(source)
        return result

    def get_configuration(self) -> dict[str, Any]:
        """Read TOBACCO's ``configuration.py`` and return values as a dict."""
        config_file = self._root / "configuration.py"
        if not config_file.is_file():
            return {"success": False, "error": "configuration.py not found"}

        config: dict[str, Any] = {}
        content = config_file.read_text()
        for line in content.strip().splitlines():
            line = line.strip()
            if line and "=" in line and not line.startswith("#"):
                key, _, value = line.partition("=")
                key = key.strip()
                value = value.strip()
                try:
                    config[key] = ast.literal_eval(value)
                except (ValueError, SyntaxError):
                    config[key] = value

        return {"success": True, "configuration": config, "file": str(config_file)}

    def set_configuration(self, key: str, value: Any) -> dict[str, Any]:
        """Set a single key in TOBACCO's ``configuration.py``."""
        config_file = self._root / "configuration.py"
        if not config_file.is_file():
            return {"success": False, "error": "configuration.py not found"}

        # Determine the text to write as the value.
        if isinstance(value, str):
            try:
                ast.literal_eval(value)
                value_text = value  # already a valid Python literal
            except (ValueError, SyntaxError):
                value_text = repr(value)
        else:
            value_text = repr(value)

        lines = config_file.read_text().splitlines(keepends=True)
        found = False
        key_pattern = re.compile(rf"^{re.escape(key)}\s*=")
        new_lines: list[str] = []
        for line in lines:
            stripped = line.strip()
            if key_pattern.match(stripped):
                new_lines.append(f"{key} = {value_text}\n")
                found = True
            else:
                new_lines.append(line)

        if not found:
            return {"success": False, "error": f"Configuration key '{key}' not found"}

        config_file.write_text("".join(new_lines))
        return {"success": True, "message": f"Set {key} = {value_text}", "file": str(config_file)}

    def status(self) -> dict[str, Any]:
        """Return overall status of the TOBACCO installation."""
        templates = self._list_cifs(self._dir_for_role("template"))
        nodes = self._list_cifs(self._dir_for_role("node"))
        edges = self._list_cifs(self._dir_for_role("edge"))
        outputs = self._collect_outputs()
        config = self.get_configuration()

        return {
            "success": True,
            "project_root": str(self._root),
            "templates_available": len(templates),
            "nodes_available": len(nodes),
            "edges_available": len(edges),
            "output_directories": len(outputs),
            "configuration": config.get("configuration", {}),
            "ready": True,
        }

    def list_outputs(self) -> dict[str, Any]:
        """List generated output CIF files organised by subdirectory."""
        outputs = self._collect_outputs()
        entries = [
            {"directory": d, "cif_files": cifs, "count": len(cifs)} for d, cifs in outputs.items()
        ]
        return {
            "success": True,
            "output_directories": len(entries),
            "outputs": entries,
            "directory": str(self._root / "output_cifs"),
        }
