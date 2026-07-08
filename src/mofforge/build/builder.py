"""Unified facade for MOF construction."""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import Any, Literal

from mofforge.build.base import BuilderBackend, BuildingBlock, BuildResult, Topology
from mofforge.build.config import BuildConfig

logger = logging.getLogger("mofforge")

_BACKENDS = ("tobacco", "pormake")


class MOFBuilder:
    """Unified interface for building MOFs from topology + building blocks."""

    def __init__(self, backend: str = "tobacco", **kwargs: Any) -> None:
        if backend not in _BACKENDS:
            raise ValueError(f"Unknown backend {backend!r}.  Choose from: {_BACKENDS}")

        self._backend_name = backend
        self._nodes: list[BuildingBlock] = []
        self._edges: list[BuildingBlock] = []

        # Load merged configuration (toml + env + kwargs)
        cfg = BuildConfig.load(**kwargs)

        if backend == "tobacco":
            from mofforge.build.tobacco_backend import TobaccoBackend

            # Resolve the data directory (auto-detected next to the installed
            # tobacco3 package when not explicitly configured).
            data_dir = cfg.resolve_tobacco_data_dir()  # raises ConfigError if unusable
            self._backend: BuilderBackend = TobaccoBackend(data_dir)  # type: ignore[assignment]

        elif backend == "pormake":
            from mofforge.build.pormake_backend import PormakeBackend

            output_dir = (
                kwargs.get("output_dir")
                or kwargs.get("pormake_output_dir")
                or cfg.pormake_output_dir
            )
            bb_dir = kwargs.get("bb_dir")
            self._backend = PormakeBackend(output_dir=output_dir, bb_dir=bb_dir)  # type: ignore[assignment]

    @property
    def backend(self) -> BuilderBackend:
        """The underlying backend instance."""
        return self._backend

    @property
    def backend_name(self) -> str:
        """Short name of the active backend."""
        return self._backend_name

    def list_topologies(self) -> list[str]:
        """Return names of all available topologies."""
        return self._backend.list_topologies()

    def describe_topology(self, name: str) -> str:
        """Return a human-readable description of a topology."""
        return self._backend.describe_topology(name)

    def add_node(
        self,
        source: str | Path,
        name: str | None = None,
        connection_points: list[int] | None = None,
    ) -> None:
        """Register a node (metal cluster / SBU) building block."""
        block = self._make_block("node", source, name, connection_points)
        result = self._backend.add_building_block(block)
        if not result.get("success", False):
            raise ValueError(
                f"Failed to register node '{block.name}': "
                f"{result.get('error', 'unknown error')}"
            )
        self._nodes.append(block)

    def add_edge(
        self,
        source: str | Path,
        name: str | None = None,
        connection_points: list[int] | None = None,
    ) -> None:
        """Register an edge (organic linker) building block."""
        block = self._make_block("edge", source, name, connection_points)
        result = self._backend.add_building_block(block)
        if not result.get("success", False):
            raise ValueError(
                f"Failed to register edge '{block.name}': "
                f"{result.get('error', 'unknown error')}"
            )
        self._edges.append(block)

    def list_nodes(self) -> list[str]:
        """List registered node building blocks."""
        return self._backend.list_building_blocks("node")

    def list_edges(self) -> list[str]:
        """List registered edge building blocks."""
        return self._backend.list_building_blocks("edge")

    def remove_nodes(self, names: list[str], dry_run: bool = True) -> dict[str, Any]:
        """Remove node building blocks."""
        result = self._backend.remove_building_blocks("node", names, dry_run=dry_run)
        if not dry_run:
            self._nodes = [n for n in self._nodes if n.name not in names]
        return result

    def remove_edges(self, names: list[str], dry_run: bool = True) -> dict[str, Any]:
        """Remove edge building blocks."""
        result = self._backend.remove_building_blocks("edge", names, dry_run=dry_run)
        if not dry_run:
            self._edges = [e for e in self._edges if e.name not in names]
        return result

    def clear_nodes(self, dry_run: bool = True) -> dict[str, Any]:
        """Remove all node building blocks."""
        result = self._backend.clear_building_blocks("node", dry_run=dry_run)
        if not dry_run:
            self._nodes.clear()
        return result

    def clear_edges(self, dry_run: bool = True) -> dict[str, Any]:
        """Remove all edge building blocks."""
        result = self._backend.clear_building_blocks("edge", dry_run=dry_run)
        if not dry_run:
            self._edges.clear()
        return result

    def build(
        self,
        topology: str,
        output_dir: str | Path = ".",
        **options: Any,
    ) -> BuildResult:
        """Build a MOF from the given topology and registered building blocks."""
        topo = Topology(name=topology)
        return self._backend.build(
            topology=topo,
            nodes=list(self._nodes),
            edges=list(self._edges),
            output_dir=Path(output_dir),
            **options,
        )

    def status(self) -> dict[str, Any]:
        """Return overall status of the active backend."""
        return self._backend.status()

    def get_configuration(self) -> dict[str, Any]:
        """Return the backend's current configuration."""
        return self._backend.get_configuration()

    def set_configuration(self, key: str, value: Any) -> dict[str, Any]:
        """Set a configuration key on the backend."""
        return self._backend.set_configuration(key, value)

    def copy_from_database(
        self,
        role: Literal["node", "edge"],
        names: list[str] | None = None,
        source: Path | None = None,
        dry_run: bool = True,
    ) -> dict[str, Any]:
        """Copy building blocks from the backend's database.

        Call with ``names=None`` to list what is available.
        """
        return self._backend.copy_from_database(
            role=role, names=names, source=source, dry_run=dry_run
        )

    @staticmethod
    def _make_block(
        role: Literal["node", "edge"],
        source: str | Path,
        name: str | None,
        connection_points: list[int] | None,
    ) -> BuildingBlock:
        """Create a :class:`BuildingBlock` with sensible defaults."""
        source_str = str(source)

        if name is None:
            p = Path(source_str)
            if p.suffix in (".cif", ".xyz", ".mol2"):
                name = p.stem
            else:
                # Assume SMILES or database name – use a prefix + hash
                # to avoid collisions between similar strings.
                digest = hashlib.sha256(source_str.encode()).hexdigest()[:8]
                name = f"{source_str[:12].replace('/', '_')}_{digest}"

        return BuildingBlock(
            name=name,
            role=role,
            source=source,
            connection_points=connection_points,
        )
