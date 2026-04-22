"""MOFBuilder: unified facade for constructing MOFs.

Provides a single, backend-agnostic interface for building MOFs from
topology + building blocks.  Delegates to :class:`TobaccoBackend` or
:class:`PormakeBackend` depending on the ``backend`` parameter.

Example::

    from mofforge.build import MOFBuilder

    # TOBACCO
    builder = MOFBuilder(backend="tobacco")
    builder.list_topologies()
    result = builder.build(topology="pcu")
    crystal = result.crystal

    # Pormake
    builder = MOFBuilder(backend="pormake")
    builder.add_node("path/to/node.xyz")
    builder.add_edge("path/to/edge.xyz")
    result = builder.build(topology="pcu", output_dir="./output")
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Literal

from mofforge.build.base import BuilderBackend, BuildingBlock, BuildResult, Topology
from mofforge.build.config import BuildConfig, ConfigError

logger = logging.getLogger("mofforge")

_BACKENDS = ("tobacco", "pormake")


class MOFBuilder:
    """Unified interface for building MOFs from topology + building blocks.

    Args:
        backend: ``"tobacco"`` or ``"pormake"``.
        **kwargs: Backend-specific options.  Notable keys:

            - ``tobacco_path``: Path to the TOBACCO installation directory
              (overrides ``mofforge.toml`` / env var).
            - ``output_dir``: Default output directory for generated CIFs.
            - ``bb_dir``: Directory for storing building-block files
              (pormake only).
            - ``pormake_output_dir``: Alias for ``output_dir`` when using
              pormake.

    Raises:
        ValueError: If *backend* is not recognised.
        ConfigError: If required configuration (e.g. TOBACCO path) is
            missing or invalid.
    """

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

            tobacco_path = kwargs.get("tobacco_path") or cfg.tobacco_path
            if tobacco_path is None:
                cfg.resolve_tobacco_path()  # raises ConfigError with instructions
            self._backend: BuilderBackend = TobaccoBackend(tobacco_path)  # type: ignore[assignment]

        elif backend == "pormake":
            from mofforge.build.pormake_backend import PormakeBackend

            output_dir = (
                kwargs.get("output_dir")
                or kwargs.get("pormake_output_dir")
                or cfg.pormake_output_dir
            )
            bb_dir = kwargs.get("bb_dir")
            self._backend = PormakeBackend(output_dir=output_dir, bb_dir=bb_dir)  # type: ignore[assignment]

    # ------------------------------------------------------------------ #
    # Properties
    # ------------------------------------------------------------------ #

    @property
    def backend(self) -> BuilderBackend:
        """The underlying backend instance."""
        return self._backend

    @property
    def backend_name(self) -> str:
        """Short name of the active backend."""
        return self._backend_name

    # ------------------------------------------------------------------ #
    # Topology
    # ------------------------------------------------------------------ #

    def list_topologies(self) -> list[str]:
        """Return names of all available topologies."""
        return self._backend.list_topologies()

    def describe_topology(self, name: str) -> str:
        """Return a human-readable description of a topology."""
        return self._backend.describe_topology(name)

    # ------------------------------------------------------------------ #
    # Building blocks
    # ------------------------------------------------------------------ #

    def add_node(
        self,
        source: str | Path,
        name: str | None = None,
        connection_points: list[int] | None = None,
    ) -> None:
        """Register a node (metal cluster / SBU) building block.

        Args:
            source: Path to a CIF or XYZ file, or a SMILES string.
            name: Optional label.  Defaults to the filename stem or
                the first 20 characters of the SMILES.
            connection_points: Atom indices for connection sites
                (required for SMILES-based blocks in pormake).
        """
        block = self._make_block("node", source, name, connection_points)
        self._nodes.append(block)
        # Also register with the backend for file-management tracking
        self._backend.add_building_block(block)

    def add_edge(
        self,
        source: str | Path,
        name: str | None = None,
        connection_points: list[int] | None = None,
    ) -> None:
        """Register an edge (organic linker) building block.

        Args:
            source: Path to a CIF or XYZ file, or a SMILES string.
            name: Optional label.
            connection_points: Atom indices for connection sites.
        """
        block = self._make_block("edge", source, name, connection_points)
        self._edges.append(block)
        self._backend.add_building_block(block)

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

    # ------------------------------------------------------------------ #
    # Build
    # ------------------------------------------------------------------ #

    def build(
        self,
        topology: str,
        output_dir: str | Path = ".",
        **options: Any,
    ) -> BuildResult:
        """Build a MOF from topology + registered building blocks.

        Args:
            topology: Topology name (e.g. ``"pcu"`` or ``"pcu.cif"``).
            output_dir: Where to write the generated CIF files.
            **options: Backend-specific options (e.g. ``parallel=True``
                for TOBACCO, ``accuracy=10`` for pormake).

        Returns:
            A :class:`BuildResult` with the outcome.
        """
        topo = Topology(name=topology)
        return self._backend.build(
            topology=topo,
            nodes=list(self._nodes),
            edges=list(self._edges),
            output_dir=Path(output_dir),
            **options,
        )

    # ------------------------------------------------------------------ #
    # Management / Configuration
    # ------------------------------------------------------------------ #

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

    # ------------------------------------------------------------------ #
    # Private helpers
    # ------------------------------------------------------------------ #

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
                # Assume SMILES or database name – use first 20 chars
                name = source_str[:20].replace("/", "_")

        return BuildingBlock(
            name=name,
            role=role,
            source=source,
            connection_points=connection_points,
        )
