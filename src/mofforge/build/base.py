"""Shared data types and backend protocol for MOF construction."""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Protocol, runtime_checkable

from mofforge.core.crystal import Crystal



@dataclass
class Topology:
    """A network topology used to construct a MOF."""

    name: str
    source: Path | str | None = None


@dataclass
class BuildingBlock:
    """A node or edge building block for MOF construction."""

    name: str
    role: Literal["node", "edge"]
    source: Path | str
    connection_points: list[int] | None = None

    @property
    def is_smiles(self) -> bool:
        """Return *True* if *source* looks like a SMILES string."""
        s = str(self.source)
        # A SMILES string won't be an existing file path and won't end
        # with a structure-file extension.
        return not s.endswith((".cif", ".xyz", ".mol2")) and not Path(s).is_file()


@dataclass
class BuildResult:
    """Outcome of a MOF build operation."""

    success: bool
    output_paths: list[Path] = field(default_factory=list)
    crystal: Crystal | None = None
    errors: list[str] = field(default_factory=list)
    elapsed_seconds: float = 0.0
    backend: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)


class Timer:
    """Minimal wall-clock timer context manager."""

    def __init__(self) -> None:
        self.elapsed: float = 0.0
        self._start: float = 0.0

    def __enter__(self) -> Timer:
        self._start = time.monotonic()
        return self

    def __exit__(self, *_: object) -> None:
        self.elapsed = time.monotonic() - self._start



@runtime_checkable
class BuilderBackend(Protocol):
    """Contract that every MOF-builder backend must satisfy.

    Backends are **not** required to inherit from this class; they
    only need to implement the same method signatures (structural
    subtyping via :class:`~typing.Protocol`).
    """

    name: str
    """Short identifier for the backend (e.g. ``"tobacco"``)."""

    def build(
        self,
        topology: Topology,
        nodes: list[BuildingBlock],
        edges: list[BuildingBlock],
        output_dir: Path,
        **options: Any,
    ) -> BuildResult:
        """Construct one or more MOF structures."""
        ...

    def list_topologies(self) -> list[str]:
        """Return names of all available topologies."""
        ...

    def describe_topology(self, name: str) -> str:
        """Return a human-readable description of a topology."""
        ...

    def list_building_blocks(self, role: Literal["node", "edge"]) -> list[str]:
        """List registered building blocks for the given role."""
        ...

    def add_building_block(self, block: BuildingBlock) -> dict[str, Any]:
        """Register (or copy) a building block into the backend."""
        ...

    def remove_building_blocks(
        self,
        role: Literal["node", "edge"],
        names: list[str],
        dry_run: bool = True,
    ) -> dict[str, Any]:
        """Remove building blocks.  Dry-run by default."""
        ...

    def clear_building_blocks(
        self,
        role: Literal["node", "edge"],
        dry_run: bool = True,
    ) -> dict[str, Any]:
        """Remove **all** building blocks for the given role."""
        ...

    def copy_from_database(
        self,
        role: Literal["node", "edge"],
        names: list[str] | None = None,
        source: Path | None = None,
        dry_run: bool = True,
    ) -> dict[str, Any]:
        """Copy building blocks from a database / source directory."""
        ...

    def get_configuration(self) -> dict[str, Any]:
        """Return current backend configuration as a dict."""
        ...

    def set_configuration(self, key: str, value: Any) -> dict[str, Any]:
        """Set a single configuration key."""
        ...

    def status(self) -> dict[str, Any]:
        """Return an overall status summary."""
        ...
