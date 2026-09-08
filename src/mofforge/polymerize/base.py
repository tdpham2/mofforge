"""Shared data types and backend protocol for amorphous POP generation.

These mirror :mod:`mofforge.build.base` so that a polymer build looks like any
other mofforge structure operation.  A porous organic polymer (CMP / PIM / HCP)
is an *amorphous* covalent network: unlike a MOF or crystalline COF it has no
periodic net, so it is produced by **simulated polymerization** (pack monomers
into a box, then iteratively bond nearby reactive sites under MD relaxation)
rather than net-based assembly.  mofforge does not implement that loop itself; it
orchestrates :mod:`pysimm` (Packmol + LAMMPS + Polymatic).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol, runtime_checkable

from mofforge.core.crystal import Crystal
from mofforge.validation import ValidationReport


@dataclass
class ReactiveSite:
    """One reactive attachment point on a monomer.

    This is the POP analogue of a MOF connection point / ``!`` anchor.  During
    simulated polymerization each site becomes a Polymatic *linker atom*: a
    heavy atom that forms one new inter-monomer covalent bond.

    Attributes
    ----------
    atom_idx:
        Index (within the monomer) of the heavy atom that forms the new bond.
    anchor_idx:
        Index of a neighbouring heavy atom that defines the outward bonding
        direction, or ``None`` when it is left to the engine.
    site_type:
        Reaction role, e.g. ``"amine"``, ``"aldehyde"``, ``"aryl_halide"``.
        Two sites bond only when :mod:`mofforge.polymerize.reactions` marks
        their types compatible.
    """

    atom_idx: int
    anchor_idx: int | None = None
    site_type: str = "generic"


@dataclass
class Monomer:
    """A monomer used to build a porous organic polymer.

    ``source`` is either a SMILES string or a path to an XYZ/CIF file (mirroring
    :class:`mofforge.build.base.BuildingBlock`).  ``functionality`` is the number
    of reactive sites: 2 gives a linear chain, >= 3 forms a cross-linked network.
    """

    name: str
    source: Path | str
    functionality: int = 2
    sites: list[ReactiveSite] = field(default_factory=list)

    @property
    def is_smiles(self) -> bool:
        """Return *True* if *source* looks like a SMILES string.

        A SMILES string won't be an existing file path and won't end with a
        structure-file extension.  (Same heuristic as ``BuildingBlock``.)
        """
        s = str(self.source)
        return not s.endswith((".cif", ".xyz", ".mol2")) and not Path(s).is_file()


@dataclass
class PopConfig:
    """Settings for a simulated-polymerization run.

    ``box_length`` (Angstrom) is derived from ``target_density`` and the packed
    monomer masses when left as ``None``.
    """

    box_length: float | None = None
    target_density: float = 0.8  # g/cm^3
    forcefield: str = "gaff2"  # gaff2 | dreiding | pcff
    target_conversion: float = 0.95  # fraction of reactive sites consumed
    n_monomers: int | None = None  # total monomers packed; None -> engine default
    equilibrate: bool = True  # run compression/decompression after bonding
    random_seed: int | None = None
    md_settings: dict[str, Any] = field(default_factory=dict)


@dataclass
class POPResult:
    """Outcome of a polymer build.

    Field-for-field identical to :class:`mofforge.build.base.BuildResult` so that
    downstream code (validation, manifests, batch processing) treats POPs like
    any other built structure.
    """

    success: bool
    output_paths: list[Path] = field(default_factory=list)
    crystal: Crystal | None = None
    errors: list[str] = field(default_factory=list)
    elapsed_seconds: float = 0.0
    backend: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)
    validation: ValidationReport | None = None


class Timer:
    """Minimal wall-clock timer context manager (copied from build.base)."""

    def __init__(self) -> None:
        self.elapsed: float = 0.0
        self._start: float = 0.0

    def __enter__(self) -> Timer:
        self._start = time.monotonic()
        return self

    def __exit__(self, *_: object) -> None:
        self.elapsed = time.monotonic() - self._start


@runtime_checkable
class PolymerizerBackend(Protocol):
    """Contract every polymerization backend must satisfy.

    Backends are **not** required to inherit from this class; they only need to
    implement the same method signatures (structural subtyping via
    :class:`~typing.Protocol`).
    """

    name: str
    """Short identifier for the backend (e.g. ``"pysimm"``)."""

    def polymerize(
        self,
        monomers: list[Monomer],
        config: PopConfig,
        output_dir: Path,
        **options: Any,
    ) -> POPResult:
        """Run simulated polymerization and return the resulting structure."""
        ...
