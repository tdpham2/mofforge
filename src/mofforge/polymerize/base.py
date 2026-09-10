"""Public inputs and results for geometric polymer-box construction."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from mofforge.core.crystal import Crystal
from mofforge.validation import ValidationReport

if TYPE_CHECKING:
    from mofforge.polymerize.state import ConstructionState


def positive(value, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a positive finite number.")
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a positive finite number.")
    return float(value)


def integer(value, name: str, minimum: int = 1) -> int:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}.")
    return value


def box_lengths(value) -> tuple[float, float, float]:
    if len(value) != 3:
        raise ValueError("box_lengths requires three orthorhombic lengths in angstroms.")
    return tuple(positive(v, "box length") for v in value)


def reject_removed(options: dict) -> None:
    replacements = {
        "target_density": "Use initial_packing_density (no universal default).",
        "box_length": "Use box_lengths=[Lx, Ly, Lz].",
        "n_monomers": "Supply count for each component.",
        "sites": "Supply explicit connectors; SMARTS detection is annotation only.",
        "forcefield": "Configure force fields and simulations in MatKit.",
        "equilibrate": "Run equilibration externally in MatKit.",
        "md_settings": "Configure MD externally in MatKit.",
        "lammps_bin": "Configure LAMMPS externally in MatKit.",
    }
    for name, message in replacements.items():
        if name in options:
            raise TypeError(f"Removed POP option {name!r}. {message}")


@dataclass(frozen=True)
class Connector:
    """Template site; atom -> anchor is its outgoing direction.

    A replaceable hydrogen is a useful anchor. PopBuilder.prepare() exposes
    template labels. Each connector is consumed by exactly one event.
    """

    label: str
    atom: str
    anchor: str
    role: str
    replaceable: tuple[str, ...] = ()

    def __post_init__(self):
        for value in (self.label, self.atom, self.anchor, self.role):
            if not isinstance(value, str) or not value:
                raise ValueError("Connector labels and roles must be nonempty strings.")
        object.__setattr__(self, "replaceable", tuple(self.replaceable))
        if self.atom == self.anchor or self.atom in self.replaceable:
            raise ValueError("A connector needs a distinct anchor and retained attachment atom.")
        if len(set(self.replaceable)) != len(self.replaceable):
            raise ValueError("Duplicate replaceable atom labels.")


@dataclass(frozen=True)
class ConnectionRule:
    """Complete pairwise edits, scoped to the two selected instances.

    delete_atoms must explicitly contain a and b lists. $replaceable expands
    the selected connector's declared labels. Other edits reference literal
    template labels, or $atom / $anchor. Bond changes contain side, atom1,
    atom2, order; charge changes contain side, atom, charge. No chemistry is
    inferred. angle_tolerance limits each outgoing vector's deviation from
    the new bond direction, in degrees.
    """

    name: str
    roles: tuple[str, str]
    bond_order: float
    bond_length: tuple[float, float]
    delete_atoms: dict[str, list[str]]
    bond_changes: tuple[dict, ...] = ()
    charge_changes: tuple[dict, ...] = ()
    minimum_cycle_size: int | None = None
    angle_tolerance: float = 30.0

    def __post_init__(self):
        object.__setattr__(self, "roles", tuple(self.roles))
        object.__setattr__(self, "bond_length", tuple(self.bond_length))
        if not self.name or len(self.roles) != 2 or not all(self.roles):
            raise ValueError("A connection rule needs a name and two roles.")
        if self.bond_order not in (1, 1.5, 2, 3) or isinstance(self.bond_order, bool):
            raise ValueError("Supported bond orders are 1, 1.5, 2, and 3.")
        if len(self.bond_length) != 2:
            raise ValueError("bond_length requires a minimum and maximum in angstroms.")
        low, high = [positive(v, "bond length") for v in self.bond_length]
        if low > high:
            raise ValueError("bond_length minimum exceeds maximum.")
        if set(self.delete_atoms) != {"a", "b"} or any(
            not isinstance(v, (list, tuple)) or any(not isinstance(x, str) for x in v)
            for v in self.delete_atoms.values()
        ):
            raise ValueError("delete_atoms must explicitly provide a and b label lists.")
        if not np.isfinite(self.angle_tolerance) or not 0 <= self.angle_tolerance <= 180:
            raise ValueError("angle_tolerance must be between 0 and 180 degrees.")
        if self.minimum_cycle_size is not None:
            integer(self.minimum_cycle_size, "minimum_cycle_size", 3)
        for edit in self.bond_changes:
            if set(edit) != {"side", "atom1", "atom2", "order"}:
                raise ValueError("Bond changes require side, atom1, atom2, order.")
            if edit["side"] not in ("a", "b") or edit["order"] not in (1, 1.5, 2, 3):
                raise ValueError("Invalid bond change side or order.")
        for edit in self.charge_changes:
            if set(edit) != {"side", "atom", "charge"}:
                raise ValueError("Charge changes require side, atom, charge.")
            if edit["side"] not in ("a", "b") or type(edit["charge"]) is not int:
                raise ValueError("Invalid charge change side or formal charge.")


@dataclass
class Monomer:
    name: str
    source: str | Path | Crystal
    count: int
    connectors: list[Connector] = field(default_factory=list)
    functionality: int | None = None
    graph: dict | None = None

    def __post_init__(self):
        integer(self.count, "count")
        self.connectors = [Connector(**c) if isinstance(c, dict) else c for c in self.connectors]
        if len({c.label for c in self.connectors}) != len(self.connectors):
            raise ValueError("Connector labels must be unique within a template.")
        if self.functionality is not None and (
            type(self.functionality) is not int or self.functionality != len(self.connectors)
        ):
            raise ValueError("functionality must equal the number of declared connectors.")
        self.functionality = len(self.connectors)


@dataclass
class POPResult:
    success: bool
    output_paths: list[Path] = field(default_factory=list)
    crystal: Crystal | None = None
    errors: list[str] = field(default_factory=list)
    elapsed_seconds: float = 0.0
    backend: str = "native"
    metadata: dict[str, Any] = field(default_factory=dict)
    validation: ValidationReport | None = None
    operation: str = "pack"
    status: str = "failed"
    state: ConstructionState | None = None
    state_path: Path | None = None

    def to_dict(self) -> dict:
        return {
            "success": self.success,
            "operation": self.operation,
            "status": self.status,
            "backend": self.backend,
            "errors": self.errors,
            "output_paths": [str(p) for p in self.output_paths],
            "state_path": str(self.state_path) if self.state_path else None,
            "state_hash": self.state.state_hash if self.state else None,
            "elapsed_seconds": self.elapsed_seconds,
            "metadata": self.metadata,
            "atoms": self.crystal.n_atoms if self.crystal is not None else None,
            "validation": self.validation.to_dict() if self.validation else None,
        }
