"""Amorphous porous organic polymer (POP) generation via simulated polymerization.

Public API mirrors :mod:`mofforge.build`.  The engine wraps :mod:`pysimm`
(Packmol + LAMMPS + Polymatic) and is imported lazily, so importing this package
never requires the optional ``pop`` dependencies.
"""

from __future__ import annotations

from mofforge.polymerize.base import (
    Monomer,
    PolymerizerBackend,
    PopConfig,
    POPResult,
    ReactiveSite,
)
from mofforge.polymerize.builder import PopBuilder

__all__ = [
    "Monomer",
    "POPResult",
    "PolymerizerBackend",
    "PopBuilder",
    "PopConfig",
    "ReactiveSite",
]
