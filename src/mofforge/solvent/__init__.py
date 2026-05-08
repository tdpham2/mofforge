"""Automatic solvent identification and removal for MOF structures."""

from mofforge.solvent.removal import (
    RemovedMolecule,
    SolventRemovalResult,
    remove_solvent,
)

__all__ = [
    "RemovedMolecule",
    "SolventRemovalResult",
    "remove_solvent",
]
