"""Curated menu of functional groups for MOF linker functionalization.

Each group is defined by *chemistry*, not geometry: a SMILES fragment plus the
index of its attachment atom (the atom that bonds to the aromatic ring in place
of the removed hydrogen).  The 3-D coordinates are generated deterministically
by RDKit at fragment-build time (see :mod:`mofforge.functionalize.generate`), so
an AI agent only ever picks a group *name* from this menu — it never authors a
SMILES string or coordinates.

Adding a new group is a one-line, offline edit to :data:`_GROUPS` by a chemist.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FunctionalGroup:
    """A substituent that can replace an aromatic hydrogen.

    Attributes
    ----------
    name:
        Short, human-facing label used by the agent (e.g. ``"NH2"``).
    smiles:
        SMILES of the group as it attaches to the ring.  The first atom in the
        SMILES is, by convention, the attachment atom unless ``attach_idx`` says
        otherwise.
    attach_idx:
        Index (into the SMILES atom order) of the atom that forms the bond to
        the aromatic ring carbon.
    description:
        One-line human description.
    """

    name: str
    smiles: str
    attach_idx: int
    description: str


# The curated palette.  Keyed by canonical name.  ``H`` is included so an agent
# can explicitly request "no substituent" (defunctionalization / reset).
_GROUPS: dict[str, FunctionalGroup] = {
    group.name: group
    for group in (
        FunctionalGroup("H", "[H]", 0, "Hydrogen (no substituent / reset)"),
        FunctionalGroup("F", "F", 0, "Fluoro"),
        FunctionalGroup("Cl", "Cl", 0, "Chloro"),
        FunctionalGroup("Br", "Br", 0, "Bromo"),
        FunctionalGroup("NH2", "N", 0, "Amino"),
        FunctionalGroup("OH", "O", 0, "Hydroxyl"),
        FunctionalGroup("CH3", "C", 0, "Methyl"),
        FunctionalGroup("NO2", "[N+](=O)[O-]", 0, "Nitro"),
        FunctionalGroup("COOH", "C(=O)O", 0, "Carboxyl"),
        FunctionalGroup("CN", "C#N", 0, "Nitrile"),
        FunctionalGroup("OCH3", "OC", 0, "Methoxy"),
        FunctionalGroup("acetamido", "NC(C)=O", 0, "Acetylamido (acetamido)"),
    )
}


def available_groups() -> list[str]:
    """Return the sorted list of functional-group names the agent may choose."""
    return sorted(_GROUPS)


def get_group(name: str) -> FunctionalGroup:
    """Look up a :class:`FunctionalGroup` by name.

    Raises
    ------
    KeyError
        If *name* is not in the curated menu; the message lists valid names.
    """
    try:
        return _GROUPS[name]
    except KeyError:
        raise KeyError(
            f"Unknown functional group {name!r}. "
            f"Available groups: {', '.join(available_groups())}"
        ) from None


def group_smiles(name: str) -> str:
    """Return the SMILES fragment for a named functional group."""
    return get_group(name).smiles
