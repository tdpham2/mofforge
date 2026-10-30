"""Optional SMARTS site annotations; no reaction recipes or compatibility inference."""

from __future__ import annotations

from dataclasses import dataclass, field

from mofforge.build.smiles_to_bb import _ensure_rdkit


@dataclass(frozen=True)
class ReactiveGroup:
    """A curated reactive functional group.

    Attributes
    ----------
    site_type:
        Short name the agent uses (e.g. ``"amine"``).
    smarts:
        SMARTS matching the group; match atom 0 is the bonding heavy atom.
    description:
        Human-readable summary.
    """

    site_type: str
    smarts: str
    description: str = ""


# Curated reactive groups.  Match atom 0 is always the atom that forms the new
# inter-monomer bond (the "linker atom" in Polymatic terms).
_GROUPS: tuple[ReactiveGroup, ...] = (
    ReactiveGroup("amine", "[NX3;H2][#6]", "primary amine (-NH2)"),
    ReactiveGroup("aldehyde", "[CX3H1](=O)[#6]", "aldehyde (-CHO)"),
    ReactiveGroup("aryl_halide", "[c][F,Cl,Br,I]", "aryl halide (Ar-X)"),
    ReactiveGroup("boronic_acid", "[#6]B(O)O", "boronic acid (-B(OH)2)"),
    ReactiveGroup("vinyl", "[CX3H1]=[CX3H2]", "terminal vinyl (-CH=CH2)"),
    ReactiveGroup("hydroxyl", "[OX2H][#6]", "hydroxyl (-OH)"),
    ReactiveGroup("nitrile", "[NX1]#[CX2][#6]", "nitrile (-C#N)"),
)

_GROUPS_BY_TYPE: dict[str, ReactiveGroup] = {g.site_type: g for g in _GROUPS}


def available_reactions() -> list[dict[str, str]]:
    """No built-in construction recipes; supply explicit ConnectionRule edits."""
    return []


def available_site_types() -> list[dict[str, str]]:
    """Return the curated reactive site types as JSON-friendly dicts."""
    return [
        {"site_type": g.site_type, "smarts": g.smarts, "description": g.description}
        for g in _GROUPS
    ]


def get_group(site_type: str) -> ReactiveGroup:
    """Return the curated :class:`ReactiveGroup` for *site_type*."""
    try:
        return _GROUPS_BY_TYPE[site_type]
    except KeyError:
        raise ValueError(
            f"Unknown reactive site type {site_type!r}. "
            f"Choose from: {', '.join(sorted(_GROUPS_BY_TYPE))}"
        ) from None


@dataclass
class DetectedSite:
    """A reactive site found on a molecule by SMARTS matching."""

    atom_idx: int
    site_type: str
    anchor_idx: int | None = None
    neighbors: list[str] = field(default_factory=list)


def detect_reactive_sites(
    smiles: str,
    site_types: list[str] | None = None,
) -> list[DetectedSite]:
    """Enumerate reactive sites on a monomer from its SMILES.

    Parameters
    ----------
    smiles:
        SMILES of the monomer.
    site_types:
        Restrict detection to these curated site types.  ``None`` scans all
        curated groups.

    Returns
    -------
    list[DetectedSite]
        Sites ordered deterministically by ``(atom_idx, site_type)``.  Each
        bonding atom is reported once (the first matching group wins), so a
        molecule with two amines yields two amine sites.
    """
    _ensure_rdkit()
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")

    groups = _GROUPS if site_types is None else [get_group(t) for t in site_types]

    seen: dict[int, DetectedSite] = {}
    for group in groups:
        pattern = Chem.MolFromSmarts(group.smarts)
        if pattern is None:  # pragma: no cover - curated SMARTS are valid
            raise ValueError(f"Invalid SMARTS for {group.site_type!r}: {group.smarts!r}")
        for match in mol.GetSubstructMatches(pattern):
            bond_atom = match[0]
            if bond_atom in seen:
                continue  # first matching group wins for a given atom
            atom = mol.GetAtomWithIdx(bond_atom)
            heavy_nbrs = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]
            anchor = heavy_nbrs[0].GetIdx() if heavy_nbrs else None
            seen[bond_atom] = DetectedSite(
                atom_idx=bond_atom,
                site_type=group.site_type,
                anchor_idx=anchor,
                neighbors=[n.GetSymbol() for n in atom.GetNeighbors()],
            )

    return sorted(seen.values(), key=lambda s: (s.atom_idx, s.site_type))
