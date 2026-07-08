"""Detect functionalizable positions on a MOF linker from its SMILES.

The functionalizable positions are the aromatic C-H hydrogens of the linker:
replacing one of these hydrogens with a functional group is exactly the
post-synthetic modification an agent wants to perform.  Metal-binding groups
(carboxylate carbons, etc.) carry no aromatic hydrogen and are therefore never
selected — the ``[cH]`` SMARTS query inherently excludes them.

Each site is tagged with a *symmetry class* (via RDKit canonical ranking that
ignores atom order but respects graph symmetry).  This disambiguates the two
meanings of an index selection for the agent:

* Two indices in **different** symmetry classes (e.g. the alpha vs beta
  positions of a naphthalene linker) target genuinely different chemical
  environments.
* Multiple indices in the **same** symmetry class on one ring define a
  *substitution pattern* (mono- vs di-substitution, ortho / meta / para).

The agent picks site *indices*; it never sees or authors coordinates.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from mofforge.build.smiles_to_bb import _ensure_rdkit

# Aromatic carbon bearing exactly one hydrogen — the substitutable positions.
_AROMATIC_CH_SMARTS = "[cH]"


@dataclass
class FunctionalizableSite:
    """One aromatic C-H position that can be functionalized.

    Attributes
    ----------
    index:
        Stable, zero-based identifier the agent uses to select this site.
        Indices are ordered by ``(symmetry_class, atom_idx)`` so they are
        deterministic across runs for a given SMILES.
    atom_idx:
        RDKit atom index of the aromatic carbon (internal; not for the agent).
    symmetry_class:
        Integer label shared by all symmetry-equivalent sites.  Equal classes
        mean the positions are chemically indistinguishable.
    element:
        Element of the ring atom (always ``"C"`` for aromatic C-H).
    ring_id:
        Index of the ring (from RDKit's SSSR) this site belongs to; ``-1`` if
        the atom is not in a detected ring.
    neighbors:
        Element symbols of the ring atom's heavy-atom neighbours, for context.
    description:
        Human-readable summary for the agent.
    """

    index: int
    atom_idx: int
    symmetry_class: int
    element: str
    ring_id: int
    neighbors: list[str] = field(default_factory=list)
    description: str = ""


def find_functionalizable_sites(linker_smiles: str) -> list[FunctionalizableSite]:
    """Enumerate the functionalizable aromatic C-H positions of a linker.

    Parameters
    ----------
    linker_smiles:
        SMILES of the organic linker (e.g. obtained from MOFid).  Metal-binding
        carboxylate groups may be present; they simply contain no aromatic C-H
        and so are never returned.

    Returns
    -------
    list[FunctionalizableSite]
        Sites ordered deterministically by ``(symmetry_class, atom_idx)``.  The
        list may be empty (e.g. an aliphatic linker with no aromatic C-H).
    """
    _ensure_rdkit()
    from rdkit import Chem

    mol = Chem.MolFromSmiles(linker_smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {linker_smiles!r}")

    # Canonical ranking with symmetry respected (breakTies=False) gives every
    # symmetry-equivalent atom the same rank -> our symmetry class.
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))

    # Map each ring atom to the ring it belongs to (first match wins).
    ring_info = mol.GetRingInfo()
    atom_ring: dict[int, int] = {}
    for ring_id, atom_ring_atoms in enumerate(ring_info.AtomRings()):
        for a in atom_ring_atoms:
            atom_ring.setdefault(a, ring_id)

    ch_pattern = Chem.MolFromSmarts(_AROMATIC_CH_SMARTS)
    matches = mol.GetSubstructMatches(ch_pattern)
    ch_atoms = sorted({m[0] for m in matches})

    sites: list[FunctionalizableSite] = []
    for atom_idx in ch_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
        sites.append(
            FunctionalizableSite(
                index=-1,  # assigned after sorting
                atom_idx=atom_idx,
                symmetry_class=int(ranks[atom_idx]),
                element=atom.GetSymbol(),
                ring_id=atom_ring.get(atom_idx, -1),
                neighbors=neighbors,
            )
        )

    # Deterministic ordering, then assign stable public indices.
    sites.sort(key=lambda s: (s.symmetry_class, s.atom_idx))

    # Re-label symmetry classes to small contiguous integers (0, 1, 2, ...) in
    # order of first appearance, so the agent sees tidy class labels.
    class_remap: dict[int, int] = {}
    for site in sites:
        if site.symmetry_class not in class_remap:
            class_remap[site.symmetry_class] = len(class_remap)

    n_classes = len(class_remap)
    for i, site in enumerate(sites):
        site.index = i
        site.symmetry_class = class_remap[site.symmetry_class]
        site.description = (
            f"aromatic C-H, symmetry class {site.symmetry_class} "
            f"of {n_classes}, ring {site.ring_id}"
        )

    return sites
