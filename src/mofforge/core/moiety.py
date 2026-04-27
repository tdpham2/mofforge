"""Fragment loading from XYZ files with anchor-atom ('!' suffix) handling."""

from __future__ import annotations

import logging
from pathlib import Path

from mofforge.core.bonding import BondingRule, infer_bonds, tagged_bonding_rules
from mofforge.core.crystal import Crystal
from mofforge.io.xyz import read_xyz
from mofforge.utils.config import config

logger = logging.getLogger("mofforge")


def anchor_indices(species: list[str], r_tag: str | None = None) -> list[int]:
    """Return indices of anchor atoms (species ending with the anchor tag)."""
    if r_tag is None:
        r_tag = config.r_tag

    indices = []
    for idx, label in enumerate(species):
        if label.endswith(r_tag):
            # Verify it's a proper tag (element + single !)
            base = label[: -len(r_tag)]
            if len(base) > 0:
                indices.append(idx)
    return indices


def untag_anchor(species: list[str], r_tag: str | None = None) -> list[str]:
    """Remove anchor tags from species labels ('H!' -> 'H', 'C!' -> 'C')."""
    if r_tag is None:
        r_tag = config.r_tag

    result = []
    for label in species:
        if label.endswith(r_tag):
            result.append(label[: -len(r_tag)])
        else:
            result.append(label)
    return result


def subtract_anchor(crystal: Crystal) -> Crystal:
    """Return a copy of the crystal with anchor atoms removed."""
    r_indices = set(anchor_indices(crystal.species))
    keep = [i for i in range(crystal.n_atoms) if i not in r_indices]
    return crystal[keep]


def fragment(
    name: str | None,
    fragment_path: str | Path | None = None,
    bonding_rules: list[BondingRule] | None = None,
    presort: bool = True,
) -> Crystal:
    """Load a molecular fragment from an XYZ file.

    Bonds are inferred (non-periodic). Atoms are sorted by bond degree
    (highest first), with anchor atoms moved to the end.
    """
    if name is None:
        return Crystal.empty(name="nothing")

    # Resolve file path
    if fragment_path is None:
        fragment_path = config.moiety_path
    if fragment_path is None:
        raise ValueError("No fragment path configured.")
    filepath = Path(fragment_path) / name

    # Read XYZ file
    species_list, cart_coords = read_xyz(filepath)

    if len(species_list) == 0:
        return Crystal.empty(name=name)

    # Create Crystal with large cubic box
    xtal = Crystal.from_xyz(species_list, cart_coords, name=name)

    # Infer bonds (non-periodic for molecular fragments)
    if bonding_rules is None:
        bonding_rules = tagged_bonding_rules()
    xtal = infer_bonds(xtal, periodic=False, bonding_rules=bonding_rules)

    # Identify anchor atoms
    r_indices = set(anchor_indices(xtal.species))
    non_r_indices = [i for i in range(xtal.n_atoms) if i not in r_indices]

    # Sort non-anchor atoms by bond degree (descending) for search efficiency
    if presort and len(non_r_indices) > 0:
        degrees = [xtal.bonds.degree(i) for i in non_r_indices]
        sorted_pairs = sorted(zip(degrees, non_r_indices, strict=True), reverse=True)
        sorted_non_r = [idx for _, idx in sorted_pairs]
    else:
        sorted_non_r = non_r_indices

    # Final ordering: sorted non-anchor atoms, then anchor atoms at the end
    r_list = sorted(r_indices)
    final_order = sorted_non_r + r_list

    # Reorder the crystal
    reordered = xtal[final_order]

    # Re-infer bonds on the reordered crystal
    reordered = infer_bonds(reordered, periodic=False, bonding_rules=bonding_rules)
    # Restore the name
    reordered.name = name

    logger.debug(
        "Loaded fragment '%s': %d atoms (%d anchor)",
        name,
        reordered.n_atoms,
        len(r_indices),
    )

    return reordered
