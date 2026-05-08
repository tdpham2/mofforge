"""Automatic solvent removal from MOF crystal structures."""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, field

import networkx as nx

from mofforge.core.bonding import infer_bonds
from mofforge.core.crystal import Crystal
from mofforge.utils.config import clean_species
from mofforge.vis.colors import METALS

logger = logging.getLogger("mofforge")


@dataclass
class RemovedMolecule:
    """Details of a single removed solvent/guest molecule."""

    atom_indices: list[int]
    formula: str
    n_atoms: int
    contains_metal: bool


@dataclass
class SolventRemovalResult:
    """Result of automatic solvent removal."""

    crystal: Crystal
    removed_molecules: list[RemovedMolecule] = field(default_factory=list)
    n_atoms_original: int = 0
    n_atoms_removed: int = 0
    n_components_removed: int = 0
    n_framework_components: int = 1

    def summary(self) -> str:
        """Return a human-readable summary of the removal."""
        lines = [
            f"Solvent Removal: {self.n_atoms_removed} atoms removed "
            f"({self.n_atoms_original} -> {self.crystal.n_atoms})",
            f"  Framework components: {self.n_framework_components}",
            f"  Molecules removed: {self.n_components_removed}",
        ]
        if self.removed_molecules:
            formula_counts = Counter(m.formula for m in self.removed_molecules)
            lines.append("  Removed species:")
            for formula, count in formula_counts.most_common():
                lines.append(f"    {count}x {formula}")
        return "\n".join(lines)


def _composition_formula(species: list[str]) -> str:
    """Generate a Hill-order chemical formula from a list of species labels.

    Hill order: C first, then H, then remaining elements alphabetically.
    """
    counts = Counter(clean_species(s) for s in species)

    elements = []
    # C first
    if "C" in counts:
        elements.append(("C", counts.pop("C")))
    # H second
    if "H" in counts:
        elements.append(("H", counts.pop("H")))
    # Remaining alphabetically
    for elem in sorted(counts):
        elements.append((elem, counts[elem]))

    parts = []
    for elem, count in elements:
        parts.append(f"{elem}{count}" if count > 1 else elem)
    return "".join(parts)


def _has_metal(species: list[str]) -> bool:
    """Check if any species in the list is a metal element."""
    return any(clean_species(s) in METALS for s in species)


def remove_solvent(
    crystal: Crystal,
    min_atoms: int = 1,
    keep_metal_containing: bool = False,
    n_framework_components: int | None = None,
    periodic: bool = True,
) -> SolventRemovalResult:
    """Automatically identify and remove uncoordinated solvent molecules.

    Analyses connected components in the periodic bond graph. The largest
    component(s) are considered framework; smaller isolated components are
    solvent/guest molecules and are removed.

    Parameters
    ----------
    crystal : Crystal
        Input crystal structure (bonds are inferred if absent).
    min_atoms : int
        Components with at least this many atoms are kept as framework
        regardless of size ranking. Default ``1`` means only the standard
        algorithm applies.
    keep_metal_containing : bool
        If ``True``, do not remove components that contain metal atoms.
    n_framework_components : int | None
        Number of largest components to keep as framework.  If ``None``
        (default), auto-detect: keep the largest component plus any
        component with at least 50 % of its size (handles interpenetrated
        MOFs).
    periodic : bool
        Whether to use periodic boundary bonds for component detection.
        Almost always ``True`` for MOF structures.

    Returns
    -------
    SolventRemovalResult
        The desolvated crystal and metadata about removed molecules.
    """
    n_original = crystal.n_atoms
    if n_original == 0:
        return SolventRemovalResult(
            crystal=crystal.copy(),
            n_atoms_original=0,
        )

    xtal = crystal.copy()

    # Infer bonds if not already present
    if xtal.n_bonds == 0:
        xtal = infer_bonds(xtal, periodic=periodic)

    # Find connected components, sorted largest-first
    components = sorted(nx.connected_components(xtal.bonds), key=len, reverse=True)

    if len(components) <= 1:
        return SolventRemovalResult(
            crystal=xtal,
            n_atoms_original=n_original,
            n_framework_components=len(components),
        )

    # Classify framework vs solvent
    largest_size = len(components[0])

    if n_framework_components is not None:
        # Explicit: keep the N largest
        n_keep = min(n_framework_components, len(components))
        framework_components = components[:n_keep]
        solvent_components = components[n_keep:]
    else:
        # Auto-detect: keep largest + anything >= 50% of its size
        framework_components = [components[0]]
        solvent_components = []
        for comp in components[1:]:
            if len(comp) >= largest_size * 0.5:
                framework_components.append(comp)
            else:
                solvent_components.append(comp)

    # Apply filters to rescue components from removal
    rescued = []
    remaining_solvent = []
    for comp in solvent_components:
        comp_species = [xtal.species[i] for i in comp]
        if keep_metal_containing and _has_metal(comp_species):
            rescued.append(comp)
            continue
        if len(comp) >= min_atoms and min_atoms > 1:
            rescued.append(comp)
            continue
        remaining_solvent.append(comp)

    framework_components.extend(rescued)
    solvent_components = remaining_solvent

    # Nothing to remove
    if not solvent_components:
        return SolventRemovalResult(
            crystal=xtal,
            n_atoms_original=n_original,
            n_framework_components=len(framework_components),
        )

    # Build list of removed molecules
    removed_molecules = []
    for comp in solvent_components:
        indices = sorted(comp)
        comp_species = [xtal.species[i] for i in indices]
        removed_molecules.append(
            RemovedMolecule(
                atom_indices=indices,
                formula=_composition_formula(comp_species),
                n_atoms=len(indices),
                contains_metal=_has_metal(comp_species),
            )
        )

    # Collect framework atom indices and extract sub-crystal
    keep_indices = sorted(
        idx for comp in framework_components for idx in comp
    )
    result_crystal = xtal[keep_indices]
    result_crystal.name = f"desolvated_{xtal.name}"

    n_removed = n_original - len(keep_indices)

    logger.info(
        "Removed %d solvent molecule(s) (%d atoms) from '%s'",
        len(removed_molecules),
        n_removed,
        xtal.name,
    )

    return SolventRemovalResult(
        crystal=result_crystal,
        removed_molecules=removed_molecules,
        n_atoms_original=n_original,
        n_atoms_removed=n_removed,
        n_components_removed=len(removed_molecules),
        n_framework_components=len(framework_components),
    )
