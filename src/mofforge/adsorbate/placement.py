"""Adsorbate placement into MOF structures.

Places adsorbate molecules at identified adsorption sites within a MOF
crystal structure.  Supports built-in molecule geometries (from
:mod:`~mofforge.adsorbate.molecules`) or custom fragments (from XYZ files
via :func:`~mofforge.core.moiety.fragment`).

The main entry point is :func:`place_adsorbate`, which:
    1. Resolves the adsorbate geometry (built-in name or Crystal fragment).
    2. Identifies adsorption sites (or uses explicitly provided sites).
    3. Applies a random 3D rotation to the adsorbate (unless single-atom).
    4. Translates the adsorbate to each target site.
    5. Combines with the host framework via ``Crystal.__add__``.
    6. Optionally validates for steric clashes.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.spatial.transform import Rotation

from mofforge.adsorbate.molecules import get_molecule
from mofforge.adsorbate.sites import AdsorptionSite, find_adsorption_sites
from mofforge.core.crystal import Crystal
from mofforge.provenance import Provenance
from mofforge.validation import validate_structure

logger = logging.getLogger("mofforge")


@dataclass
class AdsorbatePlacement:
    """Result of placing one or more adsorbate molecules into a MOF.

    Attributes:
        crystal: The combined MOF + adsorbate(s) Crystal.
        sites: The adsorption sites where adsorbates were placed.
        adsorbate_indices: For each placed adsorbate, the list of atom
            indices in the combined crystal.  ``adsorbate_indices[k]``
            corresponds to ``sites[k]``.
        n_adsorbates: Number of adsorbate molecules placed.
        adsorbate_name: Name/formula of the adsorbate.
        clashes: Number of steric clashes detected (0 if validation skipped).
    """

    crystal: Crystal
    sites: list[AdsorptionSite]
    adsorbate_indices: list[list[int]]
    n_adsorbates: int
    adsorbate_name: str
    clashes: int = 0


def place_adsorbate(
    crystal: Crystal,
    adsorbate: Crystal | str,
    site: AdsorptionSite | None = None,
    sites: list[AdsorptionSite] | None = None,
    strategy: str = "void",
    n_adsorbates: int = 1,
    orient: str = "random",
    min_intermolecular_dist: float = 3.0,
    validate: bool = True,
    clash_tolerance: float = 0.5,
    random_seed: int | None = None,
    name: str | None = None,
    **site_kwargs: Any,
) -> AdsorbatePlacement:
    """Place adsorbate molecule(s) into a MOF crystal structure.

    The adsorbate can be specified as:

    - A string name of a built-in molecule (e.g. ``"CO2"``, ``"H2O"``,
      ``"methane"``).  See :func:`~mofforge.adsorbate.molecules.available_molecules`.
    - A :class:`~mofforge.core.crystal.Crystal` object (e.g. loaded from
      an XYZ file via :func:`~mofforge.core.moiety.fragment`).

    If no sites are provided, they are auto-detected using
    :func:`~mofforge.adsorbate.sites.find_adsorption_sites`.

    Args:
        crystal: The host MOF crystal.
        adsorbate: Adsorbate molecule name (str) or Crystal fragment.
        site: A single adsorption site to place at.
        sites: Multiple adsorption sites.  Mutually exclusive with ``site``.
        strategy: Site-finding strategy if no sites given
            (``"void"``, ``"open_metal"``, ``"both"``).
        n_adsorbates: Number of adsorbate molecules to place (used only
            when auto-detecting sites; ignored if sites are explicit).
        orient: Orientation strategy:
            - ``"random"``: Apply a random 3D rotation.
            - ``"fixed"``: Keep the molecule in its default orientation.
        min_intermolecular_dist: Minimum distance between adsorbate
            centers when placing multiple molecules (Angstroms).
        validate: If True, run steric clash validation after placement.
        clash_tolerance: Tolerance for clash detection (Angstroms).
        random_seed: Random seed for reproducible orientations.
        name: Name for the resulting crystal. Defaults to
            ``"{crystal.name}_with_{adsorbate}"``.
        **site_kwargs: Extra keyword arguments passed to
            :func:`~mofforge.adsorbate.sites.find_adsorption_sites`
            (e.g. ``min_distance``, ``grid_spacing``).

    Returns:
        :class:`AdsorbatePlacement` with the combined structure and metadata.

    Raises:
        ValueError: If both ``site`` and ``sites`` are given, if the
            adsorbate name is not recognized, or if no suitable sites
            are found.
    """
    if site is not None and sites is not None:
        raise ValueError("Provide either 'site' or 'sites', not both.")

    rng = np.random.default_rng(random_seed)

    # --- Resolve adsorbate geometry ---
    if isinstance(adsorbate, str):
        ads_name = adsorbate
        species, coords = get_molecule(adsorbate)
    elif isinstance(adsorbate, Crystal):
        ads_name = adsorbate.name
        species = adsorbate.species
        coords = adsorbate.cart_coords
        # Center the adsorbate at the origin
        coords = coords - coords.mean(axis=0)
    else:
        raise TypeError(f"'adsorbate' must be a str or Crystal, got {type(adsorbate).__name__}.")

    if name is None:
        name = f"{crystal.name}_with_{ads_name}"

    is_single_atom = len(species) == 1

    # --- Resolve placement sites ---
    if site is not None:
        target_sites = [site]
    elif sites is not None:
        target_sites = list(sites)
    else:
        # Auto-detect sites
        target_sites = find_adsorption_sites(
            crystal,
            strategy=strategy,
            max_sites=n_adsorbates,
            **site_kwargs,
        )
        if not target_sites:
            raise ValueError(
                f"No adsorption sites found in '{crystal.name}' with "
                f"strategy='{strategy}'. Try lowering min_distance or "
                f"using a different strategy."
            )
        # If we found more sites than requested, take the top n_adsorbates
        if len(target_sites) > n_adsorbates:
            target_sites = target_sites[:n_adsorbates]

    # --- Filter sites by intermolecular distance ---
    if len(target_sites) > 1:
        target_sites = _filter_by_intermolecular_dist(target_sites, min_intermolecular_dist)

    # --- Place adsorbates ---
    combined = crystal.copy()
    all_indices: list[list[int]] = []

    for i, target in enumerate(target_sites):
        # Rotate adsorbate (skip for single atoms)
        if is_single_atom or orient == "fixed":
            rotated_coords = coords.copy()
        elif orient == "random":
            rot = Rotation.random(random_state=rng.integers(0, 2**31))
            rotated_coords = rot.apply(coords)
        else:
            raise ValueError(f"Unknown orient='{orient}'. Use 'random' or 'fixed'.")

        # Translate to site position (Cartesian)
        placed_coords = rotated_coords + target.cart_coords

        # Build adsorbate Crystal with the host's lattice
        ads_crystal = Crystal.from_xyz(
            species=list(species),
            cart_coords=placed_coords,
            name=f"{ads_name}_{i}",
            lattice=combined.lattice,
        )

        # Track indices
        offset = combined.n_atoms
        indices = list(range(offset, offset + ads_crystal.n_atoms))
        all_indices.append(indices)

        # Combine
        combined = combined + ads_crystal

    combined.name = name
    combined.provenance = Provenance(
        parent=crystal.name,
        operation="add_adsorbate",
        parameters={
            "adsorbate": ads_name,
            "n_adsorbates": len(target_sites),
            "strategy": strategy,
            "orient": orient,
        },
    )

    # --- Optional validation ---
    n_clashes = 0
    if validate:
        report = validate_structure(
            combined,
            check_clashes=True,
            check_bonds=False,
            check_coordination=False,
            clash_tolerance=clash_tolerance,
        )
        n_clashes = len(report.steric_clashes)
        if n_clashes > 0:
            logger.warning(
                "Placed %d adsorbate(s) in '%s' but detected %d steric clash(es). "
                "Consider adjusting min_distance or site selection.",
                len(target_sites),
                crystal.name,
                n_clashes,
            )

    logger.info(
        "Placed %d %s molecule(s) in '%s'",
        len(target_sites),
        ads_name,
        crystal.name,
    )

    return AdsorbatePlacement(
        crystal=combined,
        sites=target_sites,
        adsorbate_indices=all_indices,
        n_adsorbates=len(target_sites),
        adsorbate_name=ads_name,
        clashes=n_clashes,
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _filter_by_intermolecular_dist(
    sites: list[AdsorptionSite],
    min_dist: float,
) -> list[AdsorptionSite]:
    """Greedily filter sites so no two are closer than ``min_dist``.

    Keeps sites in order of priority (first site always kept).

    Args:
        sites: Candidate sites (pre-sorted by priority).
        min_dist: Minimum allowed distance between site centers (A).

    Returns:
        Filtered list of sites.
    """
    if min_dist <= 0 or len(sites) <= 1:
        return sites

    kept: list[AdsorptionSite] = [sites[0]]
    for candidate in sites[1:]:
        too_close = False
        for existing in kept:
            d = np.linalg.norm(candidate.cart_coords - existing.cart_coords)
            if d < min_dist:
                too_close = True
                break
        if not too_close:
            kept.append(candidate)
    return kept
