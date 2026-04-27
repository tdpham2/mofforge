"""Adsorbate placement into MOF structures."""

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
    """Result of placing one or more adsorbate molecules into a MOF."""

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
    """Place adsorbate molecule(s) into a MOF crystal structure."""
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
            raise ValueError("No adsorption sites found")
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
            raise ValueError(f"Unknown orient: {orient!r}")

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

    logger.debug(
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


def _filter_by_intermolecular_dist(
    sites: list[AdsorptionSite],
    min_dist: float,
) -> list[AdsorptionSite]:
    """Greedily filter sites so no two are closer than ``min_dist``."""
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
