"""Pattern replacement pipeline for crystal substructures."""

from __future__ import annotations

import logging
import random as pyrandom
from dataclasses import dataclass

import numpy as np

from mofforge.core.bonding import infer_bonds, remove_bonds
from mofforge.core.crystal import Crystal
from mofforge.core.moiety import anchor_indices, fragment, untag_anchor
from mofforge.replace.alignment import Alignment, apply_alignment, get_r2p_alignment
from mofforge.search.search import MatchResult, find_pattern
from mofforge.utils.periodic import is_cross_boundary, wrap_coords

logger = logging.getLogger("mofforge")

# Type alias: mapping from atom index in one context to atom index in another.
# q2r: query atom idx -> replacement atom idx
# q2p: query atom idx -> parent atom idx
# r2p: replacement atom idx -> parent atom idx
AtomMapping = dict[int, int]


@dataclass
class Installation:
    """A prepared replacement ready to be installed into the parent."""

    aligned_replacement: Crystal
    q2p: dict[int, int]
    r2p: dict[int, int]


def _find_query_in_replacement(
    match: MatchResult,
    replacement: Crystal,
) -> dict[int, int]:
    """Find the correspondence between query (unmasked) and replacement atoms."""
    # Count non-masked atoms in query
    r_indices = set(anchor_indices(match.query.species))
    nb_not_masked = sum(1 for i in range(match.query.n_atoms) if i not in r_indices)

    if nb_not_masked == 0:
        return {}

    # Extract unmasked portion of query
    unmasked_indices = [i for i in range(match.query.n_atoms) if i not in r_indices]
    query_unmasked = match.query[unmasked_indices]

    # Untag anchor species for searching
    untagged_species = untag_anchor(query_unmasked.species)
    for i in range(query_unmasked.n_atoms):
        query_unmasked.bonds.nodes[i]["species"] = untagged_species[i]

    # Search for unmasked query in replacement
    q_in_r = find_pattern(query_unmasked, replacement)

    if q_in_r.nb_locations() == 0:
        raise ValueError(
            f"Query atoms not found in replacement '{replacement.name}'."
        )

    # Take the first isomorphism.
    # query_unmasked was created via __getitem__ which renumbers indices to
    # 0..len(unmasked_indices)-1, so the keys should be exactly 0..nb_not_masked-1.
    first_isom = q_in_r.isomorphisms[0][0]
    if set(first_isom.keys()) != set(range(nb_not_masked)):
        raise ValueError(
            f"Unexpected isomorphism keys: expected {{0..{nb_not_masked - 1}}}, "
            f"got {set(first_isom.keys())}. This indicates a bug in index renumbering."
        )
    q2r = {q: first_isom[q] for q in range(nb_not_masked)}
    return q2r


def optimal_replacement(
    match: MatchResult,
    replacement: Crystal,
    q2r: dict[int, int],
    loc_id: int,
    ori_ids: list[int] | None = None,
) -> Installation:
    """Find the best orientation at a location by minimizing alignment error."""
    isomorphisms = match.isomorphisms
    parent = match.parent

    # Handle replace-with-nothing (no alignment needed, just mark atoms for deletion)
    if not q2r:
        q2p = isomorphisms[loc_id][0]
        return Installation(replacement, q2p, r2p={})

    # If ori_ids is None, evaluate all orientations to find the optimal one
    if ori_ids is None:
        n_ori = len(isomorphisms[loc_id])
        ori_ids = list(range(n_ori))

    best_alignment = Alignment(np.eye(3), np.zeros(3), np.zeros(3), float("inf"))
    best_ori = 0
    best_r2p: dict[int, int] = {}

    for ori_id in ori_ids:
        q2p = isomorphisms[loc_id][ori_id]
        # Build r2p: for each (q, r) pair in q2r, map r -> q2p[q]
        r2p = {r: q2p[q] for q, r in q2r.items()}

        try:
            test_alignment = get_r2p_alignment(replacement, parent, r2p, q2p)
        except (ValueError, AssertionError, np.linalg.LinAlgError):
            logger.debug("Alignment failed for loc=%d ori=%d, skipping", loc_id, ori_id)
            continue

        if test_alignment.error < best_alignment.error:
            best_alignment = test_alignment
            best_ori = ori_id
            best_r2p = r2p

    # Check if any alignment succeeded
    if best_alignment.error == float("inf"):
        raise ValueError(f"Alignment failed for location {loc_id}.")

    # Apply the best alignment
    aligned_rep = apply_alignment(replacement, parent, best_alignment)

    if aligned_rep.n_bonds != replacement.n_bonds:
        raise ValueError(
            f"Bond count changed during alignment: expected {replacement.n_bonds}, "
            f"got {aligned_rep.n_bonds}."
        )

    return Installation(
        aligned_replacement=aligned_rep,
        q2p=isomorphisms[loc_id][best_ori],
        r2p=best_r2p,
    )


def install_replacements(
    parent: Crystal,
    installations: list[Installation],
    name: str,
) -> Crystal:
    """Install replacement fragments into the parent crystal.

    For each installation:
        1. Add aligned replacement atoms to the child crystal.
        2. Reconstruct bonds between replacement and parent neighbor atoms.
        3. Accumulate obsolete (matched query) atoms for deletion.

    After all installations, delete obsolete atoms.
    """
    # Start with a copy of the parent (clear symmetry for combination)
    child = parent.copy()
    child.name = name

    obsolete_atoms: list[int] = []

    for installation in installations:
        rep = installation.aligned_replacement
        q2p = installation.q2p
        r2p = installation.r2p

        # Add replacement atoms to child
        if rep.n_atoms > 0:
            offset = child.n_atoms
            child = child + rep

            # Pre-compute set of query-matched parent atoms for O(1) lookup
            q2p_parent_atoms = set(q2p.values())

            # Reconstruct bonds between replacement atoms and parent neighbors
            for r_idx, p_idx in r2p.items():
                # Find neighbors of p_idx in the original parent
                if not parent.bonds.has_node(p_idx):
                    continue
                for p_nbr in parent.bonds.neighbors(p_idx):
                    # Skip if neighbor is part of the query match (will be deleted)
                    if p_nbr in q2p_parent_atoms:
                        continue

                    # Connect p_nbr to the corresponding replacement atom in child
                    r_in_child = offset + r_idx

                    if child.bonds.has_node(p_nbr) and child.bonds.has_node(r_in_child):
                        # Copy edge attributes from parent, but mark cross_boundary as unknown
                        edge_data = parent.bonds.get_edge_data(p_idx, p_nbr)
                        if edge_data is None:
                            edge_data = {}
                        child.bonds.add_edge(
                            p_nbr,
                            r_in_child,
                            distance=edge_data.get("distance", 0.0),
                            cross_boundary=None,  # will be reassessed later
                        )

        # Accumulate atoms to delete
        obsolete_atoms.extend(q2p.values())

    # Delete obsolete atoms
    unique_obsolete = sorted(set(obsolete_atoms))
    keep = [i for i in range(child.n_atoms) if i not in unique_obsolete]
    child = child[keep]

    return child


def effect_replacements(
    match: MatchResult,
    replacement: Crystal,
    configs: list[tuple[int, int]],
    name: str,
) -> Crystal:
    """Orchestrate replacement for given (location, orientation) configurations."""
    # Find query-to-replacement correspondence
    if replacement.n_atoms > 0:
        q2r = _find_query_in_replacement(match, replacement)
    else:
        q2r = {}

    # Compute optimal installations
    installations = [
        optimal_replacement(
            match,
            replacement,
            q2r,
            loc_id,
            None if ori_id is None else [ori_id],
        )
        for loc_id, ori_id in configs
    ]

    # Install into parent
    child = install_replacements(match.parent, installations, name)

    # Fix cross_boundary edge attributes for edges with None/unknown
    _fix_cross_boundary_attrs(child)

    return child


def _fix_cross_boundary_attrs(crystal: Crystal) -> None:
    """Reassess cross_boundary attributes for edges where it's unknown."""
    for u, v, data in crystal.bonds.edges(data=True):
        if data.get("cross_boundary") is None:
            dist = data.get("distance", 0.0)
            if dist > 0 and crystal.n_atoms > 0:
                cross_pb = is_cross_boundary(
                    crystal.frac_coords[u],
                    crystal.frac_coords[v],
                    crystal.lattice,
                    dist,
                )
                crystal.bonds[u][v]["cross_boundary"] = cross_pb
            else:
                crystal.bonds[u][v]["cross_boundary"] = False


def replace_pattern(
    match: MatchResult,
    replacement: Crystal | None,
    *,
    random: bool = False,
    nb_loc: int = 0,
    loc: list[int] | None = None,
    ori: list[int] | None = None,
    name: str = "new_xtal",
    verbose: bool = False,
    remove_duplicates: bool = False,
    periodic_boundaries: bool = True,
    reinfer_bonds: bool = False,
    wrap: bool = True,
    auto_supercell: bool = False,
) -> Crystal:
    """Replace substructures of match.parent matching match.query with replacement.

    Default behavior replaces at all locations with optimal orientation.

    Replacement modes:
        - Default: all locations, optimal orientation at each.
        - ``random=True``: all locations, random orientation at each.
        - ``nb_loc=N``: N random locations.
        - ``loc=[...]``: specific locations.
        - ``loc=[...], ori=[...]``: specific location+orientation pairs.
    """
    # Handle None replacement
    if replacement is None:
        replacement = fragment(None)

    if loc is None:
        loc = []
    if ori is None:
        ori = []

    n_locations = match.nb_locations()
    ori_counts = match.nb_ori_at_loc()

    # Filter out locations with zero orientations (defensive)
    valid_locs = [i for i in range(n_locations) if ori_counts[i] > 0]

    # Determine replacement mode
    if nb_loc == 0 and not loc and not ori:
        # Replace at all locations
        loc = list(valid_locs)
        nb_loc = len(loc)
        if random:
            ori = [pyrandom.randint(0, ori_counts[i] - 1) for i in loc]
            if verbose:
                logger.info("Replacing: random ori @ all %d loc", nb_loc)
        else:
            ori = [None] * nb_loc  # None = evaluate all, pick optimal
            if verbose:
                logger.info("Replacing: optimal ori @ all %d loc", nb_loc)

    elif nb_loc > 0 and not loc and not ori:
        # Random locations
        loc = pyrandom.sample(valid_locs, min(nb_loc, len(valid_locs)))
        if random:
            ori = [pyrandom.randint(0, ori_counts[i] - 1) for i in loc]
            if verbose:
                logger.info("Replacing: random ori @ %d random loc", nb_loc)
        else:
            ori = [None] * len(loc)
            if verbose:
                logger.info("Replacing: optimal ori @ %d random loc", nb_loc)

    elif loc and ori:
        # Specific locations and orientations
        if len(loc) != len(ori):
            raise ValueError("One orientation per location required")
        # Validate location and orientation indices
        for l, o in zip(loc, ori):
            if l < 0 or l >= n_locations:
                raise ValueError(f"Location index {l} out of range (0..{n_locations - 1}).")
            if o is not None and (o < 0 or o >= ori_counts[l]):
                raise ValueError(
                    f"Orientation index {o} out of range for location {l} (0..{ori_counts[l] - 1})."
                )
        nb_loc = len(loc)
        if verbose:
            logger.info("Replacing: loc=%s ori=%s", loc, ori)

    elif loc:
        # Specific locations, auto orientation
        # Validate location indices
        for l in loc:
            if l < 0 or l >= n_locations:
                raise ValueError(f"Location index {l} out of range (0..{n_locations - 1}).")
            if ori_counts[l] == 0:
                raise ValueError(f"Location {l} has no orientations available.")
        nb_loc = len(loc)
        if random:
            ori = [pyrandom.randint(0, ori_counts[i] - 1) for i in loc]
            if verbose:
                logger.info("Replacing: random ori @ loc=%s", loc)
        else:
            ori = [None] * nb_loc
            if verbose:
                logger.info("Replacing: optimal ori @ loc=%s", loc)

    # Note: pymatgen stores partial charges on Species objects, not on
    # PeriodicSite.  Charge transfer during replacement is not supported.

    # Generate configuration tuples
    configs = [(loc[i], ori[i]) for i in range(len(loc))]

    # Process replacements
    child = effect_replacements(match, replacement, configs, name)

    # Remove duplicate atoms if requested
    if remove_duplicates:
        child = _remove_duplicates(child, periodic_boundaries)

    if wrap:
        # Check if replacement is too large for the unit cell
        if np.any(np.abs(child.frac_coords) > 2.0):
            if auto_supercell:
                child = _handle_supercell(match, replacement, configs, name)
            else:
                raise ValueError("Replacement too large for unit cell.")

        # Wrap coordinates
        wrapped_frac = wrap_coords(child.frac_coords)
        child.set_frac_coords(wrapped_frac)

        # Reassess cross_boundary for all edges after wrapping
        _fix_cross_boundary_attrs(child)

    if reinfer_bonds:
        child = remove_bonds(child)
        child = infer_bonds(child, periodic_boundaries)

    return child


def _remove_duplicates(
    crystal: Crystal,
    periodic: bool = True,
    tolerance: float = 0.01,
) -> Crystal:
    """Remove duplicate atoms (same species, overlapping positions)."""
    n = crystal.n_atoms
    if n == 0:
        return crystal

    removed: set[int] = set()
    species = crystal.species

    if periodic:
        # Use pymatgen's neighbor search for efficient periodic distance queries
        all_neighbors = crystal.structure.get_all_neighbors(tolerance)
        for i, neighbors in enumerate(all_neighbors):
            if i in removed:
                continue
            for nbr in neighbors:
                j = nbr.index
                if j <= i or j in removed:
                    continue
                if species[i] == species[j]:
                    removed.add(j)
    else:
        # Use scipy KD-tree for efficient non-periodic neighbor search
        from scipy.spatial import cKDTree

        cart = crystal.cart_coords
        tree = cKDTree(cart)
        pairs = tree.query_pairs(tolerance)
        for i, j in pairs:
            if i in removed or j in removed:
                continue
            if species[i] == species[j]:
                removed.add(max(i, j))

    if removed:
        logger.debug("Removed %d duplicate atoms", len(removed))

    keep = [i for i in range(n) if i not in removed]
    return crystal[keep]


def _handle_supercell(
    match: MatchResult,
    replacement: Crystal,
    configs: list[tuple[int, int]],
    name: str,
) -> Crystal:
    """Handle replacement that's too large by creating a supercell."""
    logger.warning(
        "Replacement fragment too large for unit cell. Creating 2x2x2 supercell "
        "(auto_supercell=True). WARNING: original location/orientation selections "
        "cannot be preserved in the supercell — %d locations will be replaced "
        "with optimal orientation instead.",
        len(configs),
    )

    # Create 2x2x2 supercell of parent — rebuild Crystal from scratch so that
    # _species_labels, bonds, and structure stay in sync.
    super_struct = match.parent.structure.copy()
    super_struct.make_supercell([2, 2, 2])
    parent_super = Crystal.from_structure(super_struct, name=match.parent.name + "_super")

    # Re-infer bonds on supercell
    parent_super = infer_bonds(parent_super, periodic=True)

    # Re-search on supercell
    new_match = find_pattern(match.query, parent_super)

    # Re-derive configs: the original configs referenced locations/orientations
    # in the original match.  The supercell has different (and more) locations,
    # so we select the same number of locations with optimal orientation.
    # NOTE: orientation 0 is used because the supercell's location ordering
    # differs from the original; the original selections cannot be mapped.
    nb_loc = len(configs)
    n_super_locations = new_match.nb_locations()
    if nb_loc > n_super_locations:
        logger.warning(
            "Requested %d locations but supercell only has %d; replacing at %d locations.",
            nb_loc,
            n_super_locations,
            n_super_locations,
        )
    nb_loc = min(nb_loc, n_super_locations)
    new_configs = [(loc_id, None) for loc_id in range(nb_loc)]

    # Re-do replacement without auto_supercell to avoid recursion
    child = effect_replacements(new_match, replacement, new_configs, name)

    # Wrap
    wrapped_frac = wrap_coords(child.frac_coords)
    child.set_frac_coords(wrapped_frac)

    return child


def swap(parent: Crystal, query: Crystal, replacement: Crystal | None, **kwargs) -> Crystal:
    """Convenience function: find and replace in one call.

    Equivalent to::

        match = find_pattern(query, parent)
        child = replace_pattern(match, replacement, **kwargs)
    """
    match = find_pattern(query, parent)
    return replace_pattern(match, replacement, **kwargs)
