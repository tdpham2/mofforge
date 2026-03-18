"""Bond inference using covalent radii-based distance rules."""

from __future__ import annotations

import functools
import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

import networkx as nx
import numpy as np

from mofforge.utils.config import config
from mofforge.utils.periodic import is_cross_boundary

if TYPE_CHECKING:
    from mofforge.core.crystal import Crystal

logger = logging.getLogger("mofforge")


@dataclass
class BondingRule:
    """A rule for determining if two atoms are bonded.

    Attributes:
        species_i: First atom species (e.g. 'C').
        species_j: Second atom species (e.g. 'O').
        max_dist: Maximum bonding distance in Angstroms.
    """

    species_i: str
    species_j: str
    max_dist: float

    def matches(self, sp_a: str, sp_b: str) -> bool:
        """Check if this rule applies to a pair of species."""
        return (self.species_i == sp_a and self.species_j == sp_b) or (
            self.species_i == sp_b and self.species_j == sp_a
        )


@functools.lru_cache(maxsize=4)
def default_bonding_rules(bond_pad: float | None = None) -> tuple[BondingRule, ...]:
    """Generate bonding rules from all pairs of elements with known covalent radii.

    The maximum bonding distance is the sum of covalent radii plus a padding value
    (default 0.25 A).

    Results are cached (keyed on ``bond_pad``) to avoid regenerating ~4,000
    rules on every call.

    Args:
        bond_pad: Bond padding in Angstroms. If None, uses ``config.bond_pad``.

    Returns:
        Tuple of BondingRule objects (tuple for hashability/caching).
    """
    from mofforge.utils.config import COVALENT_RADII

    if bond_pad is None:
        bond_pad = config.bond_pad

    rules = []
    elements = sorted(COVALENT_RADII.keys())
    for i, el_i in enumerate(elements):
        for el_j in elements[i:]:
            max_dist = COVALENT_RADII[el_i] + COVALENT_RADII[el_j] + bond_pad
            rules.append(BondingRule(el_i, el_j, max_dist))
    return tuple(rules)


def tagged_bonding_rules(
    rules: list[BondingRule] | None = None,
    r_tag: str | None = None,
) -> list[BondingRule]:
    """Generate additional bonding rules for R-group tagged species.

    For each rule (A, B, dist), generates:
        (A!, B, dist), (B!, A, dist), (A!, B!, dist)

    This ensures that tagged atoms (e.g. H!) bond the same way as their
    untagged counterparts (e.g. H).

    Args:
        rules: Base rules; defaults to default_bonding_rules().
        r_tag: The R-group tag character; defaults to config.r_tag.

    Returns:
        Combined list of original + tagged rules.
    """
    if rules is None:
        rules = list(default_bonding_rules(bond_pad=config.bond_pad))
    if r_tag is None:
        r_tag = config.r_tag

    new_rules = list(rules)
    for rule in rules:
        si, sj = rule.species_i, rule.species_j
        new_rules.append(BondingRule(f"{si}{r_tag}", sj, rule.max_dist))
        new_rules.append(BondingRule(f"{sj}{r_tag}", si, rule.max_dist))
        new_rules.append(BondingRule(f"{si}{r_tag}", f"{sj}{r_tag}", rule.max_dist))
    return new_rules


def _build_rule_lookup(rules: list[BondingRule]) -> dict[tuple[str, str], float]:
    """Build a dictionary lookup from bonding rules for O(1) access.

    Keys are normalized (sorted) species pairs so that (A, B) and (B, A)
    map to the same entry.

    Args:
        rules: List of BondingRule objects.

    Returns:
        Dict mapping (species_i, species_j) -> max_dist.
    """
    lookup: dict[tuple[str, str], float] = {}
    for rule in rules:
        key = (rule.species_i, rule.species_j)
        rev_key = (rule.species_j, rule.species_i)
        # Keep the largest max_dist if there are duplicate keys
        if key not in lookup or rule.max_dist > lookup[key]:
            lookup[key] = rule.max_dist
        if rev_key not in lookup or rule.max_dist > lookup[rev_key]:
            lookup[rev_key] = rule.max_dist
    return lookup


def _get_max_bond_dist(
    species_i: str,
    species_j: str,
    rule_lookup: dict[tuple[str, str], float] | None = None,
) -> float | None:
    """Look up the maximum bonding distance for a species pair.

    Falls back to covalent radii sum + padding if no explicit rule matches.

    Args:
        species_i: First species.
        species_j: Second species.
        rule_lookup: Optional dict from _build_rule_lookup() for O(1) access.

    Returns:
        Maximum bonding distance, or None if species are unknown.
    """
    if rule_lookup is not None:
        dist = rule_lookup.get((species_i, species_j))
        if dist is not None:
            return dist

    # Fallback: compute from covalent radii
    try:
        return config.max_bond_distance(species_i, species_j)
    except ValueError:
        return None


def infer_bonds(
    crystal: Crystal,
    periodic: bool = True,
    bonding_rules: list[BondingRule] | None = None,
) -> Crystal:
    """Infer bonds for a Crystal based on interatomic distances.

    Uses covalent radii-based bonding rules. For periodic structures,
    also detects bonds across periodic boundaries.

    Args:
        crystal: The Crystal to infer bonds for.
        periodic: If True, consider periodic boundary conditions.
        bonding_rules: Optional custom bonding rules. If None, uses
            tagged default rules (including R-group species).

    Returns:
        A new Crystal with the inferred bond graph.
    """

    xtal = crystal.copy()

    if bonding_rules is None:
        bonding_rules = tagged_bonding_rules()

    # Build O(1) lookup dict from rules
    rule_lookup = _build_rule_lookup(bonding_rules)

    bonds = nx.Graph()
    n = xtal.n_atoms
    species_list = xtal.species

    # Add all nodes with species attributes
    for i in range(n):
        bonds.add_node(i, species=species_list[i])

    if n == 0:
        xtal.bonds = bonds
        return xtal

    # Compute maximum possible bonding distance for cutoff
    max_cutoff = max((rule.max_dist for rule in bonding_rules), default=0.0)
    # Add some margin
    max_cutoff += 0.5

    structure = xtal.structure
    lattice = xtal.lattice

    if periodic:
        # Use pymatgen's neighbor search with PBC
        all_neighbors = structure.get_all_neighbors(max_cutoff)
        for i, neighbors in enumerate(all_neighbors):
            sp_i = species_list[i]
            for neighbor in neighbors:
                j = neighbor.index
                if j <= i:
                    continue  # avoid duplicate edges
                dist = neighbor.nn_distance
                sp_j = species_list[j]

                max_dist = _get_max_bond_dist(sp_i, sp_j, rule_lookup)
                if max_dist is not None and dist <= max_dist:
                    cross_pb = is_cross_boundary(
                        xtal.frac_coords[i],
                        xtal.frac_coords[j],
                        lattice,
                        dist,
                    )
                    bonds.add_edge(i, j, distance=dist, cross_boundary=cross_pb)
    else:
        # Non-periodic: direct pairwise distance check
        cart = xtal.cart_coords
        for i in range(n):
            sp_i = species_list[i]
            for j in range(i + 1, n):
                sp_j = species_list[j]
                dist = float(np.linalg.norm(cart[i] - cart[j]))

                max_dist = _get_max_bond_dist(sp_i, sp_j, rule_lookup)
                if max_dist is not None and dist <= max_dist:
                    bonds.add_edge(i, j, distance=dist, cross_boundary=False)

    xtal.bonds = bonds
    logger.debug(
        "Inferred %d bonds for '%s' (periodic=%s)", bonds.number_of_edges(), xtal.name, periodic
    )
    return xtal


def remove_bonds(crystal: Crystal) -> Crystal:
    """Return a Crystal with all bonds removed.

    Node attributes (species) are preserved, but all edges are removed.

    Args:
        crystal: The Crystal to remove bonds from.

    Returns:
        A new Crystal with an empty bond graph.
    """
    xtal = crystal.copy()
    species_list = xtal.species
    new_bonds = nx.Graph()
    for i in range(xtal.n_atoms):
        new_bonds.add_node(i, species=species_list[i])
    xtal.bonds = new_bonds
    return xtal


def drop_cross_pb_bonds(graph: nx.Graph) -> nx.Graph:
    """Return a copy of the bond graph with cross-boundary bonds removed.

    Args:
        graph: A bond graph with 'cross_boundary' edge attributes.

    Returns:
        New graph with only non-cross-boundary edges.
    """
    new_graph = graph.copy()
    edges_to_remove = [
        (u, v) for u, v, d in new_graph.edges(data=True) if d.get("cross_boundary", False)
    ]
    new_graph.remove_edges_from(edges_to_remove)
    return new_graph
