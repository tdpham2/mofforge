"""Pattern matching API: MatchResult class and find_pattern function."""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field

from mofforge.core.crystal import Crystal
from mofforge.core.moiety import untag_anchor
from mofforge.search.isomorphism import find_subgraph_isomorphisms

logger = logging.getLogger("mofforge")


@dataclass
class MatchResult:
    """Stores the results of a pattern search.

    Attributes:
        parent: The crystal structure that was searched.
        query: The fragment that was searched for.
        isomorphisms: Nested list of query-to-parent atom correspondences,
            grouped by location then by orientation.
            ``isomorphisms[i_loc][i_ori]`` is a dict mapping
            query atom indices to parent atom indices.
    """

    parent: Crystal
    query: Crystal
    isomorphisms: list[list[dict[int, int]]] = field(default_factory=list)

    def nb_isomorphisms(self) -> int:
        """Total number of isomorphisms found."""
        return sum(self.nb_ori_at_loc())

    def nb_locations(self) -> int:
        """Number of unique locations (distinct sets of parent atoms)."""
        return len(self.isomorphisms)

    def nb_ori_at_loc(self) -> list[int]:
        """Number of orientations (isomorphisms) at each location."""
        return [len(loc) for loc in self.isomorphisms]

    def matched_substructures(self) -> Crystal:
        """Return a Crystal consisting of parent atoms involved in matches.

        Takes the first orientation at each location.
        """
        all_indices = []
        for loc in self.isomorphisms:
            if loc:
                all_indices.extend(loc[0].values())
        unique_indices = sorted(set(all_indices))
        return self.parent[unique_indices]

    def __repr__(self) -> str:
        return (
            f"{self.query.name} in {self.parent.name}: "
            f"{self.nb_isomorphisms()} hits in {self.nb_locations()} locations"
        )


def find_pattern(
    query: Crystal,
    parent: Crystal,
    disconnected_component: bool = False,
) -> MatchResult:
    """Search for a molecular pattern within a crystal.

    Matches are made on the basis of atomic species and chemical bonding
    networks, including bonds across unit cell periodic boundaries.

    The search fragment may contain anchor-tagged atoms (species ending
    with '!'). These tags are stripped for searching so that, e.g., 'H!'
    matches 'H' in the parent.

    Args:
        query: The fragment to search for.
        parent: The crystal structure to search in.
        disconnected_component: If True, search for exact isolated matches
            (e.g. guest molecules). Default False performs substructure search.

    Returns:
        MatchResult object containing all identified isomorphisms grouped
        by location and orientation.

    Raises:
        ValueError: If parent has no bonds.
    """
    if parent.n_atoms == 0:
        return MatchResult(parent=parent, query=query, isomorphisms=[])

    if parent.n_bonds == 0:
        raise ValueError(
            "The parent structure must have bonds. "
            "Use infer_bonds(crystal, periodic=True) to create them."
        )

    # Make a copy of query with anchor tags removed for searching
    search_query = query.copy()
    untagged_species = untag_anchor(search_query.species)

    # Rebuild the query bond graph with untagged species
    import networkx as nx

    query_graph = nx.Graph()
    for i in range(search_query.n_atoms):
        query_graph.add_node(i, species=untagged_species[i])
    for u, v, data in search_query.bonds.edges(data=True):
        query_graph.add_edge(u, v, **data)

    # Run VF2 subgraph isomorphism
    raw_isomorphisms = find_subgraph_isomorphisms(query_graph, parent.bonds, disconnected_component)

    # Group isomorphisms by location (same set of parent atoms = same location)
    location_groups: dict[tuple[int, ...], list[dict[int, int]]] = defaultdict(list)
    for isom in raw_isomorphisms:
        # The location key is the sorted tuple of parent atom indices
        parent_atoms = tuple(sorted(isom.values()))
        location_groups[parent_atoms].append(isom)

    # Convert to nested list format
    isomorphisms = list(location_groups.values())

    logger.info(
        "find_pattern '%s' in '%s': %d isomorphisms at %d locations",
        query.name,
        parent.name,
        sum(len(loc) for loc in isomorphisms),
        len(isomorphisms),
    )

    return MatchResult(parent=parent, query=query, isomorphisms=isomorphisms)
