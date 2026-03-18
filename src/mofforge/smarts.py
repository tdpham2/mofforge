"""SMARTS-like pattern matching for chemical substructure queries.

Provides a simplified SMARTS-like syntax for specifying query patterns
as strings instead of requiring XYZ file input.

Supported syntax:
    - Element symbols: C, N, O, Zn, etc.
    - Bonds: - (explicit bond separator, optional — consecutive atoms
      are implicitly bonded regardless of whether ``-`` is present)
    - Brackets: [C], [Zn] for explicit atoms
    - Wildcards: [*] or * for any atom
    - Ring closure: C1-C-C-C-C-C-1 (digit after atom marks ring closure point)

Examples:
    "Zn-O-C"        -> linear chain: Zn bonded to O bonded to C
    "[Zn]-[O]-[C]"  -> same as above
    "C1-C-C-C-C-C-1" -> 6-membered carbon ring
    "C-[*]-C"       -> two carbons bonded through any atom
"""

from __future__ import annotations

import logging
import re

import networkx as nx

from mofforge.core.crystal import Crystal
from mofforge.search.search import MatchResult

logger = logging.getLogger("mofforge")

# Regex for tokenizing SMARTS-like strings
_TOKEN_RE = re.compile(
    r"""
    \[([^\]]+)\]   |   # bracketed atom [Zn], [*]
    ([A-Z][a-z]?)  |   # element symbol
    (\*)           |   # wildcard
    (-)            |   # bond
    (\d)               # ring closure digit
    """,
    re.VERBOSE,
)


def parse_smarts(pattern: str) -> nx.Graph:
    """Parse a SMARTS-like pattern string into a NetworkX query graph.

    Args:
        pattern: SMARTS-like pattern string (e.g. "Zn-O-C").

    Returns:
        NetworkX Graph with 'species' node attributes.
        Wildcard atoms get species='*'.

    Raises:
        ValueError: If the pattern is invalid.
    """
    graph = nx.Graph()
    atom_idx = 0
    last_atom: int | None = None
    ring_opens: dict[str, int] = {}  # digit -> atom index of ring opening

    tokens = _TOKEN_RE.findall(pattern)
    if not tokens:
        raise ValueError(f"Could not parse SMARTS pattern: '{pattern}'")

    for bracketed, element, wildcard, bond, ring_digit in tokens:
        if bracketed:
            # Bracketed atom: [Zn], [*], [C]
            species = bracketed.strip()
            if species == "*":
                graph.add_node(atom_idx, species="*")
            else:
                graph.add_node(atom_idx, species=species)
            if last_atom is not None:
                graph.add_edge(last_atom, atom_idx)
            last_atom = atom_idx
            atom_idx += 1

        elif element:
            # Bare element symbol
            graph.add_node(atom_idx, species=element)
            if last_atom is not None:
                graph.add_edge(last_atom, atom_idx)
            last_atom = atom_idx
            atom_idx += 1

        elif wildcard:
            # Wildcard *
            graph.add_node(atom_idx, species="*")
            if last_atom is not None:
                graph.add_edge(last_atom, atom_idx)
            last_atom = atom_idx
            atom_idx += 1

        elif bond:
            # Explicit bond '-' — we just use this as a separator
            # The next atom will be bonded to last_atom
            pass

        elif ring_digit:
            # Ring closure
            if ring_digit in ring_opens:
                # Close the ring
                open_atom = ring_opens.pop(ring_digit)
                if last_atom is not None and last_atom != open_atom:
                    graph.add_edge(last_atom, open_atom)
            else:
                # Open a ring at the last atom
                if last_atom is not None:
                    ring_opens[ring_digit] = last_atom

    if ring_opens:
        raise ValueError(
            f"Unclosed ring(s) in pattern '{pattern}': digits {list(ring_opens.keys())}"
        )

    if graph.number_of_nodes() == 0:
        raise ValueError(f"No atoms found in pattern '{pattern}'")

    logger.debug(
        "Parsed SMARTS pattern '%s' -> %d atoms, %d bonds",
        pattern,
        graph.number_of_nodes(),
        graph.number_of_edges(),
    )

    return graph


def _wildcard_node_match(n1: dict, n2: dict) -> bool:
    """Node match function that supports wildcards."""
    s1 = n1.get("species", "")
    s2 = n2.get("species", "")
    if s1 == "*" or s2 == "*":
        return True
    return s1 == s2


def smarts_search(
    pattern: str,
    parent: Crystal,
    disconnected_component: bool = False,
) -> MatchResult:
    """Search a parent crystal using a SMARTS-like pattern string.

    This is an alternative to loading a query from an XYZ file.
    The pattern is parsed into a query graph and matched against
    the parent's bond graph.

    .. note::
        The returned ``MatchResult`` uses a dummy empty ``query`` Crystal
        and is intended for **inspection only** (e.g. counting matches,
        extracting matched atoms). It **cannot** be passed to
        ``replace_pattern()`` — use a proper XYZ-loaded query for that.

    Args:
        pattern: SMARTS-like pattern string.
        parent: The crystal structure to search in.
        disconnected_component: If True, search for exact isolated matches.

    Returns:
        MatchResult object with isomorphism results.
    """
    from collections import defaultdict

    from networkx.algorithms.isomorphism import GraphMatcher

    if parent.n_bonds == 0:
        raise ValueError(
            "The parent structure must have bonds. "
            "Use infer_bonds(crystal, periodic=True) to create them."
        )

    query_graph = parse_smarts(pattern)
    has_wildcards = any(d.get("species") == "*" for _, d in query_graph.nodes(data=True))

    node_match = (
        _wildcard_node_match
        if has_wildcards
        else lambda n1, n2: n1.get("species") == n2.get("species")
    )

    if not disconnected_component:
        gm = GraphMatcher(parent.bonds, query_graph, node_match=node_match)
        raw_isomorphisms = [{v: k for k, v in m.items()} for m in gm.subgraph_isomorphisms_iter()]
    else:
        raw_isomorphisms = []
        query_n = query_graph.number_of_nodes()
        for comp_nodes in nx.connected_components(parent.bonds):
            if len(comp_nodes) != query_n:
                continue
            comp = parent.bonds.subgraph(comp_nodes)
            gm = GraphMatcher(comp, query_graph, node_match=node_match)
            for m in gm.isomorphisms_iter():
                raw_isomorphisms.append({v: k for k, v in m.items()})

    # Group by location
    location_groups: dict[tuple[int, ...], list[dict[int, int]]] = defaultdict(list)
    for isom in raw_isomorphisms:
        parent_atoms = tuple(sorted(isom.values()))
        location_groups[parent_atoms].append(isom)

    isomorphisms = list(location_groups.values())

    # Create a dummy query Crystal for the MatchResult object
    query_crystal = Crystal.empty(name=f"smarts:{pattern}")

    return MatchResult(
        parent=parent,
        query=query_crystal,
        isomorphisms=isomorphisms,
    )
