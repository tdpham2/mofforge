"""Subgraph isomorphism using NetworkX's VF2 algorithm."""

from __future__ import annotations

import logging

import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

logger = logging.getLogger("mofforge")


def _node_match(n1: dict, n2: dict) -> bool:
    """Node compatibility check: species must match."""
    return n1.get("species") == n2.get("species")


def find_subgraph_isomorphisms(
    query_graph: nx.Graph,
    parent_graph: nx.Graph,
    disconnected_component: bool = False,
) -> list[dict[int, int]]:
    """Find all subgraph isomorphisms of a query graph within a parent graph."""
    if query_graph.number_of_nodes() == 0:
        return []

    if parent_graph.number_of_nodes() == 0:
        return []

    isomorphisms: list[dict[int, int]] = []

    if not disconnected_component:
        # Standard subgraph isomorphism search
        # NetworkX GraphMatcher(G1, G2): finds subgraphs of G1 isomorphic to G2
        gm = GraphMatcher(parent_graph, query_graph, node_match=_node_match)
        for mapping in gm.subgraph_isomorphisms_iter():
            # mapping: {parent_node -> query_node}
            # We want: {query_node -> parent_node}
            inverted = {v: k for k, v in mapping.items()}
            isomorphisms.append(inverted)
    else:
        # Disconnected component mode: find exact isomorphisms on
        # connected components that have the same number of nodes as query
        query_n = query_graph.number_of_nodes()
        for comp_nodes in nx.connected_components(parent_graph):
            if len(comp_nodes) != query_n:
                continue
            comp = parent_graph.subgraph(comp_nodes)
            gm = GraphMatcher(comp, query_graph, node_match=_node_match)
            for mapping in gm.isomorphisms_iter():
                inverted = {v: k for k, v in mapping.items()}
                isomorphisms.append(inverted)

    logger.debug(
        "Found %d subgraph isomorphisms (disconnected_component=%s)",
        len(isomorphisms),
        disconnected_component,
    )

    return isomorphisms
