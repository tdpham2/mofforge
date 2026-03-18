"""Tests for subgraph isomorphism (VF2).

Tests use synthetic labeled graphs with known isomorphism counts.
"""

import networkx as nx


def _make_graph(species: list[str], edges: list[tuple[int, int]]) -> nx.Graph:
    """Helper to create a labeled graph."""
    g = nx.Graph()
    for i, sp in enumerate(species):
        g.add_node(i, species=sp)
    for u, v in edges:
        g.add_edge(u, v)
    return g


class TestSubgraphIsomorphism:
    """Tests for find_subgraph_isomorphisms using VF2."""

    def setup_method(self):
        """Create the test graph (6 nodes, same as Julia test)."""
        self.graph = _make_graph(
            species=["A", "B", "B", "C", "A", "A"],
            edges=[(0, 1), (1, 2), (2, 3), (3, 1), (2, 4), (4, 5)],
        )

    def test_simple_edge_bc(self):
        """Subgraph s1: B-C edge. Should find 2 isomorphisms."""
        from mofforge.search.isomorphism import find_subgraph_isomorphisms

        s1 = _make_graph(["B", "C"], [(0, 1)])
        isoms = find_subgraph_isomorphisms(s1, self.graph)
        assert len(isoms) == 2
        # Verify mappings: query node 0 (B) -> parent B nodes, query node 1 (C) -> parent C node
        parent_sets = [tuple(sorted(m.values())) for m in isoms]
        assert (1, 3) in parent_sets
        assert (2, 3) in parent_sets

    def test_simple_edge_ab(self):
        """Subgraph s2: A-B edge. Should find 2 isomorphisms."""
        from mofforge.search.isomorphism import find_subgraph_isomorphisms

        s2 = _make_graph(["A", "B"], [(0, 1)])
        isoms = find_subgraph_isomorphisms(s2, self.graph)
        assert len(isoms) == 2

    def test_triangle_bbc(self):
        """Subgraph s3: B-B-C triangle. Should find 2 isomorphisms."""
        from mofforge.search.isomorphism import find_subgraph_isomorphisms

        s3 = _make_graph(["B", "B", "C"], [(0, 1), (1, 2), (2, 0)])
        isoms = find_subgraph_isomorphisms(s3, self.graph)
        assert len(isoms) == 2

    def test_triangle_no_match_missing_species(self):
        """Subgraph s4: B-D-C triangle. No D in graph, should find 0."""
        from mofforge.search.isomorphism import find_subgraph_isomorphisms

        s4 = _make_graph(["B", "D", "C"], [(0, 1), (1, 2), (2, 0)])
        isoms = find_subgraph_isomorphisms(s4, self.graph)
        assert len(isoms) == 0

    def test_star_no_match_insufficient_degree(self):
        """Subgraph s5: star with degree 4 center. No node has degree 4, should find 0."""
        from mofforge.search.isomorphism import find_subgraph_isomorphisms

        s5 = _make_graph(
            ["B", "A", "B", "C", "D"],
            [(0, 1), (0, 2), (0, 3), (0, 4)],
        )
        isoms = find_subgraph_isomorphisms(s5, self.graph)
        assert len(isoms) == 0

    def test_chain_aab(self):
        """Subgraph s7: A-A-B chain. Should find 1 isomorphism."""
        from mofforge.search.isomorphism import find_subgraph_isomorphisms

        s7 = _make_graph(["A", "A", "B"], [(0, 1), (1, 2)])
        isoms = find_subgraph_isomorphisms(s7, self.graph)
        assert len(isoms) == 1

    def test_full_graph_self_isomorphism(self):
        """Full graph should have exactly 1 subgraph isomorphism with itself."""
        from mofforge.search.isomorphism import find_subgraph_isomorphisms

        isoms = find_subgraph_isomorphisms(self.graph, self.graph)
        assert len(isoms) == 1

    def test_empty_query(self):
        """Empty query should return no isomorphisms."""
        from mofforge.search.isomorphism import find_subgraph_isomorphisms

        empty = nx.Graph()
        isoms = find_subgraph_isomorphisms(empty, self.graph)
        assert len(isoms) == 0

    def test_disconnected_component_mode(self):
        """Disconnected component mode should match isolated subgraphs exactly."""
        from mofforge.search.isomorphism import find_subgraph_isomorphisms

        # Create a parent with two disconnected components
        parent = _make_graph(
            ["A", "B", "A", "B"],
            [(0, 1), (2, 3)],
        )
        query = _make_graph(["A", "B"], [(0, 1)])

        # Subgraph mode: should find both
        isoms_sub = find_subgraph_isomorphisms(query, parent, disconnected_component=False)
        assert len(isoms_sub) == 2

        # Disconnected component mode: should also find both (each component matches)
        isoms_dc = find_subgraph_isomorphisms(query, parent, disconnected_component=True)
        assert len(isoms_dc) == 2
