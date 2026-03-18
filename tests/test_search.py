"""Tests for pattern matching."""

from tests.conftest import CRYSTAL_DIR, MOIETY_DIR


class TestPatternMatching:
    """Tests for find_pattern on real crystal structures."""

    def test_phenylene_in_irmof1(self):
        """p-phenylene in IRMOF-1: expect 96 isomorphisms at 24 locations."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.core.moiety import fragment
        from mofforge.search.search import find_pattern

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)

        query = fragment("p-phenylene.xyz", fragment_path=MOIETY_DIR)

        search = find_pattern(query, parent)

        # IRMOF-1 has 24 linker locations, each with 4 orientations
        assert search.nb_locations() == 24
        assert search.nb_isomorphisms() == 96
        # Each location has 4 orientations
        for n_ori in search.nb_ori_at_loc():
            assert n_ori == 4

    def test_phenylene_in_ti_mil_125(self):
        """p-phenylene in Ti-MIL-125: expect 48 isomorphisms at 12 locations."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.core.moiety import fragment
        from mofforge.search.search import find_pattern

        parent = Crystal.from_cif(CRYSTAL_DIR / "Ti-MIL-125.cif")
        parent = infer_bonds(parent, periodic=True)

        query = fragment("p-phenylene.xyz", fragment_path=MOIETY_DIR)

        search = find_pattern(query, parent)

        assert search.nb_locations() == 12
        assert search.nb_isomorphisms() == 48
        for n_ori in search.nb_ori_at_loc():
            assert n_ori == 4

    def test_tagged_phenylene_in_ti_mil_125(self):
        """Tagged 2-!-p-phenylene in Ti-MIL-125: same result as untagged."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.core.moiety import fragment
        from mofforge.search.search import find_pattern

        parent = Crystal.from_cif(CRYSTAL_DIR / "Ti-MIL-125.cif")
        parent = infer_bonds(parent, periodic=True)

        query = fragment("2-!-p-phenylene.xyz", fragment_path=MOIETY_DIR)

        search = find_pattern(query, parent)

        # Anchor tag shouldn't affect search results
        assert search.nb_locations() == 12
        assert search.nb_isomorphisms() == 48

    def test_search_repr(self):
        """MatchResult should have informative string representation."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.core.moiety import fragment
        from mofforge.search.search import find_pattern

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)
        query = fragment("p-phenylene.xyz", fragment_path=MOIETY_DIR)
        search = find_pattern(query, parent)

        s = repr(search)
        assert "96" in s or "hits" in s
        assert "24" in s or "locations" in s

    def test_matched_substructures(self):
        """matched_substructures should return a valid sub-crystal."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.core.moiety import fragment
        from mofforge.search.search import find_pattern

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)
        query = fragment("p-phenylene.xyz", fragment_path=MOIETY_DIR)
        search = find_pattern(query, parent)

        sub = search.matched_substructures()
        # Should have 24 locations * 10 atoms per phenylene = 240 atoms
        # (though some atoms may be shared between locations)
        assert sub.n_atoms > 0

    def test_empty_parent(self):
        """Search on empty parent should return empty results."""
        from mofforge.core.crystal import Crystal
        from mofforge.core.moiety import fragment
        from mofforge.search.search import find_pattern

        parent = Crystal.empty()
        query = fragment("p-phenylene.xyz", fragment_path=MOIETY_DIR)
        search = find_pattern(query, parent)
        assert search.nb_locations() == 0
        assert search.nb_isomorphisms() == 0
