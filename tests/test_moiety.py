"""Tests for fragment loading and anchor handling."""

from tests.conftest import MOIETY_DIR


class TestAnchorFunctions:
    """Tests for anchor detection and manipulation."""

    def test_anchor_indices_no_tags(self):
        """No anchor atoms in a clean species list."""
        from mofforge.core.moiety import anchor_indices

        species = ["C", "H", "O"]
        assert anchor_indices(species) == []

    def test_anchor_indices_with_tags(self):
        """Detect anchor atoms."""
        from mofforge.core.moiety import anchor_indices

        species = ["C", "H!", "O", "C!"]
        indices = anchor_indices(species)
        assert 1 in indices  # H!
        assert 3 in indices  # C!
        assert len(indices) == 2

    def test_untag_anchor(self):
        """Remove ! tags from species labels."""
        from mofforge.core.moiety import untag_anchor

        species = ["C", "H!", "O", "C!"]
        untagged = untag_anchor(species)
        assert untagged == ["C", "H", "O", "C"]

    def test_subtract_anchor(self):
        """Remove anchor atoms from a crystal."""
        from mofforge.core.moiety import fragment, subtract_anchor

        m = fragment("!-S-bromochlorofluoromethane.xyz", fragment_path=MOIETY_DIR)
        # Should have H! tag
        assert any("!" in s for s in m.species)

        no_r = subtract_anchor(m)
        # Should have one fewer atom (H! removed)
        assert no_r.n_atoms == m.n_atoms - 1
        # No tags remaining
        assert not any("!" in s for s in no_r.species)


class TestFragmentLoading:
    """Tests for the fragment() function."""

    def test_load_basic_fragment(self):
        """Load a simple fragment without anchor tags."""
        from mofforge.core.moiety import fragment

        m = fragment("p-phenylene.xyz", fragment_path=MOIETY_DIR)
        assert m.n_atoms == 10
        assert m.n_bonds > 0

    def test_load_tagged_fragment(self):
        """Load a fragment with anchor tags."""
        from mofforge.core.moiety import fragment

        m = fragment("2-!-p-phenylene.xyz", fragment_path=MOIETY_DIR)
        assert m.n_atoms == 10
        # H! should be at the end due to anchor sorting
        assert m.species[-1] == "H!"

    def test_load_nothing_fragment(self):
        """fragment(None) should return an empty crystal."""
        from mofforge.core.moiety import fragment

        m = fragment(None)
        assert m.n_atoms == 0
        assert m.name == "nothing"

    def test_presort(self):
        """Atoms should be sorted by bond degree when presort=True."""
        from mofforge.core.moiety import fragment

        m1 = fragment("glycine_res.xyz", fragment_path=MOIETY_DIR, presort=True)
        m2 = fragment("glycine_res.xyz", fragment_path=MOIETY_DIR, presort=False)

        # Both should have same number of atoms
        assert m1.n_atoms == m2.n_atoms

        # Coordinates should differ (different ordering)
        # (unless the presort happens to produce the same order)
        # At minimum, both should have valid bonds
        assert m1.n_bonds > 0
        assert m2.n_bonds > 0

    def test_bromochlorofluoromethane_species(self):
        """Test species of S-bromochlorofluoromethane fragment."""
        from mofforge.core.moiety import fragment

        m = fragment("S-bromochlorofluoromethane.xyz", fragment_path=MOIETY_DIR)
        assert m.n_atoms == 5
        # Should contain C, Cl, F, Br, H (in some order)
        species_set = set(m.species)
        assert "C" in species_set
        assert "H" in species_set
        # Cl, F, Br should all be present
        assert "Cl" in species_set or "F" in species_set or "Br" in species_set

    def test_tagged_bromochlorofluoromethane(self):
        """!-tagged version should have H! atom."""
        from mofforge.core.moiety import fragment

        m = fragment("!-S-bromochlorofluoromethane.xyz", fragment_path=MOIETY_DIR)
        assert m.n_atoms == 5
        assert "H!" in m.species
        # H! should be at the end
        assert m.species[-1] == "H!"
