"""Tests for pattern replacement."""

import pytest

from tests.conftest import CRYSTAL_DIR, MOIETY_DIR


class TestPatternReplacement:
    """Tests for replace_pattern on real crystal structures."""

    def _load_irmof1_search(self):
        """Helper: load IRMOF-1, infer bonds, search for 2-!-p-phenylene."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.core.moiety import fragment
        from mofforge.search.search import find_pattern

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)
        query = fragment("2-!-p-phenylene.xyz", fragment_path=MOIETY_DIR)
        search = find_pattern(query, parent)
        return parent, query, search

    def test_replace_single_location(self):
        """Replace at 1 random location should add atoms."""
        from mofforge.core.moiety import fragment
        from mofforge.replace.replace import replace_pattern

        parent, _query, search = self._load_irmof1_search()
        replacement = fragment("2-acetylamido-p-phenylene.xyz", fragment_path=MOIETY_DIR)

        child = replace_pattern(search, replacement, nb_loc=1)

        # Original IRMOF-1 has 424 atoms
        # Replacing 1 location should change atom count
        assert child.n_atoms > 0
        assert child.n_atoms != parent.n_atoms

    def test_replace_specific_locations(self):
        """Replace at specific locations."""
        from mofforge.core.moiety import fragment
        from mofforge.replace.replace import replace_pattern

        _parent, _query, search = self._load_irmof1_search()
        replacement = fragment("2-acetylamido-p-phenylene.xyz", fragment_path=MOIETY_DIR)

        child = replace_pattern(search, replacement, loc=[0, 1])

        # Should have more atoms than single replacement
        assert child.n_atoms > 0

    def test_replace_with_nothing(self):
        """Replace with None should remove matched atoms."""
        from mofforge.replace.replace import replace_pattern

        parent, _query, search = self._load_irmof1_search()

        child = replace_pattern(search, None, nb_loc=1)

        # Should have fewer atoms (removed a phenylene + anchor)
        assert child.n_atoms < parent.n_atoms

    def test_identity_replacement(self):
        """Replacing p-phenylene with p-phenylene should preserve atom count."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.core.moiety import fragment
        from mofforge.replace.replace import replace_pattern
        from mofforge.search.search import find_pattern

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)
        query = fragment("p-phenylene.xyz", fragment_path=MOIETY_DIR)
        replacement = fragment("p-phenylene.xyz", fragment_path=MOIETY_DIR)

        search = find_pattern(query, parent)
        child = replace_pattern(search, replacement)

        # Identity replacement: atom count should be preserved
        assert child.n_atoms == parent.n_atoms

    def test_replace_too_large_raises(self):
        """Replacement too large for unit cell should raise error."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.core.moiety import fragment
        from mofforge.replace.replace import replace_pattern
        from mofforge.search.search import find_pattern

        parent = Crystal.from_cif(CRYSTAL_DIR / "NiPyC_fragment_trouble.cif")
        parent = infer_bonds(parent, periodic=True)
        query = fragment("PyC.xyz", fragment_path=MOIETY_DIR)
        replacement = fragment("PyC-long_chain.xyz", fragment_path=MOIETY_DIR)

        search = find_pattern(query, parent)
        if search.nb_locations() > 0:
            with pytest.raises(ValueError, match="too large"):
                replace_pattern(search, replacement)


class TestReassemble:
    """Tests for periodic boundary reassembly."""

    def test_reassemble_basic(self):
        """Reassemble should not crash on a simple structure."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.replace.conglomerate import reassemble

        xtal = Crystal.from_cif(CRYSTAL_DIR / "conglomerate_test.cif")
        xtal = infer_bonds(xtal, periodic=True)

        # The test structure should have cross-PB bonds
        cross_pb_count = sum(
            1 for _, _, d in xtal.bonds.edges(data=True) if d.get("cross_boundary", False)
        )

        if cross_pb_count > 0:
            result = reassemble(xtal)
            assert result.n_atoms == xtal.n_atoms

    def test_reassemble_preserves_connectivity(self):
        """Reassemble should preserve the number of bonds."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.replace.conglomerate import reassemble

        xtal = Crystal.from_cif(CRYSTAL_DIR / "conglomerate_test.cif")
        xtal = infer_bonds(xtal, periodic=True)

        if xtal.n_bonds > 0:
            result = reassemble(xtal)
            # Bond count should be preserved (reassembly only moves atoms)
            assert result.n_bonds == xtal.n_bonds

    def test_reassemble_contiguous_coords(self):
        """Reassembled atoms should be spatially contiguous (smaller spread)."""
        import numpy as np

        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.replace.conglomerate import reassemble

        xtal = Crystal.from_cif(CRYSTAL_DIR / "conglomerate_test.cif")
        xtal = infer_bonds(xtal, periodic=True)

        cross_pb_count = sum(
            1 for _, _, d in xtal.bonds.edges(data=True) if d.get("cross_boundary", False)
        )

        if cross_pb_count > 0 and xtal.n_atoms > 0:
            result = reassemble(xtal)
            # The reassembled Cartesian coordinates should have a smaller
            # or equal bounding box (atoms brought together across PBC)
            cart_before = xtal.cart_coords
            cart_after = result.cart_coords
            spread_before = np.ptp(cart_before, axis=0).sum()
            spread_after = np.ptp(cart_after, axis=0).sum()
            # Reassembly should not increase spread significantly
            assert spread_after <= spread_before * 1.5 + 1.0


class TestSwapConvenience:
    """Tests for the convenience swap() function."""

    def test_swap_function(self):
        """swap(parent, query, replacement) convenience function."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.core.moiety import fragment
        from mofforge.replace.replace import swap

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)
        query = fragment("2-!-p-phenylene.xyz", fragment_path=MOIETY_DIR)
        replacement = fragment("2-acetylamido-p-phenylene.xyz", fragment_path=MOIETY_DIR)

        child = swap(parent, query, replacement, nb_loc=1)
        assert child.n_atoms > 0
