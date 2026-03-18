"""Tests for the Crystal data structure."""

import numpy as np
import pytest

from tests.conftest import CRYSTAL_DIR, MOIETY_DIR


class TestCrystal:
    """Tests for Crystal class."""

    def test_from_cif(self):
        """Load a Crystal from a CIF file."""
        from mofforge.core.crystal import Crystal

        xtal = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        assert xtal.n_atoms > 0
        assert xtal.name == "IRMOF-1"

    def test_from_xyz(self):
        """Create a Crystal from species and Cartesian coordinates."""
        from mofforge.core.crystal import Crystal

        species = ["C", "H", "O"]
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        xtal = Crystal.from_xyz(species, coords, name="test")
        assert xtal.n_atoms == 3
        assert xtal.species == ["C", "H", "O"]
        assert xtal.name == "test"

    def test_empty_crystal(self):
        """Create an empty Crystal."""
        from mofforge.core.crystal import Crystal

        xtal = Crystal.empty()
        assert xtal.n_atoms == 0
        assert xtal.n_bonds == 0
        assert xtal.name == "empty"

    def test_indexing(self):
        """Extract sub-crystal by indices."""
        from mofforge.core.crystal import Crystal

        species = ["C", "H", "O", "N"]
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        xtal = Crystal.from_xyz(species, coords)
        sub = xtal[[0, 2]]
        assert sub.n_atoms == 2
        assert sub.species == ["C", "O"]

    def test_addition(self):
        """Combine two crystals."""
        from mofforge.core.crystal import Crystal

        xtal1 = Crystal.from_xyz(["C"], np.array([[0.0, 0.0, 0.0]]), name="a")
        xtal2 = Crystal.from_xyz(["O"], np.array([[1.0, 0.0, 0.0]]), name="b")
        combined = xtal1 + xtal2
        assert combined.n_atoms == 2
        assert "C" in combined.species
        assert "O" in combined.species

    def test_copy(self):
        """Deep copy should be independent."""
        from mofforge.core.crystal import Crystal

        xtal = Crystal.from_xyz(["C", "O"], np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]))
        xtal2 = xtal.copy()
        assert xtal2.n_atoms == xtal.n_atoms
        assert xtal2 is not xtal

    def test_repr(self):
        """String representation should be informative."""
        from mofforge.core.crystal import Crystal

        xtal = Crystal.from_xyz(["C"], np.array([[0.0, 0.0, 0.0]]), name="test")
        assert "test" in repr(xtal)
        assert "n_atoms=1" in repr(xtal)

    def test_duplicate_indices_raises(self):
        """__getitem__ with duplicate indices should raise ValueError."""
        from mofforge.core.crystal import Crystal

        species = ["C", "H", "O"]
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        xtal = Crystal.from_xyz(species, coords)

        with pytest.raises(ValueError, match="Duplicate atom indices"):
            xtal[[0, 0, 1]]

    def test_set_frac_coords_wrong_shape_raises(self):
        """set_frac_coords with wrong shape should raise ValueError."""
        from mofforge.core.crystal import Crystal

        xtal = Crystal.from_xyz(["C", "O"], np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]))
        with pytest.raises(ValueError, match="Expected coordinates shape"):
            xtal.set_frac_coords(np.array([[0.0, 0.0, 0.0]]))

    def test_addition_preserves_provenance(self):
        """__add__ should preserve the left operand's provenance."""
        from mofforge.core.crystal import Crystal
        from mofforge.provenance import Provenance

        prov = Provenance(parent="test", operation="mock")
        xtal1 = Crystal.from_xyz(["C"], np.array([[0.0, 0.0, 0.0]]), name="a")
        xtal1.provenance = prov
        xtal2 = Crystal.from_xyz(["O"], np.array([[1.0, 0.0, 0.0]]), name="b")

        combined = xtal1 + xtal2
        assert combined.provenance is not None
        assert combined.provenance.parent == "test"

    def test_contains(self):
        """Crystal.__contains__ should find a substructure."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.core.moiety import fragment

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)
        query = fragment("p-phenylene.xyz", fragment_path=MOIETY_DIR)

        assert query in parent

    def test_contains_no_match(self):
        """Crystal.__contains__ should return False for non-present fragments."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.core.moiety import fragment

        # Use UiO-66 (has bonds) and search for a fragment that doesn't exist in it
        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)
        # ADC (acetylenedicarboxylate) is a different linker, not in IRMOF-1
        query = fragment("ADC.xyz", fragment_path=MOIETY_DIR)

        assert query not in parent

    def test_find(self):
        """Crystal.find() should return a MatchResult."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.core.moiety import fragment

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)
        query = fragment("p-phenylene.xyz", fragment_path=MOIETY_DIR)

        result = parent.find(query)
        assert result.nb_locations() > 0
