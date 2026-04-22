"""Tests for SVD alignment (Procrustes)."""

import numpy as np
import pytest

from mofforge.core.crystal import Crystal
from mofforge.replace.alignment import Alignment, apply_alignment, get_r2p_alignment


def _make_crystal_from_cart(species, coords, name="test"):
    """Helper: create a Crystal from Cartesian coordinates."""
    return Crystal.from_xyz(species, np.array(coords, dtype=np.float64), name=name)


class TestAlignment:
    """Tests for SVD alignment."""

    def test_identity_alignment(self):
        """Aligning a fragment to itself should give near-zero error."""
        species = ["C", "N", "O", "H", "H"]
        coords = [
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [0.0, 1.5, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.5],
        ]
        replacement = _make_crystal_from_cart(species, coords, "rep")
        parent = _make_crystal_from_cart(species, coords, "parent")

        r2p = {0: 0, 1: 1, 2: 2}
        q2p = {0: 0, 1: 1, 2: 2}

        alignment = get_r2p_alignment(replacement, parent, r2p, q2p)
        assert alignment.error < 0.01
        # Rotation should be close to identity
        np.testing.assert_allclose(alignment.rotation, np.eye(3), atol=0.1)

    def test_translated_alignment(self):
        """Aligning a translated fragment should recover the translation."""
        species = ["C", "N", "O", "H"]
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.5, 0.0, 0.0],
                [0.0, 1.5, 0.0],
                [1.0, 1.0, 1.0],
            ]
        )
        shift = np.array([10.0, 20.0, 30.0])
        shifted_coords = coords + shift

        replacement = _make_crystal_from_cart(species, coords.tolist(), "rep")
        parent = _make_crystal_from_cart(species, shifted_coords.tolist(), "parent")

        r2p = {0: 0, 1: 1, 2: 2}
        q2p = {0: 0, 1: 1, 2: 2}

        alignment = get_r2p_alignment(replacement, parent, r2p, q2p)
        assert alignment.error < 0.1

    def test_two_atom_alignment_works(self):
        """Alignment with 2 atoms should succeed (translation + partial rotation)."""
        species = ["C", "N"]
        coords = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        replacement = _make_crystal_from_cart(species, coords, "rep")
        parent = _make_crystal_from_cart(species, coords, "parent")

        r2p = {0: 0, 1: 1}
        q2p = {0: 0, 1: 1}

        alignment = get_r2p_alignment(replacement, parent, r2p, q2p)
        assert alignment.error < 0.1

    def test_too_few_alignment_points_raises(self):
        """Alignment with < 2 alignment points should raise ValueError."""
        species = ["C", "N", "O"]
        coords = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        replacement = _make_crystal_from_cart(species, coords, "rep")
        parent = _make_crystal_from_cart(species, coords, "parent")

        r2p = {0: 0}  # only 1 alignment point
        q2p = {0: 0, 1: 1, 2: 2}

        with pytest.raises(ValueError, match="at least 2 alignment points"):
            get_r2p_alignment(replacement, parent, r2p, q2p)

    def test_apply_alignment_preserves_atom_count(self):
        """apply_alignment should preserve the number of atoms."""
        species = ["C", "N", "O", "H"]
        coords = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.0, 1.5, 0.0], [1.0, 1.0, 1.0]]
        replacement = _make_crystal_from_cart(species, coords, "rep")
        parent = _make_crystal_from_cart(species, coords, "parent")

        alignment = Alignment(
            rotation=np.eye(3),
            shift_pre=np.zeros(3),
            shift_post=np.zeros(3),
            error=0.0,
        )

        result = apply_alignment(replacement, parent, alignment)
        assert result.n_atoms == replacement.n_atoms

    def test_alignment_dataclass(self):
        """Alignment dataclass should store fields correctly."""
        rot = np.eye(3)
        pre = np.array([1.0, 2.0, 3.0])
        post = np.array([4.0, 5.0, 6.0])
        a = Alignment(rotation=rot, shift_pre=pre, shift_post=post, error=0.5)
        assert a.error == 0.5
        np.testing.assert_array_equal(a.rotation, rot)
        np.testing.assert_array_equal(a.shift_pre, pre)
        np.testing.assert_array_equal(a.shift_post, post)
