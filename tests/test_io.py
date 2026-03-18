"""Tests for file I/O modules."""

import tempfile

import numpy as np
import pytest

from tests.conftest import CRYSTAL_DIR, MOIETY_DIR


class TestXYZIO:
    """Tests for XYZ file reading and writing."""

    def test_read_xyz_basic(self):
        """Read a simple moiety XYZ file."""
        from mofforge.io.xyz import read_xyz

        species, coords = read_xyz(MOIETY_DIR / "p-phenylene.xyz")
        assert len(species) == 10
        assert coords.shape == (10, 3)
        # p-phenylene has 6 C and 4 H
        assert species.count("C") == 6
        assert species.count("H") == 4

    def test_read_xyz_with_r_tags(self):
        """Read XYZ file with !-tagged R-group atoms."""
        from mofforge.io.xyz import read_xyz

        species, _coords = read_xyz(MOIETY_DIR / "2-!-p-phenylene.xyz")
        assert len(species) == 10
        # Should have one H! atom
        assert "H!" in species

    def test_read_xyz_multiple_tags(self):
        """Read XYZ file with multiple !-tagged species."""
        from mofforge.io.xyz import read_xyz

        species, _coords = read_xyz(MOIETY_DIR / "BDC.xyz")
        tagged = [s for s in species if "!" in s]
        assert len(tagged) > 0

    def test_write_and_read_roundtrip(self):
        """Write then read XYZ should preserve data."""
        from mofforge.io.xyz import read_xyz, write_xyz

        species = ["C", "H", "O"]
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])

        with tempfile.NamedTemporaryFile(suffix=".xyz", delete=False) as f:
            write_xyz(species, coords, f.name)
            sp_read, coords_read = read_xyz(f.name)

        assert sp_read == species
        np.testing.assert_allclose(coords_read, coords, atol=1e-4)

    def test_read_nonexistent_file(self):
        """Reading a nonexistent file should raise FileNotFoundError."""
        from mofforge.io.xyz import read_xyz

        with pytest.raises(FileNotFoundError):
            read_xyz("/nonexistent/file.xyz")


class TestCIFIO:
    """Tests for CIF file reading and writing."""

    def test_read_cif(self):
        """Read a CIF file."""
        from mofforge.io.cif import read_cif

        structure = read_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        assert len(structure) > 0

    def test_read_nonexistent_file(self):
        """Reading a nonexistent CIF should raise FileNotFoundError."""
        from mofforge.io.cif import read_cif

        with pytest.raises(FileNotFoundError):
            read_cif("/nonexistent/file.cif")
