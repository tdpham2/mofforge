"""Tests for SMARTS-like pattern matching."""

import pytest

from tests.conftest import CRYSTAL_DIR


class TestSMARTSParsing:
    """Tests for parse_smarts."""

    def test_simple_chain(self):
        """Parse a simple chain pattern."""
        from mofforge.smarts import parse_smarts

        g = parse_smarts("C-O-C")
        assert g.number_of_nodes() == 3
        assert g.number_of_edges() == 2
        assert g.nodes[0]["species"] == "C"
        assert g.nodes[1]["species"] == "O"
        assert g.nodes[2]["species"] == "C"

    def test_bracketed_atoms(self):
        """Parse bracketed atom notation."""
        from mofforge.smarts import parse_smarts

        g = parse_smarts("[Zn]-[O]-[C]")
        assert g.number_of_nodes() == 3
        assert g.nodes[0]["species"] == "Zn"
        assert g.nodes[1]["species"] == "O"
        assert g.nodes[2]["species"] == "C"

    def test_wildcard(self):
        """Parse wildcard atom."""
        from mofforge.smarts import parse_smarts

        g = parse_smarts("C-[*]-O")
        assert g.number_of_nodes() == 3
        assert g.nodes[1]["species"] == "*"

    def test_ring(self):
        """Parse ring closure notation."""
        from mofforge.smarts import parse_smarts

        g = parse_smarts("C1-C-C-C-C-C-1")
        assert g.number_of_nodes() == 6
        assert g.number_of_edges() == 6  # 5 chain + 1 ring closure
        # Verify ring: node 0 should be connected to node 5
        assert g.has_edge(0, 5)

    def test_empty_pattern_raises(self):
        """Empty or invalid pattern should raise ValueError."""
        from mofforge.smarts import parse_smarts

        with pytest.raises(ValueError):
            parse_smarts("")

    def test_unclosed_ring_raises(self):
        """Unclosed ring should raise ValueError."""
        from mofforge.smarts import parse_smarts

        with pytest.raises(ValueError, match="Unclosed ring"):
            parse_smarts("C1-C-C")

    def test_single_atom(self):
        """Single atom pattern."""
        from mofforge.smarts import parse_smarts

        g = parse_smarts("Zn")
        assert g.number_of_nodes() == 1
        assert g.nodes[0]["species"] == "Zn"


class TestSMARTSSearch:
    """Tests for smarts_search on real structures."""

    def test_smarts_search_basic(self):
        """Search for Zn-O pattern in IRMOF-1."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.smarts import smarts_search

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)

        result = smarts_search("[Zn]-[O]", parent)
        # IRMOF-1 has Zn-O bonds
        assert result.nb_isomorphisms() > 0
