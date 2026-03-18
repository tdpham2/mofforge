"""Tests for bond inference."""


from tests.conftest import CRYSTAL_DIR, MOIETY_DIR


class TestBondingRules:
    """Tests for bonding rule generation."""

    def test_default_rules(self):
        """Default rules should cover common element pairs."""
        from mofforge.core.bonding import default_bonding_rules

        rules = default_bonding_rules()
        assert len(rules) > 0
        # Check C-C rule exists
        cc_rules = [r for r in rules if r.matches("C", "C")]
        assert len(cc_rules) > 0
        # C-C distance should be around 1.52 + padding
        assert cc_rules[0].max_dist > 1.4
        assert cc_rules[0].max_dist < 2.0

    def test_tagged_rules(self):
        """Tagged rules should include !-suffixed species."""
        from mofforge.core.bonding import tagged_bonding_rules

        rules = tagged_bonding_rules()
        # Should have rules for H!
        h_bang_rules = [r for r in rules if "H!" in (r.species_i, r.species_j)]
        assert len(h_bang_rules) > 0

    def test_rule_matching(self):
        """BondingRule.matches should be symmetric."""
        from mofforge.core.bonding import BondingRule

        rule = BondingRule("C", "H", 1.2)
        assert rule.matches("C", "H")
        assert rule.matches("H", "C")
        assert not rule.matches("C", "O")


class TestBondInference:
    """Tests for bond inference on crystals."""

    def test_infer_bonds_moiety(self):
        """Infer bonds on a molecular fragment (non-periodic)."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.io.xyz import read_xyz

        species, coords = read_xyz(MOIETY_DIR / "p-phenylene.xyz")
        xtal = Crystal.from_xyz(species, coords, name="p-phenylene")
        xtal = infer_bonds(xtal, periodic=False)

        # p-phenylene (C6H4 with 2 extra C stubs) should have bonds
        assert xtal.n_bonds > 0
        # All C atoms should be bonded
        for i in range(xtal.n_atoms):
            if xtal.species[i] == "C":
                assert xtal.bonds.degree(i) > 0

    def test_infer_bonds_periodic(self):
        """Infer bonds on a periodic crystal structure."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal

        xtal = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        xtal = infer_bonds(xtal, periodic=True)
        assert xtal.n_bonds > 0

    def test_remove_bonds(self):
        """Remove all bonds from a crystal."""
        from mofforge.core.bonding import infer_bonds, remove_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.io.xyz import read_xyz

        species, coords = read_xyz(MOIETY_DIR / "p-phenylene.xyz")
        xtal = Crystal.from_xyz(species, coords)
        xtal = infer_bonds(xtal, periodic=False)
        assert xtal.n_bonds > 0

        xtal = remove_bonds(xtal)
        assert xtal.n_bonds == 0
        assert xtal.n_atoms > 0  # atoms preserved

    def test_drop_cross_pb_bonds(self):
        """Drop cross-boundary bonds."""
        from mofforge.core.bonding import drop_cross_pb_bonds, infer_bonds
        from mofforge.core.crystal import Crystal

        xtal = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        xtal = infer_bonds(xtal, periodic=True)

        bonds_no_pb = drop_cross_pb_bonds(xtal.bonds)
        # Should have fewer or equal bonds
        assert bonds_no_pb.number_of_edges() <= xtal.n_bonds
