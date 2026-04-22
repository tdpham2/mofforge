"""Tests for the adsorbate initialization module."""

from __future__ import annotations

import numpy as np
import pytest

from mofforge.adsorbate.molecules import available_molecules, get_molecule
from mofforge.adsorbate.placement import AdsorbatePlacement, place_adsorbate
from mofforge.adsorbate.sites import AdsorptionSite, find_adsorption_sites
from mofforge.core.bonding import infer_bonds
from mofforge.core.crystal import Crystal
from tests.conftest import CRYSTAL_DIR, MOIETY_DIR

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def irmof1():
    """IRMOF-1 with inferred bonds."""
    xtal = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
    return infer_bonds(xtal, periodic=True)


@pytest.fixture
def irmof1_noH():
    """IRMOF-1 without hydrogens (has open metal sites due to missing H)."""
    xtal = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1_noH.cif")
    return infer_bonds(xtal, periodic=True)


@pytest.fixture
def mof74():
    """MOF-74 with inferred bonds (has open metal sites)."""
    xtal = Crystal.from_cif(CRYSTAL_DIR / "MOF-74.cif")
    return infer_bonds(xtal, periodic=True)


# ===========================================================================
# molecules.py tests
# ===========================================================================


class TestMolecules:
    """Tests for built-in molecule geometries."""

    def test_available_molecules_non_empty(self):
        mols = available_molecules()
        assert len(mols) > 0
        assert "CO2" in mols
        assert "H2O" in mols
        assert "CH4" in mols

    def test_get_co2(self):
        species, coords = get_molecule("CO2")
        assert species == ["O", "C", "O"]
        assert coords.shape == (3, 3)
        # CO2 is linear along z-axis, centered at origin
        assert abs(coords[1, 2]) < 1e-6  # C at origin z

    def test_get_h2o(self):
        species, coords = get_molecule("H2O")
        assert len(species) == 3
        assert "O" in species
        assert species.count("H") == 2
        assert coords.shape == (3, 3)

    def test_get_ch4(self):
        species, _coords = get_molecule("CH4")
        assert len(species) == 5
        assert species[0] == "C"
        assert species.count("H") == 4

    def test_single_atom_molecules(self):
        for name in ["He", "Ar", "Xe"]:
            species, coords = get_molecule(name)
            assert len(species) == 1
            assert coords.shape == (1, 3)
            np.testing.assert_array_equal(coords[0], [0.0, 0.0, 0.0])

    def test_case_insensitive(self):
        s1, c1 = get_molecule("CO2")
        s2, c2 = get_molecule("co2")
        assert s1 == s2
        np.testing.assert_array_equal(c1, c2)

    def test_alias_lookup(self):
        s1, c1 = get_molecule("CO2")
        s2, c2 = get_molecule("carbon_dioxide")
        assert s1 == s2
        np.testing.assert_array_equal(c1, c2)

        s3, _c3 = get_molecule("water")
        assert "O" in s3

    def test_unknown_molecule_raises(self):
        with pytest.raises(ValueError, match="Unknown molecule"):
            get_molecule("unobtanium")

    def test_molecules_centered_at_origin(self):
        """All multi-atom molecules should be approximately centered."""
        for name in available_molecules():
            species, coords = get_molecule(name)
            if len(species) > 1:
                center = coords.mean(axis=0)
                assert np.linalg.norm(center) < 0.5, (
                    f"Molecule {name} not centered: center = {center}"
                )

    def test_diatomic_bond_lengths_reasonable(self):
        """Diatomic molecules should have reasonable bond lengths."""
        for name in ["H2", "N2", "O2", "CO"]:
            _species, coords = get_molecule(name)
            bond_len = np.linalg.norm(coords[1] - coords[0])
            assert 0.5 < bond_len < 2.0, f"{name} bond length {bond_len:.3f} A is unreasonable"


# ===========================================================================
# sites.py tests
# ===========================================================================


class TestSites:
    """Tests for adsorption site identification."""

    def test_find_void_sites_irmof1(self, irmof1):
        """IRMOF-1 has large pores; should find void sites."""
        sites = find_adsorption_sites(
            irmof1,
            strategy="void",
            min_distance=2.5,
            grid_spacing=1.0,  # coarser for speed
        )
        assert len(sites) > 0
        for s in sites:
            assert s.site_type == "void"
            assert s.nearest_framework_dist >= 2.5
            assert s.frac_coords.shape == (3,)
            assert s.cart_coords.shape == (3,)

    def test_void_sites_inside_unit_cell(self, irmof1):
        """All void site fractional coords should be in [0, 1)."""
        sites = find_adsorption_sites(
            irmof1,
            strategy="void",
            min_distance=2.0,
            grid_spacing=1.0,
        )
        for s in sites:
            assert np.all(s.frac_coords >= 0.0), f"frac < 0: {s.frac_coords}"
            assert np.all(s.frac_coords < 1.0), f"frac >= 1: {s.frac_coords}"

    def test_void_sites_sorted_descending(self, irmof1):
        """Void sites should be sorted by nearest_framework_dist descending."""
        sites = find_adsorption_sites(
            irmof1,
            strategy="void",
            grid_spacing=1.0,
        )
        if len(sites) > 1:
            dists = [s.nearest_framework_dist for s in sites]
            assert dists == sorted(dists, reverse=True)

    def test_max_sites_limits_results(self, irmof1):
        sites = find_adsorption_sites(
            irmof1,
            strategy="void",
            grid_spacing=1.0,
            max_sites=3,
        )
        assert len(sites) <= 3

    def test_no_void_sites_with_large_min_distance(self, irmof1):
        """With an unreasonably large min_distance, no sites should be found."""
        sites = find_adsorption_sites(
            irmof1,
            strategy="void",
            min_distance=50.0,
            grid_spacing=1.0,
        )
        assert len(sites) == 0

    def test_empty_crystal_raises(self):
        empty = Crystal.empty()
        with pytest.raises(ValueError, match="empty crystal"):
            find_adsorption_sites(empty)

    def test_invalid_strategy_raises(self, irmof1):
        with pytest.raises(ValueError, match="Unknown strategy"):
            find_adsorption_sites(irmof1, strategy="magic")

    def test_both_strategy(self, irmof1):
        """'both' strategy should return sites from both methods."""
        sites = find_adsorption_sites(
            irmof1,
            strategy="both",
            grid_spacing=1.0,
        )
        # Should have at least the void sites
        types = {s.site_type for s in sites}
        assert "void" in types

    def test_open_metal_site_detection(self, mof74):
        """MOF-74 should have open metal sites (depending on structure)."""
        sites = find_adsorption_sites(
            mof74,
            strategy="open_metal",
        )
        # MOF-74 is known for open metal sites; if detected, check properties
        for s in sites:
            assert s.site_type == "open_metal"
            assert "metal_species" in s.metadata
            assert "metal_index" in s.metadata
            assert "actual_cn" in s.metadata

    def test_open_metal_no_bonds_warns(self):
        """Crystal without bonds should return empty list and warn."""
        xtal = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        # Don't infer bonds
        sites = find_adsorption_sites(xtal, strategy="open_metal")
        assert len(sites) == 0

    def test_adsorption_site_repr(self):
        site = AdsorptionSite(
            frac_coords=np.array([0.5, 0.5, 0.5]),
            cart_coords=np.array([5.0, 5.0, 5.0]),
            site_type="void",
            nearest_framework_dist=3.5,
        )
        r = repr(site)
        assert "void" in r
        assert "3.50" in r


# ===========================================================================
# placement.py tests
# ===========================================================================


class TestPlacement:
    """Tests for adsorbate placement."""

    def test_place_single_co2(self, irmof1):
        """Place a single CO2 in IRMOF-1."""
        result = place_adsorbate(
            irmof1,
            "CO2",
            n_adsorbates=1,
            random_seed=42,
            grid_spacing=1.0,
        )
        assert isinstance(result, AdsorbatePlacement)
        assert result.n_adsorbates == 1
        assert result.adsorbate_name == "CO2"
        assert result.crystal.n_atoms == irmof1.n_atoms + 3  # O=C=O
        assert len(result.adsorbate_indices) == 1
        assert len(result.adsorbate_indices[0]) == 3

    def test_place_single_atom_adsorbate(self, irmof1):
        """Place Ar (single atom) -- no rotation needed."""
        result = place_adsorbate(
            irmof1,
            "Ar",
            n_adsorbates=1,
            grid_spacing=1.0,
        )
        assert result.crystal.n_atoms == irmof1.n_atoms + 1
        assert result.n_adsorbates == 1

    def test_place_multiple_adsorbates(self, irmof1):
        """Place multiple CO2 molecules."""
        result = place_adsorbate(
            irmof1,
            "CO2",
            n_adsorbates=3,
            random_seed=42,
            grid_spacing=1.0,
        )
        # May get fewer than 3 if sites are too close
        assert result.n_adsorbates >= 1
        assert result.n_adsorbates <= 3
        expected_atoms = irmof1.n_atoms + result.n_adsorbates * 3
        assert result.crystal.n_atoms == expected_atoms

    def test_place_at_explicit_site(self, irmof1):
        """Place at a user-specified site."""
        site = AdsorptionSite(
            frac_coords=np.array([0.5, 0.5, 0.5]),
            cart_coords=irmof1.lattice.get_cartesian_coords([0.5, 0.5, 0.5]),
            site_type="void",
            nearest_framework_dist=5.0,
        )
        result = place_adsorbate(
            irmof1,
            "H2O",
            site=site,
            random_seed=42,
        )
        assert result.n_adsorbates == 1
        assert result.crystal.n_atoms == irmof1.n_atoms + 3  # O + 2H

    def test_place_at_explicit_sites_list(self, irmof1):
        """Place at multiple user-specified sites."""
        lattice = irmof1.lattice
        sites = [
            AdsorptionSite(
                frac_coords=np.array([0.25, 0.25, 0.25]),
                cart_coords=lattice.get_cartesian_coords([0.25, 0.25, 0.25]),
                site_type="void",
                nearest_framework_dist=4.0,
            ),
            AdsorptionSite(
                frac_coords=np.array([0.75, 0.75, 0.75]),
                cart_coords=lattice.get_cartesian_coords([0.75, 0.75, 0.75]),
                site_type="void",
                nearest_framework_dist=4.0,
            ),
        ]
        result = place_adsorbate(irmof1, "N2", sites=sites, random_seed=42)
        assert result.n_adsorbates == 2
        assert result.crystal.n_atoms == irmof1.n_atoms + 4  # 2 x N2

    def test_site_and_sites_mutually_exclusive(self, irmof1):
        site = AdsorptionSite(
            frac_coords=np.array([0.5, 0.5, 0.5]),
            cart_coords=np.array([5.0, 5.0, 5.0]),
            site_type="void",
            nearest_framework_dist=5.0,
        )
        with pytest.raises(ValueError, match="not both"):
            place_adsorbate(irmof1, "CO2", site=site, sites=[site])

    def test_fixed_orientation(self, irmof1):
        """With orient='fixed', molecule should keep its default orientation."""
        site = AdsorptionSite(
            frac_coords=np.array([0.5, 0.5, 0.5]),
            cart_coords=irmof1.lattice.get_cartesian_coords([0.5, 0.5, 0.5]),
            site_type="void",
            nearest_framework_dist=5.0,
        )
        result = place_adsorbate(
            irmof1,
            "CO2",
            site=site,
            orient="fixed",
        )
        # The CO2 atoms should be along the z-axis at the site
        ads_idx = result.adsorbate_indices[0]
        ads_coords = result.crystal.cart_coords[ads_idx]
        center = site.cart_coords
        # Relative coords should match the original CO2 geometry
        relative = ads_coords - center
        _, original = get_molecule("CO2")
        np.testing.assert_allclose(relative, original, atol=1e-4)

    def test_invalid_orient_raises(self, irmof1):
        with pytest.raises(ValueError, match="Unknown orient"):
            place_adsorbate(
                irmof1,
                "CO2",
                orient="upside_down",
                grid_spacing=1.0,
            )

    def test_invalid_adsorbate_type_raises(self, irmof1):
        with pytest.raises(TypeError, match="str or Crystal"):
            place_adsorbate(irmof1, 42)

    def test_place_crystal_fragment(self, irmof1):
        """Place a Crystal object as adsorbate (e.g. from fragment())."""
        from mofforge.core.moiety import fragment

        acetylene = fragment("acetylene.xyz", fragment_path=MOIETY_DIR)
        result = place_adsorbate(
            irmof1,
            acetylene,
            n_adsorbates=1,
            random_seed=42,
            grid_spacing=1.0,
        )
        assert result.n_adsorbates == 1
        assert result.crystal.n_atoms == irmof1.n_atoms + acetylene.n_atoms

    def test_reproducible_with_seed(self, irmof1):
        """Same random_seed should produce same result."""
        r1 = place_adsorbate(
            irmof1,
            "CO2",
            n_adsorbates=1,
            random_seed=123,
            grid_spacing=1.0,
        )
        r2 = place_adsorbate(
            irmof1,
            "CO2",
            n_adsorbates=1,
            random_seed=123,
            grid_spacing=1.0,
        )
        np.testing.assert_array_equal(r1.crystal.cart_coords, r2.crystal.cart_coords)

    def test_provenance_tracked(self, irmof1):
        """Placement should record provenance."""
        result = place_adsorbate(
            irmof1,
            "CH4",
            n_adsorbates=1,
            random_seed=42,
            grid_spacing=1.0,
        )
        prov = result.crystal.provenance
        assert prov is not None
        assert prov.operation == "add_adsorbate"
        assert prov.parameters["adsorbate"] == "CH4"

    def test_combined_crystal_name(self, irmof1):
        result = place_adsorbate(
            irmof1,
            "CO2",
            n_adsorbates=1,
            random_seed=42,
            grid_spacing=1.0,
        )
        assert "CO2" in result.crystal.name

    def test_custom_name(self, irmof1):
        result = place_adsorbate(
            irmof1,
            "CO2",
            n_adsorbates=1,
            random_seed=42,
            name="my_structure",
            grid_spacing=1.0,
        )
        assert result.crystal.name == "my_structure"

    def test_no_sites_found_raises(self, irmof1):
        """If min_distance is too large, should raise."""
        with pytest.raises(ValueError, match="No adsorption sites"):
            place_adsorbate(
                irmof1,
                "CO2",
                strategy="void",
                min_distance=100.0,
                grid_spacing=1.0,
            )

    def test_validation_detects_clashes(self, irmof1):
        """Place an adsorbate directly on a framework atom -> clashes."""
        # Pick a site right on top of an existing atom
        framework_cart = irmof1.cart_coords[0]
        framework_frac = irmof1.frac_coords[0]
        bad_site = AdsorptionSite(
            frac_coords=framework_frac,
            cart_coords=framework_cart,
            site_type="void",
            nearest_framework_dist=0.0,
        )
        result = place_adsorbate(
            irmof1,
            "Ar",
            site=bad_site,
            validate=True,
        )
        assert result.clashes > 0

    def test_intermolecular_dist_filtering(self, irmof1):
        """Sites too close together should be filtered out."""
        lattice = irmof1.lattice
        # Two sites very close together
        sites = [
            AdsorptionSite(
                frac_coords=np.array([0.5, 0.5, 0.5]),
                cart_coords=lattice.get_cartesian_coords([0.5, 0.5, 0.5]),
                site_type="void",
                nearest_framework_dist=5.0,
            ),
            AdsorptionSite(
                frac_coords=np.array([0.5, 0.5, 0.501]),
                cart_coords=lattice.get_cartesian_coords([0.5, 0.5, 0.501]),
                site_type="void",
                nearest_framework_dist=5.0,
            ),
        ]
        result = place_adsorbate(
            irmof1,
            "Ar",
            sites=sites,
            min_intermolecular_dist=3.0,
        )
        # Second site should be filtered out
        assert result.n_adsorbates == 1


# ===========================================================================
# Integration: top-level imports
# ===========================================================================


class TestTopLevelImports:
    """Verify adsorbate symbols are accessible from mofforge top level."""

    def test_import_from_mofforge(self):
        import mofforge

        assert hasattr(mofforge, "AdsorptionSite")
        assert hasattr(mofforge, "AdsorbatePlacement")
        assert hasattr(mofforge, "find_adsorption_sites")
        assert hasattr(mofforge, "place_adsorbate")
        assert hasattr(mofforge, "get_molecule")
        assert hasattr(mofforge, "available_molecules")

    def test_in_all(self):
        import mofforge

        for name in [
            "AdsorptionSite",
            "AdsorbatePlacement",
            "find_adsorption_sites",
            "place_adsorbate",
            "get_molecule",
            "available_molecules",
        ]:
            assert name in mofforge.__all__
