"""Tests for the solvent removal module."""

import numpy as np

from tests.conftest import CRYSTAL_DIR, MOIETY_DIR


def _make_solvated_crystal():
    """Create a test crystal with artificial solvent molecules.

    Loads IRMOF-1 (clean framework), then adds two water molecules
    at positions far from the framework so they remain disconnected.
    """
    from mofforge.core.bonding import infer_bonds
    from mofforge.core.crystal import Crystal

    parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
    parent = infer_bonds(parent, periodic=True)

    # Create water molecule 1 at the large pore center (0.5, 0.5, 0.5)
    water1 = Crystal.from_xyz(
        species=["O", "H", "H"],
        cart_coords=np.array([
            [12.9, 12.9, 12.9],
            [12.9, 13.86, 12.9],
            [12.9, 12.9, 13.86],
        ]),
        name="water1",
        lattice=parent.lattice,
    )

    # Create water molecule 2 at the corner pore (0, 0, 0)
    water2 = Crystal.from_xyz(
        species=["O", "H", "H"],
        cart_coords=np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.96, 0.0],
            [0.0, 0.0, 0.96],
        ]),
        name="water2",
        lattice=parent.lattice,
    )

    solvated = parent + water1 + water2
    solvated = infer_bonds(solvated, periodic=True)
    return solvated, parent.n_atoms


class TestRemoveSolvent:
    """Tests for the core remove_solvent function."""

    def test_basic_desolvation(self):
        """Solvent molecules should be removed, framework preserved."""
        from mofforge.solvent.removal import remove_solvent

        solvated, n_framework = _make_solvated_crystal()
        result = remove_solvent(solvated)

        assert result.crystal.n_atoms == n_framework
        assert result.n_atoms_removed == 6  # 2 water * 3 atoms
        assert result.n_components_removed == 2
        assert result.n_framework_components == 1

    def test_no_solvent_unchanged(self):
        """Clean framework should return unchanged."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.solvent.removal import remove_solvent

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)

        result = remove_solvent(parent)

        assert result.crystal.n_atoms == parent.n_atoms
        assert result.n_atoms_removed == 0
        assert result.n_components_removed == 0
        assert len(result.removed_molecules) == 0

    def test_empty_crystal(self):
        """Empty crystal should not crash."""
        from mofforge.core.crystal import Crystal
        from mofforge.solvent.removal import remove_solvent

        empty = Crystal.empty()
        result = remove_solvent(empty)

        assert result.crystal.n_atoms == 0
        assert result.n_atoms_removed == 0

    def test_min_atoms_threshold(self):
        """Components with >= min_atoms should be kept."""
        from mofforge.solvent.removal import remove_solvent

        solvated, n_framework = _make_solvated_crystal()

        # Water has 3 atoms; min_atoms=3 should keep them
        result = remove_solvent(solvated, min_atoms=3)
        assert result.crystal.n_atoms == solvated.n_atoms
        assert result.n_components_removed == 0

        # min_atoms=4 should still remove water (3 < 4)
        result = remove_solvent(solvated, min_atoms=4)
        assert result.n_components_removed == 2

    def test_n_framework_components(self):
        """Explicit n_framework_components should keep that many largest."""
        from mofforge.solvent.removal import remove_solvent

        solvated, n_framework = _make_solvated_crystal()

        # Keep only 1 framework component (the default behavior)
        result = remove_solvent(solvated, n_framework_components=1)
        assert result.crystal.n_atoms == n_framework
        assert result.n_framework_components == 1

    def test_result_fields(self):
        """SolventRemovalResult should have correct field values."""
        from mofforge.solvent.removal import remove_solvent

        solvated, n_framework = _make_solvated_crystal()
        result = remove_solvent(solvated)

        assert result.n_atoms_original == solvated.n_atoms
        assert result.n_atoms_removed == solvated.n_atoms - n_framework
        assert result.n_components_removed == len(result.removed_molecules)
        assert all(m.n_atoms == 3 for m in result.removed_molecules)

    def test_removed_molecule_formula(self):
        """Removed molecules should have correct Hill-order formula."""
        from mofforge.solvent.removal import remove_solvent

        solvated, _ = _make_solvated_crystal()
        result = remove_solvent(solvated)

        formulas = {m.formula for m in result.removed_molecules}
        assert "H2O" in formulas

    def test_bonds_inferred_if_missing(self):
        """If input has no bonds, they should be inferred automatically."""
        from mofforge.core.crystal import Crystal
        from mofforge.solvent.removal import remove_solvent

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        assert parent.n_bonds == 0

        result = remove_solvent(parent)

        # Should work without error and detect a single framework component
        assert result.crystal.n_atoms == parent.n_atoms
        assert result.n_framework_components == 1

    def test_summary(self):
        """Summary method should return a readable string."""
        from mofforge.solvent.removal import remove_solvent

        solvated, _ = _make_solvated_crystal()
        result = remove_solvent(solvated)

        summary = result.summary()
        assert "Solvent Removal" in summary
        assert "atoms removed" in summary


class TestCompositionFormula:
    """Tests for the Hill-order formula generator."""

    def test_water(self):
        from mofforge.solvent.removal import _composition_formula

        assert _composition_formula(["O", "H", "H"]) == "H2O"

    def test_co2(self):
        from mofforge.solvent.removal import _composition_formula

        assert _composition_formula(["C", "O", "O"]) == "CO2"

    def test_dmf(self):
        from mofforge.solvent.removal import _composition_formula

        # DMF: C3H7NO
        species = ["C", "C", "C", "H", "H", "H", "H", "H", "H", "H", "N", "O"]
        assert _composition_formula(species) == "C3H7NO"

    def test_single_element(self):
        from mofforge.solvent.removal import _composition_formula

        assert _composition_formula(["O"]) == "O"

    def test_tagged_species(self):
        from mofforge.solvent.removal import _composition_formula

        assert _composition_formula(["O", "H!", "H!"]) == "H2O"


class TestDesolvateCommand:
    """Tests for the desolvate CLI command."""

    def test_desolvate_help(self):
        from click.testing import CliRunner

        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["desolvate", "--help"])
        assert result.exit_code == 0
        assert "solvent" in result.output.lower()

    def test_desolvate_basic(self, tmp_path):
        from click.testing import CliRunner

        from mofforge.cli import main

        cif = str(CRYSTAL_DIR / "IRMOF-1.cif")
        output = str(tmp_path / "desolvated.cif")

        runner = CliRunner()
        result = runner.invoke(main, ["desolvate", "-p", cif, "-o", output])
        assert result.exit_code == 0
        assert "Output written to" in result.output


class TestPipelineDesolvate:
    """Tests for the desolvate pipeline step."""

    def test_pipeline_desolvate(self):
        """Pipeline with desolvate step should work."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.pipeline import Pipeline

        solvated, n_framework = _make_solvated_crystal()

        child = (
            Pipeline(solvated)
            .desolvate()
            .build(name="desolvated")
        )

        assert child.n_atoms == n_framework
        assert child.name == "desolvated"

    def test_pipeline_desolvate_with_replace(self):
        """Pipeline combining replace and desolvate."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.pipeline import Pipeline

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)

        child = (
            Pipeline(parent, fragment_path=MOIETY_DIR)
            .replace(
                query="2-!-p-phenylene.xyz",
                replacement="2-acetylamido-p-phenylene.xyz",
                nb_loc=1,
            )
            .desolvate()
            .build(name="modified_and_desolvated")
        )

        assert child.n_atoms > 0
        assert child.name == "modified_and_desolvated"
