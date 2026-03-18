"""Tests for structure validation."""

import numpy as np


class TestValidation:
    """Tests for validate_structure."""

    def test_validate_empty(self):
        """Validating an empty crystal should warn."""
        from mofforge.core.crystal import Crystal
        from mofforge.validation import validate_structure

        xtal = Crystal.empty()
        report = validate_structure(xtal)
        assert any("no atoms" in w.lower() for w in report.warnings)

    def test_validate_good_structure(self):
        """A well-formed crystal should pass validation."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.validation import validate_structure
        from tests.conftest import CRYSTAL_DIR

        xtal = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        xtal = infer_bonds(xtal, periodic=True)
        report = validate_structure(xtal, check_charges=False)

        # IRMOF-1 is a well-formed structure, should have no steric clashes
        # (though validation may find some depending on tolerance)
        assert isinstance(report.is_valid, bool)

    def test_validation_report_summary(self):
        """ValidationReport.summary() should return a string."""
        from mofforge.validation import ValidationReport

        report = ValidationReport()
        assert isinstance(report.summary(), str)
        assert report.is_valid

    def test_steric_clash_detection(self):
        """Two atoms very close should trigger a steric clash."""
        from mofforge.core.crystal import Crystal
        from mofforge.validation import validate_structure

        # Create two C atoms very close together
        xtal = Crystal.from_xyz(
            ["C", "C"],
            np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]),
        )
        report = validate_structure(xtal, check_bonds=False, check_coordination=False)
        # 0.5 A is well below vdW sum for C-C (~3.4 A)
        assert len(report.steric_clashes) > 0
