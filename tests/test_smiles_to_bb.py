"""Tests for SMILES-to-building-block conversion (TOBACCO CIF and Pormake XYZ).

These tests require ``rdkit`` to be installed.  They are skipped
automatically when rdkit is unavailable.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

rdkit = pytest.importorskip("rdkit", reason="rdkit is required for these tests")

from mofforge.build.smiles_to_bb import (
    ConnectionInfo,
    detect_carboxylic_groups,
    detect_connection_points,
    smiles_to_pormake_edge_xyz,
    smiles_to_tobacco_edge_cif,
)


# ------------------------------------------------------------------ #
# detect_connection_points
# ------------------------------------------------------------------ #


class TestDetectConnectionPoints:
    """Tests for automatic connection-point detection."""

    # --- Carboxylate mode ----------------------------------------- #

    def test_bdc_carboxylate(self):
        """BDC (terephthalic acid) should be detected as carboxylate mode."""
        info = detect_connection_points("OC(=O)c1ccc(C(=O)O)cc1")
        assert info.mode == "carboxylate"
        assert len(info.carboxylate_groups) == 2
        assert len(info.connection_atom_indices) == 2
        # Each group should have valid indices
        for g in info.carboxylate_groups:
            assert g.carbon_idx >= 0
            assert g.oxy_double_idx >= 0
            assert g.oxy_single_idx >= 0
            assert g.anchor_idx >= 0
            # Anchor should be different from carboxylate atoms
            assert g.anchor_idx not in (g.carbon_idx, g.oxy_double_idx, g.oxy_single_idx)

    def test_ndc_carboxylate(self):
        """Naphthalene dicarboxylic acid should detect 2 carboxylates."""
        info = detect_connection_points("OC(=O)c1ccc2cc(C(=O)O)ccc2c1")
        assert info.mode == "carboxylate"
        assert len(info.carboxylate_groups) == 2

    def test_deprotonated_carboxylate(self):
        """Deprotonated form C(=O)[O-] should also match."""
        info = detect_connection_points("[O-]C(=O)c1ccc(C(=O)[O-])cc1")
        assert info.mode == "carboxylate"
        assert len(info.carboxylate_groups) == 2

    def test_wrong_carboxylate_count(self):
        """Should raise ValueError when carboxylate count != 2."""
        # Benzoic acid: only 1 carboxylate
        with pytest.raises(ValueError, match="Expected 2.*found 1"):
            detect_connection_points("OC(=O)c1ccccc1")

    def test_three_carboxylates(self):
        """Trimesic acid (3 COOHs) should fail for n_points=2."""
        with pytest.raises(ValueError, match="Expected 2.*found 3"):
            detect_connection_points("OC(=O)c1cc(C(=O)O)cc(C(=O)O)c1")

    # --- Direct mode ---------------------------------------------- #

    def test_biphenyl_direct(self):
        """Biphenyl should be detected as direct mode."""
        info = detect_connection_points("c1ccc(-c2ccccc2)cc1")
        assert info.mode == "direct"
        assert len(info.connection_atom_indices) == 2
        # The two connection points should be different atoms
        assert info.connection_atom_indices[0] != info.connection_atom_indices[1]

    def test_benzene_direct(self):
        """Benzene as a simple linker: direct mode, para carbons."""
        info = detect_connection_points("c1ccccc1")
        assert info.mode == "direct"
        assert len(info.connection_atom_indices) == 2

    def test_azobenzene_direct(self):
        """Azobenzene should be direct mode (no carboxylates)."""
        info = detect_connection_points("c1ccc(/N=N/c2ccccc2)cc1")
        assert info.mode == "direct"
        assert len(info.connection_atom_indices) == 2

    # --- Error cases ---------------------------------------------- #

    def test_invalid_smiles(self):
        """Invalid SMILES should raise ValueError."""
        with pytest.raises(ValueError, match="could not parse"):
            detect_connection_points("not_a_smiles_string!!!")

    def test_single_atom(self):
        """Single heavy atom should raise ValueError."""
        with pytest.raises(ValueError, match="fewer than 2"):
            detect_connection_points("[Cu]")


# ------------------------------------------------------------------ #
# ConnectionInfo dataclass
# ------------------------------------------------------------------ #


class TestConnectionInfo:
    """Tests for the ConnectionInfo dataclass."""

    def test_defaults(self):
        info = ConnectionInfo(mode="direct")
        assert info.connection_atom_indices == []
        assert info.carboxylate_groups == []
        assert info.smiles == ""


# ------------------------------------------------------------------ #
# smiles_to_tobacco_edge_cif
# ------------------------------------------------------------------ #


class TestSmilesToTobaccoEdgeCif:
    """Tests for full SMILES-to-CIF conversion."""

    def test_bdc_cif_output(self, tmp_path: Path):
        """BDC conversion should produce a valid CIF with X/Fr atoms."""
        out = tmp_path / "BDC_edge.cif"
        result = smiles_to_tobacco_edge_cif(
            "OC(=O)c1ccc(C(=O)O)cc1",
            output_path=out,
            name="BDC_edge",
        )
        assert result == out
        assert out.exists()

        content = out.read_text()

        # CIF header checks
        assert "data_BDC_edge" in content
        assert "_symmetry_space_group_name_H-M" in content
        assert "_cell_length_a" in content

        # Should have exactly 2 X atoms with Fr type
        atom_lines = [
            l for l in content.splitlines()
            if l.strip().startswith("X") and "Fr" in l
        ]
        assert len(atom_lines) == 2, f"Expected 2 Fr X atoms, got {len(atom_lines)}"

        # Should have bond entries for X atoms
        assert "X1" in content
        assert "X2" in content

        # Should have atom site and bond loops
        assert "_atom_site_label" in content
        assert "_geom_bond_atom_site_label_1" in content

    def test_biphenyl_cif_output(self, tmp_path: Path):
        """Biphenyl direct mode should produce CIF with X-labelled C atoms."""
        out = tmp_path / "biphenyl_edge.cif"
        result = smiles_to_tobacco_edge_cif(
            "c1ccc(-c2ccccc2)cc1",
            output_path=out,
            name="biphenyl_edge",
        )
        assert result == out
        assert out.exists()

        content = out.read_text()
        assert "data_biphenyl_edge" in content

        # Should have X-labelled atoms that are NOT Fr (they keep C type).
        # Atom lines have many columns (label, symbol, fx, fy, fz, ...);
        # bond lines only have 5 columns.  Filter by column count.
        x_atom_lines = [
            l for l in content.splitlines()
            if l.strip().startswith("X")
            and len(l.split()) > 5
            and "Fr" not in l
            and l.split()[1] == "C"
        ]
        # We expect exactly 2 X atoms with C type symbol
        assert len(x_atom_lines) == 2, f"Expected 2 X(C) atoms, got {len(x_atom_lines)}"

        # Should NOT have any Fr atoms
        assert "Fr" not in content

    def test_biphenyl_no_h_on_x(self, tmp_path: Path):
        """H atoms bonded to connection-point carbons should be removed."""
        out = tmp_path / "biph.cif"
        smiles_to_tobacco_edge_cif("c1ccc(-c2ccccc2)cc1", output_path=out)

        content = out.read_text()
        # Count H atoms: biphenyl has 10 H, we remove 2 (one per X)
        h_lines = [
            l for l in content.splitlines()
            if l.strip() and l.split()[0].startswith("H") and len(l.split()) > 5
        ]
        assert len(h_lines) == 8, f"Expected 8 H atoms, got {len(h_lines)}"

    def test_ndc_cif_output(self, tmp_path: Path):
        """NDC (naphthalene dicarboxylic acid) should work like BDC."""
        out = tmp_path / "NDC.cif"
        smiles_to_tobacco_edge_cif(
            "OC(=O)c1ccc2cc(C(=O)O)ccc2c1",
            output_path=out,
            name="NDC",
        )
        content = out.read_text()
        assert "data_NDC" in content

        # 2 X/Fr dummy atoms
        x_fr_lines = [
            l for l in content.splitlines()
            if l.strip().startswith("X") and "Fr" in l
        ]
        assert len(x_fr_lines) == 2

    def test_bdc_no_carboxylate_h(self, tmp_path: Path):
        """Carboxylate -OH hydrogens should be stripped from BDC."""
        out = tmp_path / "bdc.cif"
        smiles_to_tobacco_edge_cif("OC(=O)c1ccc(C(=O)O)cc1", output_path=out)

        content = out.read_text()
        # BDC: C6H4(COOH)2 has 4 ring H + 2 OH H = 6 H total.
        # After stripping carboxylate OH, we should have 4 H.
        h_lines = [
            l for l in content.splitlines()
            if l.strip() and l.split()[0].startswith("H") and len(l.split()) > 5
        ]
        assert len(h_lines) == 4, f"Expected 4 H atoms after stripping COOH, got {len(h_lines)}"

    def test_aromatic_bonds_present(self, tmp_path: Path):
        """BDC CIF should have aromatic bond entries."""
        out = tmp_path / "bdc_arom.cif"
        smiles_to_tobacco_edge_cif("OC(=O)c1ccc(C(=O)O)cc1", output_path=out)
        content = out.read_text()
        # Should have at least one aromatic bond
        assert "     A" in content, "Expected aromatic bond entries in CIF"

    def test_output_path_returned(self, tmp_path: Path):
        """Return value should be the resolved output path."""
        out = tmp_path / "sub" / "edge.cif"
        result = smiles_to_tobacco_edge_cif(
            "c1ccc(-c2ccccc2)cc1", output_path=out
        )
        assert result == out.resolve()
        assert result.exists()

    def test_invalid_smiles_raises(self, tmp_path: Path):
        """Invalid SMILES should raise ValueError."""
        with pytest.raises(ValueError):
            smiles_to_tobacco_edge_cif(
                "not_valid!!!",
                output_path=tmp_path / "bad.cif",
            )

    def test_cell_length_custom(self, tmp_path: Path):
        """Custom cell length should be reflected in the CIF."""
        out = tmp_path / "custom_cell.cif"
        smiles_to_tobacco_edge_cif(
            "c1ccc(-c2ccccc2)cc1",
            output_path=out,
            cell_length=30.0,
        )
        content = out.read_text()
        assert "30.0000" in content

    def test_bond_distances_positive(self, tmp_path: Path):
        """All bond distances should be positive."""
        out = tmp_path / "dist_check.cif"
        smiles_to_tobacco_edge_cif("OC(=O)c1ccc(C(=O)O)cc1", output_path=out)
        content = out.read_text()

        in_bond_loop = False
        for line in content.splitlines():
            if "_ccdc_geom_bond_type" in line:
                in_bond_loop = True
                continue
            if in_bond_loop and line.strip():
                parts = line.split()
                if len(parts) >= 3:
                    dist = float(parts[2])
                    assert dist > 0, f"Non-positive bond distance: {line}"


# ------------------------------------------------------------------ #
# detect_carboxylic_groups
# ------------------------------------------------------------------ #


class TestDetectCarboxylicGroups:
    """Tests for carboxylic group detection."""

    def test_bdc_carboxylic(self):
        """BDC should detect 2 COOH groups in carboxylic mode."""
        info = detect_carboxylic_groups("OC(=O)c1ccc(C(=O)O)cc1")
        assert info.mode == "carboxylic"
        assert len(info.carboxylate_groups) == 2
        assert len(info.connection_atom_indices) == 2
        for g in info.carboxylate_groups:
            assert g.anchor_idx >= 0

    def test_no_cooh_raises(self):
        """Molecule without COOH should raise ValueError."""
        with pytest.raises(ValueError, match="No COOH"):
            detect_carboxylic_groups("c1ccc(-c2ccccc2)cc1")

    def test_one_cooh_raises(self):
        """Molecule with 1 COOH should raise ValueError."""
        with pytest.raises(ValueError, match="Expected exactly 2"):
            detect_carboxylic_groups("OC(=O)c1ccccc1")

    def test_three_cooh_raises(self):
        """Trimesic acid (3 COOHs) should raise ValueError."""
        with pytest.raises(ValueError, match="Expected exactly 2"):
            detect_carboxylic_groups("OC(=O)c1cc(C(=O)O)cc(C(=O)O)c1")


# ------------------------------------------------------------------ #
# smiles_to_tobacco_edge_cif with mode="carboxylic"
# ------------------------------------------------------------------ #


class TestCarboxylicModeTobaccoCif:
    """Tests for TOBACCO CIF generation in carboxylic mode."""

    def test_bdc_carboxylic_cif(self, tmp_path: Path):
        """BDC with carboxylic mode: strip COOH, anchor becomes X."""
        out = tmp_path / "BDC_carboxylic.cif"
        result = smiles_to_tobacco_edge_cif(
            "OC(=O)c1ccc(C(=O)O)cc1",
            output_path=out,
            name="BDC_carboxylic",
            mode="carboxylic",
        )
        assert result == out
        assert out.exists()

        content = out.read_text()
        assert "data_BDC_carboxylic" in content

        # Should have exactly 2 X atoms (the anchors).
        atom_lines = [
            l for l in content.splitlines()
            if l.strip().startswith("X") and len(l.split()) > 5
        ]
        assert len(atom_lines) == 2, f"Expected 2 X atoms, got {len(atom_lines)}"

        # Should NOT have any Fr atoms (no dummy atoms added).
        assert "Fr" not in content

        # The X atoms should have type_symbol C (they are ring carbons).
        for line in atom_lines:
            parts = line.split()
            assert parts[1] == "C", f"X atom type_symbol should be C, got {parts[1]}"

    def test_bdc_carboxylic_no_carboxylate_atoms(self, tmp_path: Path):
        """Carboxylate C and O atoms should be completely stripped."""
        out = tmp_path / "bdc_strip.cif"
        smiles_to_tobacco_edge_cif(
            "OC(=O)c1ccc(C(=O)O)cc1",
            output_path=out,
            mode="carboxylic",
        )
        content = out.read_text()

        # BDC core after stripping COOH: 6 ring C (2 relabelled as X) + 4 H
        # = 4 C + 2 X + 4 H = 10 atoms total
        # Note: BDC anchor carbons are ring C with no H attached to them
        # (each has 2 ring-C neighbours and the carboxylate-C neighbour).
        atom_lines = [
            l for l in content.splitlines()
            if l.strip() and len(l.split()) > 5 and not l.strip().startswith("_")
            and not l.strip().startswith("loop")
            and not l.strip().startswith("data_")
        ]
        # Filter to only atom site lines (they have labels like C1, H2, X1)
        atom_lines = [
            l for l in atom_lines
            if l.split()[0][0].isalpha() and "." in l  # has fractional coords
            and "Uiso" in l
        ]
        assert len(atom_lines) == 10, (
            f"Expected 10 atoms (4C + 2X + 4H), got {len(atom_lines)}"
        )

        # No O atoms should remain.
        o_lines = [l for l in atom_lines if l.split()[1] == "O"]
        assert len(o_lines) == 0, f"Expected 0 O atoms, got {len(o_lines)}"

    def test_bdc_carboxylic_atom_counts(self, tmp_path: Path):
        """Check specific element counts for BDC carboxylic mode."""
        out = tmp_path / "bdc_counts.cif"
        smiles_to_tobacco_edge_cif(
            "OC(=O)c1ccc(C(=O)O)cc1",
            output_path=out,
            mode="carboxylic",
        )
        content = out.read_text()

        # Extract atom type_symbols from atom lines.
        atom_symbols = []
        for line in content.splitlines():
            parts = line.split()
            if len(parts) >= 9 and "Uiso" in line:
                atom_symbols.append(parts[1])

        from collections import Counter
        counts = Counter(atom_symbols)
        assert counts["C"] == 6, f"Expected 6 C-type atoms, got {counts.get('C', 0)}"
        assert counts["H"] == 4, f"Expected 4 H atoms, got {counts.get('H', 0)}"
        assert counts.get("O", 0) == 0, f"Expected 0 O atoms, got {counts.get('O', 0)}"
        assert counts.get("Fr", 0) == 0, f"Expected 0 Fr atoms, got {counts.get('Fr', 0)}"

    def test_ndc_carboxylic_cif(self, tmp_path: Path):
        """NDC in carboxylic mode should also work."""
        out = tmp_path / "NDC_carboxylic.cif"
        smiles_to_tobacco_edge_cif(
            "OC(=O)c1ccc2cc(C(=O)O)ccc2c1",
            output_path=out,
            name="NDC_carboxylic",
            mode="carboxylic",
        )
        content = out.read_text()
        assert "data_NDC_carboxylic" in content

        # 2 X atoms, no Fr.
        x_lines = [
            l for l in content.splitlines()
            if l.strip().startswith("X") and len(l.split()) > 5
        ]
        assert len(x_lines) == 2
        assert "Fr" not in content

    def test_bpdc_carboxylic_cif(self, tmp_path: Path):
        """BPDC in carboxylic mode."""
        out = tmp_path / "BPDC_carboxylic.cif"
        smiles_to_tobacco_edge_cif(
            "OC(=O)c1ccc(-c2ccc(C(=O)O)cc2)cc1",
            output_path=out,
            name="BPDC_carboxylic",
            mode="carboxylic",
        )
        content = out.read_text()
        assert "data_BPDC_carboxylic" in content

        x_lines = [
            l for l in content.splitlines()
            if l.strip().startswith("X") and len(l.split()) > 5
        ]
        assert len(x_lines) == 2
        assert "Fr" not in content
        assert "O" not in [l.split()[1] for l in content.splitlines()
                           if l.strip() and len(l.split()) > 5 and "Uiso" in l]

    def test_carboxylic_no_cooh_raises(self, tmp_path: Path):
        """Molecule without COOH should raise ValueError in carboxylic mode."""
        with pytest.raises(ValueError, match="No COOH"):
            smiles_to_tobacco_edge_cif(
                "c1ccc(-c2ccccc2)cc1",
                output_path=tmp_path / "bad.cif",
                mode="carboxylic",
            )

    def test_carboxylic_bond_distances_positive(self, tmp_path: Path):
        """All bond distances should be positive in carboxylic mode."""
        out = tmp_path / "dist_check.cif"
        smiles_to_tobacco_edge_cif(
            "OC(=O)c1ccc(C(=O)O)cc1",
            output_path=out,
            mode="carboxylic",
        )
        content = out.read_text()

        in_bond_loop = False
        for line in content.splitlines():
            if "_ccdc_geom_bond_type" in line:
                in_bond_loop = True
                continue
            if in_bond_loop and line.strip():
                parts = line.split()
                if len(parts) >= 3:
                    dist = float(parts[2])
                    assert dist > 0, f"Non-positive bond distance: {line}"

    def test_auto_mode_is_default(self, tmp_path: Path):
        """Default mode='auto' should behave like the original (carboxylate)."""
        out_auto = tmp_path / "auto.cif"
        out_explicit = tmp_path / "explicit.cif"

        smiles_to_tobacco_edge_cif(
            "OC(=O)c1ccc(C(=O)O)cc1",
            output_path=out_auto,
            name="auto_test",
        )
        smiles_to_tobacco_edge_cif(
            "OC(=O)c1ccc(C(=O)O)cc1",
            output_path=out_explicit,
            name="auto_test",
            mode="auto",
        )
        # Both should produce a carboxylate-mode CIF with Fr atoms.
        for out in (out_auto, out_explicit):
            content = out.read_text()
            assert "Fr" in content


# ------------------------------------------------------------------ #
# Helper to parse extended XYZ
# ------------------------------------------------------------------ #


def _parse_pormake_xyz(content: str):
    """Parse a Pormake-format extended XYZ string.

    Returns (atoms, bonds, x_indices_from_comment) where:
    - atoms: list of (symbol, x, y, z) tuples
    - bonds: list of (i, j, bond_type) tuples
    - x_indices_from_comment: list of ints from the comment line
    """
    lines = content.strip().splitlines()
    n_atoms = int(lines[0].strip())
    comment = lines[1].strip()
    x_indices = [int(x) for x in comment.split()] if comment else []

    atoms = []
    for line in lines[2 : 2 + n_atoms]:
        parts = line.split()
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))

    bonds = []
    for line in lines[2 + n_atoms :]:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        bonds.append((int(parts[0]), int(parts[1]), parts[2]))

    return atoms, bonds, x_indices


# ------------------------------------------------------------------ #
# smiles_to_pormake_edge_xyz
# ------------------------------------------------------------------ #


class TestSmilesToPormakeEdgeXyz:
    """Tests for full SMILES-to-Pormake XYZ conversion."""

    def test_bdc_xyz_output(self, tmp_path: Path):
        """BDC carboxylate mode: should strip COO groups and add X atoms."""
        out = tmp_path / "BDC.xyz"
        result = smiles_to_pormake_edge_xyz(
            "OC(=O)c1ccc(C(=O)O)cc1",
            output_path=out,
        )
        assert result == out
        assert out.exists()

        atoms, bonds, x_indices = _parse_pormake_xyz(out.read_text())

        # Exactly 2 X atoms
        x_atoms = [(i, a) for i, a in enumerate(atoms) if a[0] == "X"]
        assert len(x_atoms) == 2, f"Expected 2 X atoms, got {len(x_atoms)}"

        # Comment line should list the X atom indices
        assert len(x_indices) == 2
        for xi in x_indices:
            assert atoms[xi][0] == "X"

        # No carboxylate C or O should remain.
        # BDC core is benzene: 6 C + 4 H = 10 atoms + 2 X = 12 total.
        symbols = [a[0] for a in atoms]
        assert "O" not in symbols, "Carboxylate oxygens should be stripped"
        assert symbols.count("C") == 6
        assert symbols.count("H") == 4
        assert symbols.count("X") == 2
        assert len(atoms) == 12

    def test_biphenyl_xyz_output(self, tmp_path: Path):
        """Biphenyl direct mode: should add X atoms at terminal positions."""
        out = tmp_path / "biphenyl.xyz"
        result = smiles_to_pormake_edge_xyz(
            "c1ccc(-c2ccccc2)cc1",
            output_path=out,
        )
        assert result == out
        assert out.exists()

        atoms, bonds, x_indices = _parse_pormake_xyz(out.read_text())

        x_atoms = [(i, a) for i, a in enumerate(atoms) if a[0] == "X"]
        assert len(x_atoms) == 2

        # Biphenyl: 12 C + 10 H originally; 2 H removed, 2 X added = 22
        # Actually: 12 C + 8 H (2 removed) + 2 X = 22
        symbols = [a[0] for a in atoms]
        assert symbols.count("C") == 12
        assert symbols.count("H") == 8
        assert symbols.count("X") == 2

    def test_x_atom_distance(self, tmp_path: Path):
        """X atoms should be ~0.75 A from their bonded atom."""
        out = tmp_path / "biph_dist.xyz"
        smiles_to_pormake_edge_xyz("c1ccc(-c2ccccc2)cc1", output_path=out)

        atoms, bonds, x_indices = _parse_pormake_xyz(out.read_text())

        for xi in x_indices:
            x_pos = np.array(atoms[xi][1:4])
            # Find the atom bonded to this X
            bonded_idx = None
            for i, j, bt in bonds:
                if i == xi:
                    bonded_idx = j
                    break
                if j == xi:
                    bonded_idx = i
                    break
            assert bonded_idx is not None, f"X atom {xi} has no bond"
            bonded_pos = np.array(atoms[bonded_idx][1:4])
            dist = np.linalg.norm(x_pos - bonded_pos)
            assert abs(dist - 0.75) < 0.01, f"X-atom distance {dist:.4f} != 0.75"

    def test_x_bond_type_single(self, tmp_path: Path):
        """Bonds from X to the connection atom should be single (S)."""
        out = tmp_path / "bond_type.xyz"
        smiles_to_pormake_edge_xyz("OC(=O)c1ccc(C(=O)O)cc1", output_path=out)

        atoms, bonds, x_indices = _parse_pormake_xyz(out.read_text())
        x_set = set(x_indices)
        for i, j, bt in bonds:
            if i in x_set or j in x_set:
                assert bt == "S", f"X bond should be 'S', got '{bt}'"

    def test_ndc_xyz_output(self, tmp_path: Path):
        """NDC carboxylate mode should also work."""
        out = tmp_path / "NDC.xyz"
        smiles_to_pormake_edge_xyz(
            "OC(=O)c1ccc2cc(C(=O)O)ccc2c1",
            output_path=out,
        )
        atoms, bonds, x_indices = _parse_pormake_xyz(out.read_text())
        x_atoms = [a for a in atoms if a[0] == "X"]
        assert len(x_atoms) == 2

    def test_aromatic_bonds_in_xyz(self, tmp_path: Path):
        """Aromatic bonds should be labelled 'A' in the bond table."""
        out = tmp_path / "arom.xyz"
        smiles_to_pormake_edge_xyz("c1ccc(-c2ccccc2)cc1", output_path=out)
        atoms, bonds, _ = _parse_pormake_xyz(out.read_text())
        aromatic = [b for b in bonds if b[2] == "A"]
        assert len(aromatic) > 0, "Expected aromatic bonds in biphenyl"

    def test_atom_count_line(self, tmp_path: Path):
        """First line should be the total atom count."""
        out = tmp_path / "count.xyz"
        smiles_to_pormake_edge_xyz("c1ccc(-c2ccccc2)cc1", output_path=out)
        content = out.read_text()
        first_line = content.splitlines()[0].strip()
        atoms, _, _ = _parse_pormake_xyz(content)
        assert int(first_line) == len(atoms)

    def test_output_creates_parent_dirs(self, tmp_path: Path):
        """Should create parent directories if they don't exist."""
        out = tmp_path / "deep" / "nested" / "edge.xyz"
        result = smiles_to_pormake_edge_xyz(
            "c1ccc(-c2ccccc2)cc1", output_path=out,
        )
        assert result == out.resolve()
        assert result.exists()

    def test_invalid_smiles_raises(self, tmp_path: Path):
        """Invalid SMILES should raise ValueError."""
        with pytest.raises(ValueError):
            smiles_to_pormake_edge_xyz(
                "not_valid!!!",
                output_path=tmp_path / "bad.xyz",
            )

    def test_pormake_loads_generated_xyz(self, tmp_path: Path):
        """Pormake should be able to load the generated XYZ file."""
        pm = pytest.importorskip("pormake", reason="pormake required")

        out = tmp_path / "biph_pm.xyz"
        smiles_to_pormake_edge_xyz("c1ccc(-c2ccccc2)cc1", output_path=out)

        bb = pm.BuildingBlock(str(out))
        assert len(bb.connection_point_indices) == 2

    def test_bdc_pormake_loads(self, tmp_path: Path):
        """Pormake should load BDC carboxylate edge correctly."""
        pm = pytest.importorskip("pormake", reason="pormake required")

        out = tmp_path / "bdc_pm.xyz"
        smiles_to_pormake_edge_xyz("OC(=O)c1ccc(C(=O)O)cc1", output_path=out)

        bb = pm.BuildingBlock(str(out))
        assert len(bb.connection_point_indices) == 2
