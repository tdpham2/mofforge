"""Tests for the visualization / rendering subsystem."""

import tempfile

import numpy as np
import pytest

from tests.conftest import CRYSTAL_DIR, MOIETY_DIR


class TestColors:
    """Tests for vis.colors module."""

    def test_jmol_colors_populated(self):
        """JMOL_COLORS should contain common elements."""
        from mofforge.vis.colors import JMOL_COLORS

        assert "C" in JMOL_COLORS
        assert "Zn" in JMOL_COLORS
        assert "H" in JMOL_COLORS
        assert "O" in JMOL_COLORS
        assert len(JMOL_COLORS) > 80

    def test_metals_set(self):
        """METALS should contain common MOF metals."""
        from mofforge.vis.colors import METALS

        assert "Zn" in METALS
        assert "Cu" in METALS
        assert "Fe" in METALS
        assert "Zr" in METALS
        # Non-metals should not be present
        assert "C" not in METALS
        assert "H" not in METALS
        assert "O" not in METALS

    def test_metals_is_frozenset(self):
        """METALS should be immutable."""
        from mofforge.vis.colors import METALS

        assert isinstance(METALS, frozenset)

    def test_get_element_color_known(self):
        """get_element_color should return correct colors for known elements."""
        from mofforge.vis.colors import get_element_color

        assert get_element_color("C") == "#909090"
        assert get_element_color("O") == "#FF0D0D"
        assert get_element_color("Zn") == "#7D80B0"

    def test_get_element_color_unknown(self):
        """get_element_color should return default color for unknown elements."""
        from mofforge.vis.colors import DEFAULT_COLOR, get_element_color

        assert get_element_color("Xx") == DEFAULT_COLOR

    def test_get_element_color_strips_r_tag(self):
        """get_element_color should strip R-group '!' tag."""
        from mofforge.vis.colors import get_element_color

        assert get_element_color("H!") == get_element_color("H")
        assert get_element_color("C!") == get_element_color("C")


class TestAtomLabels:
    """Tests for generate_atom_labels function."""

    def test_sequential_labels(self):
        """Sequential labels should be Element+global_index."""
        from mofforge.vis.render import generate_atom_labels

        species = ["C", "N", "O", "C", "H"]
        labels = generate_atom_labels(species, mode="sequential")
        assert labels == ["C1", "N2", "O3", "C4", "H5"]

    def test_per_element_labels(self):
        """Per-element labels should count each element separately."""
        from mofforge.vis.render import generate_atom_labels

        species = ["C", "N", "O", "C", "H"]
        labels = generate_atom_labels(species, mode="per_element")
        assert labels == ["C1", "N1", "O1", "C2", "H1"]

    def test_none_labels(self):
        """Mode 'none' should return an empty list."""
        from mofforge.vis.render import generate_atom_labels

        labels = generate_atom_labels(["C", "O"], mode="none")
        assert labels == []

    def test_invalid_mode(self):
        """Invalid label mode should raise ValueError."""
        from mofforge.vis.render import generate_atom_labels

        with pytest.raises(ValueError, match="Unknown label mode"):
            generate_atom_labels(["C"], mode="invalid")

    def test_labels_strip_r_tags(self):
        """Labels should use clean element symbols, stripping '!' tags."""
        from mofforge.vis.render import generate_atom_labels

        species = ["C", "H!", "O"]
        labels = generate_atom_labels(species, mode="sequential")
        assert labels == ["C1", "H2", "O3"]

    def test_empty_species_list(self):
        """Empty species list should return empty labels."""
        from mofforge.vis.render import generate_atom_labels

        labels = generate_atom_labels([], mode="sequential")
        assert labels == []


class TestBuildHtml:
    """Tests for build_html function."""

    def _make_crystal(self):
        """Create a simple test Crystal."""
        from mofforge.core.crystal import Crystal

        species = ["C", "O", "H", "H"]
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.2, 0.0, 0.0],
                [-0.5, 0.9, 0.0],
                [-0.5, -0.9, 0.0],
            ]
        )
        return Crystal.from_xyz(species, coords, name="test_mol")

    def test_build_html_returns_string(self):
        """build_html should return a non-empty HTML string."""
        from mofforge.vis.render import build_html

        crystal = self._make_crystal()
        html = build_html(crystal)
        assert isinstance(html, str)
        assert len(html) > 100

    def test_build_html_contains_3dmol_script(self):
        """HTML should include the 3Dmol.js CDN script."""
        from mofforge.vis.render import build_html

        crystal = self._make_crystal()
        html = build_html(crystal)
        assert "3Dmol" in html
        assert "3dmol.org/build/3Dmol-min.js" in html

    def test_build_html_contains_viewer_div(self):
        """HTML should include the viewer div."""
        from mofforge.vis.render import build_html

        crystal = self._make_crystal()
        html = build_html(crystal)
        assert 'id="viewer"' in html

    def test_build_html_contains_ready_signal(self):
        """HTML should set the __3dmol_ready flag."""
        from mofforge.vis.render import build_html

        crystal = self._make_crystal()
        html = build_html(crystal)
        assert "__3dmol_ready" in html

    def test_build_html_with_labels(self):
        """HTML should include atom labels when label_mode is not 'none'."""
        from mofforge.vis.render import build_html

        crystal = self._make_crystal()
        html = build_html(crystal, label_mode="sequential")
        assert "C1" in html
        assert "O2" in html

    def test_build_html_without_labels(self):
        """HTML should not include atom index labels when label_mode is 'none'."""
        from mofforge.vis.render import build_html

        crystal = self._make_crystal()
        html = build_html(crystal, label_mode="none", show_formula=False)
        assert "addLabel" not in html

    def test_build_html_with_formula(self):
        """HTML should include chemical formula when show_formula=True."""
        from mofforge.vis.render import build_html

        crystal = self._make_crystal()
        html = build_html(crystal, show_formula=True)
        # pymatgen's reduced formula for CH2O
        assert "fontColor: 'yellow'" in html

    def test_build_html_without_formula(self):
        """HTML should not include formula overlay when show_formula=False."""
        from mofforge.vis.render import build_html

        crystal = self._make_crystal()
        html = build_html(crystal, show_formula=False, label_mode="none")
        assert "fontColor: 'yellow'" not in html

    def test_build_html_dimensions(self):
        """HTML should use specified width and height."""
        from mofforge.vis.render import build_html

        crystal = self._make_crystal()
        html = build_html(crystal, width=1024, height=768)
        assert "1024px" in html
        assert "768px" in html

    def test_build_html_ball_stick(self):
        """ball_stick representation should include both sphere and stick styles."""
        from mofforge.vis.render import build_html

        crystal = self._make_crystal()
        html = build_html(crystal, representation="ball_stick")
        assert "sphere" in html
        assert "stick" in html

    def test_build_html_stick_only(self):
        """stick representation should only include stick style."""
        from mofforge.vis.render import build_html

        crystal = self._make_crystal()
        html = build_html(crystal, representation="stick")
        assert "stick" in html

    def test_build_html_sphere_only(self):
        """sphere representation should only include sphere style."""
        from mofforge.vis.render import build_html

        crystal = self._make_crystal()
        html = build_html(crystal, representation="sphere")
        assert "sphere" in html

    def test_build_html_with_rotation(self):
        """HTML should include rotation commands when rotate is set."""
        from mofforge.vis.render import build_html

        crystal = self._make_crystal()
        html = build_html(crystal, rotate=(45, 30, 0))
        assert "viewer.rotate(45, 'x')" in html
        assert "viewer.rotate(30, 'y')" in html

    def test_build_html_bg_color(self):
        """HTML should use the specified background color."""
        from mofforge.vis.render import build_html

        crystal = self._make_crystal()
        html = build_html(crystal, bg_color="black")
        assert '"black"' in html


class TestBuildHtmlUnitCell:
    """Tests for unit cell rendering in build_html."""

    def _make_periodic_crystal(self):
        """Create a Crystal from a CIF file (periodic)."""
        from mofforge.core.crystal import Crystal

        return Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")

    def test_unit_cell_disabled_by_default(self):
        """Unit cell should not be drawn by default."""
        from mofforge.vis.render import build_html

        crystal = self._make_periodic_crystal()
        html = build_html(crystal, show_unit_cell=False, label_mode="none")
        assert "addCylinder" not in html

    def test_unit_cell_enabled(self):
        """Unit cell edges should appear when show_unit_cell=True."""
        from mofforge.vis.render import build_html

        crystal = self._make_periodic_crystal()
        html = build_html(crystal, show_unit_cell=True, label_mode="none")
        assert "addCylinder" in html
        # 12 edges in a parallelepiped
        assert html.count("addCylinder") == 12


class TestBuildHtmlWithMetals:
    """Tests for metal atom handling in build_html."""

    def _make_metal_crystal(self):
        """Create a Crystal with a metal atom."""
        from mofforge.core.crystal import Crystal

        species = ["Zn", "O", "C", "C"]
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [3.2, 0.0, 0.0],
                [4.4, 0.0, 0.0],
            ]
        )
        return Crystal.from_xyz(species, coords, name="zn_test")

    def test_metal_atoms_get_larger_spheres(self):
        """Metal atoms should have larger sphere scale in ball_stick mode."""
        from mofforge.vis.render import build_html

        crystal = self._make_metal_crystal()
        html = build_html(crystal, representation="ball_stick")
        assert "elem: 'Zn'" in html
        assert "scale: 0.75" in html


class TestCrystalToXyz:
    """Tests for internal Crystal to XYZ string conversion."""

    def test_xyz_string_format(self):
        """XYZ string should be valid format with atom count and coords."""
        from mofforge.core.crystal import Crystal
        from mofforge.vis.render import _crystal_to_xyz_string

        species = ["C", "O"]
        coords = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]])
        crystal = Crystal.from_xyz(species, coords, name="co")

        xyz = _crystal_to_xyz_string(crystal)
        lines = xyz.strip().split("\n")
        assert lines[0].strip() == "2"  # atom count
        assert lines[1].strip() == "co"  # name
        assert "C" in lines[2]
        assert "O" in lines[3]

    def test_xyz_string_strips_r_tags(self):
        """XYZ string should use clean species (no '!' tags)."""
        from mofforge.core.crystal import Crystal
        from mofforge.vis.render import _crystal_to_xyz_string

        species = ["C", "H!"]
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        crystal = Crystal.from_xyz(species, coords, name="test")

        xyz = _crystal_to_xyz_string(crystal)
        assert "H!" not in xyz
        # H should appear (without the !)
        lines = xyz.strip().split("\n")
        assert lines[3].strip().startswith("H")


class TestFileLoading:
    """Tests for _load_file helper."""

    def test_load_cif(self):
        """Should load a CIF file into a Crystal."""
        from mofforge.vis.render import _load_file

        crystal = _load_file(str(CRYSTAL_DIR / "IRMOF-1.cif"))
        assert crystal.n_atoms > 0
        assert crystal.name == "IRMOF-1"

    def test_load_xyz(self):
        """Should load an XYZ file into a Crystal."""
        from mofforge.vis.render import _load_file

        crystal = _load_file(str(MOIETY_DIR / "p-phenylene.xyz"))
        assert crystal.n_atoms == 10

    def test_load_unsupported_format(self):
        """Should raise ValueError for unsupported formats."""
        from mofforge.vis.render import _load_file

        with pytest.raises(ValueError, match="Unsupported file format"):
            _load_file("structure.pdb")


class TestRenderCLI:
    """Tests for the render CLI command."""

    def test_render_help(self):
        """Render command should show help."""
        from click.testing import CliRunner

        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["render", "--help"])
        assert result.exit_code == 0
        assert "input" in result.output.lower()
        assert "output" in result.output.lower()

    def test_render_missing_input(self):
        """Render command should fail without --input."""
        from click.testing import CliRunner

        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["render"])
        assert result.exit_code != 0


class TestRenderIntegration:
    """Integration tests for render_to_png (requires Playwright + Chromium).

    These tests are marked with ``pytest.mark.integration`` and will be
    skipped if Playwright is not installed or Chromium is not available.
    """

    @pytest.fixture(autouse=True)
    def _check_playwright(self):
        """Skip tests if Playwright is not available."""
        try:
            from playwright.sync_api import sync_playwright

            with sync_playwright() as p:
                browser = p.chromium.launch(headless=True)
                browser.close()
        except Exception:
            pytest.skip("Playwright with Chromium not available")

    def test_render_xyz_to_png(self):
        """Should render an XYZ file to PNG."""
        from mofforge.vis.render import render_file_to_png

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            output = render_file_to_png(
                input_file=str(MOIETY_DIR / "p-phenylene.xyz"),
                output_file=f.name,
                label_mode="sequential",
            )
        assert output.endswith(".png")
        # Check that the file is non-empty (valid PNG)
        import os

        assert os.path.getsize(output) > 1000  # PNG should be at least 1KB

    def test_render_cif_to_png(self):
        """Should render a CIF file to PNG with unit cell."""
        from mofforge.vis.render import render_file_to_png

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            output = render_file_to_png(
                input_file=str(CRYSTAL_DIR / "IRMOF-1.cif"),
                output_file=f.name,
                label_mode="none",
                show_unit_cell=True,
            )
        assert output.endswith(".png")
        import os

        assert os.path.getsize(output) > 1000

    def test_render_crystal_object(self):
        """Should render a Crystal object directly to PNG."""
        from mofforge.core.crystal import Crystal
        from mofforge.vis.render import render_to_png

        species = ["C", "O", "H", "H"]
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.2, 0.0, 0.0],
                [-0.5, 0.9, 0.0],
                [-0.5, -0.9, 0.0],
            ]
        )
        crystal = Crystal.from_xyz(species, coords, name="test")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
            output = render_to_png(crystal, output_file=f.name)

        assert output.endswith(".png")
        import os

        assert os.path.getsize(output) > 1000
