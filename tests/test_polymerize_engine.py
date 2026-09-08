"""Unit tests for M1: monomer prep, box sizing, and System->Crystal conversion.

These require rdkit but NOT pysimm / Packmol / LAMMPS: the pysimm ``System`` is
replaced with a lightweight fake so the conversion logic is testable in CI.  The
real end-to-end run is covered by ``test_polymerize_integration.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest

from mofforge.polymerize import Monomer, PopBuilder
from mofforge.polymerize.convert import system_to_crystal
from mofforge.polymerize.engine import PysimmBackend
from mofforge.polymerize.monomer import box_length_for_density, prepare_monomer

# ---------------------------------------------------------------------------
# Monomer preparation (needs rdkit)
# ---------------------------------------------------------------------------


def test_prepare_smiles_monomer_detects_sites(tmp_path):
    pytest.importorskip("rdkit")
    prep = prepare_monomer(Monomer(name="diamine", source="NCCN"), tmp_path, random_seed=42)
    assert [s.site_type for s in prep.sites] == ["amine", "amine"]
    assert prep.mol_path.is_file()
    assert prep.n_atoms == 12  # NCCN + 8 H
    assert prep.molar_mass == pytest.approx(60.1, abs=0.2)


def test_prepare_monomer_no_sites_raises(tmp_path):
    pytest.importorskip("rdkit")
    # Methane has no curated reactive group.
    with pytest.raises(ValueError, match="No reactive sites"):
        prepare_monomer(Monomer(name="methane", source="C"), tmp_path)


def test_prepare_monomer_explicit_sites(tmp_path):
    pytest.importorskip("rdkit")
    from mofforge.polymerize import ReactiveSite

    mono = Monomer(name="x", source="c1ccccc1", sites=[ReactiveSite(atom_idx=0, site_type="vinyl")])
    prep = prepare_monomer(mono, tmp_path, random_seed=1)
    assert [s.site_type for s in prep.sites] == ["vinyl"]


def test_prepare_from_xyz_with_anchor_tags(tmp_path):
    pytest.importorskip("rdkit")
    xyz = tmp_path / "frag.xyz"
    xyz.write_text(
        "3\ntest\nC!   0.0 0.0 0.0\nC    1.5 0.0 0.0\nC!   3.0 0.0 0.0\n"
    )
    prep = prepare_monomer(Monomer(name="frag", source=str(xyz)), tmp_path)
    assert [s.atom_idx for s in prep.sites] == [0, 2]
    # The re-emitted XYZ must have the '!' tags stripped.
    assert "!" not in prep.mol_path.read_text()


# ---------------------------------------------------------------------------
# Box sizing
# ---------------------------------------------------------------------------


def test_box_length_scales_with_density(tmp_path):
    pytest.importorskip("rdkit")
    prep = prepare_monomer(Monomer(name="diamine", source="NCCN"), tmp_path, random_seed=42)
    low = box_length_for_density([prep], [10], 0.5)
    high = box_length_for_density([prep], [10], 1.5)
    # Higher density -> smaller box.
    assert high < low


def test_box_length_bad_density(tmp_path):
    pytest.importorskip("rdkit")
    prep = prepare_monomer(Monomer(name="diamine", source="NCCN"), tmp_path, random_seed=42)
    with pytest.raises(ValueError, match="target_density must be positive"):
        box_length_for_density([prep], [10], 0.0)


def test_monomer_counts_split():
    assert PysimmBackend._monomer_counts(2, 20) == [10, 10]
    assert PysimmBackend._monomer_counts(3, 50) == [17, 17, 16]
    assert PysimmBackend._monomer_counts(1, None) == [20]


# ---------------------------------------------------------------------------
# System -> Crystal conversion (fake pysimm System)
# ---------------------------------------------------------------------------


@dataclass
class _FakeType:
    elem: str


@dataclass
class _FakeParticle:
    tag: int
    x: float
    y: float
    z: float
    elem: str

    @property
    def type(self):
        return _FakeType(self.elem)

    def get_chem_element(self):
        return self.elem


@dataclass
class _FakeBond:
    a: _FakeParticle
    b: _FakeParticle


class _FakeDim:
    xlo, xhi = 0.0, 10.0
    ylo, yhi = 0.0, 12.0
    zlo, zhi = 0.0, 14.0
    dx, dy, dz = 10.0, 12.0, 14.0


class _FakeSystem:
    def __init__(self, particles, bonds):
        self.particles = particles
        self.bonds = bonds
        self.dim = _FakeDim()


def test_system_to_crystal_roundtrips_atoms_and_bonds():
    p1 = _FakeParticle(1, 1.0, 1.0, 1.0, "C")
    p2 = _FakeParticle(2, 2.5, 1.0, 1.0, "N")
    p3 = _FakeParticle(3, 4.0, 1.0, 1.0, "C")
    system = _FakeSystem([p1, p2, p3], [_FakeBond(p1, p2), _FakeBond(p2, p3)])

    crystal = system_to_crystal(system, name="pop_test")

    assert crystal.n_atoms == 3
    assert crystal.species == ["C", "N", "C"]
    assert crystal.n_bonds == 2
    # Orthorhombic box from the fake dim.
    lengths = crystal.lattice.abc
    assert lengths == pytest.approx((10.0, 12.0, 14.0))
    assert np.allclose(crystal.cart_coords[1], [2.5, 1.0, 1.0])


def test_system_to_crystal_empty_raises():
    with pytest.raises(ValueError, match="no particles"):
        system_to_crystal(_FakeSystem([], []))


def test_system_to_crystal_bad_box_raises():
    class _BadDim(_FakeDim):
        dx = 0.0

    sys = _FakeSystem([_FakeParticle(1, 0, 0, 0, "C")], [])
    sys.dim = _BadDim()
    with pytest.raises(ValueError, match="non-positive box dimension"):
        system_to_crystal(sys)


# ---------------------------------------------------------------------------
# Backend graceful failure without binaries
# ---------------------------------------------------------------------------


def test_build_without_binaries_reports_error(tmp_path, monkeypatch):
    pytest.importorskip("rdkit")
    monkeypatch.setenv("PATH", "")
    monkeypatch.delenv("MOFFORGE_PACKMOL_BIN", raising=False)
    monkeypatch.delenv("MOFFORGE_LAMMPS_BIN", raising=False)
    b = PopBuilder()
    b.add_monomer("NCCN", name="diamine")
    b.add_monomer("O=Cc1ccc(C=O)cc1", name="dial")
    result = b.build(output_dir=str(tmp_path), random_seed=42)
    assert result.success is False
    assert result.backend == "pysimm"
    assert any("not found" in e for e in result.errors)
