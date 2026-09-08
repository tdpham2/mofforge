"""Unit tests for the amorphous POP subsystem scaffolding (M0).

These cover the dataclasses, curated reaction/site tables, RDKit reactive-site
detection, the PopBuilder facade, and external-binary resolution.  No external
binaries (Packmol / LAMMPS) or pysimm are required; binary resolution is tested
by pointing the environment at a dummy executable.
"""

from __future__ import annotations

import stat

import pytest

from mofforge.polymerize import Monomer, PopBuilder, PopConfig, POPResult, ReactiveSite, reactions
from mofforge.polymerize import config as popconfig

# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


def test_monomer_is_smiles():
    assert Monomer(name="m", source="NCCN").is_smiles
    assert not Monomer(name="m", source="linker.xyz").is_smiles
    assert not Monomer(name="m", source="linker.cif").is_smiles


def test_popconfig_defaults():
    cfg = PopConfig()
    assert cfg.target_density == pytest.approx(0.8)
    assert cfg.forcefield == "gaff2"
    assert cfg.equilibrate is True
    assert cfg.box_length is None


def test_popresult_defaults():
    r = POPResult(success=False)
    assert r.output_paths == []
    assert r.errors == []
    assert r.crystal is None


def test_reactive_site_fields():
    site = ReactiveSite(atom_idx=0, anchor_idx=1, site_type="amine")
    assert site.atom_idx == 0
    assert site.site_type == "amine"


# ---------------------------------------------------------------------------
# Curated reactions / site types
# ---------------------------------------------------------------------------


def test_available_site_types():
    types = {g["site_type"] for g in reactions.available_site_types()}
    assert {"amine", "aldehyde", "aryl_halide", "vinyl"}.issubset(types)


def test_available_reactions():
    names = {r["reaction"] for r in reactions.available_reactions()}
    assert any("imine" in n for n in names)


def test_compatibility_table():
    assert reactions.are_compatible("amine", "aldehyde")
    assert reactions.are_compatible("aldehyde", "amine")  # symmetric
    assert reactions.are_compatible("aryl_halide", "aryl_halide")  # homo-coupling
    assert not reactions.are_compatible("amine", "amine")
    assert reactions.reaction_name("amine", "aldehyde") is not None
    assert reactions.reaction_name("amine", "vinyl") is None


def test_get_group_unknown_raises():
    with pytest.raises(ValueError, match="Unknown reactive site type"):
        reactions.get_group("nonsense")


# ---------------------------------------------------------------------------
# Reactive-site detection (needs rdkit)
# ---------------------------------------------------------------------------


def test_detect_sites_diamine():
    pytest.importorskip("rdkit")
    sites = reactions.detect_reactive_sites("NCCN")
    assert [s.site_type for s in sites] == ["amine", "amine"]


def test_detect_sites_trialdehyde():
    pytest.importorskip("rdkit")
    sites = reactions.detect_reactive_sites("O=Cc1cc(C=O)cc(C=O)c1")
    assert [s.site_type for s in sites] == ["aldehyde", "aldehyde", "aldehyde"]


def test_detect_sites_restricted_type():
    pytest.importorskip("rdkit")
    # Only look for aldehydes on a molecule that also has an amine.
    sites = reactions.detect_reactive_sites("NCc1ccc(C=O)cc1", site_types=["aldehyde"])
    assert [s.site_type for s in sites] == ["aldehyde"]


def test_detect_sites_bad_smiles():
    pytest.importorskip("rdkit")
    with pytest.raises(ValueError, match="could not parse SMILES"):
        reactions.detect_reactive_sites("not-a-smiles((")


# ---------------------------------------------------------------------------
# PopBuilder facade
# ---------------------------------------------------------------------------


def test_builder_unknown_backend():
    with pytest.raises(ValueError, match="Unknown backend"):
        PopBuilder(backend="nope")


def test_builder_add_and_list_monomers():
    b = PopBuilder()
    b.add_monomer("NCCN", name="diamine")
    b.add_monomer("O=Cc1ccc(C=O)cc1", name="dial")
    assert b.list_monomers() == ["diamine", "dial"]


def test_builder_build_without_monomers_fails():
    result = PopBuilder().build()
    assert result.success is False
    assert "No monomers" in result.errors[0]


def test_builder_rejects_unknown_option():
    b = PopBuilder()
    b.add_monomer("NCCN")
    with pytest.raises(TypeError, match="Unknown build option"):
        b.build(bogus_option=1)


def test_builder_make_config_passes_known_fields():
    cfg = PopBuilder._make_config({"target_density": 1.1, "forcefield": "pcff"})
    assert cfg.target_density == pytest.approx(1.1)
    assert cfg.forcefield == "pcff"


# ---------------------------------------------------------------------------
# External-binary resolution
# ---------------------------------------------------------------------------


@pytest.fixture
def fake_binary(tmp_path):
    """Create an executable dummy file and return its path."""
    exe = tmp_path / "packmol"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    return exe


def test_resolve_packmol_via_env(fake_binary, monkeypatch):
    monkeypatch.setenv("MOFFORGE_PACKMOL_BIN", str(fake_binary))
    cfg = popconfig.PopBuildConfig.load()
    assert cfg.resolve_packmol_binary() == fake_binary.resolve()


def test_resolve_packmol_via_kwarg(fake_binary):
    cfg = popconfig.PopBuildConfig.load(packmol_bin=str(fake_binary))
    assert cfg.resolve_packmol_binary() == fake_binary.resolve()


def test_resolve_missing_binary_raises(monkeypatch):
    # Ensure nothing is discoverable: clear env and blank PATH.
    monkeypatch.delenv("MOFFORGE_PACKMOL_BIN", raising=False)
    monkeypatch.delenv("MOFFORGE_LAMMPS_BIN", raising=False)
    monkeypatch.setenv("PATH", "")
    cfg = popconfig.PopBuildConfig.load()
    with pytest.raises(popconfig.ConfigError, match="Packmol executable not found"):
        cfg.resolve_packmol_binary()
    with pytest.raises(popconfig.ConfigError, match="LAMMPS executable not found"):
        cfg.resolve_lammps_binary()


def test_doctor_reports_all_tools(monkeypatch):
    monkeypatch.setenv("PATH", "")
    monkeypatch.delenv("MOFFORGE_PACKMOL_BIN", raising=False)
    monkeypatch.delenv("MOFFORGE_LAMMPS_BIN", raising=False)
    report = popconfig.doctor()
    assert set(report) == {"pysimm", "packmol", "lammps"}
    assert report["packmol"]["available"] is False
