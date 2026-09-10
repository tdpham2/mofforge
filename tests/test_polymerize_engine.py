"""Molecular input preparation and Packmol failure handling."""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pytest

from mofforge.polymerize import PopBuilder, engine
from mofforge.polymerize.config import ConfigError, PopBuildConfig, doctor, probe_packmol


@pytest.fixture
def builder():
    pytest.importorskip("rdkit")
    b = PopBuilder()
    b.add_monomer("CC", count=2)
    return b


def test_smiles_preparation_counts_mass_labels_and_seed(builder):
    a, b = builder.prepare(random_seed=42)[0], builder.prepare(random_seed=42)[0]
    assert a.n_atoms == 8 and a.molar_mass == pytest.approx(30.07)
    assert np.array_equal(a.coordinates, b.coordinates)
    assert "h:a0:1" in [atom.label for atom in a.atoms]
    assert a.provenance["embedding_seed"] == b.provenance["embedding_seed"]


def test_cleanup_must_converge(builder, monkeypatch):
    from rdkit.Chem import AllChem

    monkeypatch.setattr(AllChem, "UFFOptimizeMolecule", lambda *args, **kwargs: 1)
    with pytest.raises(ValueError, match="did not converge"):
        builder.prepare(random_seed=42)


def test_coordinate_input_preserves_geometry_and_hydrogen_count(tmp_path):
    Chem = pytest.importorskip("rdkit.Chem")
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    path = tmp_path / "ethane.mol"
    Chem.MolToMolFile(mol, str(path))
    original = path.read_bytes()
    b = PopBuilder()
    b.add_monomer(path, count=2, name="../unsafe")
    prepared = b.prepare()[0]
    assert prepared.n_atoms == 8
    assert np.allclose(prepared.coordinates, mol.GetConformer().GetPositions(), atol=0.0001)
    assert path.read_bytes() == original
    Chem.MolToMolFile(Chem.RemoveHs(mol), str(path))
    with pytest.raises(ValueError, match="explicit hydrogen"):
        b.prepare()


def test_xyz_needs_explicit_graph_and_preserves_source(tmp_path):
    pytest.importorskip("rdkit")
    path = tmp_path / "hydrogen.xyz"
    original = "2\nsource\nH 0 0 0\nH 0.74 0 0\n"
    path.write_text(original)
    b = PopBuilder()
    b.add_monomer(path, count=1)
    with pytest.raises(ValueError, match="graph"):
        b.prepare()
    b = PopBuilder()
    b.add_monomer(
        path,
        count=1,
        graph={
            "atoms": [{"label": "h0", "species": "H"}, {"label": "h1", "species": "H"}],
            "bonds": [{"atom1": "h0", "atom2": "h1", "order": 1}],
        },
    )
    assert b.prepare()[0].n_atoms == 2
    assert path.read_text() == original


def test_crystal_oligomer_input_unwraps_finite_topology(builder):
    from mofforge.polymerize.provision import instantiate

    template = builder.prepare(random_seed=42)
    state = instantiate(builder._monomers, template, (10, 10, 10), one_copy=True).wrap()
    b = PopBuilder()
    b.add_monomer(state.crystal, count=3)
    prepared = b.prepare()[0]
    assert prepared.n_atoms == 8
    assert prepared.molar_mass == pytest.approx(template[0].molar_mass)
    assert np.allclose(
        np.sort(prepared.coordinates[:, 0] - prepared.coordinates[0, 0]),
        np.sort(template[0].coordinates[:, 0] - template[0].coordinates[0, 0]),
    )


@pytest.mark.parametrize(
    "options",
    [
        {},
        {"box_lengths": [10, 10, 10], "initial_packing_density": 0.3},
        {"box_lengths": [10, -1, 10]},
        {"initial_packing_density": float("nan")},
        {"box_lengths": [10, 10, 10], "random_seed": True},
    ],
)
def test_bad_packing_options_fail_before_execution(builder, monkeypatch, tmp_path, options):
    monkeypatch.setattr(engine, "probe_packmol", lambda *_: pytest.fail("engine must not run"))
    with pytest.raises(ValueError):
        builder.pack(output_dir=tmp_path, **options)


def mock_engine(monkeypatch, runner):
    monkeypatch.setattr(PopBuildConfig, "resolve_packmol_binary", lambda _: Path("/fake/packmol"))
    monkeypatch.setattr(engine, "probe_packmol", lambda *_: ("21.2.3", "Version 21.2.3"))
    monkeypatch.setattr(engine.subprocess, "run", runner)


def test_timeout_preserves_logs_and_input(builder, monkeypatch, tmp_path):
    def timeout(*args, **kwargs):
        assert kwargs["stdin"].seekable()
        raise subprocess.TimeoutExpired(
            args[0], 0.1, output=b"partial stdout", stderr=b"partial stderr"
        )

    mock_engine(monkeypatch, timeout)
    result = builder.pack(output_dir=tmp_path, box_lengths=[10, 10, 10], timeout=0.1)
    assert not result.success and result.status == "failed"
    run = result.output_paths[0]
    assert (run / "stdout.log").read_text() == "partial stdout"
    assert (run / "stderr.log").read_text() == "partial stderr"
    assert (run / "packmol.inp").is_file()
    assert (run / "result.json").is_file()
    assert "timed out" in result.errors[0]
    assert result.metadata["random_seed"] is not None
    assert result.metadata["requested_counts"] == [2]


def test_malformed_zero_exit_rejected(builder, monkeypatch, tmp_path):
    def malformed(*args, **kwargs):
        (kwargs["cwd"] / "packed.xyz").write_text("1\nwrong\nHe 0 0 0\n")
        return subprocess.CompletedProcess(args[0], 0, "pretend success", "")

    mock_engine(monkeypatch, malformed)
    result = builder.pack(output_dir=tmp_path, box_lengths=[10, 10, 10])
    assert not result.success
    assert "atom count" in result.errors[0]
    assert not list(tmp_path.rglob("manifest.json"))


def test_configured_invalid_binary_does_not_fall_back(monkeypatch):
    monkeypatch.setenv("MOFFORGE_PACKMOL_BIN", "/missing/packmol")
    with pytest.raises(ConfigError, match="Configured Packmol"):
        PopBuildConfig().resolve_packmol_binary()


def test_omitted_builder_override_preserves_toml_configuration(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("MOFFORGE_LAMMPS_BIN", raising=False)
    (tmp_path / "mofforge.toml").write_text('[backends.pop]\npackmol_bin = "/configured/packmol"\n')
    assert PopBuildConfig.load(packmol_bin=None).packmol_bin == "/configured/packmol"
    assert PopBuildConfig.load(packmol_bin="/explicit/packmol").packmol_bin == "/explicit/packmol"


def test_version_feature_gate(monkeypatch):
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *a, **kw: subprocess.CompletedProcess(a[0], 0, "Version 20.14.0", ""),
    )
    with pytest.raises(ConfigError, match=r"20\.15\.0"):
        probe_packmol("/fake")


def test_doctor_reports_rdkit_and_packmol(monkeypatch):
    monkeypatch.setenv("MOFFORGE_PACKMOL_BIN", "/missing/packmol")
    monkeypatch.delenv("MOFFORGE_LAMMPS_BIN", raising=False)
    report = doctor()
    assert set(report) == {"rdkit", "packmol"}
    assert not report["packmol"]["available"]


@pytest.mark.parametrize(
    "bad",
    [
        {"target_conversion": 1.2},
        {"target_conversion": "0.5"},
        {"candidate_attempt_budget": 0},
        {"min_nonbonded_distance": -1},
        {"allow_cycles": "false"},
        {"allow_cycles": True},
    ],
)
def test_bad_connection_settings_fail_before_packing(builder, monkeypatch, bad):
    from mofforge.polymerize import ConnectionRule

    monkeypatch.setattr(engine, "pack", lambda *a, **kw: pytest.fail("must not pack"))
    options = dict(target_conversion=1, candidate_attempt_budget=1, min_nonbonded_distance=1.5)
    options.update(bad)
    with pytest.raises(ValueError):
        builder.build(
            rules=[ConnectionRule("r", ("a", "b"), 1, (1.4, 1.6), {"a": [], "b": []})],
            box_lengths=[20, 20, 20],
            **options,
        )


def test_crystal_preserves_explicit_bond_orders_and_isotope_masses():
    pytest.importorskip("rdkit")
    from dataclasses import replace

    from mofforge.polymerize.provision import instantiate

    builder = PopBuilder()
    builder.add_monomer("[2H]c1ccccc1", count=1)
    original = instantiate(builder._monomers, builder.prepare(random_seed=42), (20, 20, 20))
    # Explicit Kekule orders must survive RDKit's aromaticity perception.
    aromatic = [i for i, b in enumerate(original.bonds) if b.order == 1.5]
    for n, i in enumerate(aromatic):
        original.bonds[i] = replace(original.bonds[i], order=1 if n % 2 else 2)
    builder = PopBuilder()
    builder.add_monomer(original.crystal, count=1)
    prepared = builder.prepare()[0]
    assert prepared.molar_mass == pytest.approx(original.statistics["mass_g_per_mol"])
    assert [b.order for b in prepared.bonds] == [b.order for b in original.bonds]
