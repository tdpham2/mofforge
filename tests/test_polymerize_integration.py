"""Native execution evidence with real Packmol, required in the dedicated CI job."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from mofforge.polymerize import ConnectionRule, ConstructionState, PopBuilder, connect
from mofforge.polymerize.config import PopBuildConfig, probe_packmol
from tests.test_polymerize_connections import aryl_connectors

pytestmark = pytest.mark.integration


@pytest.fixture(autouse=True)
def real_dependencies():
    required = os.environ.get("MOFFORGE_RUN_POP_TESTS") == "1"
    try:
        import rdkit  # noqa: F401

        binary = PopBuildConfig.load().resolve_packmol_binary()
        version, _ = probe_packmol(binary)
        if required:
            assert version == "21.2.3", f"Required pinned Packmol 21.2.3; found {version}"
    except (ImportError, ValueError, OSError) as exc:
        if required:
            pytest.fail(f"Required native POP integration dependency missing: {exc}")
        pytest.skip(str(exc))


def test_real_single_species_and_same_seed_coordinates(tmp_path):
    results = []
    for i in range(2):
        builder = PopBuilder()
        builder.add_monomer("CC", count=3)
        result = builder.pack(
            output_dir=tmp_path / str(i), box_lengths=[15, 16, 17], random_seed=42
        )
        assert result.success, result.errors
        assert result.metadata["composition"] == {"C": 6, "H": 18}
        results.append(result)
    assert np.array_equal(results[0].state.coordinates, results[1].state.coordinates)
    assert results[0].state.state_hash == results[1].state.state_hash


def test_real_counted_mixture_and_density(tmp_path):
    builder = PopBuilder()
    builder.add_monomer("O", count=3, name="same")
    builder.add_monomer("CCO", count=2, name="same")
    result = builder.pack(
        output_dir=tmp_path,
        initial_packing_density=0.05,
        random_seed=42,
        material_state="solvent-containing",
    )
    assert result.success, result.errors
    assert result.metadata["actual_counts"] == [3, 2]
    assert result.metadata["composition"] == {"C": 4, "H": 18, "O": 5}
    assert result.metadata["current_density"] == pytest.approx(0.05)
    run = result.output_paths[0]
    assert (run / "t0000.xyz").exists() and (run / "t0001.xyz").exists()


def test_real_finite_oligomer_chain_input(tmp_path):
    builder = PopBuilder()
    builder.add_monomer("CCCCCCCC", count=1)
    from mofforge.polymerize.provision import instantiate

    template = instantiate(builder._monomers, builder.prepare(random_seed=42), (30, 30, 30))
    builder = PopBuilder()
    builder.add_monomer(template.crystal, count=3)
    result = builder.pack(output_dir=tmp_path, box_lengths=[30, 30, 30], random_seed=42)
    assert result.success, result.errors
    assert len(result.state.bonds) == 3 * len(template.bonds)
    assert result.metadata["composition"] == {"C": 24, "H": 54}


def test_real_build_save_update_resume(tmp_path):
    builder = PopBuilder()
    builder.add_monomer("c1ccccc1", count=3, connectors=aryl_connectors())
    rule = ConnectionRule(
        "aryl_fixture",
        ("aryl", "aryl"),
        1,
        (1.45, 1.60),
        {"a": ["$replaceable"], "b": ["$replaceable"]},
    )
    first = builder.build(
        output_dir=tmp_path / "first",
        rules=[rule],
        box_lengths=[25, 25, 25],
        random_seed=42,
        target_conversion=2 / 3,
        candidate_attempt_budget=1,
        min_nonbonded_distance=1.5,
    )
    assert first.status == "partial" and first.metadata["new_bonds"] == 1
    state = ConstructionState.load(first.state_path)
    atoms = [
        {
            "id": a.id,
            "species": a.species,
            "coordinates": (pos + np.array([0.1, 0.2, 0.3])).tolist(),
        }
        for a, pos in zip(state.atoms, state.coordinates, strict=True)
    ]
    updated = state.update_geometry(parent_hash=state.state_hash, atoms=atoms, unwrapped=True)
    final = connect(
        updated,
        [rule],
        target_conversion=2 / 3,
        candidate_attempt_budget=30,
        min_nonbonded_distance=1.5,
        output_dir=tmp_path / "resumed",
    )
    assert final.success, final.errors
    assert final.metadata["new_bonds"] == 2
    assert final.metadata["composition"] == {"C": 18, "H": 14}
    assert final.validation.is_valid
    assert ConstructionState.load(final.state_path).state_hash == final.state.state_hash
    assert {a.id for a in final.state.atoms} <= {a.id for a in first.state.atoms}


def test_real_branch(tmp_path):
    builder = PopBuilder()
    builder.add_monomer("c1ccccc1", count=1, connectors=aryl_connectors((0, 2, 4), "center"))
    builder.add_monomer("c1ccccc1", count=3, connectors=aryl_connectors((0,), "leaf"))
    rule = ConnectionRule(
        "branch_fixture",
        ("center", "leaf"),
        1,
        (1.45, 1.60),
        {"a": ["$replaceable"], "b": ["$replaceable"]},
    )
    result = builder.build(
        output_dir=tmp_path,
        rules=[rule],
        box_lengths=[35, 35, 35],
        random_seed=42,
        target_conversion=1,
        candidate_attempt_budget=30,
        min_nonbonded_distance=1.5,
    )
    assert result.success, result.errors
    assert result.metadata["new_bonds"] == 3
    assert result.metadata["composition"] == {"C": 24, "H": 18}
    assert result.metadata["component_sizes"] == [42]


def test_source_file_never_overwritten(tmp_path):
    path = tmp_path / "source.xyz"
    source = "2\nsource!\nH 0 0 0\nH 0.74 0 0\n"
    path.write_text(source)
    graph = {
        "atoms": [{"label": "h0", "species": "H"}, {"label": "h1", "species": "H"}],
        "bonds": [{"atom1": "h0", "atom2": "h1", "order": 1}],
    }
    builder = PopBuilder()
    builder.add_monomer(path, count=2, name="source", graph=graph)
    builder.add_monomer(path, count=1, name="../source", graph=graph)
    result = builder.pack(output_dir=tmp_path, box_lengths=[10, 10, 10], random_seed=42)
    assert result.success, result.errors
    assert path.read_text() == source


def test_impossible_packing_preserves_attempt(tmp_path):
    builder = PopBuilder()
    builder.add_monomer("C", count=8)
    result = builder.pack(
        output_dir=tmp_path, box_lengths=[1, 1, 1], minimum_distance=2, random_seed=42, timeout=0.5
    )
    assert not result.success
    assert result.output_paths
    assert Path(result.output_paths[0], "result.json").exists()
