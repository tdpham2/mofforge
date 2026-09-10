"""Input contracts and engine-independent native state invariants."""

from __future__ import annotations

import subprocess
import sys
from dataclasses import replace

import numpy as np
import pytest

from mofforge.polymerize import Atom, Bond, ConnectionRule, ConstructionState, PopBuilder
from mofforge.polymerize.geometry import validate_box
from mofforge.polymerize.state import finite_cycle_too_small


def two_atoms(*, image=(0, 0, 0), bonded=True):
    atoms = [Atom("a", "t", "m", "a", "C", 12.011, 0), Atom("b", "t", "m", "b", "C", 12.011, 0)]
    return ConstructionState(
        atoms,
        [[0.1, 1, 1], [1.6, 1, 1]],
        (10, 10, 10),
        [Bond("a", "b", 1, image)] if bonded else [],
    )


@pytest.mark.parametrize("count", [0, -1, True, 1.5, "2"])
def test_explicit_positive_counts(count):
    with pytest.raises(ValueError, match="count"):
        PopBuilder().add_monomer("C", count=count)


def test_count_is_required():
    with pytest.raises(TypeError, match="count"):
        PopBuilder().add_monomer("C")


def test_removed_options_are_actionable():
    builder = PopBuilder()
    for key in ("forcefield", "equilibrate", "md_settings", "target_density", "n_monomers"):
        with pytest.raises(TypeError, match="Removed POP option"):
            builder.pack(**{key: 1})
    with pytest.raises(ValueError, match="native"):
        PopBuilder(backend="pysimm")
    with pytest.raises(ValueError, match="rules"):
        builder.build()


def test_functionality_must_match_connectors():
    with pytest.raises(ValueError, match="functionality"):
        PopBuilder().add_monomer("C", count=1, functionality=2)


def test_rule_requires_complete_deletion_declaration():
    with pytest.raises(ValueError, match="delete_atoms"):
        ConnectionRule("r", ("a", "b"), 1, (1.4, 1.6), {"a": []})
    with pytest.raises(ValueError, match="maximum"):
        ConnectionRule("r", ("a", "b"), 1, (1.6, 1.4), {"a": [], "b": []})


def test_native_roundtrip_and_topology_hash(tmp_path):
    state = two_atoms()
    path = state.save(tmp_path)
    loaded = ConstructionState.load(path)
    assert loaded.state_hash == state.state_hash
    assert loaded.to_dict()["components"] == state.component_records
    assert loaded.crystal.structure.site_properties["component_id"] == ["a", "a"]
    with pytest.raises(ValueError, match=r"state\.json"):
        ConstructionState.load(path.parent / "nonexistent.json")
    altered = state.copy()
    altered.bonds[0] = replace(altered.bonds[0], order=2)
    assert altered.state_hash != state.state_hash
    assert altered.crystal.periodic_bonds[0].order == 2
    assert state.crystal.periodic_bonds[0].order == 1
    (path.parent / "box.xyz").write_text("changed")
    with pytest.raises(ValueError, match="hash mismatch"):
        ConstructionState.load(path)


def test_missing_and_interrupted_bundle_rejected(tmp_path, monkeypatch):
    state = two_atoms()
    path = state.save(tmp_path)
    (path.parent / "box.cif").unlink()
    with pytest.raises(OSError):
        ConstructionState.load(path)
    from mofforge.core.crystal import Crystal

    monkeypatch.setattr(
        Crystal, "write_cif", lambda *_: (_ for _ in ()).throw(OSError("interrupted"))
    )
    with pytest.raises(OSError, match="interrupted"):
        state.save(tmp_path / "interrupted")
    bundle = next((tmp_path / "interrupted").iterdir())
    assert not (bundle / "manifest.json").exists()
    with pytest.raises(OSError):
        ConstructionState.load(bundle)


def test_native_versions_and_statistics_are_checked():
    data = two_atoms().to_dict()
    data["schema_version"] = 99
    with pytest.raises(ValueError, match="schema_version"):
        ConstructionState.from_dict(data)
    data = two_atoms().to_dict()
    data["statistics"]["mass_g_per_mol"] += 1
    with pytest.raises(ValueError, match="statistics"):
        ConstructionState.from_dict(data)
    data = two_atoms().to_dict()
    data["components"][0]["id"] = "changed"
    with pytest.raises(ValueError, match="component identities"):
        ConstructionState.from_dict(data)


def test_native_load_without_optional_dependencies(tmp_path):
    path = two_atoms().save(tmp_path)
    program = """
import importlib.abc, sys
class BlockOptional(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'rdkit', 'packmol', 'pysimm', 'mcp'}:
            raise ModuleNotFoundError(fullname)
sys.meta_path.insert(0, BlockOptional())
from mofforge.polymerize import ConstructionState
s = ConstructionState.load(sys.argv[1])
print(s.statistics['mass_g_per_mol'])
"""
    result = subprocess.run(
        [sys.executable, "-c", program, str(path)], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    assert "24.022" in result.stdout


def test_wrap_preserves_bond_images_orders_and_atom_metadata():
    state = two_atoms()
    state.coordinates[1, 0] += 10
    state.bonds = [replace(state.bonds[0], image=(-1, 0, 0), order=2)]
    wrapped = state.wrap()
    assert wrapped.bonds[0].image == (0, 0, 0)
    crystal = state.crystal
    for modified in (crystal.copy(), crystal.wrap(), crystal[[1, 0]], crystal + crystal):
        assert all(b.order == 2 for b in modified.periodic_bonds)
        assert set(modified.structure.site_properties) >= {
            "atom_id",
            "formal_charge",
            "template_label",
        }
    assert np.linalg.norm(wrapped.crystal.bond_vectors(0)[0]) == pytest.approx(1.5)


def test_periodic_components_distinguish_finite_and_percolating():
    state = two_atoms(image=(1, 0, 0))
    assert state.component_info[0][2] == 0  # A boundary-crossing tree is finite.
    state.bonds.append(Bond("a", "b", 1, (0, 0, 0)))
    assert state.component_info[0][2] == 1
    state.bonds.extend([Bond("a", "b", 1, (0, 1, 0)), Bond("a", "b", 1, (0, 0, 1))])
    assert state.component_info[0][2] == 3
    state.check()
    with pytest.raises(ValueError, match="Duplicate periodic"):
        state.bonds.append(Bond("b", "a", 1, (-1, 0, 0)))
        state.check()


def test_minimum_cycle_size_uses_atom_images():
    state = two_atoms()
    assert finite_cycle_too_small(state.atoms, state.bonds, "a", "b", (0, 0, 0), 3)
    assert not finite_cycle_too_small(state.atoms, state.bonds, "a", "b", (1, 0, 0), 3)


def test_construction_rejects_warning_only_contacts():
    state = two_atoms(bonded=False)
    report, _ = validate_box(state, separation=2.0)
    assert not report.is_valid
    assert report.close_contacts


def test_own_periodic_images_and_packing_tolerance():
    state = ConstructionState(
        [Atom("a", "t", "m", "a", "He", 4.0026)], [[0, 0, 0]], (1.98, 10, 10), []
    )
    assert not validate_box(state, separation=2.0, packing=True)[0].is_valid
    state.box_lengths = (1.995, 10, 10)
    assert validate_box(state, separation=2.0, packing=True)[0].is_valid


def test_geometry_mapping_and_cell_updates():
    state = two_atoms()
    state.metadata.update(
        initial_packing_density=0.3, validation={"old": True}, derived_results={"x": 1}
    )
    records = [
        {"id": a.id, "species": a.species, "coordinates": pos.tolist()}
        for a, pos in zip(state.atoms, state.coordinates, strict=True)
    ]
    updated = state.update_geometry(
        parent_hash=state.state_hash,
        atoms=list(reversed(records)),
        box_lengths=[20, 10, 10],
        unwrapped=True,
    )
    assert updated.statistics["current_density"] == pytest.approx(
        state.statistics["current_density"] / 2
    )
    assert updated.metadata["initial_packing_density"] == 0.3
    assert "validation" not in updated.metadata and "derived_results" not in updated.metadata
    assert updated.metadata["geometry_requires_validation"]
    assert updated.bonds == state.bonds
    for malformed in (records[:1], [records[0], records[0]]):
        with pytest.raises(ValueError, match="every atom"):
            state.update_geometry(parent_hash=state.state_hash, atoms=malformed, unwrapped=True)
    with pytest.raises(ValueError, match="parent_hash"):
        state.update_geometry(parent_hash="wrong", atoms=records, unwrapped=True)
    with pytest.raises(ValueError, match="wrapping offsets"):
        state.update_geometry(parent_hash=state.state_hash, atoms=records)
    records[0]["species"] = "N"
    with pytest.raises(ValueError, match="species"):
        state.update_geometry(parent_hash=state.state_hash, atoms=records, unwrapped=True)


def test_explicit_wrapping_offsets_preserve_physical_bonds():
    state = two_atoms()
    coords = state.coordinates + np.array([9, 0, 0])
    records = [
        {
            "id": atom.id,
            "species": atom.species,
            "coordinates": (pos % state.box_lengths).tolist(),
            "wrapping_offset": np.floor(pos / state.box_lengths).astype(int).tolist(),
        }
        for atom, pos in zip(state.atoms, coords, strict=True)
    ]
    updated = state.update_geometry(parent_hash=state.state_hash, atoms=records)
    assert updated.bonds[0].image == (1, 0, 0)
    assert np.linalg.norm(updated.crystal.bond_vectors(0)[0]) == pytest.approx(1.5)


def test_site_annotations_do_not_claim_reaction_recipes():
    from mofforge.polymerize.reactions import available_reactions, available_site_types

    assert available_reactions() == []
    assert any(s["site_type"] == "amine" for s in available_site_types())
