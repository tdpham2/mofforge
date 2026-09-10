"""Native connection transactions, valence, branching, and periodic closures."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from mofforge.polymerize import ConnectionRule, Connector, PopBuilder, connect
from mofforge.polymerize.connect import _align
from mofforge.polymerize.provision import instantiate


@pytest.fixture
def aryl_rule():
    return ConnectionRule(
        "aryl_fixture",
        ("aryl", "aryl"),
        1,
        (1.45, 1.60),
        {"a": ["$replaceable"], "b": ["$replaceable"]},
        minimum_cycle_size=6,
    )


def aryl_connectors(indices=(0, 3), role="aryl"):
    return [Connector(f"s{i}", f"a{i}", f"h:a{i}:1", role, (f"h:a{i}:1",)) for i in indices]


def prepared_chain(count=3):
    pytest.importorskip("rdkit")
    b = PopBuilder()
    b.add_monomer("c1ccccc1", count=count, connectors=aryl_connectors())
    state = instantiate(b._monomers, b.prepare(random_seed=42), (30, 30, 30))
    for n in range(count):
        state.coordinates[n * 12 : (n + 1) * 12] += [7 * n + 5, 10, 10]
    return state


def test_chain_and_exhaustion_preserve_formula_and_connectors(aryl_rule, tmp_path):
    state = prepared_chain()
    before = state.state_hash
    result = connect(
        state,
        [aryl_rule],
        target_conversion=1,
        candidate_attempt_budget=30,
        min_nonbonded_distance=1.5,
        random_seed=42,
        output_dir=tmp_path,
    )
    assert state.state_hash == before  # Transactional input remains untouched.
    assert result.status == "partial" and not result.success
    assert result.state.statistics["composition"] == {"C": 18, "H": 14}
    assert result.state.statistics["new_bonds"] == 2
    assert result.state.conversion == pytest.approx(2 / 3)
    assert result.validation.is_valid
    assert result.state_path.is_file()
    chain_builder = PopBuilder()
    chain_builder.add_monomer(result.state.crystal, count=2)
    template = chain_builder.prepare()[0]
    assert template.n_atoms == 32
    assert template.molar_mass == pytest.approx(result.state.statistics["mass_g_per_mol"])
    assert len({a.label for a in template.atoms}) == 32
    assert sum(len(e["removed_atoms"]) for e in result.state.events if e["committed"]) == 4
    repeated = connect(
        result.state,
        [aryl_rule],
        target_conversion=1,
        candidate_attempt_budget=1,
        min_nonbonded_distance=1.5,
    )
    assert repeated.state.statistics["new_bonds"] == 2  # No site reuse.


def test_budget_limits_committed_events_and_can_resume(aryl_rule):
    first = connect(
        prepared_chain(),
        [aryl_rule],
        target_conversion=2 / 3,
        candidate_attempt_budget=1,
        min_nonbonded_distance=1.5,
        random_seed=42,
    )
    assert first.status == "partial"
    assert first.metadata["candidate_attempts"] == 1
    second = connect(
        first.state,
        [aryl_rule],
        target_conversion=2 / 3,
        candidate_attempt_budget=10,
        min_nonbonded_distance=1.5,
    )
    assert second.success
    assert second.state.conversion == pytest.approx(2 / 3)


def test_incomplete_edits_fail_valence_without_mutating_state(aryl_rule):
    bad = replace(aryl_rule, delete_atoms={"a": [], "b": []})
    state = prepared_chain(2)
    result = connect(
        state, [bad], target_conversion=0.5, candidate_attempt_budget=1, min_nonbonded_distance=1.5
    )
    assert not result.success
    assert result.state.conversion == 0
    assert len(result.state.atoms) == len(state.atoms)
    assert any("valence" in reason for e in result.state.events for reason in e.get("reasons", []))


def test_aromatic_double_connection_requires_internal_edits(aryl_rule):
    result = connect(
        prepared_chain(2),
        [replace(aryl_rule, bond_order=2)],
        target_conversion=0.5,
        candidate_attempt_budget=1,
        min_nonbonded_distance=1.5,
    )
    assert not result.success and result.state.conversion == 0
    assert any("valence" in reason for e in result.state.events for reason in e.get("reasons", []))


def test_rule_unknown_labels_rejected_before_trials(aryl_rule):
    bad = replace(aryl_rule, delete_atoms={"a": ["nonexistent"], "b": []})
    with pytest.raises(ValueError, match="Unknown template label"):
        connect(
            prepared_chain(),
            [bad],
            target_conversion=1,
            candidate_attempt_budget=1,
            min_nonbonded_distance=1.5,
        )


def test_invalidating_untouched_connector_is_rejected(aryl_rule):
    state = prepared_chain(2)
    # An extra declared site intentionally shares a leaving atom with the selected site.
    extra = replace(state.connectors[0], id="extra", label="extra", role="unused")
    state.connectors.append(extra)
    state.initial_connector_count += 1
    # Both candidate attachment sites on instance 0 share the protected atom.
    for site in state.connectors[:2]:
        site.replaceable = list(dict.fromkeys([*site.replaceable, extra.anchor]))
    result = connect(
        state,
        [aryl_rule],
        target_conversion=1,
        candidate_attempt_budget=8,
        min_nonbonded_distance=1.5,
    )
    assert result.state.conversion == 0
    assert any(
        "untouched connector" in reason
        for event in result.state.events
        for reason in event.get("reasons", [])
    )
    assert result.state.atoms == state.atoms


def test_explicit_bond_order_and_charge_edits_are_recorded():
    pytest.importorskip("rdkit")
    builder = PopBuilder()
    builder.add_monomer(
        "C=O",
        count=2,
        connectors=[Connector("s", "a0", "h:a0:1", "carbonyl", ("h:a0:1",))],
    )
    state = instantiate(builder._monomers, builder.prepare(random_seed=42), (25, 25, 25))
    state.coordinates[:4] += [5, 10, 10]
    state.coordinates[4:] += [15, 10, 10]
    rule = ConnectionRule(
        "explicit_edit_fixture",
        ("carbonyl", "carbonyl"),
        1,
        (1.45, 1.6),
        {"a": ["$replaceable"], "b": ["$replaceable"]},
        bond_changes=tuple(
            {"side": side, "atom1": "a0", "atom2": "a1", "order": 1} for side in ("a", "b")
        ),
        charge_changes=tuple(
            {"side": side, "atom": atom, "charge": charge}
            for side in ("a", "b")
            for atom, charge in (("a0", 1), ("a1", -1))
        ),
    )
    result = connect(
        state,
        [rule],
        target_conversion=1,
        candidate_attempt_budget=1,
        min_nonbonded_distance=1.0,
        random_seed=42,
    )
    assert result.success, result.errors
    assert result.state.statistics["composition"] == {"C": 2, "H": 2, "O": 2}
    assert all(b.order == 1 for b in result.state.bonds)
    assert [a.formal_charge for a in result.state.atoms if a.species == "O"] == [-1, -1]
    event = result.state.events[-1]
    assert len(event["changed_bonds"]) == 2 and len(event["charge_changes"]) == 4
    assert all(e["before"]["order"] == 2 for e in event["changed_bonds"])
    assert all(a.formal_charge == 0 for a in state.atoms)


def test_periodic_closure_is_opt_in_and_uses_winding_topology(aryl_rule):
    state = prepared_chain(1)
    idx = state.index
    left, right = state.connectors
    displacement = state.coordinates[idx[right.atom]] - state.coordinates[idx[left.atom]]
    distance = np.linalg.norm(displacement)
    rotation = _align(displacement / distance, np.array([1.0, 0.0, 0.0]))
    state.coordinates = (state.coordinates - state.coordinates[idx[left.atom]]) @ rotation.T
    state.coordinates += [0.7625, 10, 10]
    state.box_lengths = (distance + 1.525, 20, 20)
    state = state.wrap()
    blocked = connect(
        state,
        [aryl_rule],
        target_conversion=1,
        candidate_attempt_budget=5,
        min_nonbonded_distance=1.0,
    )
    assert not blocked.success and blocked.state.conversion == 0
    result = connect(
        state,
        [aryl_rule],
        target_conversion=1,
        candidate_attempt_budget=5,
        min_nonbonded_distance=1.0,
        allow_cycles=True,
    )
    assert result.success, result.errors
    assert result.state.statistics["composition"] == {"C": 6, "H": 4}
    assert result.state.statistics["periodic_connectivity_ranks"] == [1]
    assert any(b.image != (0, 0, 0) for b in result.state.bonds)
    assert result.state.events[-1]["new_bond"]["image"] != (0, 0, 0)


def test_cycle_minimum_is_required(aryl_rule):
    with pytest.raises(ValueError, match="minimum_cycle_size"):
        connect(
            prepared_chain(),
            [replace(aryl_rule, minimum_cycle_size=None)],
            target_conversion=1,
            candidate_attempt_budget=2,
            min_nonbonded_distance=1.5,
            allow_cycles=True,
        )


def test_failed_bundle_write_retains_valid_connection_result(aryl_rule, monkeypatch, tmp_path):
    from mofforge.core.crystal import Crystal

    monkeypatch.setattr(
        Crystal, "write_cif", lambda *_: (_ for _ in ()).throw(OSError("interrupted"))
    )
    result = connect(
        prepared_chain(2),
        [aryl_rule],
        target_conversion=0.5,
        candidate_attempt_budget=1,
        min_nonbonded_distance=1.5,
        output_dir=tmp_path,
    )
    assert not result.success and result.status == "failed"
    assert result.state.conversion == 0.5 and result.validation.is_valid
    assert "interrupted" in result.errors[0]
    assert result.state_path is None
    assert result.output_paths
    assert not list(tmp_path.rglob("manifest.json"))
