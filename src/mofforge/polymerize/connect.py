"""Bounded, transactional geometric connection construction (no dynamics)."""

from __future__ import annotations

import copy
import time
from dataclasses import asdict, replace

import numpy as np
from scipy.spatial.transform import Rotation

from mofforge.polymerize.base import ConnectionRule, POPResult, integer, positive
from mofforge.polymerize.engine import persist_result
from mofforge.polymerize.geometry import inspect_box, validate_box
from mofforge.polymerize.monomer import check_valence
from mofforge.polymerize.state import Bond, adjacency, finite_cycle_too_small
from mofforge.provenance import derive_seed, effective_seed


def rules_from_input(rules):
    if not rules:
        raise ValueError("Explicit connection rules are required; use pack() for packing alone.")
    result = [ConnectionRule(**r) if isinstance(r, dict) else copy.deepcopy(r) for r in rules]
    if any(not isinstance(r, ConnectionRule) for r in result):
        raise ValueError("rules must contain ConnectionRule objects or their JSON definitions.")
    if len({r.name for r in result}) != len(result):
        raise ValueError("Connection rule names must be unique.")
    return result


def check_options(rules, target_conversion, candidate_attempt_budget, separation, allow_cycles):
    """Validate construction settings before any external packing work."""
    integer(candidate_attempt_budget, "candidate_attempt_budget")
    positive(separation, "min_nonbonded_distance")
    if positive(target_conversion, "target_conversion") > 1:
        raise ValueError("target_conversion must be in (0, 1].")
    if type(allow_cycles) is not bool:
        raise ValueError("allow_cycles must be boolean.")
    if allow_cycles and any(r.minimum_cycle_size is None for r in rules):
        raise ValueError("allow_cycles requires minimum_cycle_size on every rule.")


def _resolve(state, site, label):
    if label == "$atom":
        return site.atom
    if label == "$anchor":
        return site.anchor
    for atom in state.atoms:
        if atom.instance_id == site.instance_id and atom.label == label:
            return atom.id
    raise ValueError(f"Unknown template label {label!r} for connector {site.id}.")


def _deleted(state, site, labels):
    result = []
    for label in labels:
        result.extend(
            site.replaceable if label == "$replaceable" else [_resolve(state, site, label)]
        )
    if len(set(result)) != len(result) or site.atom in result:
        raise ValueError("Deletion lists must be unique and retain the attachment atom.")
    return set(result)


def check_rules(state, rules):
    roles = {s.role for s in state.connectors}
    for rule in rules:
        if not set(rule.roles) <= roles:
            raise ValueError(f"Rule {rule.name!r} refers to undeclared connector roles.")
        for side, role in zip(("a", "b"), rule.roles, strict=True):
            for site in state.connectors:
                if site.role != role or site.consumed:
                    continue
                deleted = _deleted(state, site, rule.delete_atoms[side])
                for edit in rule.bond_changes:
                    if edit["side"] == side:
                        i = _resolve(state, site, edit["atom1"])
                        j = _resolve(state, site, edit["atom2"])
                        matches = [b for b in state.bonds if {b.i, b.j} == {i, j}]
                        if len(matches) != 1 or i in deleted or j in deleted:
                            raise ValueError(
                                "Bond-order edit must identify one retained existing bond."
                            )
                for edit in rule.charge_changes:
                    if edit["side"] == side and _resolve(state, site, edit["atom"]) in deleted:
                        raise ValueError("Cannot change formal charge on a deleted atom.")


def persistent_violations(state, rules, separation):
    """Geometry violations that no connection trial can remove.

    Trials delete only declared atoms, exempt new 1-2/1-3 pairs (each involving
    a connector atom), and move finite components rigidly. A violation persists
    when it avoids those atoms and is rigid: within one finite component's own
    image, within one percolating component, or between percolating components.
    """
    mutable = set()
    for site in state.connectors:
        if site.consumed:
            continue
        mutable.add(site.atom)
        for rule in rules:
            for side, role in zip(("a", "b"), rule.roles, strict=True):
                if role == site.role:
                    mutable |= _deleted(state, site, rule.delete_atoms[side])
    owner = {atom: comp for comp in state.component_info for atom in comp[0]}
    result = []
    for i, j, image, message in inspect_box(state, separation=separation)[2]:
        if {i, j} & mutable:
            continue
        first, second = owner[i], owner[j]
        if first is second:
            internal = tuple(int(x) for x in first[1][j] - first[1][i]) == tuple(image)
            if first[2] or internal:
                result.append(message)
        elif first[2] and second[2]:
            result.append(message)
    return result


def persistent_message(violations, separation):
    return (
        f"{len(violations)} geometry violations cannot be removed by any connection at "
        f"min_nonbonded_distance={separation} angstrom. Only bonded and 1-3 pairs are exempt, "
        "including within a molecule; lower the distance or fix the geometry. "
        + "; ".join(violations[:3])
    )


def _outward(state, site, graph=None):
    graph = adjacency(state.atoms, state.bonds) if graph is None else graph
    images = [image for neighbor, image in graph[site.atom] if neighbor == site.anchor]
    if len(images) != 1:
        raise ValueError("Orientation anchor must be a uniquely bonded neighbor.")
    idx = state.index
    vector = (
        state.coordinates[idx[site.anchor]]
        + images[0] * state.box_lengths
        - state.coordinates[idx[site.atom]]
    )
    norm = np.linalg.norm(vector)
    if norm < 1e-8:
        raise ValueError("Degenerate connector orientation.")
    return vector / norm


def _align(source, target):
    axis = np.cross(source, target)
    sine = np.linalg.norm(axis)
    cosine = np.clip(np.dot(source, target), -1, 1)
    if sine < 1e-10:
        if cosine > 0:
            return np.eye(3)
        basis = np.eye(3)[np.argmin(np.abs(source))]
        axis = np.cross(source, basis)
        return Rotation.from_rotvec(np.pi * axis / np.linalg.norm(axis)).as_matrix()
    return Rotation.from_rotvec(np.arctan2(sine, cosine) * axis / sine).as_matrix()


def _placement_error(state, a, b, rule, vector, graph=None):
    """Return why endpoints at this separation vector cannot bond, or None."""
    distance = np.linalg.norm(vector)
    if not rule.bond_length[0] - 1e-8 <= distance <= rule.bond_length[1] + 1e-8:
        return "Attachment distance outside the rule's bond_length interval."
    for site, target in ((a, vector), (b, -vector)):
        direction = _outward(state, site, graph)
        angle = np.degrees(np.arccos(np.clip(np.dot(direction, target / distance), -1, 1)))
        if angle > rule.angle_tolerance + 1e-6:
            return "Attachment orientation outside the rule's angle_tolerance."
    return None


def _trial(state, a, b, rule, image, component_a, component_b, torsion, separation, graph=None):
    trial = state.copy(history=False)
    idx = trial.index
    direction_a, direction_b = _outward(state, a, graph), _outward(state, b, graph)
    same = component_a is component_b
    moving = None
    if not same:
        if component_b[2] == 0 and (
            component_a[2] != 0 or len(component_b[0]) <= len(component_a[0])
        ):
            moving, moving_site, fixed_site = component_b, b, a
            moving_vector, fixed_vector = direction_b, direction_a
        elif component_a[2] == 0:
            moving, moving_site, fixed_site = component_a, a, b
            moving_vector, fixed_vector = direction_a, direction_b
        if moving is not None:
            members, potentials, _ = moving
            ids = sorted(members)
            positions = np.array(
                [trial.coordinates[idx[i]] + potentials[i] * trial.box_lengths for i in ids]
            )
            origin = (
                trial.coordinates[idx[moving_site.atom]]
                + potentials[moving_site.atom] * trial.box_lengths
            )
            rotation = Rotation.from_rotvec(torsion * fixed_vector).as_matrix() @ _align(
                moving_vector, -fixed_vector
            )
            target = trial.coordinates[idx[fixed_site.atom]] + fixed_vector * np.mean(
                rule.bond_length
            )
            positions = (positions - origin) @ rotation.T + target
            for atom_id, position in zip(ids, positions, strict=True):
                trial.coordinates[idx[atom_id]] = position
            trial.bonds = [
                replace(bond, image=(0, 0, 0)) if bond.i in members else bond
                for bond in trial.bonds
            ]
            image = (0, 0, 0)
    vector = (
        trial.coordinates[idx[b.atom]]
        + np.array(image) * trial.box_lengths
        - trial.coordinates[idx[a.atom]]
    )
    # Moving a component resets its bond images, so only a fixed trial reuses the graph.
    error = _placement_error(trial, a, b, rule, vector, graph if moving is None else None)
    if error:
        raise ValueError(error)
    new_bond = Bond(a.atom, b.atom, rule.bond_order, tuple(image), rule.bond_length)
    if new_bond.key in {bond.key for bond in trial.bonds}:
        raise ValueError("Duplicate periodic bond.")
    if same and finite_cycle_too_small(
        trial.atoms, trial.bonds, a.atom, b.atom, image, rule.minimum_cycle_size
    ):
        raise ValueError("Closure would create a finite ring below minimum_cycle_size.")
    sides = {"a": a, "b": b}
    deleted = _deleted(trial, a, rule.delete_atoms["a"]) | _deleted(
        trial, b, rule.delete_atoms["b"]
    )
    if deleted & {a.atom, b.atom}:
        raise ValueError("Deletion must retain both selected attachment atoms.")
    for site in trial.connectors:
        if site.id in (a.id, b.id):
            site.consumed = True
        elif not site.consumed and deleted & {site.atom, site.anchor, *site.replaceable}:
            raise ValueError("Deletion would invalidate an untouched connector.")
    changed = []
    for edit in rule.bond_changes:
        side = sides[edit["side"]]
        i, j = _resolve(trial, side, edit["atom1"]), _resolve(trial, side, edit["atom2"])
        if {i, j} & deleted:
            raise ValueError("Cannot change a bond removed by the other side's edits.")
        for n, bond in enumerate(trial.bonds):
            if {bond.i, bond.j} == {i, j}:
                trial.bonds[n] = replace(bond, order=edit["order"])
                changed.append({"before": asdict(bond), "after": asdict(trial.bonds[n])})
    charges = []
    for edit in rule.charge_changes:
        atom_id = _resolve(trial, sides[edit["side"]], edit["atom"])
        if atom_id in deleted:
            raise ValueError("Cannot change charge on an atom removed by the other side's edits.")
        atom = trial.atoms[idx[atom_id]]
        trial.atoms[idx[atom_id]] = replace(atom, formal_charge=edit["charge"])
        charges.append({"atom": atom_id, "before": atom.formal_charge, "after": edit["charge"]})
    removed_bonds = [asdict(bond) for bond in trial.bonds if {bond.i, bond.j} & deleted]
    removed_atoms = [
        {**asdict(atom), "coordinates": trial.coordinates[n].tolist()}
        for n, atom in enumerate(trial.atoms)
        if atom.id in deleted
    ]
    keep = [n for n, atom in enumerate(trial.atoms) if atom.id not in deleted]
    trial.atoms = [trial.atoms[n] for n in keep]
    trial.coordinates = trial.coordinates[keep]
    trial.bonds = [bond for bond in trial.bonds if not {bond.i, bond.j} & deleted] + [new_bond]
    trial = trial.wrap()
    check_valence(trial.atoms, trial.bonds)
    report, _ = validate_box(trial, separation=separation)
    if not report.is_valid:
        raise ValueError("; ".join(report.errors[:3]))
    committed_bond = trial.bonds[-1]
    trial.events.append(
        {
            "committed": True,
            "rule": rule.name,
            "sites": [a.id, b.id],
            "new_bond": asdict(committed_bond),
            "changed_bonds": changed,
            "removed_bonds": removed_bonds,
            "removed_atoms": removed_atoms,
            "charge_changes": charges,
        }
    )
    return trial


def connect(
    state,
    rules,
    *,
    target_conversion,
    candidate_attempt_budget,
    min_nonbonded_distance,
    allow_cycles=False,
    random_seed=None,
    output_dir=None,
):
    rules = rules_from_input(rules)
    check_options(
        rules, target_conversion, candidate_attempt_budget, min_nonbonded_distance, allow_cycles
    )
    separation = float(min_nonbonded_distance)
    if type(random_seed) is bool:
        raise ValueError("random_seed must be an integer.")
    seed = effective_seed(state.metadata.get("random_seed") if random_seed is None else random_seed)
    state.check()
    check_rules(state, rules)
    started = time.monotonic()
    current = state.copy()
    # Validate each retained trial after deletions. Leaving atoms can clash at
    # otherwise admissible closure endpoints; they must not veto that trial.
    check_valence(current.atoms, current.bonds)
    attempts = 0
    reason = "No compatible candidates remain."
    persistent = persistent_violations(current, rules, separation)
    if persistent:
        reason = persistent_message(persistent, separation)
    # Pairs rejected while both components were unchanged are retried last.
    rejected = {}
    while (
        not persistent
        and current.conversion + 1e-12 < target_conversion
        and attempts < candidate_attempt_budget
    ):
        idx = current.index
        comps = current.component_info
        owner = {atom: comp for comp in comps for atom in comp[0]}
        identity = {id(comp): (min(comp[0]), len(comp[0])) for comp in comps}
        graph = adjacency(current.atoms, current.bonds)
        active = sorted((s for s in current.connectors if not s.consumed), key=lambda s: s.id)
        candidates = []
        excluded = 0
        for n, first in enumerate(active):
            for second in active[n + 1 :]:
                same = owner[first.atom] is owner[second.atom]
                if same and not allow_cycles:
                    continue
                for rule in rules:
                    if (first.role, second.role) == rule.roles:
                        a, b = first, second
                    elif (second.role, first.role) == rule.roles:
                        a, b = second, first
                    else:
                        continue
                    delta = current.coordinates[idx[b.atom]] - current.coordinates[idx[a.atom]]
                    image = -np.rint(delta / current.box_lengths).astype(int)
                    vector = delta + image * current.box_lengths
                    # Endpoints that cannot move must already fit; screening them
                    # here keeps impossible closures from spending the budget.
                    if same or (owner[a.atom][2] and owner[b.atom][2]):
                        try:
                            misplaced = _placement_error(current, a, b, rule, vector, graph)
                        except ValueError:
                            misplaced = True
                        if misplaced:
                            excluded += 1
                            continue
                    signature = (identity[id(owner[a.atom])], identity[id(owner[b.atom])])
                    candidates.append(
                        (
                            rejected.get((a.id, b.id, rule.name)) == signature,
                            float(np.linalg.norm(vector)),
                            a.id,
                            b.id,
                            rule.name,
                            a,
                            b,
                            rule,
                            tuple(int(x) for x in image),
                            signature,
                        )
                    )
        committed = False
        tried = 0
        for *_, a, b, rule, image, signature in sorted(candidates, key=lambda item: item[:5]):
            if attempts >= candidate_attempt_budget:
                break
            attempts += 1
            tried += 1
            rng = np.random.default_rng(
                derive_seed(seed, [len(current.events), a.id, b.id, rule.name])
            )
            torsions = [0.0, *rng.uniform(0, 2 * np.pi, 11)]
            if owner[a.atom] is owner[b.atom] or (owner[a.atom][2] and owner[b.atom][2]):
                torsions = [0.0]
            rejections = []
            for torsion in torsions:
                try:
                    trial = _trial(
                        current,
                        a,
                        b,
                        rule,
                        image,
                        owner[a.atom],
                        owner[b.atom],
                        torsion,
                        separation,
                        graph,
                    )
                except ValueError as exc:
                    rejections.append(str(exc))
                    continue
                current = trial
                committed = True
                break
            if committed:
                break
            rejected[(a.id, b.id, rule.name)] = signature
            current.events.append(
                {
                    "committed": False,
                    "rule": rule.name,
                    "sites": [a.id, b.id],
                    "trials": len(torsions),
                    "reasons": list(dict.fromkeys(rejections)),
                }
            )
        if not committed:
            if tried:
                reason = f"All {tried} remaining candidates were rejected; see uncommitted events."
            elif excluded:
                reason = (
                    f"No movable candidates remain; {excluded} pairs that cannot move lie "
                    "outside bond_length or angle_tolerance."
                )
            break
    complete = current.conversion + 1e-12 >= target_conversion
    if not complete and attempts >= candidate_attempt_budget:
        reason = "Candidate-attempt budget reached."
    report, checks = validate_box(current, separation=separation)
    current.metadata.update({"construction_separation": separation, "validation": report.to_dict()})
    current.metadata.pop("geometry_requires_validation", None)
    current.provenance.append(
        {
            "operation": "connect",
            "parent_hash": state.state_hash,
            "random_seed": seed,
            "rules": [asdict(r) for r in rules],
            "target_conversion": target_conversion,
            "attempts": attempts,
            "candidate_attempt_budget": candidate_attempt_budget,
            "allow_cycles": allow_cycles,
        }
    )
    result = POPResult(
        complete and report.is_valid,
        crystal=current.crystal,
        state=current,
        operation="connect",
        status=("completed" if complete else "partial") if report.is_valid else "failed",
        validation=report,
        errors=([] if complete else [reason]) + report.errors,
        metadata={
            **current.metadata,
            **current.statistics,
            **checks,
            "target_conversion": target_conversion,
            "candidate_attempts": attempts,
        },
    )
    if output_dir is not None:
        persist_result(result, output_dir)
    result.elapsed_seconds = time.monotonic() - started
    return result
