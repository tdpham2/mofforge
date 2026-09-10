"""Strict geometric checks for packing and transactional connection trials."""

from __future__ import annotations

import numpy as np

from mofforge.polymerize.base import positive
from mofforge.polymerize.state import adjacency
from mofforge.utils.config import config
from mofforge.validation import ValidationReport


def validate_box(state, *, separation: float, packing: bool = False):
    positive(separation, "separation")
    state.check()
    report = ValidationReport(checks_performed=["geometry", "overlaps", "periodic_separation"])
    index = state.index
    for bond in state.bonds:
        distance = np.linalg.norm(
            state.coordinates[index[bond.j]]
            + np.array(bond.image) * state.box_lengths
            - state.coordinates[index[bond.i]]
        )
        if (
            bond.length_range
            and not bond.length_range[0] - 1e-8 <= distance <= bond.length_range[1] + 1e-8
        ):
            report.errors.append(
                f"Bond {bond.i}--{bond.j} length {distance:.6f} outside {bond.length_range}."
            )
    graph = adjacency(state.atoms, state.bonds)
    exempt = {}
    offsets = {}
    for _members, potentials, rank in state.component_info:
        if packing and rank:
            raise ValueError("Packing input contains a periodically percolating component.")
        offsets.update(potentials)
    if not packing:
        for atom in state.atoms:
            pairs = set()
            for j, image in graph[atom.id]:
                pairs.add((j, tuple(image)))
                for k, shift in graph[j]:
                    pairs.add((k, tuple(image + shift)))
            exempt[atom.id] = pairs
    radii = {
        a.species: (config.max_bond_distance(a.species, a.species) - config.bond_pad) / 2
        for a in state.atoms
    }
    cutoff = max(separation, 2 * max(radii.values()))
    structure = state.crystal.structure
    centers, neighbors, images, distances = structure.get_neighbor_list(cutoff, exclude_self=False)
    tolerance = 0.01 if packing else 1e-8
    minimum = None
    for raw_i, raw_j, raw_image, distance in zip(
        centers, neighbors, images, distances, strict=True
    ):
        i, j = int(raw_i), int(raw_j)
        image = tuple(int(x) for x in raw_image)
        if j < i or (i == j and image <= (0, 0, 0)):
            continue
        a, b = state.atoms[i], state.atoms[j]
        if distance < 0.6 * (radii[a.species] + radii[b.species]):
            report.steric_clashes.append((i, j, float(distance)))
            report.errors.append(f"Severe overlap: {a.id}--{b.id} image {image}.")
            continue
        if packing:
            internal = a.instance_id == b.instance_id and image == tuple(
                offsets[b.id] - offsets[a.id]
            )
        else:
            internal = (b.id, image) in exempt[a.id]
        if internal:
            continue
        minimum = float(distance) if minimum is None else min(minimum, float(distance))
        if distance < separation - tolerance:
            report.close_contacts.append((i, j, float(distance)))
            report.errors.append(
                f"Nonbonded separation {a.id}--{b.id} image {image}: "
                f"{distance:.6f} < {separation} angstrom."
            )
    return report, {
        "requested_separation": separation,
        "tolerance_angstrom": tolerance,
        "minimum_contact_within_cutoff": minimum,
        "separation_passed": report.is_valid,
    }
