"""Replicate prepared molecular topology without inventing intermolecular bonds."""

from dataclasses import replace

from mofforge.polymerize.state import ConstructionState, Site


def instantiate(monomers, prepared, lengths, *, one_copy=False):
    atoms, coordinates, bonds, sites = [], [], [], []
    for monomer, template in zip(monomers, prepared, strict=True):
        for n in range(1 if one_copy else monomer.count):
            instance = f"{template.template_id}:m{n:06d}"
            mapping = {a.id: f"{instance}:{a.label}" for a in template.atoms}
            atoms.extend(replace(a, id=mapping[a.id], instance_id=instance) for a in template.atoms)
            coordinates.extend(template.coordinates)
            bonds.extend(replace(b, i=mapping[b.i], j=mapping[b.j]) for b in template.bonds)
            sites.extend(
                Site(
                    f"{instance}:site:{s.label}",
                    instance,
                    s.label,
                    mapping[s.atom],
                    mapping[s.anchor],
                    s.role,
                    [mapping[a] for a in s.replaceable],
                )
                for s in monomer.connectors
            )
    return ConstructionState(atoms, coordinates, lengths, bonds, sites, len(sites))
