"""Conversion between a pysimm ``System`` and a mofforge ``Crystal``.

The polymerized network comes back from pysimm as a ``System`` (particles +
bonds in a periodic simulation box).  This module lifts it into a mofforge
``Crystal`` — a cubic P1 cell with a NetworkX bond graph — so the rest of
mofforge (validation, provenance, CIF output, rendering) can treat a POP like
any other structure.

pysimm is imported lazily inside the functions, so importing this module never
requires the optional ``pop`` dependency.
"""

from __future__ import annotations

import logging

import networkx as nx
import numpy as np
from pymatgen.core import Lattice, Structure

from mofforge.core.crystal import Crystal

logger = logging.getLogger("mofforge")


def system_to_crystal(system, name: str = "pop") -> Crystal:
    """Convert a pysimm ``System`` into a mofforge ``Crystal``.

    Uses the system's simulation-box dimensions as an orthorhombic P1 lattice,
    the particle coordinates as Cartesian positions, and the system's bonds as
    the connectivity graph.  Particle indices in pysimm are 1-based; the returned
    ``Crystal`` uses 0-based indices throughout.
    """
    particles = list(system.particles)
    if not particles:
        raise ValueError("pysimm System has no particles to convert.")

    # 1-based pysimm particle tags -> 0-based contiguous indices.
    tag_to_idx: dict[int, int] = {}
    species: list[str] = []
    coords = np.empty((len(particles), 3), dtype=np.float64)
    for idx, p in enumerate(particles):
        tag_to_idx[p.tag] = idx
        species.append(_element_of(p))
        coords[idx] = (float(p.x), float(p.y), float(p.z))

    lattice = _lattice_from_dim(system.dim)

    structure = Structure(
        lattice,
        species,
        coords,
        coords_are_cartesian=True,
    )

    bonds = nx.Graph()
    for i, sp in enumerate(species):
        bonds.add_node(i, species=sp)
    for b in system.bonds:
        ai = tag_to_idx.get(b.a.tag)
        bj = tag_to_idx.get(b.b.tag)
        if ai is None or bj is None:  # pragma: no cover - defensive
            continue
        bonds.add_edge(ai, bj)

    crystal = Crystal(name=name, structure=structure, bonds=bonds, species_labels=species)
    return crystal


def _element_of(particle) -> str:
    """Best-effort element symbol for a pysimm particle."""
    getter = getattr(particle, "get_chem_element", None)
    if callable(getter):
        try:
            elem = getter()
            if elem:
                return str(elem)
        except Exception:  # pragma: no cover - fall through to attributes
            logger.debug("get_chem_element failed for particle", exc_info=True)
    if getattr(particle, "elem", None):
        return str(particle.elem)
    ptype = getattr(particle, "type", None)
    if ptype is not None and getattr(ptype, "elem", None):
        return str(ptype.elem)
    raise ValueError(f"Cannot determine element for pysimm particle tag={particle.tag}.")


def _lattice_from_dim(dim) -> Lattice:
    """Build an orthorhombic :class:`Lattice` from a pysimm ``Dimension``."""
    dx = float(getattr(dim, "dx", dim.xhi - dim.xlo))
    dy = float(getattr(dim, "dy", dim.yhi - dim.ylo))
    dz = float(getattr(dim, "dz", dim.zhi - dim.zlo))
    if min(dx, dy, dz) <= 0:
        raise ValueError(f"pysimm System has a non-positive box dimension: ({dx}, {dy}, {dz}).")
    return Lattice.from_parameters(dx, dy, dz, 90.0, 90.0, 90.0)
