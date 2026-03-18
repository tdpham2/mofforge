"""Periodic boundary reassembly for cross-PB substructures.

When a matched substructure spans periodic boundaries, its atoms may be
scattered across the unit cell. Reassembly shifts the disconnected
components so they form a contiguous cluster, which is necessary for
correct SVD alignment.
"""

from __future__ import annotations

import logging

import networkx as nx

from mofforge.core.bonding import drop_cross_pb_bonds
from mofforge.core.crystal import Crystal
from mofforge.utils.periodic import nearest_image

logger = logging.getLogger("mofforge")


def reassemble(crystal: Crystal) -> Crystal:
    """Reassemble a crystal substructure split across periodic boundaries.

    Algorithm:
        1. Check connected components before removing cross-PB bonds.
           If >1 components, warn and return (substructure not connected).
        2. Remove cross-PB bonds to identify components split by boundaries.
        3. If only 1 component, already reassembled; return.
        4. Use the largest component as the reference.
        5. Iteratively shift unshifted components:
            - Find cross-PB edges connecting shifted and unshifted atoms.
            - Compute displacement: dx = xf[shifted] - xf[unshifted]
            - Compute nearest image: n_dx = nearest_image(dx)
            - Shift all atoms in the unshifted component by (dx - n_dx).
        6. Repeat until all components are shifted.

    Args:
        crystal: The Crystal substructure to reassemble.

    Returns:
        A new Crystal with shifted coordinates forming a contiguous cluster.
    """
    xtal = crystal.copy()

    if xtal.n_atoms == 0:
        return xtal

    # Check if the structure is connected before removing cross-PB bonds
    if not nx.is_connected(xtal.bonds) and xtal.bonds.number_of_edges() > 0:
        n_comp = nx.number_connected_components(xtal.bonds)
        logger.warning(
            "# connected components in parent substructure = %d > 1. "
            "Assuming the substructure does not cross the periodic boundary...",
            n_comp,
        )
        return xtal

    # Remove cross-PB bonds to find components split across boundaries
    bonds_no_pb = drop_cross_pb_bonds(xtal.bonds)
    conn_comps = list(nx.connected_components(bonds_no_pb))

    # If only one component, already reassembled
    if len(conn_comps) == 1:
        return xtal

    # Track which components have been shifted
    comp_shifted = [False] * len(conn_comps)
    # Reference component = largest one
    ref_comp_id = max(range(len(conn_comps)), key=lambda i: len(conn_comps[i]))
    comp_shifted[ref_comp_id] = True

    # Map each atom to its component index
    atom_to_comp: dict[int, int] = {}
    for comp_id, comp_nodes in enumerate(conn_comps):
        for node in comp_nodes:
            atom_to_comp[node] = comp_id

    # Get mutable fractional coordinates
    frac_coords = xtal.frac_coords  # (N, 3)

    def is_shifted(atom: int) -> bool:
        return comp_shifted[atom_to_comp[atom]]

    # Iterate until all components are shifted
    max_iterations = len(conn_comps) * 2  # safety bound
    iteration = 0
    while not all(comp_shifted):
        iteration += 1
        if iteration > max_iterations:
            logger.warning(
                "Reassemble: max iterations reached. Some components may not be shifted."
            )
            break

        # Loop over cross-PB edges in the original bond graph
        for u, v, data in xtal.bonds.edges(data=True):
            if not data.get("cross_boundary", False):
                continue

            # Check if one atom is shifted and the other is not
            u_shifted = is_shifted(u)
            v_shifted = is_shifted(v)

            if u_shifted and not v_shifted:
                p_ref, p_unshifted = u, v
            elif v_shifted and not u_shifted:
                p_ref, p_unshifted = v, u
            else:
                continue  # both shifted or both unshifted

            # Component to shift
            comp_id = atom_to_comp[p_unshifted]

            # Compute displacement
            dx = frac_coords[p_ref] - frac_coords[p_unshifted]
            n_dx = nearest_image(dx)
            shift = dx - n_dx

            # Shift all atoms in this component
            for atom_idx in conn_comps[comp_id]:
                frac_coords[atom_idx] += shift

            comp_shifted[comp_id] = True

    # Update coordinates in the crystal
    xtal.set_frac_coords(frac_coords)

    return xtal
