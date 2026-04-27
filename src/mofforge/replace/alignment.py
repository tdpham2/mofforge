"""SVD-based Procrustes alignment for substructure replacement."""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.linalg import svd

from mofforge.core.crystal import Crystal
from mofforge.replace.conglomerate import reassemble

logger = logging.getLogger("mofforge")


@dataclass
class Alignment:
    """The transformation is applied as: X_new = rotation @ (X + shift_pre) + shift_post."""

    rotation: np.ndarray
    shift_pre: np.ndarray
    shift_post: np.ndarray
    error: float


def get_r2p_alignment(
    replacement: Crystal,
    parent: Crystal,
    r2p: dict[int, int],
    q2p: dict[int, int],
) -> Alignment:
    """Compute the optimal rigid-body alignment of replacement onto parent."""
    n_align = len(r2p)
    if n_align < 2:
        raise ValueError(
            f"Need at least 2 alignment points, got {n_align}."
        )
    if replacement.n_atoms == 0 or parent.n_atoms == 0:
        raise ValueError("Parent and replacement must each have at least 1 atom.")

    # -- Replacement coordinates (atoms involved in alignment) --
    r_indices = [r for r, p in r2p.items()]
    r_sub = replacement[r_indices]
    X_r = r_sub.cart_coords.T  # (3, n_align)
    x_r_center = X_r.mean(axis=1, keepdims=True)
    X_r_centered = X_r - x_r_center

    # -- Parent coordinates (atoms involved in alignment) --
    # Extract parent substructure corresponding to the query match
    p_indices_from_q = [p for q, p in sorted(q2p.items())]
    parent_sub = parent[p_indices_from_q]
    parent_sub = reassemble(parent_sub)

    # Build map: parent atom index -> index in parent_sub
    p2ps = {p: i for i, p in enumerate(p_indices_from_q)}

    # Extract the subset of parent_sub atoms that correspond to alignment atoms
    ps_indices = [p2ps[p] for r, p in r2p.items()]
    parent_align_sub = parent_sub[ps_indices]

    # Get Cartesian coords using parent's lattice
    X_p = parent.lattice.get_cartesian_coords(parent_align_sub.frac_coords).T  # (3, n_align)
    x_p_center = X_p.mean(axis=1, keepdims=True)
    X_p_centered = X_p - x_p_center

    # -- Solve orthogonal Procrustes via SVD --
    H = X_r_centered @ X_p_centered.T  # (3, 3)
    U, _S, Vt = svd(H)
    # Correct for possible reflection: ensure det(R) = +1 (proper rotation)
    d = np.linalg.det(Vt.T @ U.T)
    correction = np.diag([1.0, 1.0, 1.0 if d >= 0 else -1.0])
    rotation = Vt.T @ correction @ U.T

    # Alignment error
    error = float(np.linalg.norm(rotation @ X_r_centered - X_p_centered))

    return Alignment(
        rotation=rotation,
        shift_pre=-x_r_center.flatten(),
        shift_post=x_p_center.flatten(),
        error=error,
    )


def apply_alignment(
    replacement: Crystal,
    parent: Crystal,
    alignment: Alignment,
) -> Crystal:
    """Apply the computed alignment to place replacement into parent's coordinate system."""
    # Get replacement in Cartesian
    cart = replacement.cart_coords.T  # (3, N)

    # Apply transformation: X_new = R @ (X + shift_pre) + shift_post
    cart_aligned = (
        alignment.rotation @ (cart + alignment.shift_pre[:, np.newaxis])
        + alignment.shift_post[:, np.newaxis]
    )

    # Convert to fractional coords in parent's lattice
    frac_aligned = parent.lattice.get_fractional_coords(cart_aligned.T)

    # Build new Crystal with parent's lattice
    # Use clean species for pymatgen, keep original labels separately
    from mofforge.core.crystal import _clean_species

    clean_species = [_clean_species(s) for s in replacement.species]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from pymatgen.core import Structure

        new_structure = Structure(
            parent.lattice,
            clean_species,
            frac_aligned,
        )

    aligned = Crystal(
        name=replacement.name,
        structure=new_structure,
        bonds=replacement.bonds.copy(),
        species_labels=replacement.species,
    )

    return aligned
