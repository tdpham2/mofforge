"""Adsorption site identification in MOF structures."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from mofforge.core.crystal import Crystal
from mofforge.utils.config import clean_species as _clean_species, config
from mofforge.validation import EXPECTED_COORDINATION

logger = logging.getLogger("mofforge")

_SITE_TYPES = ("void", "open_metal")


@dataclass
class AdsorptionSite:
    """A candidate site for placing an adsorbate molecule."""

    frac_coords: np.ndarray
    cart_coords: np.ndarray
    site_type: str
    nearest_framework_dist: float
    metadata: dict[str, Any] = field(default_factory=dict)

    def __repr__(self) -> str:
        return (
            f"AdsorptionSite(type={self.site_type!r}, "
            f"dist={self.nearest_framework_dist:.2f} A, "
            f"frac={np.array2string(self.frac_coords, precision=3)})"
        )


def find_adsorption_sites(
    crystal: Crystal,
    strategy: str = "void",
    min_distance: float = 2.5,
    grid_spacing: float = 0.5,
    cluster_tolerance: float = 1.0,
    max_sites: int | None = None,
) -> list[AdsorptionSite]:
    """Identify candidate adsorption sites in a MOF structure."""
    if crystal.n_atoms == 0:
        raise ValueError("Cannot find adsorption sites in an empty crystal.")

    valid = {"void", "open_metal", "both"}
    if strategy not in valid:
        raise ValueError(f"Unknown strategy '{strategy}'. Choose from {valid}.")

    sites: list[AdsorptionSite] = []

    if strategy in ("void", "both"):
        void_sites = _find_void_sites(
            crystal,
            min_distance=min_distance,
            grid_spacing=grid_spacing,
            cluster_tolerance=cluster_tolerance,
        )
        sites.extend(void_sites)

    if strategy in ("open_metal", "both"):
        oms_sites = _find_open_metal_sites(crystal)
        sites.extend(oms_sites)

    if max_sites is not None and len(sites) > max_sites:
        sites = sites[:max_sites]

    logger.debug(
        "Found %d adsorption site(s) in '%s' (strategy=%s)",
        len(sites),
        crystal.name,
        strategy,
    )
    return sites


def _find_void_sites(
    crystal: Crystal,
    min_distance: float = 2.5,
    grid_spacing: float = 0.5,
    cluster_tolerance: float = 1.0,
) -> list[AdsorptionSite]:
    """Find pore-center sites using a 3D grid sampling approach."""
    lattice = crystal.lattice
    frac = crystal.frac_coords

    abc = lattice.abc  # (a, b, c) lengths
    n_a = max(2, int(np.ceil(abc[0] / grid_spacing)))
    n_b = max(2, int(np.ceil(abc[1] / grid_spacing)))
    n_c = max(2, int(np.ceil(abc[2] / grid_spacing)))

    fa = np.linspace(0.0, 1.0, n_a, endpoint=False)
    fb = np.linspace(0.0, 1.0, n_b, endpoint=False)
    fc = np.linspace(0.0, 1.0, n_c, endpoint=False)

    grid_frac = np.array(np.meshgrid(fa, fb, fc, indexing="ij")).reshape(3, -1).T  # (M, 3)
    grid_cart = lattice.get_cartesian_coords(grid_frac)

    # Replicate framework atoms in a 3x3x3 supercell for periodic boundaries
    shifts = np.array(
        [[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
    )  # (27, 3)
    supercell_frac = (frac[np.newaxis, :, :] + shifts[:, np.newaxis, :]).reshape(-1, 3)
    supercell_cart = lattice.get_cartesian_coords(supercell_frac)

    tree = cKDTree(supercell_cart)
    min_dists, _ = tree.query(grid_cart)

    mask = min_dists >= min_distance
    if not np.any(mask):
        logger.debug("No void sites found with min_distance=%.2f A", min_distance)
        return []

    void_frac = grid_frac[mask]
    void_cart = grid_cart[mask]
    void_dists = min_dists[mask]

    sites = _cluster_void_points(void_frac, void_cart, void_dists, lattice, cluster_tolerance)

    sites.sort(key=lambda s: s.nearest_framework_dist, reverse=True)

    logger.debug("Found %d void site(s)", len(sites))
    return sites


def _cluster_void_points(
    frac_coords: np.ndarray,
    cart_coords: np.ndarray,
    min_dists: np.ndarray,
    lattice,
    tolerance: float,
) -> list[AdsorptionSite]:
    """Cluster nearby void points and pick the best representative per cluster."""
    if len(frac_coords) == 0:
        return []

    # Sort by min_dist descending (best points first)
    order = np.argsort(-min_dists)
    frac_sorted = frac_coords[order]
    cart_sorted = cart_coords[order]
    dist_sorted = min_dists[order]

    centers_cart: list[np.ndarray] = []
    centers_frac: list[np.ndarray] = []
    centers_dist: list[float] = []

    for i in range(len(frac_sorted)):
        pt = cart_sorted[i]
        merged = False
        for _j, center in enumerate(centers_cart):
            if np.linalg.norm(pt - center) < tolerance:
                merged = True
                break
        if not merged:
            centers_cart.append(pt)
            centers_frac.append(frac_sorted[i])
            centers_dist.append(float(dist_sorted[i]))

    sites = []
    for fc, cc, d in zip(centers_frac, centers_cart, centers_dist, strict=True):
        sites.append(
            AdsorptionSite(
                frac_coords=fc,
                cart_coords=cc,
                site_type="void",
                nearest_framework_dist=d,
                metadata={"estimated_pore_radius": d},
            )
        )
    return sites


def _find_open_metal_sites(crystal: Crystal) -> list[AdsorptionSite]:
    """Find under-coordinated metal centers and compute placement vectors.

    For each metal whose bond-graph degree is below the expected minimum
    coordination number, the site is placed opposite the average bond
    direction at a distance typical for metal-adsorbate interactions.
    """
    if crystal.n_bonds == 0:
        logger.warning(
            "Crystal '%s' has no bonds; open-metal-site detection "
            "requires inferred bonds. Call infer_bonds() first.",
            crystal.name,
        )
        return []

    species = crystal.species
    cart = crystal.cart_coords
    lattice = crystal.lattice

    sites: list[AdsorptionSite] = []

    for i in range(crystal.n_atoms):
        elem = _clean_species(species[i])
        if elem not in EXPECTED_COORDINATION:
            continue

        cn_min, cn_max = EXPECTED_COORDINATION[elem]
        actual_cn = crystal.bonds.degree(i)

        if actual_cn >= cn_min:
            continue  # not under-coordinated

        # Compute average bond direction vector from this metal
        neighbors = list(crystal.bonds.neighbors(i))
        if not neighbors:
            continue

        metal_cart = cart[i]

        # Get neighbor Cartesian positions (handling PBC via lattice)
        neighbor_vecs = []
        for j in neighbors:
            # Use minimum-image convention
            diff_frac = crystal.frac_coords[j] - crystal.frac_coords[i]
            # Wrap to [-0.5, 0.5)
            diff_frac -= np.round(diff_frac)
            diff_cart = lattice.get_cartesian_coords(diff_frac)
            neighbor_vecs.append(diff_cart)

        avg_bond_dir = np.mean(neighbor_vecs, axis=0)
        avg_bond_dir_norm = np.linalg.norm(avg_bond_dir)

        if avg_bond_dir_norm < 1e-6:
            # Bonds are symmetric around the metal; pick an arbitrary
            # direction perpendicular to the first bond.
            v = neighbor_vecs[0]
            # Find a non-parallel vector
            ref = np.array([1.0, 0.0, 0.0])
            if abs(np.dot(v / np.linalg.norm(v), ref)) > 0.9:
                ref = np.array([0.0, 1.0, 0.0])
            open_dir = np.cross(v, ref)
            open_dir = open_dir / np.linalg.norm(open_dir)
        else:
            # Open site is opposite to average bond direction
            open_dir = -avg_bond_dir / avg_bond_dir_norm

        # Place at a typical metal-adsorbate distance (vdW radius of metal
        # + ~1.5 A for the adsorbate interaction distance)
        try:
            metal_vdw = config.get_vdw_radius(elem)
        except ValueError:
            metal_vdw = 1.5
        placement_dist = metal_vdw + 0.5

        site_cart = metal_cart + open_dir * placement_dist
        site_frac = lattice.get_fractional_coords(site_cart)
        # Wrap to [0, 1)
        site_frac = site_frac % 1.0
        site_cart = lattice.get_cartesian_coords(site_frac)

        # Compute actual min dist to framework
        # (should be ~placement_dist from the metal, but check all atoms)
        diffs = cart - site_cart
        dists_to_framework = np.linalg.norm(diffs, axis=1)
        min_dist = float(np.min(dists_to_framework))

        sites.append(
            AdsorptionSite(
                frac_coords=site_frac,
                cart_coords=site_cart,
                site_type="open_metal",
                nearest_framework_dist=min_dist,
                metadata={
                    "metal_index": i,
                    "metal_species": elem,
                    "actual_cn": actual_cn,
                    "expected_cn_range": (cn_min, cn_max),
                    "missing_cn": cn_min - actual_cn,
                    "placement_distance": placement_dist,
                },
            )
        )

    # Sort by missing CN (most under-coordinated first)
    sites.sort(key=lambda s: s.metadata.get("missing_cn", 0), reverse=True)

    logger.debug("Found %d open-metal site(s)", len(sites))
    return sites
