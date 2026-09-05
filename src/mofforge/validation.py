"""Post-replacement structure validation."""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field

import numpy as np

from mofforge.core.crystal import Crystal
from mofforge.utils.config import config

logger = logging.getLogger("mofforge")

# Expected coordination numbers for common metal centers in MOFs
EXPECTED_COORDINATION: dict[str, tuple[int, int]] = {
    "Zn": (4, 6),
    "Cu": (4, 6),
    "Fe": (4, 6),
    "Co": (4, 6),
    "Ni": (4, 6),
    "Mn": (4, 6),
    "Cr": (4, 6),
    "Ti": (4, 6),
    "Zr": (6, 8),
    "Hf": (6, 8),
    "Al": (4, 6),
    "In": (4, 6),
    "V": (4, 6),
    "Mo": (4, 6),
    "W": (4, 6),
    "Pd": (4, 4),
    "Pt": (4, 4),
    "Ag": (2, 4),
    "Au": (2, 4),
}


@dataclass
class ValidationReport:
    """Results of structure validation."""

    steric_clashes: list[tuple[int, int, float]] = field(default_factory=list)
    unusual_bonds: list[tuple[int, int, float, float]] = field(default_factory=list)
    coordination_issues: list[tuple[int, str, int, tuple[int, int]]] = field(default_factory=list)
    charge_balance: float | None = None
    warnings: list[str] = field(default_factory=list)

    close_contacts: list[tuple[int, int, float]] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    checks_performed: list[str] = field(default_factory=list)
    checks_skipped: dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict:
        """JSON-compatible result shared by all workflow interfaces."""
        return {**asdict(self), "is_valid": self.is_valid}

    @property
    def is_valid(self) -> bool:
        """True when mandatory geometry checks completed without errors."""
        return (
            not self.errors
            and not self.steric_clashes
            and {"geometry", "overlaps"}.issubset(self.checks_performed)
        )

    def summary(self) -> str:
        """Return a human-readable summary of the validation results."""
        lines = [f"Validation Report (valid={self.is_valid}):"]

        if self.steric_clashes:
            lines.append(f"  Steric clashes: {len(self.steric_clashes)}")
            for i, j, d in self.steric_clashes[:5]:
                lines.append(f"    atoms {i}-{j}: {d:.3f} A")
            if len(self.steric_clashes) > 5:
                lines.append(f"    ... and {len(self.steric_clashes) - 5} more")

        if self.unusual_bonds:
            lines.append(f"  Unusual bonds: {len(self.unusual_bonds)}")
            for i, j, actual, expected in self.unusual_bonds[:5]:
                lines.append(f"    atoms {i}-{j}: {actual:.3f} A (expected ~{expected:.3f} A)")
            if len(self.unusual_bonds) > 5:
                lines.append(f"    ... and {len(self.unusual_bonds) - 5} more")

        if self.coordination_issues:
            lines.append(f"  Coordination issues: {len(self.coordination_issues)}")
            for i, sp, cn, expected in self.coordination_issues:
                lines.append(f"    atom {i} ({sp}): CN={cn} (expected {expected})")

        if self.charge_balance is not None:
            lines.append(f"  Charge balance: {self.charge_balance:.3f}")

        for error in self.errors:
            lines.append(f"  ERROR: {error}")
        for check, reason in self.checks_skipped.items():
            lines.append(f"  NOT CHECKED: {check}: {reason}")
        for w in self.warnings:
            lines.append(f"  WARNING: {w}")

        if self.is_valid and not self.unusual_bonds and not self.warnings:
            lines.append("  All checks passed.")

        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.summary()


def validate_structure(
    crystal: Crystal,
    check_clashes: bool = True,
    check_bonds: bool = True,
    check_coordination: bool = True,
    check_charges: bool = False,
    clash_tolerance: float = 0.5,
    bond_tolerance: float = 0.3,
) -> ValidationReport:
    """Validate a crystal structure after modification."""
    if not np.isfinite([clash_tolerance, bond_tolerance]).all():
        raise ValueError("Validation tolerances must be finite.")
    if clash_tolerance < 0 or not 0 <= bond_tolerance < 1:
        raise ValueError("clash_tolerance must be nonnegative and bond_tolerance in [0, 1).")
    report = ValidationReport()
    requested = {
        "clashes": check_clashes,
        "bonds": check_bonds,
        "coordination": check_coordination,
        "charges": check_charges,
    }
    for check, enabled in requested.items():
        if not enabled:
            report.checks_skipped[check] = "Disabled by caller."
    if crystal.n_atoms == 0:
        report.errors.append("Crystal has no atoms.")
        report.warnings.append("Crystal has no atoms.")
        report.checks_skipped["overlaps"] = "No atoms."
        return report
    if (
        not np.isfinite(crystal.frac_coords).all()
        or not np.isfinite(crystal.lattice.matrix).all()
        or abs(np.linalg.det(crystal.lattice.matrix)) < 1e-8
    ):
        report.errors.append("Structure must have finite coordinates and a nonsingular lattice.")
        return report
    report.checks_performed.append("geometry")
    if not crystal.structure.is_ordered:
        report.checks_skipped["overlaps"] = (
            "Resolve partial occupancies before geometric validation."
        )
        report.warnings.append("Disordered structures require an explicit occupancy model.")
        return report

    # Never mutate the caller's graph. Infer it when chemical checks require it.
    from mofforge.core.bonding import infer_bonds

    checked = crystal.copy()
    if (
        checked.n_bonds == 0
        and checked.periodic_bonds is None
        and (check_bonds or check_coordination)
    ):
        checked = infer_bonds(checked)
    checked.refresh_bond_geometry()
    _check_contacts(checked, report, check_clashes, clash_tolerance, bond_tolerance)
    if check_bonds:
        _check_bond_lengths(checked, report, bond_tolerance)
        report.checks_performed.append("bonds")
    if check_coordination:
        _check_coordination(checked, report)
        report.checks_performed.append("coordination")
    if check_charges:
        _check_charges(checked, report)
    return report


def _check_contacts(crystal, report, check_clashes, clash_tolerance, bond_tolerance):
    """Hard overlaps are checked even for bonded or same-molecule atoms."""
    covalent, vdw = {}, {}
    for species in set(crystal.species):
        try:
            covalent[species] = (config.max_bond_distance(species, species) - config.bond_pad) / 2
            vdw[species] = config.get_vdw_radius(species)
        except ValueError:
            report.warnings.append(f"No radius available for {species}.")
    if len(covalent) != len(set(crystal.species)):
        report.checks_skipped["overlaps"] = "Missing covalent radii."
    else:
        report.checks_performed.append("overlaps")
    if check_clashes:
        if len(vdw) != len(set(crystal.species)):
            report.checks_skipped["clashes"] = "Missing van der Waals radii."
        else:
            report.checks_performed.append("clashes")
    cutoff = max(2 * max(covalent.values(), default=0), 2 * max(vdw.values(), default=0))
    if cutoff <= 0:
        return
    centers, neighbors, images, distances = crystal.structure.get_neighbor_list(
        cutoff, exclude_self=False
    )
    species = crystal.species
    molecules = crystal.structure.site_properties.get("mofforge_molecule", [None] * crystal.n_atoms)
    adjacency = {i: set(crystal.bonds.neighbors(i)) for i in range(crystal.n_atoms)}
    bond_images = None
    if crystal.periodic_bonds is not None:
        bond_images = {(b.i, b.j, b.image) for b in crystal.periodic_bonds}
    for raw_i, raw_j, raw_image, raw_dist in zip(
        centers, neighbors, images, distances, strict=True
    ):
        i, j = int(raw_i), int(raw_j)
        image = tuple(int(x) for x in raw_image)
        if j < i or (i == j and image <= (0, 0, 0)):
            continue
        distance = float(raw_dist)
        a, b = species[i], species[j]
        if (
            a in covalent
            and b in covalent
            and distance < (1 - bond_tolerance) * (covalent[a] + covalent[b])
        ):
            report.steric_clashes.append((i, j, distance))
            report.errors.append(
                f"Severe overlap: atoms {i}-{j} at {distance:.3f} A (image {image})."
            )
            continue
        if not check_clashes or a not in vdw or b not in vdw:
            continue
        if bond_images is not None:
            bonded = (i, j, image) in bond_images
        else:
            bonded = crystal.bonds.has_edge(i, j) and image == (0, 0, 0)
        # Molecular internal contacts are not host/guest clashes. Exact overlap
        # checks above still apply, including across periodic images.
        same_molecule = (
            molecules[i] is not None and molecules[i] == molecules[j] and image == (0, 0, 0)
        )
        near_bond = bool(adjacency[i] & adjacency[j]) and i != j
        if bonded or same_molecule or near_bond:
            continue
        if distance < vdw[a] + vdw[b] - clash_tolerance:
            report.close_contacts.append((i, j, distance))
            report.warnings.append(
                f"Nonbonded contact: atoms {i}-{j} at {distance:.3f} A (image {image})."
            )


def _check_bond_lengths(
    crystal: Crystal,
    report: ValidationReport,
    tolerance: float,
) -> None:
    """Check for bonds with unusual lengths."""
    species = crystal.species

    for u, v, data in crystal.bonds.edges(data=True):
        dist = data.get("distance", 0.0)
        if dist <= 0:
            continue

        try:
            expected = config.max_bond_distance(species[u], species[v]) - config.bond_pad
        except ValueError:
            continue

        deviation = abs(dist - expected) / expected if expected > 0 else 0
        if deviation > tolerance:
            report.unusual_bonds.append((u, v, dist, expected))
            report.warnings.append(f"Unusual bond length for atoms {u}-{v}: {dist:.3f} A.")


def _check_coordination(crystal: Crystal, report: ValidationReport) -> None:
    """Check metal coordination numbers against expected ranges."""
    species = crystal.species

    for i in range(crystal.n_atoms):
        sp = species[i].removesuffix(config.r_tag)
        if sp not in EXPECTED_COORDINATION:
            continue

        cn = crystal.coordination_number(i)
        expected_range = EXPECTED_COORDINATION[sp]

        if cn < expected_range[0] or cn > expected_range[1]:
            report.coordination_issues.append((i, sp, cn, expected_range))
            report.warnings.append(
                f"Unusual coordination for atom {i} ({sp}): {cn}; expected {expected_range}."
            )


def _check_charges(crystal: Crystal, report: ValidationReport) -> None:
    """Check charge balance (placeholder — requires oxidation state info)."""
    # pymatgen structures may have oxidation states.
    # Use site.species (Composition) rather than deprecated site.specie.
    try:
        total_charge = 0.0
        for site in crystal.structure:
            # site.species is a Composition; iterate its Species objects
            for sp, amt in site.species.items():
                if not hasattr(sp, "oxi_state"):
                    report.checks_skipped["charges"] = "Oxidation states are missing."
                    report.warnings.append(
                        "Charge balance is unknown: oxidation states are missing."
                    )
                    return
                total_charge += sp.oxi_state * amt
        report.charge_balance = total_charge
        report.checks_performed.append("charges")
        if abs(total_charge) > 0.1:
            report.warnings.append(f"Non-zero net charge: {total_charge:.3f}")
    except (AttributeError, TypeError):
        report.checks_skipped["charges"] = "Could not compute charge balance."
        report.warnings.append("Could not compute charge balance.")
