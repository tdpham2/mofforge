"""Post-replacement structure validation.

Checks for steric clashes, unusual bond lengths, coordination geometry
issues, and charge balance after crystal modification.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

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
    """Results of structure validation.

    Attributes:
        steric_clashes: List of (atom_i, atom_j, distance) tuples for
            atoms closer than expected.
        unusual_bonds: List of (atom_i, atom_j, actual_dist, expected_dist)
            for bonds with unusual lengths.
        coordination_issues: List of (atom_i, species, coord_number, expected_range)
            for atoms with unexpected coordination numbers.
        charge_balance: Net charge (None if not computed).
        warnings: List of warning messages.
        is_valid: True if no serious issues found.
    """

    steric_clashes: list[tuple[int, int, float]] = field(default_factory=list)
    unusual_bonds: list[tuple[int, int, float, float]] = field(default_factory=list)
    coordination_issues: list[tuple[int, str, int, tuple[int, int]]] = field(default_factory=list)
    charge_balance: float | None = None
    warnings: list[str] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        """True if no steric clashes or coordination issues found."""
        return len(self.steric_clashes) == 0 and len(self.coordination_issues) == 0

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
    """Validate a crystal structure after modification.

    Args:
        crystal: The Crystal to validate.
        check_clashes: If True, check for steric clashes.
        check_bonds: If True, check for unusual bond lengths.
        check_coordination: If True, check metal coordination numbers.
        check_charges: If True, check charge balance.
        clash_tolerance: Distance below sum of vdW radii to flag (Angstroms).
        bond_tolerance: Fractional deviation from expected bond length to flag.

    Returns:
        ValidationReport with all findings.
    """
    report = ValidationReport()

    if crystal.n_atoms == 0:
        report.warnings.append("Crystal has no atoms.")
        return report

    if check_clashes:
        _check_steric_clashes(crystal, report, clash_tolerance)

    if check_bonds:
        _check_bond_lengths(crystal, report, bond_tolerance)

    if check_coordination:
        _check_coordination(crystal, report)

    if check_charges:
        _check_charges(crystal, report)

    logger.info("Validation: %s", "PASSED" if report.is_valid else "ISSUES FOUND")
    return report


def _check_steric_clashes(
    crystal: Crystal,
    report: ValidationReport,
    tolerance: float,
) -> None:
    """Check for atoms that are too close together."""
    species = crystal.species

    # Use pymatgen's neighbor finder for efficiency
    max_vdw = 4.0  # reasonable max vdW diameter
    all_neighbors = crystal.structure.get_all_neighbors(max_vdw)

    for i, neighbors in enumerate(all_neighbors):
        for neighbor in neighbors:
            j = neighbor.index
            if j <= i:
                continue
            dist = neighbor.nn_distance

            # Skip bonded atoms
            if crystal.bonds.has_edge(i, j):
                continue

            try:
                vdw_sum = config.get_vdw_radius(species[i]) + config.get_vdw_radius(species[j])
            except ValueError:
                continue

            if dist < vdw_sum - tolerance:
                report.steric_clashes.append((i, j, dist))


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


def _check_coordination(crystal: Crystal, report: ValidationReport) -> None:
    """Check metal coordination numbers against expected ranges."""
    species = crystal.species

    for i in range(crystal.n_atoms):
        sp = species[i].removesuffix(config.r_tag)
        if sp not in EXPECTED_COORDINATION:
            continue

        cn = crystal.bonds.degree(i)
        expected_range = EXPECTED_COORDINATION[sp]

        if cn < expected_range[0] or cn > expected_range[1]:
            report.coordination_issues.append((i, sp, cn, expected_range))


def _check_charges(crystal: Crystal, report: ValidationReport) -> None:
    """Check charge balance (placeholder — requires oxidation state info)."""
    # pymatgen structures may have oxidation states.
    # Use site.species (Composition) rather than deprecated site.specie.
    try:
        total_charge = 0.0
        for site in crystal.structure:
            # site.species is a Composition; iterate its Species objects
            for sp, amt in site.species.items():
                total_charge += getattr(sp, "oxi_state", 0) * amt
        report.charge_balance = total_charge
        if abs(total_charge) > 0.1:
            report.warnings.append(f"Non-zero net charge: {total_charge:.3f}")
    except (AttributeError, TypeError):
        report.warnings.append("Could not compute charge balance.")
