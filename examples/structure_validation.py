#!/usr/bin/env python3
"""Example 10: Structure Validation

Demonstrates post-modification structure validation to catch problems
like steric clashes, unusual bond lengths, and coordination geometry
issues.

Usage:
    python structure_validation.py
    python structure_validation.py path/to/structure.cif
"""

import argparse
from pathlib import Path

from mofforge import Crystal, infer_bonds, validate_structure

SCRIPT_DIR = Path(__file__).parent
CRYSTAL_DIR = SCRIPT_DIR / "data" / "crystals"


def run_validation(crystal_path: str):
    """Validate a crystal structure and report findings."""

    # -------------------------------------------------------------------------
    # Load and prepare the structure
    # -------------------------------------------------------------------------
    print(f"Loading structure: {crystal_path}")
    xtal = Crystal.from_cif(crystal_path)
    xtal = infer_bonds(xtal, periodic=True)

    print(f"  Name:  {xtal.name}")
    print(f"  Atoms: {xtal.n_atoms}")
    print(f"  Bonds: {xtal.n_bonds}")

    # Count species
    from collections import Counter

    species_counts = Counter(xtal.species)
    print(f"  Species: {dict(species_counts)}")
    print()

    # -------------------------------------------------------------------------
    # Run full validation
    # -------------------------------------------------------------------------
    print("Running validation...")
    print("=" * 60)

    report = validate_structure(
        xtal,
        check_clashes=True,
        check_bonds=True,
        check_coordination=True,
        check_charges=False,  # Most MOF CIFs don't have oxidation states
        clash_tolerance=0.5,  # Flag atoms closer than vdW_sum - 0.5 A
        bond_tolerance=0.3,  # Flag bonds deviating >30% from expected
    )

    # -------------------------------------------------------------------------
    # Display results
    # -------------------------------------------------------------------------
    print(report.summary())
    print()

    # -------------------------------------------------------------------------
    # Detailed breakdown
    # -------------------------------------------------------------------------
    if report.steric_clashes:
        print("STERIC CLASHES (atoms too close):")
        for i, j, dist in report.steric_clashes[:10]:
            sp_i = xtal.species[i]
            sp_j = xtal.species[j]
            print(f"  {sp_i}[{i}] -- {sp_j}[{j}]: {dist:.3f} A")
        if len(report.steric_clashes) > 10:
            print(f"  ... and {len(report.steric_clashes) - 10} more")
        print()

    if report.unusual_bonds:
        print("UNUSUAL BOND LENGTHS:")
        for i, j, actual, expected in report.unusual_bonds[:10]:
            sp_i = xtal.species[i]
            sp_j = xtal.species[j]
            deviation = (actual - expected) / expected * 100
            print(
                f"  {sp_i}[{i}] -- {sp_j}[{j}]: {actual:.3f} A "
                f"(expected ~{expected:.3f} A, {deviation:+.1f}%)"
            )
        if len(report.unusual_bonds) > 10:
            print(f"  ... and {len(report.unusual_bonds) - 10} more")
        print()

    if report.coordination_issues:
        print("COORDINATION GEOMETRY ISSUES:")
        for i, sp, cn, expected_range in report.coordination_issues:
            print(f"  {sp}[{i}]: CN={cn} (expected {expected_range[0]}-{expected_range[1]})")
        print()

    # -------------------------------------------------------------------------
    # Overall verdict
    # -------------------------------------------------------------------------
    if report.is_valid:
        print("VERDICT: Structure looks good.")
    else:
        print("VERDICT: Issues found. Review warnings above.")

    return report


def run_all_validations():
    """Validate all crystal structures in the test data directory."""

    cif_files = sorted(CRYSTAL_DIR.glob("*.cif"))
    print(f"Found {len(cif_files)} CIF files to validate\n")

    results = {}
    for cif in cif_files:
        print(f"\n{'#' * 60}")
        print(f"# {cif.stem}")
        print(f"{'#' * 60}")
        try:
            report = run_validation(str(cif))
            results[cif.stem] = report
        except Exception as e:
            print(f"  ERROR: {e}")
            results[cif.stem] = None

    # Summary table
    print(f"\n\n{'=' * 60}")
    print("SUMMARY")
    print(f"{'=' * 60}")
    for name, report in results.items():
        if report is None:
            status = "ERROR"
        elif report.is_valid:
            status = "OK"
        else:
            issues = []
            if report.steric_clashes:
                issues.append(f"{len(report.steric_clashes)} clashes")
            if report.coordination_issues:
                issues.append(f"{len(report.coordination_issues)} coord")
            status = f"ISSUES: {', '.join(issues)}"
        print(f"  {name:35s} {status}")


def main():
    parser = argparse.ArgumentParser(description="Validate crystal structures")
    parser.add_argument(
        "structure",
        nargs="?",
        default=None,
        help="CIF file to validate (omit to validate all test data)",
    )
    args = parser.parse_args()

    if args.structure:
        run_validation(args.structure)
    else:
        run_all_validations()


if __name__ == "__main__":
    main()
