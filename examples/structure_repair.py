#!/usr/bin/env python3
"""Example 4: Structure Repair — Correct Missing Hydrogen Atoms

Demonstrates repairing a crystal structure that has missing hydrogen atoms
on its linkers — a common artifact of X-ray crystallography.

The query represents the bare carbon ring (no H), and the replacement is
the complete ring with H atoms properly positioned. This is a case where
the replacement is a superset of the query (no '!' masking needed).

Input files:
    crystals/IRMOF-1_noH.cif           - Parent MOF with missing H atoms
    moieties/1,4-C-phenylene_noH.xyz    - Query (bare phenylene, no H)
    moieties/1,4-C-phenylene.xyz        - Replacement (phenylene with H)

Usage:
    python structure_repair.py
    python structure_repair.py --output simulation_ready.cif
"""

import argparse
from pathlib import Path

from mofforge import Crystal, infer_bonds, fragment, find_pattern, replace_pattern

SCRIPT_DIR = Path(__file__).parent
CRYSTAL_DIR = SCRIPT_DIR / "data" / "crystals"
MOIETY_DIR = SCRIPT_DIR / "data" / "moieties"


def run_correct_missing_h(output: str, fragment_path: str):
    """Repair missing hydrogen atoms in IRMOF-1."""

    # -------------------------------------------------------------------------
    # Step 1: Load the corrupted parent (H atoms missing)
    # -------------------------------------------------------------------------
    print("Loading parent with missing H: IRMOF-1_noH.cif")
    parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1_noH.cif")
    parent = infer_bonds(parent, periodic=True)
    print(f"  Atoms: {parent.n_atoms}")
    print(f"  Bonds: {parent.n_bonds}")

    # Count H atoms in the parent
    h_count_before = sum(1 for s in parent.species if s == "H")
    print(f"  Hydrogen atoms: {h_count_before}")

    # -------------------------------------------------------------------------
    # Step 2: Load query — the bare phenylene ring WITHOUT hydrogen atoms
    #
    # Note: no '!' tags are used here. The query is simply the fragment
    # as it appears in the damaged structure.
    # -------------------------------------------------------------------------
    print(f"\nLoading query: 1,4-C-phenylene_noH.xyz")
    query = fragment("1,4-C-phenylene_noH.xyz", fragment_path=fragment_path)
    print(f"  Atoms: {query.n_atoms}, Species: {query.species}")

    # -------------------------------------------------------------------------
    # Step 3: Load replacement — the COMPLETE phenylene ring WITH hydrogens
    #
    # The replacement has more atoms than the query. The extra atoms (H)
    # will be added to the parent structure at each matched location.
    # -------------------------------------------------------------------------
    print(f"\nLoading replacement: 1,4-C-phenylene.xyz")
    replacement = fragment("1,4-C-phenylene.xyz", fragment_path=fragment_path)
    print(f"  Atoms: {replacement.n_atoms}, Species: {replacement.species}")
    print(f"  -> {replacement.n_atoms - query.n_atoms} extra atoms (H) will be added per location")

    # -------------------------------------------------------------------------
    # Step 4: Search and replace at ALL locations (default mode)
    # -------------------------------------------------------------------------
    print("\nSearching...")
    search = find_pattern(query, parent)
    print(f"  Found {search.nb_locations()} locations")

    print("Replacing at all locations with optimal orientation...")
    child = replace_pattern(
        search,
        replacement,
        name="IRMOF-1_repaired",
    )

    h_count_after = sum(1 for s in child.species if s == "H")
    print(f"\n  Parent atoms:  {parent.n_atoms} (H: {h_count_before})")
    print(f"  Child atoms:   {child.n_atoms} (H: {h_count_after})")
    print(f"  H atoms added: {h_count_after - h_count_before}")

    # -------------------------------------------------------------------------
    # Step 5: Write output
    # -------------------------------------------------------------------------
    child.write_cif(output)
    print(f"\nOutput written to: {output}")

    return child


def main():
    parser = argparse.ArgumentParser(description="Repair missing H atoms in a crystal structure")
    parser.add_argument("--output", default="simulation_ready_IRMOF-1.cif")
    parser.add_argument("--fragment-path", default=str(MOIETY_DIR))
    args = parser.parse_args()

    run_correct_missing_h(args.output, args.fragment_path)


if __name__ == "__main__":
    main()
