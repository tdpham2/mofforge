#!/usr/bin/env python3
"""Cleanup and Repair — Correct Disorder and Remove Guest Molecules

Demonstrates a two-step structure repair:
    Step 1: Fix disordered (multi-conformation) pyridyl rings in linkers
    Step 2: Remove guest acetylene molecules from pores

This showcases multi-step modification and the `disconnected_component`
mode for finding isolated guest molecules.

Input files:
    data/crystals/SIFSIX-2-Cu-i.cif      - Parent MOF with disorder + guests
    data/moieties/disordered_ligand!.xyz  - Query (disordered ring, masked atoms)
    data/moieties/4-pyridyl.xyz           - Replacement (single-conformation ring)
    data/moieties/acetylene.xyz           - Guest molecule to remove

Usage:
    python repair/cleanup_and_repair.py
"""

from pathlib import Path

from mofforge import (
    Crystal,
    infer_bonds,
    fragment,
    find_pattern,
    replace_pattern,
)

SCRIPT_DIR = Path(__file__).parent
EXAMPLES_DIR = SCRIPT_DIR.parent
CRYSTAL_DIR = EXAMPLES_DIR / "data" / "crystals"
MOIETY_DIR = EXAMPLES_DIR / "data" / "moieties"


def run_cleanup_and_repair():
    """Fix disorder and remove guest molecules from SIFSIX-2-Cu-i."""

    # -------------------------------------------------------------------------
    # Step 1: Load parent crystal
    #
    # This structure has:
    #   - Disordered pyridyl rings (two conformations overlapping)
    #   - Acetylene guest molecules trapped in the pores
    # -------------------------------------------------------------------------
    print("Loading parent: SIFSIX-2-Cu-i.cif")
    parent = Crystal.from_cif(CRYSTAL_DIR / "SIFSIX-2-Cu-i.cif")
    parent = infer_bonds(parent, periodic=True)
    print(f"  Atoms: {parent.n_atoms}, Bonds: {parent.n_bonds}")

    # =========================================================================
    # STEP A: Fix disordered ligands
    #
    # The query contains '!'-tagged atoms representing the alternate
    # conformation to be removed. The replacement is a single clean ring.
    # =========================================================================
    print("\n--- STEP A: Fix Disordered Ligands ---")

    print("Loading query: disordered_ligand!.xyz")
    query_disorder = fragment("disordered_ligand!.xyz", fragment_path=str(MOIETY_DIR))
    r_group = [s for s in query_disorder.species if "!" in s]
    print(f"  Atoms: {query_disorder.n_atoms}")
    print(f"  R-group atoms (to remove): {len(r_group)} atoms -> {set(r_group)}")

    print("Loading replacement: 4-pyridyl.xyz")
    replacement_ring = fragment("4-pyridyl.xyz", fragment_path=str(MOIETY_DIR))
    print(f"  Atoms: {replacement_ring.n_atoms}")

    print("Searching for disordered ligands...")
    search_disorder = find_pattern(query_disorder, parent)
    print(f"  Found {search_disorder.nb_locations()} disordered location(s)")

    if search_disorder.nb_locations() > 0:
        print("Replacing disordered ligands with ordered ones...")
        child = replace_pattern(
            search_disorder,
            replacement_ring,
            name="fixed_disorder",
        )
        print(f"  Atoms after fix: {child.n_atoms} (was {parent.n_atoms})")
    else:
        print("  No disordered ligands found, skipping...")
        child = parent.copy()

    # Re-infer bonds on the corrected structure for next step
    child = infer_bonds(child, periodic=True)

    # =========================================================================
    # STEP B: Remove guest molecules (acetylene)
    #
    # Key: use disconnected_component=True to match only ISOLATED molecules,
    # not acetylene-like fragments within the framework itself.
    #
    # Replace with None to delete the matched molecules.
    # =========================================================================
    print("\n--- STEP B: Remove Guest Molecules ---")

    print("Loading guest query: acetylene.xyz")
    query_guest = fragment("acetylene.xyz", fragment_path=str(MOIETY_DIR))
    print(f"  Atoms: {query_guest.n_atoms}, Species: {query_guest.species}")

    print("Searching for isolated guest molecules (disconnected_component=True)...")
    search_guests = find_pattern(
        query_guest,
        child,
        disconnected_component=True,  # Only match isolated molecules!
    )
    print(f"  Found {search_guests.nb_locations()} guest molecule(s)")

    if search_guests.nb_locations() > 0:
        print("Removing guest molecules (replacing with nothing)...")
        child = replace_pattern(
            search_guests,
            None,  # Replace with nothing = delete
            name="cleaned",
        )
        atoms_removed = parent.n_atoms - child.n_atoms
        print(f"  Atoms after removal: {child.n_atoms}")
    else:
        print("  No guest molecules found.")

    # -------------------------------------------------------------------------
    # Output
    # -------------------------------------------------------------------------
    output = "simulation_ready_SIFSIX-2-Cu-i.cif"
    child.write_cif(output)
    print(f"\n{'=' * 60}")
    print(f"DONE: {parent.n_atoms} atoms -> {child.n_atoms} atoms")
    print(f"Output: {output}")

    return child


if __name__ == "__main__":
    run_cleanup_and_repair()
