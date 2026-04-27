#!/usr/bin/env python3
"""Selective Modification

Demonstrates all five replacement modes available when substituting
functional groups onto MOF linkers.

Modes:
    1. Default: all locations, optimal orientation
    2. N random locations, optimal orientation
    3. Specific locations, optimal orientation
    4. Specific locations AND orientations
    5. Random orientations (skip alignment optimization)

Input files:
    data/crystals/IRMOF-1.cif                - Parent MOF
    data/moieties/2-!-p-phenylene.xyz        - Query (phenylene with H! tag)
    data/moieties/2-nitro-p-phenylene.xyz    - Replacement (phenylene with nitro group)

Usage:
    python modify/selective_modification.py
"""

from pathlib import Path

from mofforge import Crystal, infer_bonds, fragment, find_pattern, replace_pattern

SCRIPT_DIR = Path(__file__).parent
EXAMPLES_DIR = SCRIPT_DIR.parent
CRYSTAL_DIR = EXAMPLES_DIR / "data" / "crystals"
MOIETY_DIR = EXAMPLES_DIR / "data" / "moieties"


def run_replacement_modes():
    """Demonstrate all five replacement modes."""

    # -------------------------------------------------------------------------
    # Setup: load parent, query, replacement and run search
    # -------------------------------------------------------------------------
    print("Loading IRMOF-1...")
    parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
    parent = infer_bonds(parent, periodic=True)
    print(f"  Parent: {parent.n_atoms} atoms, {parent.n_bonds} bonds")

    query = fragment("2-!-p-phenylene.xyz", fragment_path=str(MOIETY_DIR))
    replacement = fragment("2-nitro-p-phenylene.xyz", fragment_path=str(MOIETY_DIR))

    search = find_pattern(query, parent)
    print(f"  Search: {search.nb_isomorphisms()} isomorphisms at {search.nb_locations()} locations")
    print(f"  Orientations per location: {search.nb_ori_at_loc()[:3]}...")
    print()

    # =========================================================================
    # Mode 1: Default — all locations, optimal orientation
    #
    # Replaces at every matched location. At each location, the orientation
    # with the lowest alignment error (best geometric fit) is chosen.
    # =========================================================================
    print("MODE 1: All locations, optimal orientation (default)")
    child1 = replace_pattern(search, replacement, name="mode1_all_optimal")
    print(f"  Result: {child1.n_atoms} atoms")
    print()

    # =========================================================================
    # Mode 2: N random locations, optimal orientation
    #
    # Randomly selects `nb_loc` locations from the available matches.
    # At each selected location, uses optimal orientation.
    # =========================================================================
    print("MODE 2: 8 random locations, optimal orientation")
    child2 = replace_pattern(
        search,
        replacement,
        nb_loc=8,
        name="mode2_8random_optimal",
    )
    print(f"  Result: {child2.n_atoms} atoms")
    print()

    # =========================================================================
    # Mode 3: Specific locations, optimal orientation
    #
    # Replace only at the listed location indices.
    # Location indices are 0-based (unlike Julia's 1-based).
    # =========================================================================
    print("MODE 3: Specific locations [0, 5, 10, 15], optimal orientation")
    child3 = replace_pattern(
        search,
        replacement,
        loc=[0, 5, 10, 15],
        name="mode3_specific_optimal",
    )
    print(f"  Result: {child3.n_atoms} atoms")
    print()

    # =========================================================================
    # Mode 4: Specific locations AND orientations
    #
    # Each location gets a specific orientation index.
    # ori=0 means "use optimal orientation for this location".
    # ori=1,2,3 select specific orientations.
    # =========================================================================
    print("MODE 4: Specific locations [0, 1, 2, 3] with orientations [0, 1, 2, 3]")
    child4 = replace_pattern(
        search,
        replacement,
        loc=[0, 1, 2, 3],
        ori=[0, 1, 2, 3],
        name="mode4_specific_loc_ori",
    )
    print(f"  Result: {child4.n_atoms} atoms")
    print()

    # =========================================================================
    # Mode 5: Random orientations (skip alignment optimization)
    #
    # Replaces at all (or nb_loc) locations, but instead of choosing the
    # orientation with the best geometric fit, picks a random one.
    # Useful for generating diverse hypothetical structures.
    # =========================================================================
    print("MODE 5: All locations, random orientations")
    child5 = replace_pattern(
        search,
        replacement,
        random=True,
        name="mode5_all_random",
    )
    print(f"  Result: {child5.n_atoms} atoms")
    print()

    # =========================================================================
    # Summary
    # =========================================================================
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Parent atoms:            {parent.n_atoms}")
    print(f"  Mode 1 (all, optimal):   {child1.n_atoms}")
    print(f"  Mode 2 (8 rand, optim):  {child2.n_atoms}")
    print(f"  Mode 3 (4 spec, optim):  {child3.n_atoms}")
    print(f"  Mode 4 (spec loc+ori):   {child4.n_atoms}")
    print(f"  Mode 5 (all, random):    {child5.n_atoms}")

    child1.write_cif("child1.cif")

    return {
        "parent": parent,
        "mode1": child1,
        "mode2": child2,
        "mode3": child3,
        "mode4": child4,
        "mode5": child5,
    }


if __name__ == "__main__":
    run_replacement_modes()
