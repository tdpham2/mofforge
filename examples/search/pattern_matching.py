#!/usr/bin/env python3
"""Pattern Matching

Demonstrates finding a chemical substructure (p-phenylene linker fragment)
within a MOF crystal structure (IRMOF-1).

This is the most basic operation: search only, no replacement.

Input files:
    data/crystals/IRMOF-1.cif        - Parent MOF crystal structure
    data/moieties/p-phenylene.xyz    - Query fragment (p-phenylene ring)

Usage:
    python search/pattern_matching.py
    python search/pattern_matching.py --crystal my_mof.cif --query my_fragment.xyz
"""

import argparse
from pathlib import Path

from mofforge import Crystal, infer_bonds, fragment, find_pattern

# Default paths (relative to this script, data is at the examples root)
SCRIPT_DIR = Path(__file__).parent
EXAMPLES_DIR = SCRIPT_DIR.parent
CRYSTAL_DIR = EXAMPLES_DIR / "data" / "crystals"
MOIETY_DIR = EXAMPLES_DIR / "data" / "moieties"


def run_search(crystal_path: str, query_name: str, fragment_path: str):
    """Run a substructure search and print results."""

    # -------------------------------------------------------------------------
    # Step 1: Load the parent crystal structure
    # -------------------------------------------------------------------------
    print(f"Loading parent crystal: {crystal_path}")
    parent = Crystal.from_cif(crystal_path)
    print(f"  Atoms: {parent.n_atoms}")
    print(
        f"  Lattice: a={parent.lattice.a:.2f}, b={parent.lattice.b:.2f}, c={parent.lattice.c:.2f} A"
    )

    # -------------------------------------------------------------------------
    # Step 2: Infer bonds (including across periodic boundaries)
    # -------------------------------------------------------------------------
    print("Inferring bonds (periodic=True)...")
    parent = infer_bonds(parent, periodic=True)
    print(f"  Bonds: {parent.n_bonds}")

    # -------------------------------------------------------------------------
    # Step 3: Load the query fragment (fragment to search for)
    # -------------------------------------------------------------------------
    print(f"\nLoading query fragment: {query_name}")
    query = fragment(query_name, fragment_path=fragment_path)
    print(f"  Atoms: {query.n_atoms}")
    print(f"  Species: {query.species}")
    print(f"  Bonds: {query.n_bonds}")

    # -------------------------------------------------------------------------
    # Step 4: Perform pattern search
    # -------------------------------------------------------------------------
    print("\nSearching for pattern...")
    search = find_pattern(query, parent)

    # -------------------------------------------------------------------------
    # Step 5: Examine results
    # -------------------------------------------------------------------------
    print(f"\n{'=' * 60}")
    print(f"SEARCH RESULTS: {search}")
    print(f"{'=' * 60}")
    print(f"  Total isomorphisms: {search.nb_isomorphisms()}")
    print(f"  Unique locations:   {search.nb_locations()}")
    print(f"  Orientations/loc:   {search.nb_ori_at_loc()}")

    # Show the first isomorphism mapping
    if search.nb_locations() > 0:
        first_isom = search.isomorphisms[0][0]
        print(f"\n  First isomorphism (query_atom -> parent_atom):")
        for q_idx, p_idx in sorted(first_isom.items()):
            q_sp = query.species[q_idx]
            p_sp = parent.species[p_idx]
            print(f"    query[{q_idx}] ({q_sp}) -> parent[{p_idx}] ({p_sp})")

    # Extract the matched substructures
    hits = search.matched_substructures()
    print(f"\n  Atoms involved in matches: {hits.n_atoms}")

    return search


def main():
    parser = argparse.ArgumentParser(description="Pattern matching in a crystal")
    parser.add_argument(
        "--crystal", default=str(CRYSTAL_DIR / "IRMOF-1.cif"), help="Path to parent CIF file"
    )
    parser.add_argument("--query", default="p-phenylene.xyz", help="Query fragment XYZ filename")
    parser.add_argument(
        "--fragment-path", default=str(MOIETY_DIR), help="Directory containing fragment XYZ files"
    )
    args = parser.parse_args()

    run_search(args.crystal, args.query, args.fragment_path)


if __name__ == "__main__":
    main()
