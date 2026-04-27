#!/usr/bin/env python3
"""SMARTS-like Pattern Matching

Demonstrates searching for substructures using string-based chemical
patterns instead of XYZ file input.

Supported syntax:
    - Element symbols: C, N, O, Zn, etc.
    - Bonds: '-' (connector between atoms)
    - Brackets: [Zn], [Cu] for explicit atoms
    - Wildcards: [*] for any atom
    - Ring closure: C1-C-C-C-C-C-1 (digit marks ring start/end)

Usage:
    python search/string_pattern_search.py
"""

from pathlib import Path

from mofforge import Crystal, infer_bonds
from mofforge.smarts import parse_smarts, smarts_search

SCRIPT_DIR = Path(__file__).parent
EXAMPLES_DIR = SCRIPT_DIR.parent
CRYSTAL_DIR = EXAMPLES_DIR / "data" / "crystals"


def run_smarts_examples():
    """Demonstrate various SMARTS pattern searches."""

    # -------------------------------------------------------------------------
    # Load a MOF structure
    # -------------------------------------------------------------------------
    print("Loading IRMOF-1...")
    parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
    parent = infer_bonds(parent, periodic=True)
    print(f"  Atoms: {parent.n_atoms}, Bonds: {parent.n_bonds}")
    print()

    # =========================================================================
    # Pattern 1: Simple metal-ligand bond
    # =========================================================================
    pattern = "[Zn]-[O]"
    print(f"Pattern: '{pattern}'  (Zn-O bonds)")
    graph = parse_smarts(pattern)
    print(f"  Parsed: {graph.number_of_nodes()} atoms, {graph.number_of_edges()} bonds")

    result = smarts_search(pattern, parent)
    print(
        f"  Matches: {result.nb_isomorphisms()} isomorphisms at {result.nb_locations()} locations"
    )
    print()

    # =========================================================================
    # Pattern 2: Metal-ligand-carbon chain
    # =========================================================================
    pattern = "[Zn]-[O]-C"
    print(f"Pattern: '{pattern}'  (Zn-O-C chain)")
    result = smarts_search(pattern, parent)
    print(
        f"  Matches: {result.nb_isomorphisms()} isomorphisms at {result.nb_locations()} locations"
    )
    print()

    # =========================================================================
    # Pattern 3: Carboxylate group (O-C-O)
    # =========================================================================
    pattern = "O-C-O"
    print(f"Pattern: '{pattern}'  (carboxylate O-C-O)")
    result = smarts_search(pattern, parent)
    print(
        f"  Matches: {result.nb_isomorphisms()} isomorphisms at {result.nb_locations()} locations"
    )
    print()

    # =========================================================================
    # Pattern 4: Wildcard — any atom bonded to Zn
    # =========================================================================
    pattern = "[Zn]-[*]"
    print(f"Pattern: '{pattern}'  (anything bonded to Zn)")
    result = smarts_search(pattern, parent)
    print(
        f"  Matches: {result.nb_isomorphisms()} isomorphisms at {result.nb_locations()} locations"
    )
    print()

    # =========================================================================
    # Pattern 5: 6-membered carbon ring (benzene)
    # =========================================================================
    pattern = "C1-C-C-C-C-C-1"
    print(f"Pattern: '{pattern}'  (6-membered C ring)")
    result = smarts_search(pattern, parent)
    print(
        f"  Matches: {result.nb_isomorphisms()} isomorphisms at {result.nb_locations()} locations"
    )
    print()

    # =========================================================================
    # Pattern 6: C-H bond
    # =========================================================================
    pattern = "C-H"
    print(f"Pattern: '{pattern}'  (C-H bonds)")
    result = smarts_search(pattern, parent)
    print(
        f"  Matches: {result.nb_isomorphisms()} isomorphisms at {result.nb_locations()} locations"
    )
    print()

    # =========================================================================
    # Try on a different MOF
    # =========================================================================
    print("=" * 60)
    print("Loading Ti-MIL-125...")
    parent2 = Crystal.from_cif(CRYSTAL_DIR / "Ti-MIL-125.cif")
    parent2 = infer_bonds(parent2, periodic=True)
    print(f"  Atoms: {parent2.n_atoms}")

    pattern = "[Ti]-[O]"
    print(f"\nPattern: '{pattern}'  (Ti-O bonds in Ti-MIL-125)")
    result = smarts_search(pattern, parent2)
    print(
        f"  Matches: {result.nb_isomorphisms()} isomorphisms at {result.nb_locations()} locations"
    )


if __name__ == "__main__":
    run_smarts_examples()
