#!/usr/bin/env python3
"""Example 2: Linker Functionalization

Demonstrates generating a hypothetical mixed-linker MOF by decorating
selected BDC linkers in IRMOF-1 with acetylamido functional groups.

The query uses a '!'-tagged hydrogen to mark the replacement site.
The replacement fragment has the acetylamido group in place of that hydrogen.

Input files:
    crystals/IRMOF-1.cif                    - Parent MOF
    moieties/2-!-p-phenylene.xyz            - Query (phenylene with H! tag)
    moieties/2-acetylamido-p-phenylene.xyz   - Replacement (phenylene with acetylamido)

Usage:
    python linker_functionalization.py
    python linker_functionalization.py --nb-loc 6 --output my_mof.cif
"""

import argparse
from pathlib import Path

from mofforge import Crystal, infer_bonds, fragment, find_pattern, replace_pattern

SCRIPT_DIR = Path(__file__).parent
CRYSTAL_DIR = SCRIPT_DIR / "data" / "crystals"
MOIETY_DIR = SCRIPT_DIR / "data" / "moieties"


def run_hypothetical_mof(crystal_path: str, nb_loc: int, output: str, fragment_path: str):
    """Build a hypothetical functionalized MOF."""

    # -------------------------------------------------------------------------
    # Step 1: Load parent crystal and infer bonds
    # -------------------------------------------------------------------------
    print(f"Loading parent: {crystal_path}")
    parent = Crystal.from_cif(crystal_path)
    parent = infer_bonds(parent, periodic=True)
    print(f"  Atoms: {parent.n_atoms}, Bonds: {parent.n_bonds}")

    # -------------------------------------------------------------------------
    # Step 2: Load query fragment (phenylene with masked H)
    #
    # The '!' suffix on the H atom in 2-!-p-phenylene.xyz marks it as an
    # R-group atom. During replacement, this H and its corresponding parent
    # atom will be removed and the acetylamido group installed in its place.
    # -------------------------------------------------------------------------
    print("\nLoading query: 2-!-p-phenylene.xyz")
    query = fragment("2-!-p-phenylene.xyz", fragment_path=fragment_path)
    print(f"  Atoms: {query.n_atoms}")
    print(f"  Species: {query.species}")
    r_group = [s for s in query.species if "!" in s]
    print(f"  R-group atoms (will be replaced): {r_group}")

    # -------------------------------------------------------------------------
    # Step 3: Load replacement fragment (phenylene with acetylamido group)
    # -------------------------------------------------------------------------
    print("\nLoading replacement: 2-acetylamido-p-phenylene.xyz")
    replacement = fragment("2-acetylamido-p-phenylene.xyz", fragment_path=fragment_path)
    print(f"  Atoms: {replacement.n_atoms}")
    print(f"  Species: {replacement.species}")

    # -------------------------------------------------------------------------
    # Step 4: Search for query in parent
    # -------------------------------------------------------------------------
    print("\nSearching...")
    search = find_pattern(query, parent)
    print(f"  Found {search.nb_isomorphisms()} isomorphisms at {search.nb_locations()} locations")

    # -------------------------------------------------------------------------
    # Step 5: Replace at selected locations
    #
    # nb_loc=6 means: randomly select 6 of the 24 BDC linker locations.
    # Omitting nb_loc would replace ALL locations.
    # Use loc=[1,5,10] to specify exact locations.
    # -------------------------------------------------------------------------
    print(f"\nReplacing at {nb_loc} random location(s)...")
    child = replace_pattern(
        search,
        replacement,
        nb_loc=nb_loc,
        name="acetylamido_IRMOF-1",
    )
    print(f"  Child atoms: {child.n_atoms} (parent had {parent.n_atoms})")
    print(f"  Atoms added per replacement: {(child.n_atoms - parent.n_atoms) / nb_loc:.0f}")

    # -------------------------------------------------------------------------
    # Step 6: Write output
    # -------------------------------------------------------------------------
    child.write_cif(output)
    print(f"\nOutput written to: {output}")


def main():
    parser = argparse.ArgumentParser(description="Generate a hypothetical functionalized MOF")
    parser.add_argument("--crystal", default=str(CRYSTAL_DIR / "IRMOF-1.cif"))
    parser.add_argument(
        "--nb-loc", type=int, default=6, help="Number of random locations to functionalize"
    )
    parser.add_argument("--output", default="acetylamido_IRMOF-1.cif")
    parser.add_argument("--fragment-path", default=str(MOIETY_DIR))
    args = parser.parse_args()

    run_hypothetical_mof(args.crystal, args.nb_loc, args.output, args.fragment_path)


if __name__ == "__main__":
    main()
