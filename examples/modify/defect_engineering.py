#!/usr/bin/env python3
"""Introduce Missing-Linker Defects

Demonstrates engineering missing-linker defects into a MOF by:
    1. Searching for BDC linkers in UiO-66
    2. Replacing selected linkers with formate caps
    3. This effectively removes the phenylene ring core while
       capping the exposed carboxylate groups with formate ions

Input files:
    data/crystals/UiO-66.cif         - Parent MOF
    data/moieties/BDC.xyz            - Query (full BDC linker with !-tagged core)
    data/moieties/formate_caps.xyz   - Replacement (pair of formate ions)

Usage:
    python modify/defect_engineering.py
    python modify/defect_engineering.py --loc 0 5 --output defected.cif
"""

import argparse
from pathlib import Path

from mofforge import Crystal, infer_bonds, fragment, find_pattern, replace_pattern

SCRIPT_DIR = Path(__file__).parent
EXAMPLES_DIR = SCRIPT_DIR.parent
CRYSTAL_DIR = EXAMPLES_DIR / "data" / "crystals"
MOIETY_DIR = EXAMPLES_DIR / "data" / "moieties"


def run_missing_linker_defect(locations: list[int], output: str, fragment_path: str):
    """Introduce missing-linker defects into UiO-66."""

    # -------------------------------------------------------------------------
    # Step 1: Load parent and infer bonds
    # -------------------------------------------------------------------------
    print("Loading parent: UiO-66.cif")
    parent = Crystal.from_cif(CRYSTAL_DIR / "UiO-66.cif")
    parent = infer_bonds(parent, periodic=True)
    print(f"  Atoms: {parent.n_atoms}, Bonds: {parent.n_bonds}")

    # -------------------------------------------------------------------------
    # Step 2: Load query — BDC linker
    #
    # BDC.xyz has '!'-tagged atoms on the p-phenylene core and its hydrogens.
    # These are the atoms that will be removed during replacement.
    # The unmasked atoms (carboxylate oxygens and their carbons) are used
    # for alignment with the formate caps.
    # -------------------------------------------------------------------------
    print("\nLoading query: BDC.xyz")
    query = fragment("BDC.xyz", fragment_path=fragment_path)
    print(f"  Total atoms: {query.n_atoms}")
    r_group = [s for s in query.species if "!" in s]
    unmasked = [s for s in query.species if "!" not in s]
    print(f"  Unmasked atoms (kept for alignment): {len(unmasked)} -> {unmasked}")
    print(f"  R-group atoms (to be removed): {len(r_group)} -> {r_group}")

    # -------------------------------------------------------------------------
    # Step 3: Load replacement — formate caps
    #
    # The formate caps provide two formate ions that fill the gap left by
    # removing the phenylene core of BDC. They bond to the same metal nodes.
    # -------------------------------------------------------------------------
    print("\nLoading replacement: formate_caps.xyz")
    replacement = fragment("formate_caps.xyz", fragment_path=fragment_path)
    print(f"  Atoms: {replacement.n_atoms}, Species: {replacement.species}")

    # -------------------------------------------------------------------------
    # Step 4: Search and examine locations
    # -------------------------------------------------------------------------
    print("\nSearching for BDC linkers...")
    search = find_pattern(query, parent)
    print(f"  Found {search.nb_locations()} BDC linker locations")
    print(f"  Available location indices: 0..{search.nb_locations() - 1}")

    # -------------------------------------------------------------------------
    # Step 5: Replace at specific locations to create defects
    #
    # Choosing specific locations allows strategic placement of defects,
    # e.g. creating a connected channel through the MOF.
    # -------------------------------------------------------------------------
    print(f"\nIntroducing defects at locations: {locations}")
    child = replace_pattern(
        search,
        replacement,
        loc=locations,
        name="defected_UiO-66",
    )

    atoms_removed = parent.n_atoms - child.n_atoms
    print(f"  Parent atoms: {parent.n_atoms}")
    print(f"  Child atoms:  {child.n_atoms}")
    print(f"  Net change:   {atoms_removed:+d} atoms ({len(locations)} defects introduced)")

    # -------------------------------------------------------------------------
    # Step 6: Write output
    # -------------------------------------------------------------------------
    child.write_cif(output)
    print(f"\nOutput written to: {output}")

    return child


def main():
    parser = argparse.ArgumentParser(description="Introduce missing-linker defects into a MOF")
    parser.add_argument(
        "--loc",
        nargs="+",
        type=int,
        default=[2, 8],
        help="Location indices for defect introduction",
    )
    parser.add_argument("--output", default="defected_UiO-66.cif")
    parser.add_argument("--fragment-path", default=str(MOIETY_DIR))
    args = parser.parse_args()

    run_missing_linker_defect(args.loc, args.output, args.fragment_path)


if __name__ == "__main__":
    main()
