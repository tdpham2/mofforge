#!/usr/bin/env python3
"""Example 7: Symmetry-Aware Replacement

Demonstrates substructure replacement on a crystal structure while
examining the symmetry properties, and optionally building a supercell.

In mofforge we use pymatgen's SpacegroupAnalyzer for symmetry analysis
and Structure.make_supercell for replication.

Input files:
    crystals/NiPyC_fragment_trouble.cif  - Parent crystal (NiPyC)
    moieties/PyC.xyz                     - Query (PyC linker with H! tag)
    moieties/PyC-CH3.xyz                 - Replacement (PyC with methyl group)

Usage:
    python symmetry_analysis.py
"""

from pathlib import Path

from mofforge import Crystal, infer_bonds, fragment, find_pattern, replace_pattern

SCRIPT_DIR = Path(__file__).parent
CRYSTAL_DIR = SCRIPT_DIR / "data" / "crystals"
MOIETY_DIR = SCRIPT_DIR / "data" / "moieties"


def run_symmetry_example():
    """Demonstrate replacement with symmetry analysis."""

    # -------------------------------------------------------------------------
    # Step 1: Load parent crystal and examine symmetry
    # -------------------------------------------------------------------------
    print("Loading parent: NiPyC_fragment_trouble.cif")
    parent = Crystal.from_cif(CRYSTAL_DIR / "NiPyC_fragment_trouble.cif")
    print(f"  Atoms: {parent.n_atoms}")
    print(
        f"  Lattice: a={parent.lattice.a:.3f}, b={parent.lattice.b:.3f}, c={parent.lattice.c:.3f} A"
    )
    print(
        f"  Angles: alpha={parent.lattice.alpha:.1f}, "
        f"beta={parent.lattice.beta:.1f}, gamma={parent.lattice.gamma:.1f}"
    )

    # Analyze symmetry with pymatgen
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    sga = SpacegroupAnalyzer(parent.structure)
    spacegroup = sga.get_space_group_symbol()
    sg_number = sga.get_space_group_number()
    print(f"  Space group: {spacegroup} (#{sg_number})")

    # Symmetry operations
    sym_ops = sga.get_symmetry_operations()
    print(f"  Symmetry operations: {len(sym_ops)}")

    # -------------------------------------------------------------------------
    # Step 2: Infer bonds
    # -------------------------------------------------------------------------
    print("\nInferring bonds...")
    parent = infer_bonds(parent, periodic=True)
    print(f"  Bonds: {parent.n_bonds}")

    # -------------------------------------------------------------------------
    # Step 3: Load query and replacement fragments
    # -------------------------------------------------------------------------
    print("\nLoading query: PyC.xyz")
    query = fragment("PyC.xyz", fragment_path=str(MOIETY_DIR))
    r_atoms = [s for s in query.species if "!" in s]
    print(f"  Atoms: {query.n_atoms}, R-group: {len(r_atoms)} ({set(r_atoms)})")

    print("Loading replacement: PyC-CH3.xyz")
    replacement = fragment("PyC-CH3.xyz", fragment_path=str(MOIETY_DIR))
    print(f"  Atoms: {replacement.n_atoms}")

    # -------------------------------------------------------------------------
    # Step 4: Search
    # -------------------------------------------------------------------------
    print("\nSearching for PyC linkers...")
    search = find_pattern(query, parent)
    print(f"  Found {search.nb_isomorphisms()} isomorphisms at {search.nb_locations()} locations")

    if search.nb_locations() == 0:
        print("  No matches found. This may be because the CIF structure is")
        print("  in a non-P1 representation. Try a different input CIF.")
        return

    # -------------------------------------------------------------------------
    # Step 5: Replace at 1 location
    # -------------------------------------------------------------------------
    print(f"\nReplacing at 1 location...")
    child = replace_pattern(
        search,
        replacement,
        nb_loc=1,
        name="NiPyC_CH3",
    )
    print(f"  Parent atoms: {parent.n_atoms}")
    print(f"  Child atoms:  {child.n_atoms}")

    # -------------------------------------------------------------------------
    # Step 6: Analyze child symmetry
    # -------------------------------------------------------------------------
    print("\nChild symmetry analysis:")
    try:
        sga_child = SpacegroupAnalyzer(child.structure)
        print(
            f"  Space group: {sga_child.get_space_group_symbol()} "
            f"(#{sga_child.get_space_group_number()})"
        )
    except Exception as e:
        print(f"  Could not determine space group: {e}")

    # -------------------------------------------------------------------------
    # Step 7: Make supercell (2x2x2)
    # -------------------------------------------------------------------------
    print("\nBuilding 2x2x2 supercell...")
    supercell_struct = child.structure.copy()
    supercell_struct.make_supercell([2, 2, 2])
    supercell = Crystal.from_structure(supercell_struct, name="NiPyC_CH3_supercell")
    supercell = infer_bonds(supercell, periodic=True)
    print(f"  Supercell atoms: {supercell.n_atoms}")
    print(f"  Supercell bonds: {supercell.n_bonds}")

    # -------------------------------------------------------------------------
    # Step 8: Write outputs
    # -------------------------------------------------------------------------
    child.write_cif("NiPyC_CH3.cif")
    supercell.write_cif("NiPyC_CH3_supercell.cif")
    print(f"\nOutputs: NiPyC_CH3.cif, NiPyC_CH3_supercell.cif")

    return child, supercell


if __name__ == "__main__":
    run_symmetry_example()
