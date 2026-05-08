#!/usr/bin/env python3
"""Desolvate MOF structures — remove all uncoordinated solvent molecules.

Usage:
    # Single CIF file
    python scripts/desolvate.py structure.cif

    # Multiple CIF files
    python scripts/desolvate.py *.cif

    # With options
    python scripts/desolvate.py -o output_dir/ --keep-metals --verbose structures/*.cif

    # Using a text file listing CIF paths (one per line)
    python scripts/desolvate.py --from-list cif_list.txt
"""

from __future__ import annotations

import argparse
import sys
import warnings
from collections import Counter
from pathlib import Path

from mofforge import Crystal, infer_bonds
from mofforge.solvent.removal import remove_solvent


def _load_cif(cif_path: Path) -> Crystal:
    """Load a CIF file with relaxed parsing for CSD-style CIFs."""
    from pymatgen.io.cif import CifParser

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parser = CifParser(str(cif_path), occupancy_tolerance=100)
        structures = parser.parse_structures(primitive=False)

    if not structures:
        raise ValueError(f"No structures found in CIF file: {cif_path}")

    return Crystal.from_structure(structures[0], name=cif_path.stem)


def process_one(cif_path: Path, output_dir: Path, args: argparse.Namespace) -> bool:
    """Desolvate a single CIF file. Returns True on success."""
    try:
        xtal = _load_cif(cif_path)
        xtal = infer_bonds(xtal, periodic=True)

        result = remove_solvent(
            xtal,
            min_atoms=args.min_atoms,
            keep_metal_containing=args.keep_metals,
            n_framework_components=args.n_frameworks,
        )

        out_path = output_dir / f"{cif_path.stem}_desolvated.cif"
        result.crystal.write_cif(out_path)

        print(f"  {cif_path.name}: {xtal.n_atoms} -> {result.crystal.n_atoms} atoms "
              f"({result.n_components_removed} molecules removed)")

        if args.verbose and result.removed_molecules:
            formula_counts = Counter(m.formula for m in result.removed_molecules)
            for formula, count in formula_counts.most_common():
                print(f"    {count}x {formula}")

        return True

    except Exception as e:
        print(f"  {cif_path.name}: FAILED — {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Remove uncoordinated solvent molecules from MOF CIF files."
    )
    parser.add_argument(
        "cif_files",
        nargs="*",
        help="CIF file(s) to desolvate.",
    )
    parser.add_argument(
        "--from-list",
        metavar="FILE",
        help="Text file with one CIF path per line.",
    )
    parser.add_argument(
        "-o", "--output-dir",
        default=".",
        help="Output directory for desolvated CIF files (default: current dir).",
    )
    parser.add_argument(
        "--min-atoms",
        type=int,
        default=1,
        help="Keep components with at least this many atoms (default: 1).",
    )
    parser.add_argument(
        "--keep-metals",
        action="store_true",
        help="Do not remove components containing metal atoms.",
    )
    parser.add_argument(
        "--n-frameworks",
        type=int,
        default=None,
        help="Number of framework components to keep (default: auto-detect).",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print detailed removal info.",
    )

    args = parser.parse_args()

    # Collect CIF paths
    cif_paths: list[Path] = []
    for f in args.cif_files:
        p = Path(f)
        if p.is_file():
            cif_paths.append(p)
        else:
            print(f"Warning: {f} not found, skipping.", file=sys.stderr)

    if args.from_list:
        list_file = Path(args.from_list)
        if not list_file.is_file():
            print(f"Error: list file {args.from_list} not found.", file=sys.stderr)
            sys.exit(1)
        for line in list_file.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                p = Path(line)
                if p.is_file():
                    cif_paths.append(p)
                else:
                    print(f"Warning: {line} not found, skipping.", file=sys.stderr)

    if not cif_paths:
        print("No CIF files provided. Use -h for help.", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Desolvating {len(cif_paths)} structure(s) -> {output_dir}/")
    successes = sum(process_one(p, output_dir, args) for p in cif_paths)
    failures = len(cif_paths) - successes

    print(f"\nDone: {successes} succeeded, {failures} failed.")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
