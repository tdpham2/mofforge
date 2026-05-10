#!/usr/bin/env python3
"""Extract descriptive metadata for all building blocks in the pormake database.

Produces a structured output file (JSONL, JSON, or CSV) with enriched
information for each building block: chemical formula, connection points,
element composition, metal content, and SMILES string (when possible).

Usage:
    python scripts/describe_bbs.py
    python scripts/describe_bbs.py --format csv -o bbs.csv
    python scripts/describe_bbs.py --limit 50
    python scripts/describe_bbs.py --names N409 E41 N10
    python scripts/describe_bbs.py --no-smiles   # skip SMILES generation

Requirements:
    This script requires a separate Python environment with compatible
    package versions.  pormake 0.2.x depends on pymatgen<2024 which
    requires numpy<2.  See scripts/requirements.txt for details.

    Install into a dedicated environment:
        pip install -r scripts/requirements.txt

    Or manually:
        pip install pormake 'pymatgen>=2023.8.10,<2024.0.0' 'numpy<2'
        pip install rdkit   # optional, for SMILES generation
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

_COMPAT_MSG = (
    "Error: Failed to import pormake.\n\n"
    "These scripts require a Python environment with compatible packages.\n"
    "Install with:\n"
    "  pip install -r scripts/requirements.txt\n\n"
    "The key constraint is numpy<2 (pormake's pymatgen<2024 dependency\n"
    "uses APIs removed in numpy 2.0).  The mofforge .venv may have\n"
    "incompatible versions -- use a separate environment for these scripts."
)


def get_database():
    """Return a pormake Database instance.

    Catches import failures caused by version conflicts (e.g. numpy>=2
    with pymatgen<2024) and prints an actionable error message.
    """
    try:
        import pormake as pm
    except (ImportError, AttributeError) as exc:
        print(f"{_COMPAT_MSG}\n\nUnderlying error: {exc}", file=sys.stderr)
        sys.exit(1)

    return pm.Database()


# ------------------------------------------------------------------ #
# SMILES conversion
# ------------------------------------------------------------------ #

_RDKIT_AVAILABLE = None


def _check_rdkit() -> bool:
    global _RDKIT_AVAILABLE
    if _RDKIT_AVAILABLE is None:
        try:
            from rdkit import Chem  # noqa: F401

            _RDKIT_AVAILABLE = True
        except ImportError:
            _RDKIT_AVAILABLE = False
    return _RDKIT_AVAILABLE


def bb_to_smiles(bb) -> tuple[str | None, str | None]:
    """Convert a pormake BuildingBlock to a SMILES string.

    Connection points (X atoms) are mapped to dummy atoms (*).

    Returns:
        (smiles, error) -- one of them is always None.
    """
    if not _check_rdkit():
        return None, "rdkit not installed"

    from rdkit import Chem
    from rdkit import RDLogger

    # Suppress noisy RDKit warnings during conversion.
    RDLogger.logger().setLevel(RDLogger.ERROR)

    try:
        rw = Chem.RWMol()

        # Map bb atom index -> rdkit atom index.
        idx_map = {}
        symbols = list(bb.atoms.symbols)
        for i, sym in enumerate(symbols):
            if sym == "X":
                atom = Chem.Atom(0)  # dummy atom -> '*'
            else:
                atom = Chem.Atom(sym)
            idx_map[i] = rw.AddAtom(atom)

        # Bond type mapping.
        bond_type_map = {
            "S": Chem.BondType.SINGLE,
            "D": Chem.BondType.DOUBLE,
            "T": Chem.BondType.TRIPLE,
            "A": Chem.BondType.AROMATIC,
        }

        # Add bonds with deduplication (some XYZ files list bonds twice).
        seen_bonds = set()
        if bb.bonds is not None:
            for (i, j), bt in zip(bb.bonds, bb.bond_types):
                pair = (min(int(i), int(j)), max(int(i), int(j)))
                if pair in seen_bonds:
                    continue
                seen_bonds.add(pair)
                rw.AddBond(
                    int(idx_map[i]),
                    int(idx_map[j]),
                    bond_type_map.get(bt, Chem.BondType.SINGLE),
                )

        mol = rw.GetMol()

        # Try full sanitization first.
        try:
            Chem.SanitizeMol(mol)
            return Chem.MolToSmiles(mol), None
        except Exception:
            pass

        # Fallback: skip aromaticity perception, output kekulized SMILES.
        try:
            Chem.SanitizeMol(
                mol,
                Chem.SanitizeFlags.SANITIZE_ALL
                ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
            )
            return Chem.MolToSmiles(mol, kekuleSmiles=True), None
        except Exception as exc:
            return None, str(exc)

    except Exception as exc:
        return None, str(exc)


# ------------------------------------------------------------------ #
# Building block description
# ------------------------------------------------------------------ #


def describe_one_bb(db, name: str, compute_smiles: bool = True) -> dict | None:
    """Extract metadata for a single building block.

    Returns a dict on success or None on failure.
    """
    try:
        bb = db.get_bb(name)
    except Exception as exc:
        print(f"  SKIP {name}: {exc}", file=sys.stderr)
        return None

    # Element list excluding connection point markers.
    all_symbols = list(bb.atoms.symbols)
    elements = sorted({s for s in all_symbols if s != "X"})

    # Chemical formula (includes X for connection points).
    formula_full = bb.atoms.get_chemical_formula()

    # Formula without X: recount excluding X atoms.
    from collections import Counter

    real_symbols = [s for s in all_symbols if s != "X"]
    elem_counts = Counter(real_symbols)
    # Build Hill-system formula (C first, H second, then alphabetical).
    formula_parts = []
    for el in ["C", "H"]:
        if el in elem_counts:
            cnt = elem_counts.pop(el)
            formula_parts.append(f"{el}{cnt}" if cnt > 1 else el)
    for el in sorted(elem_counts):
        cnt = elem_counts[el]
        formula_parts.append(f"{el}{cnt}" if cnt > 1 else el)
    formula_clean = "".join(formula_parts)

    # SMILES.
    smiles = None
    smiles_error = None
    if compute_smiles:
        smiles, smiles_error = bb_to_smiles(bb)

    return {
        "name": bb.name,
        "role": "node" if bb.is_node else "edge",
        "formula": formula_full,
        "formula_clean": formula_clean,
        "n_atoms": int(bb.n_atoms),
        "n_connection_points": int(bb.n_connection_points),
        "connection_point_indices": bb.connection_point_indices.tolist(),
        "has_metal": bool(bb.has_metal),
        "elements": elements,
        "smiles": smiles,
        "smiles_error": smiles_error,
    }


# ------------------------------------------------------------------ #
# Writers
# ------------------------------------------------------------------ #


def write_jsonl(records: list[dict], path: Path) -> None:
    with path.open("w") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")


def write_json(records: list[dict], path: Path) -> None:
    with path.open("w") as f:
        json.dump(records, f, indent=2, ensure_ascii=False)


def write_csv(records: list[dict], path: Path) -> None:
    """Write a flat CSV (list fields are JSON-encoded)."""
    if not records:
        return

    fields = [
        "name",
        "role",
        "formula",
        "formula_clean",
        "n_atoms",
        "n_connection_points",
        "connection_point_indices",
        "has_metal",
        "elements",
        "smiles",
        "smiles_error",
    ]

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for rec in records:
            row = {}
            for k in fields:
                v = rec.get(k)
                if isinstance(v, (list, dict)):
                    row[k] = json.dumps(v)
                else:
                    row[k] = v
            writer.writerow(row)


# ------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------ #


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract descriptive metadata for pormake building blocks.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output file path. Default: scripts/output/bbs.<format>",
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=["jsonl", "json", "csv"],
        default="jsonl",
        help="Output format (default: jsonl).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N building blocks (for testing).",
    )
    parser.add_argument(
        "--names",
        nargs="+",
        default=None,
        help="Process only these specific building block names.",
    )
    parser.add_argument(
        "--no-smiles",
        action="store_true",
        help="Skip SMILES generation (faster, no rdkit dependency).",
    )
    args = parser.parse_args()

    # Resolve output path.
    if args.output is None:
        out_dir = Path(__file__).resolve().parent / "output"
        out_dir.mkdir(parents=True, exist_ok=True)
        args.output = out_dir / f"bbs.{args.format}"

    print("Loading pormake database...", file=sys.stderr)
    db = get_database()

    # Determine which BBs to process.
    if args.names:
        names = args.names
    else:
        names = sorted(db.bb_list)
        if args.limit:
            names = names[: args.limit]

    total = len(names)
    compute_smiles = not args.no_smiles

    if compute_smiles and not _check_rdkit():
        print(
            "WARNING: rdkit not installed -- SMILES will be skipped. "
            "Install with: pip install rdkit",
            file=sys.stderr,
        )
        compute_smiles = False

    print(
        f"Processing {total} building blocks "
        f"(SMILES={'on' if compute_smiles else 'off'})...",
        file=sys.stderr,
    )

    records = []
    skipped = 0
    smiles_ok = 0
    smiles_fail = 0
    t_start = time.time()

    for i, name in enumerate(names, start=1):
        rec = describe_one_bb(db, name, compute_smiles=compute_smiles)
        if rec is not None:
            records.append(rec)
            if compute_smiles:
                if rec["smiles"] is not None:
                    smiles_ok += 1
                else:
                    smiles_fail += 1
        else:
            skipped += 1

        # Progress every 100 items or at the end.
        if i % 100 == 0 or i == total:
            elapsed = time.time() - t_start
            print(
                f"\r  [{i}/{total}] {len(records)} ok, {skipped} skipped ({elapsed:.1f}s)",
                end="",
                file=sys.stderr,
            )

    print(file=sys.stderr)  # newline after progress

    # Write output.
    writers = {"jsonl": write_jsonl, "json": write_json, "csv": write_csv}
    writers[args.format](records, args.output)

    elapsed_total = time.time() - t_start
    summary = (
        f"Done. Wrote {len(records)} building blocks to {args.output} "
        f"({skipped} skipped, {elapsed_total:.1f}s total)."
    )
    if compute_smiles:
        summary += (
            f"\n  SMILES: {smiles_ok} succeeded, {smiles_fail} failed "
            f"({smiles_ok / max(smiles_ok + smiles_fail, 1) * 100:.0f}% success rate)."
        )
    print(summary, file=sys.stderr)


if __name__ == "__main__":
    main()
