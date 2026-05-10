#!/usr/bin/env python3
"""Standalone demo / validation script for SMILES-to-building-block conversion.

Converts a curated set of MOF linker molecules to both TOBACCO CIF and
Pormake XYZ building-block formats, printing detailed diagnostics for
each step.

Output files are written to::

    scripts/output/tobacco/<name>.cif
    scripts/output/pormake/<name>.xyz

Usage::

    python scripts/test_smiles_to_bb_demo.py

Requires: ``rdkit``, ``numpy``, ``networkx`` (same as smiles_to_bb.py).
Optional: ``pormake`` -- if installed, the script also verifies that
pormake can load each generated XYZ file.
"""

from __future__ import annotations

import sys
import traceback
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Ensure the package is importable when running from the scripts/ dir.
_repo_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_repo_root / "src"))

from mofforge.build.smiles_to_bb import (
    detect_carboxylic_groups,
    detect_connection_points,
    smiles_to_pormake_edge_xyz,
    smiles_to_tobacco_edge_cif,
)

# ------------------------------------------------------------------ #
# Molecule registry
# ------------------------------------------------------------------ #

MOLECULES = [
    {
        "name": "BDC",
        "smiles": "OC(=O)c1ccc(C(=O)O)cc1",
        "description": "Terephthalic acid (1,4-benzenedicarboxylic acid)",
    },
    {
        "name": "NDC",
        "smiles": "OC(=O)c1ccc2cc(C(=O)O)ccc2c1",
        "description": "2,6-Naphthalenedicarboxylic acid",
    },
    {
        "name": "BPDC",
        "smiles": "OC(=O)c1ccc(-c2ccc(C(=O)O)cc2)cc1",
        "description": "Biphenyl-4,4'-dicarboxylic acid",
    },
    {
        "name": "stilbene_dc",
        "smiles": "OC(=O)c1ccc(/C=C/c2ccc(C(=O)O)cc2)cc1",
        "description": "4,4'-Stilbenedicarboxylic acid (flexible linker)",
    },
    {
        "name": "biphenyl",
        "smiles": "c1ccc(-c2ccccc2)cc1",
        "description": "Biphenyl (direct mode, no functional groups)",
    },
    {
        "name": "azobenzene",
        "smiles": "c1ccc(/N=N/c2ccccc2)cc1",
        "description": "Azobenzene (direct mode, N=N bridge)",
    },
    {
        "name": "DABCO",
        "smiles": "C1CN2CCN1CC2",
        "description": "1,4-Diazabicyclo[2.2.2]octane (aliphatic pillar)",
    },
]


# ------------------------------------------------------------------ #
# Lightweight file parsers
# ------------------------------------------------------------------ #


def _parse_pormake_xyz(content: str):
    """Parse a Pormake-format extended XYZ string.

    Returns (atoms, bonds, x_indices) where:
    - atoms: list of (symbol, x, y, z)
    - bonds: list of (i, j, bond_type)
    - x_indices: list of int from the comment line
    """
    lines = content.strip().splitlines()
    n_atoms = int(lines[0].strip())
    comment = lines[1].strip()
    x_indices = [int(x) for x in comment.split()] if comment else []

    atoms = []
    for line in lines[2 : 2 + n_atoms]:
        parts = line.split()
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))

    bonds = []
    for line in lines[2 + n_atoms :]:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        bonds.append((int(parts[0]), int(parts[1]), parts[2]))

    return atoms, bonds, x_indices


def _parse_tobacco_cif(content: str):
    """Lightweight parse of a TOBACCO CIF to extract atom/bond counts.

    Returns a dict with keys:
        atom_count, bond_count, x_atom_count, has_fr,
        bond_types (dict of type -> count),
        atom_labels (list of label strings).
    """
    lines = content.splitlines()
    atoms: list[str] = []
    bonds: list[str] = []
    bond_types: dict[str, int] = {}
    has_fr = False

    # Two-pass approach: first find loop boundaries, then parse data.
    # Each loop starts with "loop_" followed by "_key" lines, then data.

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line == "loop_":
            # Collect key lines
            j = i + 1
            keys = []
            while j < len(lines) and lines[j].strip().startswith("_"):
                keys.append(lines[j].strip())
                j += 1

            is_atom_loop = "_atom_site_label" in keys
            is_bond_loop = any(k == "_geom_bond_atom_site_label_1" for k in keys)

            # Read data lines until next loop_, blank line, or key line
            while j < len(lines):
                dline = lines[j].strip()
                if not dline or dline == "loop_" or dline.startswith("_"):
                    break
                parts = dline.split()

                if is_atom_loop and len(parts) >= 2:
                    atoms.append(parts[0])
                    if parts[1] == "Fr":
                        has_fr = True
                elif is_bond_loop and len(parts) >= 5:
                    bt = parts[4]
                    bond_types[bt] = bond_types.get(bt, 0) + 1
                    bonds.append(f"{parts[0]}-{parts[1]}")

                j += 1

            # Continue from where we left off (don't skip the boundary
            # line -- it may be the next "loop_").
            i = j
            continue

        i += 1

    x_atom_count = sum(1 for a in atoms if a.startswith("X"))

    return {
        "atom_count": len(atoms),
        "bond_count": len(bonds),
        "x_atom_count": x_atom_count,
        "has_fr": has_fr,
        "bond_types": bond_types,
        "atom_labels": atoms,
    }


# ------------------------------------------------------------------ #
# Result container
# ------------------------------------------------------------------ #


@dataclass
class MoleculeResult:
    name: str
    description: str
    smiles: str
    mode: str = ""
    canonical_smiles: str = ""
    tobacco_atoms: int = 0
    tobacco_bonds: int = 0
    pormake_atoms: int = 0
    pormake_bonds: int = 0
    x_count_tobacco: int = 0
    x_count_pormake: int = 0
    pormake_loadable: str = "N/A"
    error: str = ""


# ------------------------------------------------------------------ #
# Pretty printers
# ------------------------------------------------------------------ #

_SEP = "-" * 72


def _print_header(name: str, description: str, smiles: str):
    print(f"\n{'=' * 72}")
    print(f"  {name}  --  {description}")
    print(f"  SMILES: {smiles}")
    print(f"{'=' * 72}")


def _print_detection(info):
    print(f"\n  Detection")
    print(f"  {_SEP}")
    print(f"  Mode:             {info.mode}")
    print(f"  Canonical SMILES: {info.smiles}")
    print(f"  Connection atoms: {info.connection_atom_indices}")
    if info.carboxylate_groups:
        for i, g in enumerate(info.carboxylate_groups):
            print(
                f"    Carboxylate {i + 1}: C={g.carbon_idx}, "
                f"=O={g.oxy_double_idx}, -O(H)={g.oxy_single_idx}, "
                f"anchor={g.anchor_idx}"
            )


def _print_tobacco_summary(cif_info: dict, cif_path: Path):
    print(f"\n  TOBACCO CIF")
    print(f"  {_SEP}")
    print(f"  Output:      {cif_path}")
    print(f"  Atoms:       {cif_info['atom_count']}")
    print(f"  Bonds:       {cif_info['bond_count']}")
    print(f"  X atoms:     {cif_info['x_atom_count']}")
    print(f"  Has Fr:      {cif_info['has_fr']}")
    bt = cif_info["bond_types"]
    print(f"  Bond types:  {', '.join(f'{k}={v}' for k, v in sorted(bt.items()))}")


def _print_pormake_summary(
    atoms: list,
    bonds: list,
    x_indices: list[int],
    xyz_path: Path,
):
    print(f"\n  Pormake XYZ")
    print(f"  {_SEP}")
    print(f"  Output:      {xyz_path}")
    print(f"  Atoms:       {len(atoms)}")
    print(f"  Bonds:       {len(bonds)}")
    print(f"  X indices:   {x_indices}")

    # Element composition
    from collections import Counter

    comp = Counter(a[0] for a in atoms)
    print(f"  Composition: {', '.join(f'{el}={n}' for el, n in sorted(comp.items()))}")

    # Bond type distribution
    bt_dist = Counter(b[2] for b in bonds)
    print(f"  Bond types:  {', '.join(f'{k}={v}' for k, v in sorted(bt_dist.items()))}")

    # X-atom distances to bonded atom
    for xi in x_indices:
        x_pos = np.array(atoms[xi][1:4])
        bonded_idx = None
        for bi, bj, _ in bonds:
            if bi == xi:
                bonded_idx = bj
                break
            if bj == xi:
                bonded_idx = bi
                break
        if bonded_idx is not None:
            bonded_pos = np.array(atoms[bonded_idx][1:4])
            dist = float(np.linalg.norm(x_pos - bonded_pos))
            bonded_sym = atoms[bonded_idx][0]
            print(f"  X{xi} -> {bonded_sym}{bonded_idx} distance: {dist:.4f} A")
        else:
            print(f"  X{xi}: WARNING -- no bond found!")


def _print_pormake_load(xyz_path: Path) -> str:
    """Try loading the XYZ with pormake. Returns status string."""
    try:
        import pormake as pm
    except ImportError:
        print("  Pormake:     (not installed -- skipping load test)")
        return "N/A"

    try:
        bb = pm.BuildingBlock(str(xyz_path))
        n_cp = len(bb.connection_point_indices)
        print(f"  Pormake:     loaded OK, {n_cp} connection point(s)")
        return f"OK ({n_cp} cp)"
    except Exception as exc:
        print(f"  Pormake:     FAILED to load -- {exc}")
        return f"FAIL: {exc}"


# ------------------------------------------------------------------ #
# Summary table
# ------------------------------------------------------------------ #


def _print_summary_table(results: list[MoleculeResult]):
    print(f"\n\n{'=' * 100}")
    print("  SUMMARY")
    print(f"{'=' * 100}")

    hdr = (
        f"{'Name':<18s} {'Mode':<13s} "
        f"{'TOB atoms':>9s} {'TOB bonds':>9s} "
        f"{'PM atoms':>8s} {'PM bonds':>8s} "
        f"{'X(TOB)':>6s} {'X(PM)':>5s} "
        f"{'PM Load':<15s} {'Status':<6s}"
    )
    print(f"  {hdr}")
    print(f"  {'-' * len(hdr)}")

    for r in results:
        status = "OK" if not r.error else "FAIL"
        row = (
            f"{r.name:<18s} {r.mode:<13s} "
            f"{r.tobacco_atoms:>9d} {r.tobacco_bonds:>9d} "
            f"{r.pormake_atoms:>8d} {r.pormake_bonds:>8d} "
            f"{r.x_count_tobacco:>6d} {r.x_count_pormake:>5d} "
            f"{r.pormake_loadable:<15s} {status:<6s}"
        )
        print(f"  {row}")
        if r.error:
            print(f"    ERROR: {r.error}")

    print()

    # Count pass / fail
    passed = sum(1 for r in results if not r.error)
    failed = len(results) - passed
    print(f"  {passed}/{len(results)} molecules converted successfully", end="")
    if failed:
        print(f", {failed} FAILED")
    else:
        print()
    print()


# ------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------ #


def _has_carboxylate(smiles: str) -> bool:
    """Return True if the SMILES contains at least one COOH/COO- group."""
    try:
        info = detect_connection_points(smiles)
        return info.mode == "carboxylate"
    except ValueError:
        return False


def main():
    output_root = Path(__file__).resolve().parent / "output"
    tobacco_dir = output_root / "tobacco"
    tobacco_carboxylic_dir = output_root / "tobacco_carboxylic"
    pormake_dir = output_root / "pormake"
    tobacco_dir.mkdir(parents=True, exist_ok=True)
    tobacco_carboxylic_dir.mkdir(parents=True, exist_ok=True)
    pormake_dir.mkdir(parents=True, exist_ok=True)

    print("SMILES-to-Building-Block Conversion Demo")
    print(f"Output directory: {output_root}")
    print(f"Molecules to convert: {len(MOLECULES)}")

    results: list[MoleculeResult] = []

    for mol in MOLECULES:
        name = mol["name"]
        smiles = mol["smiles"]
        description = mol["description"]

        res = MoleculeResult(
            name=name,
            description=description,
            smiles=smiles,
        )

        _print_header(name, description, smiles)

        try:
            # --- 1. Detection ----------------------------------------- #
            info = detect_connection_points(smiles)
            res.mode = info.mode
            res.canonical_smiles = info.smiles
            _print_detection(info)

            # --- 2. TOBACCO CIF (auto / carboxylate mode) ------------- #
            cif_path = tobacco_dir / f"{name}.cif"
            smiles_to_tobacco_edge_cif(
                smiles,
                output_path=cif_path,
                name=f"{name}_edge",
            )

            cif_content = cif_path.read_text()
            cif_info = _parse_tobacco_cif(cif_content)
            res.tobacco_atoms = cif_info["atom_count"]
            res.tobacco_bonds = cif_info["bond_count"]
            res.x_count_tobacco = cif_info["x_atom_count"]
            _print_tobacco_summary(cif_info, cif_path)

            # --- 2b. TOBACCO CIF (carboxylic mode) -------------------- #
            #     Only for molecules with exactly 2 COOH/COO- groups.
            if _has_carboxylate(smiles):
                print(f"\n  TOBACCO CIF (carboxylic mode)")
                print(f"  {_SEP}")

                carboxylic_info = detect_carboxylic_groups(smiles)
                print(f"  Mode:             {carboxylic_info.mode}")
                print(f"  Connection atoms: {carboxylic_info.connection_atom_indices}")
                for i, g in enumerate(carboxylic_info.carboxylate_groups):
                    print(
                        f"    COOH group {i + 1}: C={g.carbon_idx}, "
                        f"=O={g.oxy_double_idx}, -O(H)={g.oxy_single_idx}, "
                        f"anchor={g.anchor_idx}"
                    )

                cif_carboxylic_path = tobacco_carboxylic_dir / f"{name}.cif"
                smiles_to_tobacco_edge_cif(
                    smiles,
                    output_path=cif_carboxylic_path,
                    name=f"{name}_carboxylic_edge",
                    mode="carboxylic",
                )

                cif_carboxylic_content = cif_carboxylic_path.read_text()
                cif_carboxylic_info = _parse_tobacco_cif(cif_carboxylic_content)
                _print_tobacco_summary(cif_carboxylic_info, cif_carboxylic_path)

                # Show the differences.
                print(f"\n  Comparison (carboxylate vs carboxylic)")
                print(f"  {_SEP}")
                print(f"  {'':20s} {'carboxylate':>12s} {'carboxylic':>12s}")
                print(f"  {'Atoms':<20s} {cif_info['atom_count']:>12d} {cif_carboxylic_info['atom_count']:>12d}")
                print(f"  {'Bonds':<20s} {cif_info['bond_count']:>12d} {cif_carboxylic_info['bond_count']:>12d}")
                print(f"  {'X atoms':<20s} {cif_info['x_atom_count']:>12d} {cif_carboxylic_info['x_atom_count']:>12d}")
                print(f"  {'Has Fr':<20s} {str(cif_info['has_fr']):>12s} {str(cif_carboxylic_info['has_fr']):>12s}")

            # --- 3. Pormake XYZ --------------------------------------- #
            xyz_path = pormake_dir / f"{name}.xyz"
            smiles_to_pormake_edge_xyz(
                smiles,
                output_path=xyz_path,
            )

            xyz_content = xyz_path.read_text()
            atoms, bonds, x_indices = _parse_pormake_xyz(xyz_content)
            res.pormake_atoms = len(atoms)
            res.pormake_bonds = len(bonds)
            res.x_count_pormake = len(x_indices)
            _print_pormake_summary(atoms, bonds, x_indices, xyz_path)

            # --- 4. Pormake load test --------------------------------- #
            res.pormake_loadable = _print_pormake_load(xyz_path)

        except Exception as exc:
            res.error = str(exc)
            print(f"\n  *** ERROR: {exc} ***")
            traceback.print_exc()

        results.append(res)

    # --- Final summary table ---------------------------------------------- #
    _print_summary_table(results)

    # Exit with non-zero status if any molecule failed.
    if any(r.error for r in results):
        sys.exit(1)


if __name__ == "__main__":
    main()
