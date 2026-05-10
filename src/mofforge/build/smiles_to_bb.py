"""Convert a SMILES string to a TOBACCO or Pormake edge building block.

Detects connection points via three strategies: carboxylate-terminated,
direct aromatic, and carboxylic-stripped (``mode="carboxylic"``).
Requires ``rdkit``.
"""

from __future__ import annotations

import datetime
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np

logger = logging.getLogger("mofforge")

_rdkit_loaded = False


def _ensure_rdkit():
    global _rdkit_loaded  # noqa: PLW0603
    if _rdkit_loaded:
        return
    try:
        import rdkit  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "rdkit is required for SMILES-to-building-block conversion.  "
            "Install it with:  pip install rdkit"
        ) from exc
    _rdkit_loaded = True


# SMARTS for a carboxylate group:  C(=O)[O,OH]
# Matches both protonated (COOH) and deprotonated (COO-) forms.
_CARBOXYLATE_SMARTS = "[CX3](=[OX1])[OX1,OX2H1]"

# Default cubic cell length (Angstrom) for the non-periodic CIF.
_CELL_LENGTH = 40.0



@dataclass
class CarboxylateGroup:
    """Indices (in the RDKit mol) for one carboxylate group."""

    carbon_idx: int
    oxy_double_idx: int
    oxy_single_idx: int
    anchor_idx: int = -1


@dataclass
class ConnectionInfo:
    """Result of automatic connection-point detection."""

    mode: Literal["carboxylate", "direct", "carboxylic"]
    connection_atom_indices: list[int] = field(default_factory=list)
    carboxylate_groups: list[CarboxylateGroup] = field(default_factory=list)
    smiles: str = ""


def detect_connection_points(smiles: str, n_points: int = 2) -> ConnectionInfo:
    """Analyse a SMILES string and identify edge connection points."""
    _ensure_rdkit()
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")

    canonical = Chem.MolToSmiles(mol)

    # --- Try carboxylate detection first ---------------------------------- #
    pattern = Chem.MolFromSmarts(_CARBOXYLATE_SMARTS)
    matches = mol.GetSubstructMatches(pattern)

    if matches:
        if len(matches) != n_points:
            raise ValueError(
                f"Expected {n_points} carboxylate groups but found "
                f"{len(matches)} in SMILES: {smiles!r}"
            )

        groups: list[CarboxylateGroup] = []
        for match in matches:
            c_idx, oxy_dbl, oxy_sgl = match[0], match[1], match[2]
            # Find the anchor atom: the heavy-atom neighbour of the
            # carboxylate carbon that is NOT one of the carboxylate
            # oxygens.
            c_atom = mol.GetAtomWithIdx(c_idx)
            anchor = -1
            for nbr in c_atom.GetNeighbors():
                if nbr.GetIdx() not in (oxy_dbl, oxy_sgl):
                    anchor = nbr.GetIdx()
                    break
            groups.append(
                CarboxylateGroup(
                    carbon_idx=c_idx,
                    oxy_double_idx=oxy_dbl,
                    oxy_single_idx=oxy_sgl,
                    anchor_idx=anchor,
                )
            )

        return ConnectionInfo(
            mode="carboxylate",
            connection_atom_indices=[g.anchor_idx for g in groups],
            carboxylate_groups=groups,
            smiles=canonical,
        )

    # --- Fallback: direct (no carboxylate) mode --------------------------- #
    return _detect_direct_connection_points(mol, canonical, n_points)


def _detect_direct_connection_points(
    mol,  # rdkit.Chem.Mol
    canonical_smiles: str,
    n_points: int,
) -> ConnectionInfo:
    """Identify terminal backbone atoms for direct connection."""
    import networkx as nx

    # Build heavy-atom graph
    G = nx.Graph()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1:  # skip H
            G.add_node(atom.GetIdx())

    for bond in mol.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if mol.GetAtomWithIdx(a).GetAtomicNum() != 1 and mol.GetAtomWithIdx(b).GetAtomicNum() != 1:
            G.add_edge(a, b)

    if len(G) < 2:
        raise ValueError("Molecule has fewer than 2 heavy atoms; cannot identify connection points")

    # Find the pair with the longest shortest-path distance (diameter).
    # For efficiency with small molecules, just compute eccentricities.
    if not nx.is_connected(G):
        raise ValueError(
            "Heavy-atom graph is disconnected; cannot identify connection "
            f"points for SMILES: {canonical_smiles!r}"
        )

    # Find diameter endpoints via double BFS (faster than all-pairs).
    # Pick an arbitrary start, find the farthest node, then from that
    # node find the farthest again – the two endpoints of the longest
    # shortest path.
    start = next(iter(G.nodes()))
    lengths_from_start = nx.single_source_shortest_path_length(G, start)
    u = max(lengths_from_start, key=lengths_from_start.get)
    lengths_from_u = nx.single_source_shortest_path_length(G, u)
    v = max(lengths_from_u, key=lengths_from_u.get)

    if n_points != 2:
        raise ValueError(
            f"Direct mode currently supports n_points=2 only, got {n_points}"
        )

    return ConnectionInfo(
        mode="direct",
        connection_atom_indices=[u, v],
        smiles=canonical_smiles,
    )


def detect_carboxylic_groups(smiles: str) -> ConnectionInfo:
    """Detect COOH / COO- groups and return connection info in *carboxylic* mode."""
    _ensure_rdkit()
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")

    canonical = Chem.MolToSmiles(mol)

    pattern = Chem.MolFromSmarts(_CARBOXYLATE_SMARTS)
    matches = mol.GetSubstructMatches(pattern)

    if not matches:
        raise ValueError(
            f"No COOH/COO- groups found in SMILES: {smiles!r}.  "
            f"The carboxylic mode requires exactly 2 carboxylate groups."
        )

    if len(matches) != 2:
        raise ValueError(
            f"Expected exactly 2 COOH/COO- groups but found "
            f"{len(matches)} in SMILES: {smiles!r}"
        )

    groups: list[CarboxylateGroup] = []
    for match in matches:
        c_idx, oxy_dbl, oxy_sgl = match[0], match[1], match[2]
        c_atom = mol.GetAtomWithIdx(c_idx)
        anchor = -1
        for nbr in c_atom.GetNeighbors():
            if nbr.GetIdx() not in (oxy_dbl, oxy_sgl):
                anchor = nbr.GetIdx()
                break
        groups.append(
            CarboxylateGroup(
                carbon_idx=c_idx,
                oxy_double_idx=oxy_dbl,
                oxy_single_idx=oxy_sgl,
                anchor_idx=anchor,
            )
        )

    return ConnectionInfo(
        mode="carboxylic",
        connection_atom_indices=[g.anchor_idx for g in groups],
        carboxylate_groups=groups,
        smiles=canonical,
    )


def smiles_to_tobacco_edge_cif(
    smiles: str,
    output_path: str | Path,
    name: str = "edge",
    cell_length: float = _CELL_LENGTH,
    uff_max_iters: int = 2000,
    mode: Literal["auto", "carboxylate", "direct", "carboxylic"] = "auto",
) -> Path:
    """Convert a SMILES string to a TOBACCO-format edge building-block CIF.

    Parses SMILES, generates 3-D coordinates, detects connection points,
    and writes a TOBACCO-format CIF.
    """
    _ensure_rdkit()
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdmolops

    output_path = Path(output_path).resolve()

    # ---- 1. Parse and generate 3-D coordinates --------------------------- #
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")

    # Detect connection points on the *implicit-H* molecule (the same
    # representation the SMILES encodes) so that atom indices are stable
    # before we add explicit Hs.
    if mode == "carboxylic":
        conn_info = detect_carboxylic_groups(smiles)
    elif mode in ("carboxylate", "direct"):
        conn_info = detect_connection_points(smiles)
        if conn_info.mode != mode:
            raise ValueError(
                f"Requested mode {mode!r} but detection found "
                f"{conn_info.mode!r} for SMILES: {smiles!r}"
            )
    else:
        # "auto" — the original behaviour
        conn_info = detect_connection_points(smiles)

    # Add explicit hydrogens (needed for 3-D embedding).
    mol_h = Chem.AddHs(mol)

    # Map original atom indices to the mol_h indices.  AddHs preserves
    # the order of existing atoms and appends Hs, so indices for heavy
    # atoms are unchanged.

    # Embed 3-D coordinates
    embed_result = AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())
    if embed_result != 0:
        # Retry with random coordinates
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        embed_result = AllChem.EmbedMolecule(mol_h, params)
        if embed_result != 0:
            raise ValueError(
                f"RDKit could not generate 3-D coordinates for: {smiles!r}"
            )

    # UFF optimisation
    try:
        AllChem.UFFOptimizeMolecule(mol_h, maxIters=uff_max_iters)
    except Exception:
        logger.warning("UFF optimisation failed for %r; using embedded coords", smiles)

    conf = mol_h.GetConformer()

    # ---- 2. Build atom/bond lists depending on mode ---------------------- #
    if conn_info.mode == "carboxylic":
        atoms, bonds = _build_carboxylic_edge(mol_h, conf, conn_info, cell_length)
    elif conn_info.mode == "carboxylate":
        atoms, bonds = _build_carboxylate_edge(mol_h, conf, conn_info, cell_length)
    else:
        atoms, bonds = _build_direct_edge(mol_h, conf, conn_info, cell_length)

    # ---- 3. Write CIF --------------------------------------------------- #
    _write_tobacco_cif(atoms, bonds, cell_length, name, output_path)

    logger.debug(
        "Wrote TOBACCO edge CIF: %s  (mode=%s, %d atoms, %d bonds)",
        output_path,
        conn_info.mode,
        len(atoms),
        len(bonds),
    )
    return output_path


# Small data carriers for CIF writing

@dataclass
class _CifAtom:
    label: str          # e.g. "C1", "X3", "Fr7"
    type_symbol: str    # element symbol written in _atom_site_type_symbol
    fx: float           # fractional x
    fy: float           # fractional y
    fz: float           # fractional z
    charge: float = 0.0


@dataclass
class _CifBond:
    label1: str
    label2: str
    distance: float
    bond_type: str  # "S", "D", "A"


def _bond_type_char(bond) -> str:
    """Map an RDKit bond type to the TOBACCO CIF bond-type character."""
    from rdkit import Chem

    bt = bond.GetBondType()
    if bt == Chem.rdchem.BondType.AROMATIC:
        return "A"
    if bt == Chem.rdchem.BondType.DOUBLE:
        return "D"
    if bt == Chem.rdchem.BondType.TRIPLE:
        return "T"
    return "S"


def _make_label(symbol: str, counter: dict[str, int]) -> str:
    """Return a unique atom label like ``C1``, ``O2``, etc."""
    counter.setdefault(symbol, 0)
    counter[symbol] += 1
    return f"{symbol}{counter[symbol]}"


def _frac_coords(pos: np.ndarray, cell_length: float) -> tuple[float, float, float]:
    """Convert Cartesian position to fractional coordinates in a cubic cell."""
    return (
        float(pos[0] / cell_length),
        float(pos[1] / cell_length),
        float(pos[2] / cell_length),
    )


def _center_positions(positions: np.ndarray, cell_length: float) -> np.ndarray:
    """Translate positions so the centroid sits at the cell centre."""
    centroid = positions.mean(axis=0)
    target = np.array([cell_length / 2, cell_length / 2, cell_length / 2])
    return positions + (target - centroid)


def _build_direct_edge(
    mol_h,       # RDKit Mol with explicit Hs
    conf,        # RDKit Conformer
    conn_info: ConnectionInfo,
    cell_length: float,
) -> tuple[list[_CifAtom], list[_CifBond]]:
    """Build atom/bond lists for a *direct* (non-carboxylate) edge."""
    conn_set = set(conn_info.connection_atom_indices)

    # Identify H atoms to remove: those bonded to a connection-point atom.
    h_to_remove: set[int] = set()
    for idx in conn_set:
        atom = mol_h.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 1:
                h_to_remove.add(nbr.GetIdx())

    # Collect surviving atom indices and their positions.
    surviving = [
        i for i in range(mol_h.GetNumAtoms()) if i not in h_to_remove
    ]

    # Gather Cartesian positions for surviving atoms.
    positions = np.array([list(conf.GetAtomPosition(i)) for i in surviving])
    positions = _center_positions(positions, cell_length)

    # Build index maps: old_idx -> position index, label, type_symbol
    old_to_pos: dict[int, int] = {old: pos_i for pos_i, old in enumerate(surviving)}
    counter: dict[str, int] = {}
    idx_to_label: dict[int, str] = {}
    idx_to_sym: dict[int, str] = {}

    atoms: list[_CifAtom] = []
    for pos_i, old_idx in enumerate(surviving):
        atom = mol_h.GetAtomWithIdx(old_idx)
        sym = atom.GetSymbol()
        if old_idx in conn_set:
            label = _make_label("X", counter)
            type_sym = sym  # keep original element as type_symbol
        else:
            label = _make_label(sym, counter)
            type_sym = sym
        idx_to_label[old_idx] = label
        idx_to_sym[old_idx] = type_sym
        fx, fy, fz = _frac_coords(positions[pos_i], cell_length)
        atoms.append(_CifAtom(label=label, type_symbol=type_sym, fx=fx, fy=fy, fz=fz))

    # Build bond list (only for surviving atoms).
    surviving_set = set(surviving)
    bonds: list[_CifBond] = []
    for bond in mol_h.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a in surviving_set and b in surviving_set:
            pa = positions[old_to_pos[a]]
            pb = positions[old_to_pos[b]]
            dist = float(np.linalg.norm(pa - pb))
            bt = _bond_type_char(bond)
            bonds.append(_CifBond(idx_to_label[a], idx_to_label[b], dist, bt))

    return atoms, bonds


def _build_carboxylate_edge(
    mol_h,       # RDKit Mol with explicit Hs
    conf,        # RDKit Conformer
    conn_info: ConnectionInfo,
    cell_length: float,
) -> tuple[list[_CifAtom], list[_CifBond]]:
    """Build atom/bond lists for a *carboxylate*-terminated edge."""
    from rdkit import Chem

    # Collect all atom indices belonging to carboxylate groups.
    carboxylate_atom_set: set[int] = set()
    for g in conn_info.carboxylate_groups:
        carboxylate_atom_set.update([g.carbon_idx, g.oxy_double_idx, g.oxy_single_idx])

    # Identify H atoms to remove: H bonded to any carboxylate oxygen.
    h_to_remove: set[int] = set()
    for g in conn_info.carboxylate_groups:
        for oxy_idx in (g.oxy_double_idx, g.oxy_single_idx):
            oxy_atom = mol_h.GetAtomWithIdx(oxy_idx)
            for nbr in oxy_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 1:
                    h_to_remove.add(nbr.GetIdx())

    # Surviving atoms (original molecule atoms minus removed Hs).
    surviving = [
        i for i in range(mol_h.GetNumAtoms()) if i not in h_to_remove
    ]

    # Gather positions for surviving atoms.
    positions = np.array([list(conf.GetAtomPosition(i)) for i in surviving])
    positions = _center_positions(positions, cell_length)

    # Build index map and labels for surviving atoms.
    old_to_pos: dict[int, int] = {old: pos_i for pos_i, old in enumerate(surviving)}
    counter: dict[str, int] = {}
    idx_to_label: dict[int, str] = {}
    atoms: list[_CifAtom] = []
    for pos_i, old_idx in enumerate(surviving):
        atom = mol_h.GetAtomWithIdx(old_idx)
        sym = atom.GetSymbol()
        label = _make_label(sym, counter)
        idx_to_label[old_idx] = label
        fx, fy, fz = _frac_coords(positions[pos_i], cell_length)
        atoms.append(_CifAtom(label=label, type_symbol=sym, fx=fx, fy=fy, fz=fz))

    # Build bond list for surviving atoms.
    surviving_set = set(surviving)
    bonds: list[_CifBond] = []
    for bond in mol_h.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a in surviving_set and b in surviving_set:
            pa = positions[old_to_pos[a]]
            pb = positions[old_to_pos[b]]
            dist = float(np.linalg.norm(pa - pb))
            bt = _bond_type_char(bond)
            bonds.append(_CifBond(idx_to_label[a], idx_to_label[b], dist, bt))

    # Add X dummy atoms (Fr) at carboxylate centroids.
    for g in conn_info.carboxylate_groups:
        pos_o1 = positions[old_to_pos[g.oxy_double_idx]]
        pos_o2 = positions[old_to_pos[g.oxy_single_idx]]
        centroid = (pos_o1 + pos_o2) / 2.0

        x_label = _make_label("X", counter)
        fx, fy, fz = _frac_coords(centroid, cell_length)
        atoms.append(_CifAtom(label=x_label, type_symbol="Fr", fx=fx, fy=fy, fz=fz))

        # Bonds from X to both oxygens.
        d1 = float(np.linalg.norm(centroid - pos_o1))
        d2 = float(np.linalg.norm(centroid - pos_o2))
        bonds.append(_CifBond(x_label, idx_to_label[g.oxy_double_idx], d1, "S"))
        bonds.append(_CifBond(x_label, idx_to_label[g.oxy_single_idx], d2, "S"))

    return atoms, bonds


def _build_carboxylic_edge(
    mol_h,       # RDKit Mol with explicit Hs
    conf,        # RDKit Conformer
    conn_info: ConnectionInfo,
    cell_length: float,
) -> tuple[list[_CifAtom], list[_CifBond]]:
    """Build atom/bond lists for a *carboxylic*-stripped edge."""
    # Collect all atoms to remove:
    #   - carboxylate C, both O atoms
    #   - H atoms bonded to carboxylate oxygens
    atoms_to_remove: set[int] = set()
    anchor_set: set[int] = set()
    for g in conn_info.carboxylate_groups:
        atoms_to_remove.add(g.carbon_idx)
        atoms_to_remove.add(g.oxy_double_idx)
        atoms_to_remove.add(g.oxy_single_idx)
        anchor_set.add(g.anchor_idx)
        # Remove H on carboxylate oxygens
        for oxy_idx in (g.oxy_double_idx, g.oxy_single_idx):
            oxy_atom = mol_h.GetAtomWithIdx(oxy_idx)
            for nbr in oxy_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 1:
                    atoms_to_remove.add(nbr.GetIdx())

    # Also remove H atoms bonded to anchor atoms (since the anchor
    # becomes a connection point, its H is no longer needed).
    for anchor_idx in anchor_set:
        anchor_atom = mol_h.GetAtomWithIdx(anchor_idx)
        for nbr in anchor_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 1:
                atoms_to_remove.add(nbr.GetIdx())

    # Surviving atom indices.
    surviving = [
        i for i in range(mol_h.GetNumAtoms()) if i not in atoms_to_remove
    ]

    # Gather positions for surviving atoms, centred in the cell.
    positions = np.array([list(conf.GetAtomPosition(i)) for i in surviving])
    positions = _center_positions(positions, cell_length)

    # Build index map and labels.
    old_to_pos: dict[int, int] = {old: pos_i for pos_i, old in enumerate(surviving)}
    counter: dict[str, int] = {}
    idx_to_label: dict[int, str] = {}
    atoms: list[_CifAtom] = []
    for pos_i, old_idx in enumerate(surviving):
        atom = mol_h.GetAtomWithIdx(old_idx)
        sym = atom.GetSymbol()
        if old_idx in anchor_set:
            # Relabel the anchor atom as X, keep original element as type_symbol
            label = _make_label("X", counter)
            type_sym = sym
        else:
            label = _make_label(sym, counter)
            type_sym = sym
        idx_to_label[old_idx] = label
        fx, fy, fz = _frac_coords(positions[pos_i], cell_length)
        atoms.append(_CifAtom(label=label, type_symbol=type_sym, fx=fx, fy=fy, fz=fz))

    # Build bond list for surviving atoms.
    surviving_set = set(surviving)
    bonds: list[_CifBond] = []
    for bond in mol_h.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a in surviving_set and b in surviving_set:
            pa = positions[old_to_pos[a]]
            pb = positions[old_to_pos[b]]
            dist = float(np.linalg.norm(pa - pb))
            bt = _bond_type_char(bond)
            bonds.append(_CifBond(idx_to_label[a], idx_to_label[b], dist, bt))

    return atoms, bonds


def _write_tobacco_cif(
    atoms: list[_CifAtom],
    bonds: list[_CifBond],
    cell_length: float,
    name: str,
    path: Path,
) -> None:
    """Write a TOBACCO-format CIF file."""
    lines: list[str] = []

    # Header
    lines.append(f"data_{name}")
    lines.append(f"_audit_creation_date              {datetime.date.today().isoformat()}")
    lines.append("_audit_creation_method            'mofforge'")
    lines.append("_symmetry_space_group_name_H-M    'P1'")
    lines.append("_symmetry_Int_Tables_number       1")
    lines.append("_symmetry_cell_setting            triclinic")
    lines.append("loop_")
    lines.append("_symmetry_equiv_pos_as_xyz")
    lines.append("  x,y,z")

    # Cell parameters
    cl = f"{cell_length:.4f}"
    lines.append(f"_cell_length_a                    {cl}")
    lines.append(f"_cell_length_b                    {cl}")
    lines.append(f"_cell_length_c                    {cl}")
    lines.append("_cell_angle_alpha                 90.0000")
    lines.append("_cell_angle_beta                  90.0000")
    lines.append("_cell_angle_gamma                 90.0000")

    # Atom site loop
    lines.append("loop_")
    lines.append("_atom_site_label")
    lines.append("_atom_site_type_symbol")
    lines.append("_atom_site_fract_x")
    lines.append("_atom_site_fract_y")
    lines.append("_atom_site_fract_z")
    lines.append("_atom_site_U_iso_or_equiv")
    lines.append("_atom_site_adp_type")
    lines.append("_atom_site_occupancy")
    lines.append("_atom_site_charge")

    for a in atoms:
        lines.append(
            f"{a.label:<10s} {a.type_symbol:<3s} "
            f"{a.fx:>12.7f} {a.fy:>12.7f} {a.fz:>12.7f}"
            f"   0.00000  Uiso   1.00      {a.charge:.5f}"
        )

    # Bond loop
    lines.append("loop_")
    lines.append("_geom_bond_atom_site_label_1")
    lines.append("_geom_bond_atom_site_label_2")
    lines.append("_geom_bond_distance")
    lines.append("_geom_bond_site_symmetry_2")
    lines.append("_ccdc_geom_bond_type")

    for b in bonds:
        lines.append(
            f"{b.label1:<10s} {b.label2:<10s} "
            f"{b.distance:.5f}    .     {b.bond_type}"
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


# Distance (Angstrom) from the connection-point atom to the X dummy atom.
# Pormake's built-in database consistently uses 0.75 A.
_PORMAKE_X_DISTANCE = 0.75


def smiles_to_pormake_edge_xyz(
    smiles: str,
    output_path: str | Path,
    uff_max_iters: int = 2000,
    mode: Literal["auto", "carboxylate", "direct", "carboxylic"] = "auto",
) -> Path:
    """Convert a SMILES string to a Pormake-format edge building-block XYZ.

    Parses SMILES, generates 3-D coordinates, detects connection points,
    and writes a Pormake-format extended XYZ.
    """
    _ensure_rdkit()
    from rdkit import Chem
    from rdkit.Chem import AllChem

    output_path = Path(output_path).resolve()

    # ---- 1. Parse and generate 3-D coordinates --------------------------- #
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")

    if mode == "carboxylic":
        conn_info = detect_carboxylic_groups(smiles)
    elif mode in ("carboxylate", "direct"):
        conn_info = detect_connection_points(smiles)
        if conn_info.mode != mode:
            raise ValueError(
                f"Requested mode {mode!r} but detection found "
                f"{conn_info.mode!r} for SMILES: {smiles!r}"
            )
    else:
        # "auto" — the original behaviour
        conn_info = detect_connection_points(smiles)

    mol_h = Chem.AddHs(mol)

    embed_result = AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())
    if embed_result != 0:
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        embed_result = AllChem.EmbedMolecule(mol_h, params)
        if embed_result != 0:
            raise ValueError(
                f"RDKit could not generate 3-D coordinates for: {smiles!r}"
            )

    try:
        AllChem.UFFOptimizeMolecule(mol_h, maxIters=uff_max_iters)
    except Exception:
        logger.warning("UFF optimisation failed for %r; using embedded coords", smiles)

    conf = mol_h.GetConformer()

    # ---- 2. Build atom/bond lists depending on mode ---------------------- #
    if conn_info.mode == "carboxylic":
        xyz_atoms, xyz_bonds, x_indices = _build_pormake_carboxylic_edge(
            mol_h, conf, conn_info,
        )
    elif conn_info.mode == "carboxylate":
        xyz_atoms, xyz_bonds, x_indices = _build_pormake_carboxylate_edge(
            mol_h, conf, conn_info,
        )
    else:
        xyz_atoms, xyz_bonds, x_indices = _build_pormake_direct_edge(
            mol_h, conf, conn_info,
        )

    # ---- 3. Write XYZ ---------------------------------------------------- #
    _write_pormake_xyz(xyz_atoms, xyz_bonds, x_indices, output_path)

    logger.debug(
        "Wrote Pormake edge XYZ: %s  (mode=%s, %d atoms, %d bonds)",
        output_path,
        conn_info.mode,
        len(xyz_atoms),
        len(xyz_bonds),
    )
    return output_path



@dataclass
class _XyzAtom:
    symbol: str  # element symbol or "X" for connection point
    x: float
    y: float
    z: float


@dataclass
class _XyzBond:
    idx_i: int  # 0-indexed
    idx_j: int
    bond_type: str  # "S", "D", "T", "A"


def _outward_direction(
    conn_pos: np.ndarray,
    anchor_pos: np.ndarray,
) -> np.ndarray:
    """Unit vector pointing from *anchor_pos* through *conn_pos* and beyond.

    If the two positions coincide (degenerate case), returns a default
    direction along +x.
    """
    vec = conn_pos - anchor_pos
    norm = float(np.linalg.norm(vec))
    if norm < 1e-8:
        return np.array([1.0, 0.0, 0.0])
    return vec / norm


def _build_pormake_direct_edge(
    mol_h,
    conf,
    conn_info: ConnectionInfo,
) -> tuple[list[_XyzAtom], list[_XyzBond], list[int]]:
    """Build Pormake atom/bond lists for a *direct* (non-carboxylate) edge."""
    conn_set = set(conn_info.connection_atom_indices)

    # Identify H atoms to remove.
    h_to_remove: set[int] = set()
    for idx in conn_set:
        atom = mol_h.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 1:
                h_to_remove.add(nbr.GetIdx())

    # Surviving atom indices.
    surviving = [i for i in range(mol_h.GetNumAtoms()) if i not in h_to_remove]

    # Re-index: old_idx -> new_idx
    old_to_new: dict[int, int] = {old: new for new, old in enumerate(surviving)}

    # Build Cartesian positions for surviving atoms, centered on origin.
    positions = np.array([list(conf.GetAtomPosition(i)) for i in surviving])
    centroid = positions.mean(axis=0)
    positions -= centroid

    # Build atom list.
    atoms: list[_XyzAtom] = []
    for pos_i, old_idx in enumerate(surviving):
        sym = mol_h.GetAtomWithIdx(old_idx).GetSymbol()
        p = positions[pos_i]
        atoms.append(_XyzAtom(symbol=sym, x=p[0], y=p[1], z=p[2]))

    # Build bond list (surviving atoms only).
    bonds: list[_XyzBond] = []
    surviving_set = set(surviving)
    for bond in mol_h.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a in surviving_set and b in surviving_set:
            bt = _bond_type_char(bond)
            bonds.append(_XyzBond(old_to_new[a], old_to_new[b], bt))

    # Place X dummy atoms.
    x_indices: list[int] = []
    for conn_idx in conn_info.connection_atom_indices:
        conn_new = old_to_new[conn_idx]
        conn_pos = positions[surviving.index(conn_idx)]

        # Find the heavy-atom neighbour to define the outward direction.
        atom = mol_h.GetAtomWithIdx(conn_idx)
        heavy_nbr_pos = None
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 1 and nbr.GetIdx() in surviving_set:
                nbr_new_i = surviving.index(nbr.GetIdx())
                heavy_nbr_pos = positions[nbr_new_i]
                break

        if heavy_nbr_pos is None:
            # Fallback: direction from centroid
            direction = _outward_direction(conn_pos, np.zeros(3))
        else:
            direction = _outward_direction(conn_pos, heavy_nbr_pos)

        x_pos = conn_pos + direction * _PORMAKE_X_DISTANCE
        x_new_idx = len(atoms)
        atoms.append(_XyzAtom(symbol="X", x=x_pos[0], y=x_pos[1], z=x_pos[2]))
        bonds.append(_XyzBond(x_new_idx, conn_new, "S"))
        x_indices.append(x_new_idx)

    return atoms, bonds, x_indices


def _build_pormake_carboxylate_edge(
    mol_h,
    conf,
    conn_info: ConnectionInfo,
) -> tuple[list[_XyzAtom], list[_XyzBond], list[int]]:
    """Build Pormake atom/bond lists for a *carboxylate*-terminated edge."""
    # Collect all atoms to remove: carboxylate C, both O atoms, and
    # any H atoms bonded to those oxygens.
    atoms_to_remove: set[int] = set()
    for g in conn_info.carboxylate_groups:
        atoms_to_remove.add(g.carbon_idx)
        atoms_to_remove.add(g.oxy_double_idx)
        atoms_to_remove.add(g.oxy_single_idx)
        # Remove H on carboxylate oxygens
        for oxy_idx in (g.oxy_double_idx, g.oxy_single_idx):
            oxy_atom = mol_h.GetAtomWithIdx(oxy_idx)
            for nbr in oxy_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 1:
                    atoms_to_remove.add(nbr.GetIdx())

    # Surviving atom indices.
    surviving = [i for i in range(mol_h.GetNumAtoms()) if i not in atoms_to_remove]
    old_to_new: dict[int, int] = {old: new for new, old in enumerate(surviving)}

    # Build Cartesian positions, centered on origin.
    positions = np.array([list(conf.GetAtomPosition(i)) for i in surviving])
    centroid = positions.mean(axis=0)
    positions -= centroid

    # All positions (including removed atoms) for direction computation.
    all_positions = np.array([list(conf.GetAtomPosition(i)) for i in range(mol_h.GetNumAtoms())])
    all_positions -= centroid

    # Build atom list.
    atoms: list[_XyzAtom] = []
    for pos_i, old_idx in enumerate(surviving):
        sym = mol_h.GetAtomWithIdx(old_idx).GetSymbol()
        p = positions[pos_i]
        atoms.append(_XyzAtom(symbol=sym, x=p[0], y=p[1], z=p[2]))

    # Build bond list (surviving atoms only).
    bonds: list[_XyzBond] = []
    surviving_set = set(surviving)
    for bond in mol_h.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a in surviving_set and b in surviving_set:
            bt = _bond_type_char(bond)
            bonds.append(_XyzBond(old_to_new[a], old_to_new[b], bt))

    # Place X dummy atoms at each stripped carboxylate position.
    x_indices: list[int] = []
    for g in conn_info.carboxylate_groups:
        anchor_idx = g.anchor_idx
        anchor_new = old_to_new[anchor_idx]
        anchor_pos = positions[surviving.index(anchor_idx)]

        # Direction: anchor -> carboxylate carbon (from original coords).
        carboxylate_c_pos = all_positions[g.carbon_idx]
        anchor_orig_pos = all_positions[anchor_idx]
        direction = _outward_direction(carboxylate_c_pos, anchor_orig_pos)

        x_pos = anchor_pos + direction * _PORMAKE_X_DISTANCE
        x_new_idx = len(atoms)
        atoms.append(_XyzAtom(symbol="X", x=x_pos[0], y=x_pos[1], z=x_pos[2]))
        bonds.append(_XyzBond(x_new_idx, anchor_new, "S"))
        x_indices.append(x_new_idx)

    return atoms, bonds, x_indices


def _build_pormake_carboxylic_edge(
    mol_h,
    conf,
    conn_info: ConnectionInfo,
) -> tuple[list[_XyzAtom], list[_XyzBond], list[int]]:
    """Build Pormake atom/bond lists for a *carboxylic*-stripped edge.

    Strips the entire COOH group and places X dummy atoms at the anchor
    atoms, matching the TOBACCO carboxylic mode logic.
    """
    # Collect all atoms to remove: carboxylate C, both O atoms, and
    # H atoms bonded to those oxygens or to anchor atoms.
    atoms_to_remove: set[int] = set()
    anchor_set: set[int] = set()
    for g in conn_info.carboxylate_groups:
        atoms_to_remove.add(g.carbon_idx)
        atoms_to_remove.add(g.oxy_double_idx)
        atoms_to_remove.add(g.oxy_single_idx)
        anchor_set.add(g.anchor_idx)
        for oxy_idx in (g.oxy_double_idx, g.oxy_single_idx):
            oxy_atom = mol_h.GetAtomWithIdx(oxy_idx)
            for nbr in oxy_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 1:
                    atoms_to_remove.add(nbr.GetIdx())

    # Remove H atoms bonded to anchor atoms.
    for anchor_idx in anchor_set:
        anchor_atom = mol_h.GetAtomWithIdx(anchor_idx)
        for nbr in anchor_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 1:
                atoms_to_remove.add(nbr.GetIdx())

    # Surviving atom indices.
    surviving = [i for i in range(mol_h.GetNumAtoms()) if i not in atoms_to_remove]
    old_to_new: dict[int, int] = {old: new for new, old in enumerate(surviving)}

    # Build Cartesian positions, centered on origin.
    positions = np.array([list(conf.GetAtomPosition(i)) for i in surviving])
    centroid = positions.mean(axis=0)
    positions -= centroid

    # All positions (including removed atoms) for direction computation.
    all_positions = np.array([list(conf.GetAtomPosition(i)) for i in range(mol_h.GetNumAtoms())])
    all_positions -= centroid

    # Build atom list.
    atoms: list[_XyzAtom] = []
    for pos_i, old_idx in enumerate(surviving):
        sym = mol_h.GetAtomWithIdx(old_idx).GetSymbol()
        p = positions[pos_i]
        atoms.append(_XyzAtom(symbol=sym, x=p[0], y=p[1], z=p[2]))

    # Build bond list (surviving atoms only).
    bonds: list[_XyzBond] = []
    surviving_set = set(surviving)
    for bond in mol_h.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a in surviving_set and b in surviving_set:
            bt = _bond_type_char(bond)
            bonds.append(_XyzBond(old_to_new[a], old_to_new[b], bt))

    # Place X dummy atoms at each stripped carboxylate position.
    x_indices: list[int] = []
    for g in conn_info.carboxylate_groups:
        anchor_idx = g.anchor_idx
        anchor_new = old_to_new[anchor_idx]
        anchor_pos = positions[surviving.index(anchor_idx)]

        # Direction: anchor -> carboxylate carbon (from original coords).
        carboxylate_c_pos = all_positions[g.carbon_idx]
        anchor_orig_pos = all_positions[anchor_idx]
        direction = _outward_direction(carboxylate_c_pos, anchor_orig_pos)

        x_pos = anchor_pos + direction * _PORMAKE_X_DISTANCE
        x_new_idx = len(atoms)
        atoms.append(_XyzAtom(symbol="X", x=x_pos[0], y=x_pos[1], z=x_pos[2]))
        bonds.append(_XyzBond(x_new_idx, anchor_new, "S"))
        x_indices.append(x_new_idx)

    return atoms, bonds, x_indices


def _write_pormake_xyz(
    atoms: list[_XyzAtom],
    bonds: list[_XyzBond],
    x_indices: list[int],
    path: Path,
) -> None:
    """Write a Pormake-format extended XYZ file."""
    lines: list[str] = []

    # Line 1: atom count
    lines.append(str(len(atoms)))

    # Line 2: comment = space-separated X-atom indices
    lines.append("   ".join(str(i) for i in x_indices))

    # Atom block
    for a in atoms:
        lines.append(f"{a.symbol:<4s} {a.x:>10.4f} {a.y:>10.4f} {a.z:>10.4f}")

    # Bond block
    for b in bonds:
        lines.append(f"{b.idx_i:>4d} {b.idx_j:>4d} {b.bond_type}")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
