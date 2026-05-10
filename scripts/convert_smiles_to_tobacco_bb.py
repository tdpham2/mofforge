"""Convert SMILES strings to TOBACCO-format edge building-block CIF files.

This script is a thin CLI wrapper around
:func:`mofforge.build.smiles_to_bb.smiles_to_tobacco_edge_cif`.
Connection points are detected automatically:

* **Carboxylate-terminated linkers** (``COO`` / ``COOH`` groups) --
  an ``X`` dummy atom (element ``Fr``) is placed at the centroid of
  each pair of carboxylate oxygens.

* **Direct aromatic linkers** (no carboxylate groups) -- the two
  terminal backbone atoms are relabelled with an ``X`` prefix and
  their hydrogen atoms are removed.

Usage::

    python convert_smiles_to_tobacco_bb.py "OC(=O)c1ccc(C(=O)O)cc1" -o BDC_edge.cif
    python convert_smiles_to_tobacco_bb.py "c1ccc(-c2ccccc2)cc1" -o biphenyl_edge.cif

Requires: ``rdkit`` (``pip install rdkit``).

Legacy behaviour
~~~~~~~~~~~~~~~~
The old semi-empirical workflow (OpenBabel + ASE + manual connection
sites) is preserved below in the ``legacy_*`` functions for reference.
Use ``--legacy`` to run the old code path (requires ``obabel``, ``ase``,
``networkx``).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert a SMILES string to a TOBACCO edge building-block CIF.",
    )
    parser.add_argument("smiles", help="SMILES string for the organic linker")
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output CIF path (default: <name>.cif)",
    )
    parser.add_argument(
        "-n", "--name",
        default="edge",
        help="CIF data-block name (default: 'edge')",
    )
    parser.add_argument(
        "--cell-length",
        type=float,
        default=40.0,
        help="Cubic cell side length in Angstrom (default: 40)",
    )
    parser.add_argument(
        "--legacy",
        action="store_true",
        help="Use the legacy (semi-empirical) conversion pipeline",
    )
    args = parser.parse_args()

    if args.legacy:
        print(
            "Legacy mode requires manual editing of connection_sites in "
            "this script.  See the legacy_* functions at the bottom of "
            "the file.",
            file=sys.stderr,
        )
        sys.exit(1)

    from mofforge.build.smiles_to_bb import (
        detect_connection_points,
        smiles_to_tobacco_edge_cif,
    )

    output = Path(args.output) if args.output else Path(f"{args.name}.cif")

    # Show detection info first
    info = detect_connection_points(args.smiles)
    print(f"SMILES:           {args.smiles}")
    print(f"Canonical:        {info.smiles}")
    print(f"Detection mode:   {info.mode}")
    print(f"Connection atoms: {info.connection_atom_indices}")
    if info.carboxylate_groups:
        for i, g in enumerate(info.carboxylate_groups):
            print(
                f"  Carboxylate {i + 1}: C={g.carbon_idx}, "
                f"O={g.oxy_double_idx}, O(H)={g.oxy_single_idx}, "
                f"anchor={g.anchor_idx}"
            )

    result = smiles_to_tobacco_edge_cif(
        args.smiles,
        output_path=output,
        name=args.name,
        cell_length=args.cell_length,
    )
    print(f"Wrote: {result}")


# ------------------------------------------------------------------ #
# Legacy functions (kept for reference)
# ------------------------------------------------------------------ #
# The code below is the original semi-empirical pipeline that required
# OpenBabel, ASE, and manual specification of connection sites.
# It is NOT used by the new automatic pipeline above.

import re

PT = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
    'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr',
    'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
    'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
    'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'Hf',
    'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
    'Po', 'At', 'Rn', 'Fr', 'Ra', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm',
    'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th',
    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md',
    'No', 'Lr', 'FG', 'X',
]


def _nn(string):
    return re.sub('[^a-zA-Z]', '', string)


def _isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def _iscoord(line):
    return _nn(line[0]) in PT and line[1] in PT and all(map(_isfloat, line[2:5]))


def _isbond(line):
    return _nn(line[0]) in PT and _nn(line[1]) in PT and '.' in line[3] and len(line) == 5


def legacy_update_cif_with_connection_site(cif_in, connection_sites, cif_out):
    """Replace connection-site atoms with X labels and remove bonded H atoms.

    ``connection_sites`` must include a trailing space, e.g. ``['C1 ', 'C11 ']``.
    """
    import subprocess

    with open(cif_in, 'r') as f:
        data = [i.strip() for i in f.readlines() if i.strip()]

    X_sites = ['X' + str(i[1:]) for i in connection_sites]
    rem_index = []
    H_site = []

    for site in connection_sites:
        for k, line in enumerate(data):
            line = line.split()
            if _isbond(line) and site[:-1] in line:
                if 'H' in line[0]:
                    rem_index.append(k)
                    H_site.append(line[0])
                if 'H' in line[1]:
                    rem_index.append(k)
                    H_site.append(line[1])

    for k, line in enumerate(data):
        line = line.split()
        for i in H_site:
            if _iscoord(line) and i in line:
                rem_index.append(k)

    for index in sorted(set(rem_index), reverse=True):
        del data[index]

    with open(cif_out, 'w') as f:
        for line in data:
            f.write(line + '\n')

    for x, c in zip(X_sites, connection_sites):
        subprocess.call("sed -i 's|{}|{}|g' {}".format(c, x, cif_out), shell=True)


def legacy_update_aromatic_bond(cif_in, cif_out):
    """Relabel single bonds in rings as aromatic using NetworkX cycle detection."""
    import networkx as nx

    with open(cif_in, 'r') as f:
        data = [i.strip() for i in f.readlines() if i.strip()]

    nodes = []
    edges = []
    for line in data:
        parts = line.split()
        if _iscoord(parts):
            nodes.append(parts[0])
        if _isbond(parts):
            edges.append((parts[0], parts[1]))

    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    basis_cycles = nx.cycle_basis(G)

    new_file = []
    for line_ in data:
        parts = line_.split()
        if _isbond(parts):
            new_line = line_
            for cycle in basis_cycles:
                if parts[0] in cycle and parts[1] in cycle and parts[-1] == 'S':
                    new_line = '    '.join([parts[0], parts[1], parts[2], parts[3], 'A'])
                    break
            new_file.append(new_line)
        else:
            new_file.append(line_)

    with open(cif_out, 'w') as f:
        for line in new_file:
            f.write(line + '\n')


if __name__ == "__main__":
    main()
