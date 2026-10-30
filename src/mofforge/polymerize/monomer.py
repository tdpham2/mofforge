"""Explicit molecular preparation, with RDKit imported only when needed."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from pymatgen.core import Element

from mofforge.core.crystal import Crystal
from mofforge.io.xyz import read_xyz
from mofforge.polymerize.base import Monomer, positive
from mofforge.polymerize.state import Atom, Bond, components
from mofforge.provenance import effective_seed, file_hash


def chemistry():
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise ImportError("Molecular preparation needs RDKit; install mofforge[pop].") from exc
    return Chem


@dataclass
class PreparedMonomer:
    name: str
    template_id: str
    atoms: list[Atom]
    coordinates: np.ndarray
    bonds: list[Bond]
    provenance: dict

    @property
    def n_atoms(self):
        return len(self.atoms)

    @property
    def species(self):
        return [a.species for a in self.atoms]

    @property
    def molar_mass(self):
        return sum(a.mass for a in self.atoms)

    def to_dict(self):
        return {
            "name": self.name,
            "template_id": self.template_id,
            "atoms": [asdict(a) for a in self.atoms],
            "coordinates": self.coordinates.tolist(),
            "bonds": [asdict(b) for b in self.bonds],
            "provenance": self.provenance,
        }


def _graph_molecule(species, graph):
    Chem = chemistry()
    if not isinstance(graph, dict) or set(graph) != {"atoms", "bonds"}:
        raise ValueError(
            "Molecular graph requires atoms and bonds lists with explicit labels/orders."
        )
    atoms = graph["atoms"]
    if len(atoms) != len(species) or [a["species"] for a in atoms] != species:
        raise ValueError("Molecular graph species/order must map the coordinate atoms exactly.")
    labels = [a["label"] for a in atoms]
    if len(set(labels)) != len(labels) or not all(isinstance(s, str) and s for s in labels):
        raise ValueError("Template atom labels must be unique nonempty strings.")
    mol = Chem.RWMol()
    for record in atoms:
        atom = Chem.Atom(record["species"])
        charge = record.get("formal_charge")
        if charge is not None:
            if type(charge) is not int:
                raise ValueError("Formal charge must be an integer.")
            atom.SetFormalCharge(charge)
        if record.get("isotope"):
            atom.SetIsotope(record["isotope"])
        mol.AddAtom(atom)
    lookup = {label: i for i, label in enumerate(labels)}
    types = {
        1: Chem.BondType.SINGLE,
        1.5: Chem.BondType.AROMATIC,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }
    for bond in graph["bonds"]:
        if set(bond) != {"atom1", "atom2", "order"} or bond["order"] not in types:
            raise ValueError("Each molecular bond requires atom1, atom2 and an explicit order.")
        mol.AddBond(lookup[bond["atom1"]], lookup[bond["atom2"]], types[bond["order"]])
    mol = mol.GetMol()
    Chem.SanitizeMol(mol)
    return mol, labels, [a.get("formal_charge") for a in atoms]


def prepare_monomer(monomer: Monomer, *, template_id="t0000", random_seed=None):
    Chem = chemistry()
    from rdkit.Chem import AllChem

    source = monomer.source
    provenance = {}
    supplied_bonds = None
    supplied_masses = None
    supplied = isinstance(source, (Crystal, Path))
    path = None
    if isinstance(source, (str, Path)):
        value = str(source)
        # Extension testing does not interpret SMILES ring/bond punctuation as a filename.
        is_file = Path(value).is_file() if len(value) < 240 else False
        supplied = (
            supplied or is_file or value.lower().endswith((".xyz", ".mol", ".sdf", ".cif", ".mol2"))
        )
        if supplied:
            path = Path(value)
    if isinstance(source, Crystal):
        if not source.structure.is_ordered:
            raise ValueError("Molecular Crystal inputs must have fully occupied atomic sites.")
        labels = source.structure.site_properties.get(
            "template_label", [f"a{i}" for i in range(source.n_atoms)]
        )
        # A constructed finite chain contains repeated original template labels.
        # Its full stable atom IDs provide unique labels for the new chain template.
        if len(set(labels)) != len(labels):
            labels = source.structure.site_properties.get(
                "atom_id", [f"a{i}" for i in range(source.n_atoms)]
            )
        charges = source.structure.site_properties.get("formal_charge", [None] * source.n_atoms)
        records = [
            {"label": label, "species": sp, "formal_charge": ch}
            for label, sp, ch in zip(labels, source.species, charges, strict=True)
        ]
        raw_bonds = source.periodic_bonds
        if raw_bonds is None:
            graph_bonds = [
                Bond(labels[i], labels[j], data.get("order"))
                for i, j, data in source.bonds.edges(data=True)
            ]
        else:
            graph_bonds = [Bond(labels[b.i], labels[b.j], b.order, b.image) for b in raw_bonds]
        dummy_atoms = [
            Atom(label, template_id, template_id, label, sp, float(Element(sp).atomic_mass), charge)
            for label, sp, charge in zip(labels, source.species, charges, strict=True)
        ]
        comps = components(dummy_atoms, graph_bonds)
        if len(comps) != 1 or comps[0][2] != 0:
            raise ValueError("A molecular Crystal must contain one finite connected component.")
        offsets = comps[0][1]
        coords = (
            source.cart_coords
            + np.array([offsets[label] for label in labels]) @ source.lattice.matrix
        )
        graph = {
            "atoms": records,
            "bonds": [{"atom1": b.i, "atom2": b.j, "order": b.order} for b in graph_bonds],
        }
        mol, labels, charges = _graph_molecule(source.species, graph)
        supplied_bonds = graph["bonds"]
        supplied_masses = source.structure.site_properties.get("atom_mass")
        provenance = {"source_type": "Crystal", "name": source.name}
    elif supplied:
        if path is None or not path.is_file():
            raise ValueError(f"Molecular source file not found: {path}")
        provenance = {"source_type": path.suffix.lower(), "source_sha256": file_hash(path)}
        if path.suffix.lower() == ".xyz":
            species, coords = read_xyz(path)
            mol, labels, charges = _graph_molecule(species, monomer.graph)
            supplied_bonds = monomer.graph["bonds"]
        elif path.suffix.lower() in (".mol", ".sdf"):
            if path.suffix.lower() == ".mol":
                mol = Chem.MolFromMolFile(str(path), removeHs=False)
            else:
                molecules = list(Chem.SDMolSupplier(str(path), removeHs=False))
                if len(molecules) != 1:
                    raise ValueError("SDF input must contain exactly one molecule.")
                mol = molecules[0]
            if mol is None or not mol.GetNumConformers():
                raise ValueError(
                    "MOL/SDF input requires a valid molecule and supplied coordinates."
                )
            coords = np.array(mol.GetConformer().GetPositions())
            labels = _labels(mol)
            charges = [a.GetFormalCharge() for a in mol.GetAtoms()]
        else:
            raise ValueError("Use SMILES, MOL/SDF, a finite Crystal, or XYZ with a mapped graph.")
    else:
        params = Chem.SmilesParserParams()
        params.removeHs = False
        mol = Chem.MolFromSmiles(str(source), params)
        if mol is None:
            raise ValueError(f"RDKit could not parse SMILES: {source!r}")
        labels = _labels(mol)
        original_count = mol.GetNumAtoms()
        mol = Chem.AddHs(mol)
        hydrogen_counts = {}
        for atom in list(mol.GetAtoms())[original_count:]:
            parent = labels[atom.GetNeighbors()[0].GetIdx()]
            hydrogen_counts[parent] = hydrogen_counts.get(parent, 0) + 1
            labels.append(f"h:{parent}:{hydrogen_counts[parent]}")
        seed = effective_seed(random_seed)
        embed = AllChem.ETKDGv3()
        embed.randomSeed = seed
        if AllChem.EmbedMolecule(mol, embed) != 0:
            embed.useRandomCoords = True
            if AllChem.EmbedMolecule(mol, embed) != 0:
                raise ValueError("RDKit conformer embedding failed after two attempts.")
        if not AllChem.UFFHasAllMoleculeParams(mol):
            raise ValueError("SMILES geometry cleanup lacks UFF parameters; supply coordinates.")
        if AllChem.UFFOptimizeMolecule(mol, maxIters=2000) != 0:
            raise ValueError("Bounded SMILES geometry cleanup did not converge.")
        coords = np.array(mol.GetConformer().GetPositions())
        charges = [a.GetFormalCharge() for a in mol.GetAtoms()]
        provenance = {"source_type": "smiles", "smiles": str(source), "embedding_seed": seed}
    if mol.GetNumAtoms() == 0 or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Each component template must be one nonempty connected molecule.")
    if any(a.GetNumImplicitHs() or a.GetNumExplicitHs() for a in mol.GetAtoms()):
        raise ValueError(
            "Coordinate inputs require explicit hydrogen atoms and coordinates; "
            "prepare a complete molecule or use SMILES."
        )
    if coords.shape != (mol.GetNumAtoms(), 3) or not np.isfinite(coords).all():
        raise ValueError("Molecular coordinates must be finite and map every atom.")
    if len(labels) != len(set(labels)):
        raise ValueError("Duplicate template atom labels (including atom maps).")
    masses = (
        [float(atom.GetMass()) for atom in mol.GetAtoms()]
        if supplied_masses is None
        else [positive(mass, "atom mass") for mass in supplied_masses]
    )
    atoms = [
        Atom(label, template_id, template_id, label, atom.GetSymbol(), mass, charge)
        for atom, label, charge, mass in zip(mol.GetAtoms(), labels, charges, masses, strict=True)
    ]
    bonds = []
    by_label = {label: i for i, label in enumerate(labels)}
    records = (
        supplied_bonds
        if supplied_bonds is not None
        else [
            {
                "atom1": labels[b.GetBeginAtomIdx()],
                "atom2": labels[b.GetEndAtomIdx()],
                "order": b.GetBondTypeAsDouble(),
            }
            for b in mol.GetBonds()
        ]
    )
    for record in records:
        i, j = by_label[record["atom1"]], by_label[record["atom2"]]
        distance = float(np.linalg.norm(coords[j] - coords[i]))
        if distance < 0.1:
            raise ValueError("Supplied molecular geometry has overlapping bonded atoms.")
        bonds.append(
            Bond(
                labels[i],
                labels[j],
                record["order"],
                length_range=(distance * 0.65, distance * 1.35),
            )
        )
    for connector in monomer.connectors:
        if not {connector.atom, connector.anchor, *connector.replaceable} <= set(labels):
            raise ValueError(f"Connector {connector.label!r} references unknown template labels.")
        if (
            np.linalg.norm(coords[by_label[connector.anchor]] - coords[by_label[connector.atom]])
            < 1e-6
        ):
            raise ValueError("Connector orientation anchor coincides with attachment atom.")
        if sum({bond.i, bond.j} == {connector.atom, connector.anchor} for bond in bonds) != 1:
            raise ValueError("Connector orientation anchor must be a uniquely bonded neighbor.")
    return PreparedMonomer(monomer.name, template_id, atoms, coords, bonds, provenance)


def _labels(mol):
    return [
        f"map:{a.GetAtomMapNum()}" if a.GetAtomMapNum() else f"a{a.GetIdx()}"
        for a in mol.GetAtoms()
    ]


def check_valence(atoms, bonds):
    """Check local valence using all image neighbors, without a quotient RDKit graph."""
    Chem = chemistry()
    types = {
        1: Chem.BondType.SINGLE,
        1.5: Chem.BondType.AROMATIC,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }
    incident = {a.id: [] for a in atoms}
    for bond in bonds:
        incident[bond.i].append(bond.order)
        incident[bond.j].append(bond.order)
    for atom in atoms:
        orders = incident[atom.id]
        # RDKit can round aromatic valence down on a disconnected star graph.
        # Neutral aromatic carbon has one pi bond across its aromatic neighbors;
        # three aromatic neighbors are allowed at a fused-ring junction.
        if atom.species == "C" and 1.5 in orders:
            sigma_and_pi = sum(1 if order == 1.5 else order for order in orders) + 1
            if len(orders) > 3 or (not atom.formal_charge and sigma_and_pi > 4):
                raise ValueError(f"Invalid local valence at aromatic atom {atom.id}.")
        local = Chem.RWMol()
        center = Chem.Atom(atom.species)
        center.SetFormalCharge(atom.formal_charge or 0)
        center.SetNoImplicit(True)
        local.AddAtom(center)
        for order in orders:
            # Only the central atom's property cache is validated.
            neighbor = local.AddAtom(Chem.Atom("C"))
            local.AddBond(0, neighbor, types[order])
        try:
            local.GetAtomWithIdx(0).UpdatePropertyCache(strict=True)
        except Exception as exc:
            raise ValueError(f"Invalid local valence at atom {atom.id}.") from exc
