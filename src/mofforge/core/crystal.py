"""Crystal data structure wrapping pymatgen Structure with a NetworkX bond graph."""

from __future__ import annotations

import copy
import logging
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import networkx as nx
import numpy as np
from pymatgen.core import Lattice, Structure

from mofforge.io.cif import read_cif, write_cif
from mofforge.io.xyz import write_xyz
from mofforge.utils.config import clean_species as _clean_species
from mofforge.utils.periodic import wrap_coords

if TYPE_CHECKING:
    from mofforge.provenance import Provenance
    from mofforge.search.search import MatchResult

logger = logging.getLogger("mofforge")


@dataclass(frozen=True)
class PeriodicBond:
    """One undirected bond from atom i to an image of atom j.

    Each physical bond is stored once; a self-image bond contributes two
    neighbors. ``Crystal.bonds`` remains a simple graph for pattern matching.
    """

    i: int
    j: int
    image: tuple[int, int, int]
    distance: float


class Crystal:
    """A crystal structure with an associated bonding graph.

    Wraps a pymatgen Structure (lattice + species + fractional coordinates)
    and adds a NetworkX Graph for the bonding network.

    The ``_species_labels`` list holds the original species labels (which
    may include R-group tags like 'H!' or 'C!'). The pymatgen Structure
    only stores clean element symbols.

    """

    def __init__(
        self,
        name: str,
        structure: Structure,
        bonds: nx.Graph | None = None,
        provenance: Provenance | None = None,
        species_labels: list[str] | None = None,
        periodic_bonds: list[PeriodicBond] | None = None,
    ):
        self.name = name
        self.structure = structure
        self.bonds = bonds if bonds is not None else nx.Graph()
        self.provenance = provenance
        self.periodic_bonds = periodic_bonds

        # If species_labels provided, use them; otherwise derive from structure
        if species_labels is not None:
            self._species_labels = list(species_labels)
        else:
            self._species_labels = [
                site.specie.symbol if site.is_ordered else site.species_string
                for site in self.structure
            ]

        if len(self._species_labels) != self.n_atoms:
            raise ValueError("Species labels must have one entry per atom.")

        # Ensure bond graph has nodes for all atoms with species attributes
        if self.bonds.number_of_nodes() == 0 and self.n_atoms > 0:
            for i in range(self.n_atoms):
                self.bonds.add_node(i, species=self._species_labels[i])

    @property
    def n_atoms(self) -> int:
        """Number of atoms in the crystal."""
        return len(self.structure)

    @property
    def species(self) -> list[str]:
        """List of species labels (may include R-group tags like 'H!', 'C!')."""
        return list(self._species_labels)

    @property
    def frac_coords(self) -> np.ndarray:
        """Fractional coordinates, shape (N, 3)."""
        return self.structure.frac_coords.copy().reshape((-1, 3))

    @property
    def cart_coords(self) -> np.ndarray:
        """Cartesian coordinates, shape (N, 3) in Angstroms."""
        return self.structure.cart_coords.copy().reshape((-1, 3))

    @property
    def lattice(self) -> Lattice:
        """The pymatgen Lattice of this crystal."""
        return self.structure.lattice

    @property
    def n_bonds(self) -> int:
        """Number of bonds (edges in the bond graph)."""
        return self.bonds.number_of_edges()

    def coordination_number(self, index: int) -> int:
        """Count neighbors, including distinct periodic images of an atom."""
        if self.periodic_bonds is None:
            return self.bonds.degree(index)
        return sum(int(b.i == index) + int(b.j == index) for b in self.periodic_bonds)

    def bond_vectors(self, index: int) -> list[np.ndarray]:
        """Cartesian vectors from an atom to its bonded neighbors."""
        frac = self.frac_coords
        if self.periodic_bonds is None:
            from mofforge.utils.periodic import nearest_image

            return [
                self.to_cart(nearest_image(frac[j] - frac[index], self.lattice))
                for j in self.bonds.neighbors(index)
            ]
        vectors = []
        for bond in self.periodic_bonds:
            vector = self.to_cart(frac[bond.j] + bond.image - frac[bond.i])
            if bond.i == index:
                vectors.append(vector)
            if bond.j == index:
                vectors.append(-vector)
        return vectors

    def refresh_bond_geometry(self) -> None:
        """Refresh distances and image metadata without changing connectivity."""
        frac = self.frac_coords
        by_pair: dict[tuple[int, int], list[PeriodicBond]] = {}
        if self.periodic_bonds is not None:
            updated = []
            for bond in self.periodic_bonds:
                if bond.i != bond.j and not self.bonds.has_edge(bond.i, bond.j):
                    continue
                distance = float(
                    np.linalg.norm(self.to_cart(frac[bond.j] + bond.image - frac[bond.i]))
                )
                refreshed = PeriodicBond(bond.i, bond.j, bond.image, distance)
                updated.append(refreshed)
                by_pair.setdefault((bond.i, bond.j), []).append(refreshed)
            self.periodic_bonds = updated
        for u, v, data in self.bonds.edges(data=True):
            i, j = sorted((u, v))
            candidates = by_pair.get((i, j))
            if candidates:
                bond = min(candidates, key=lambda b: b.distance)
                distance, image = bond.distance, bond.image
            elif self.periodic_bonds is not None or data.get("cross_boundary"):
                distance, raw_image = self.lattice.get_distance_and_image(frac[i], frac[j])
                image = tuple(int(x) for x in raw_image)
                if self.periodic_bonds is not None:
                    self.periodic_bonds.append(PeriodicBond(i, j, image, float(distance)))
            else:
                distance = np.linalg.norm(self.to_cart(frac[j] - frac[i]))
                image = (0, 0, 0)
            data.update(distance=float(distance), image=image, cross_boundary=any(image))

    @classmethod
    def from_cif(cls, filepath: str | Path, name: str | None = None) -> Crystal:
        """Load a Crystal from a CIF file."""
        filepath = Path(filepath)
        if name is None:
            name = filepath.stem
        structure = read_cif(filepath)
        from mofforge.provenance import Provenance, file_hash, software_versions

        return cls(
            name=name,
            structure=structure,
            provenance=Provenance(
                parent=str(filepath),
                operation="load",
                input_hash=file_hash(filepath),
                parameters={"input_path": str(filepath)},
                software_versions=dict(software_versions()),
            ),
        )

    @classmethod
    def from_structure(
        cls,
        structure: Structure,
        name: str = "crystal",
    ) -> Crystal:
        """Create a Crystal from an existing pymatgen Structure."""
        return cls(name=name, structure=structure)

    @classmethod
    def from_xyz(
        cls,
        species: list[str],
        cart_coords: np.ndarray,
        name: str = "moiety",
        lattice: Lattice | None = None,
    ) -> Crystal:
        """Create a Crystal from species labels and Cartesian coordinates.

        Species labels may include R-group tags (e.g. 'H!', 'C!').
        These are stored separately; pymatgen only sees clean element symbols.
        """
        if lattice is None:
            lattice = Lattice.cubic(100.0)
        if len(species) == 0:
            structure = Structure(lattice, [], [])
            return cls(name=name, structure=structure, species_labels=[])

        # Clean species for pymatgen (strip ! tags)
        clean_species = [_clean_species(s) for s in species]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            structure = Structure(
                lattice,
                clean_species,
                cart_coords,
                coords_are_cartesian=True,
            )
        return cls(name=name, structure=structure, species_labels=species)

    @classmethod
    def empty(cls, name: str = "empty") -> Crystal:
        """Create an empty Crystal with no atoms."""
        lattice = Lattice.cubic(100.0)
        structure = Structure(lattice, [], [])
        return cls(name=name, structure=structure, species_labels=[])

    def __getitem__(self, indices: list[int] | np.ndarray) -> Crystal:
        """Extract a sub-crystal containing only the specified atom indices.

        Bonds between selected atoms are preserved. Node IDs are
        renumbered to 0..len(indices)-1.
        """
        indices = list(indices)
        if len(indices) != len(set(indices)):
            raise ValueError("Duplicate atom indices.")

        new_labels = [self._species_labels[i] for i in indices]
        new_structure = self.structure.copy()
        # Structure.from_sites preserves partial occupancy and site properties.
        if indices:
            new_structure = Structure.from_sites(
                [self.structure[i] for i in indices],
                properties=copy.deepcopy(self.structure.properties),
            )
        else:
            new_structure.remove_sites(range(self.n_atoms))

        # Build new bond graph with renumbered nodes
        old_to_new = {old: new for new, old in enumerate(indices)}
        new_bonds = nx.Graph()
        for new_idx, old_idx in enumerate(indices):
            new_bonds.add_node(new_idx, species=self._species_labels[old_idx])
        for u, v, data in self.bonds.edges(data=True):
            if u in old_to_new and v in old_to_new:
                new_bonds.add_edge(old_to_new[u], old_to_new[v], **copy.deepcopy(data))

        periodic_bonds = None
        if self.periodic_bonds is not None:
            periodic_bonds = []
            for b in self.periodic_bonds:
                if b.i in old_to_new and b.j in old_to_new:
                    i, j = old_to_new[b.i], old_to_new[b.j]
                    image = b.image
                    if i > j:
                        i, j = j, i
                        image = tuple(-x for x in image)
                    periodic_bonds.append(PeriodicBond(i, j, image, b.distance))

        return Crystal(
            name=f"subset_{self.name}",
            structure=new_structure,
            bonds=new_bonds,
            species_labels=new_labels,
            provenance=copy.deepcopy(self.provenance),
            periodic_bonds=periodic_bonds,
        )

    def __contains__(self, query: Crystal) -> bool:
        """Check whether *query* is a substructure of this crystal."""
        from mofforge.search.search import find_pattern

        result = find_pattern(query, self)
        return result.nb_locations() > 0

    def find(self, query: Crystal) -> MatchResult:
        """Search for *query* as a substructure and return the full result."""
        from mofforge.search.search import find_pattern

        return find_pattern(query, self)

    def __add__(self, other: Crystal) -> Crystal:
        """Combine two crystals into one (add atoms from other into self)."""
        if other.n_atoms == 0:
            return self.copy()
        if self.n_atoms == 0:
            return other.copy()

        # Combine species labels
        combined_labels = self._species_labels + other._species_labels

        # Convert other's frac coords to self's lattice
        other_cart = other.cart_coords
        other_frac_in_self = self.lattice.get_fractional_coords(other_cart)
        combined_frac = np.vstack([self.frac_coords, other_frac_in_self])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            new_structure = Structure(
                self.lattice,
                [site.species for site in self.structure]
                + [site.species for site in other.structure],
                combined_frac,
                site_properties={
                    key: copy.deepcopy(
                        self.structure.site_properties.get(key, [None] * self.n_atoms)
                    )
                    + copy.deepcopy(
                        other.structure.site_properties.get(key, [None] * other.n_atoms)
                    )
                    for key in self.structure.site_properties.keys()
                    | other.structure.site_properties.keys()
                },
                properties=copy.deepcopy(self.structure.properties),
            )

        # Merge bond graphs
        offset = self.n_atoms
        new_bonds = copy.deepcopy(self.bonds)
        for node, data in other.bonds.nodes(data=True):
            new_bonds.add_node(node + offset, **copy.deepcopy(data))
        for u, v, data in other.bonds.edges(data=True):
            new_bonds.add_edge(u + offset, v + offset, **copy.deepcopy(data))

        combined_provenance = copy.deepcopy(
            self.provenance if self.provenance is not None else other.provenance
        )
        if self.provenance is not None and other.provenance is not None:
            combined_provenance = combined_provenance.chain(copy.deepcopy(other.provenance))

        periodic_bonds = None
        if self.periodic_bonds is not None or other.periodic_bonds is not None:
            periodic_bonds = list(self.periodic_bonds or [])
            if other.periodic_bonds:
                if not np.allclose(self.lattice.matrix, other.lattice.matrix):
                    raise ValueError("Cannot combine periodic bond images from different lattices.")
                periodic_bonds.extend(
                    PeriodicBond(b.i + offset, b.j + offset, b.image, b.distance)
                    for b in other.periodic_bonds
                )

        combined = Crystal(
            name=self.name,
            structure=new_structure,
            bonds=new_bonds,
            provenance=combined_provenance,
            species_labels=combined_labels,
            periodic_bonds=periodic_bonds,
        )
        combined.refresh_bond_geometry()
        return combined

    def to_cart(self, frac: np.ndarray) -> np.ndarray:
        """Convert fractional coordinates to Cartesian."""
        return self.lattice.get_cartesian_coords(frac)

    def to_frac(self, cart: np.ndarray) -> np.ndarray:
        """Convert Cartesian coordinates to fractional."""
        return self.lattice.get_fractional_coords(cart)

    def wrap(self) -> Crystal:
        """Return a new Crystal with fractional coordinates wrapped to [0, 1)."""
        new_xtal = self.copy()
        wrapped = wrap_coords(new_xtal.frac_coords)
        new_xtal.set_frac_coords(wrapped)
        return new_xtal

    def set_frac_coords(self, new_coords: np.ndarray) -> None:
        """Update fractional coordinates in-place."""
        new_coords = np.asarray(new_coords, dtype=float)
        if new_coords.shape != (self.n_atoms, 3):
            raise ValueError(
                f"Expected coordinates shape ({self.n_atoms}, 3), got {new_coords.shape}."
            )
        if not np.isfinite(new_coords).all():
            raise ValueError("Coordinates must be finite.")
        shifts = self.frac_coords - new_coords
        # Account for exact cell translations (wrap/reassembly), while retaining
        # bond images for actual geometric displacements.
        shifts = np.where(np.isclose(shifts, np.round(shifts)), np.round(shifts), 0).astype(int)
        if self.periodic_bonds is not None:
            self.periodic_bonds = [
                PeriodicBond(
                    b.i,
                    b.j,
                    tuple(int(x) for x in np.array(b.image) + shifts[b.j] - shifts[b.i]),
                    b.distance,
                )
                for b in self.periodic_bonds
            ]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.structure = Structure(
                self.lattice,
                [site.species for site in self.structure],
                new_coords,
                site_properties=copy.deepcopy(self.structure.site_properties),
                properties=copy.deepcopy(self.structure.properties),
            )
        self.refresh_bond_geometry()

    def set_cart_coords(self, new_coords: np.ndarray) -> None:
        """Update Cartesian coordinates in-place (converted to fractional)."""
        frac = self.to_frac(new_coords)
        self.set_frac_coords(frac)

    def write_cif(self, filepath: str | Path) -> None:
        """Write this crystal to a CIF file."""
        write_cif(self.structure, filepath)
        if self.provenance is not None:
            from mofforge.provenance import write_manifest

            write_manifest(self, filepath)

    def write_xyz(self, filepath: str | Path, comment: str = "") -> None:
        """Write this crystal's atoms to an XYZ file (Cartesian coordinates)."""
        write_xyz(self._species_labels, self.cart_coords, filepath, comment)
        if self.provenance is not None:
            from mofforge.provenance import write_manifest

            write_manifest(self, filepath)

    def copy(self) -> Crystal:
        """Return a deep copy of this Crystal."""
        return Crystal(
            name=self.name,
            structure=self.structure.copy(),
            bonds=copy.deepcopy(self.bonds),
            provenance=copy.deepcopy(self.provenance),
            species_labels=list(self._species_labels),
            periodic_bonds=copy.deepcopy(self.periodic_bonds),
        )

    def __repr__(self) -> str:
        return f"Crystal('{self.name}', n_atoms={self.n_atoms}, n_bonds={self.n_bonds})"

    def __len__(self) -> int:
        return self.n_atoms
