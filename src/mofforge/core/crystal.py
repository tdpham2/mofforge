"""Crystal data structure wrapping pymatgen Structure with a NetworkX bond graph."""

from __future__ import annotations

import copy
import logging
import re
import warnings
from pathlib import Path
from typing import TYPE_CHECKING

import networkx as nx
import numpy as np
from pymatgen.core import Lattice, Structure

from mofforge.io.cif import read_cif, write_cif
from mofforge.io.xyz import write_xyz
from mofforge.utils.periodic import wrap_coords

if TYPE_CHECKING:
    from mofforge.provenance import Provenance
    from mofforge.search.search import MatchResult

logger = logging.getLogger("mofforge")

# Regex to strip R-group tag and oxidation state info, keeping only the element symbol
_ELEMENT_RE = re.compile(r"^([A-Z][a-z]?)")


def _clean_species(species: str) -> str:
    """Strip R-group tag (!) from a species label for pymatgen compatibility.

    'H!' -> 'H', 'C!' -> 'C', 'Zn' -> 'Zn'.
    """
    m = _ELEMENT_RE.match(species)
    if m:
        return m.group(1)
    return species


class Crystal:
    """A crystal structure with an associated bonding graph.

    Wraps a pymatgen Structure (lattice + species + fractional coordinates)
    and adds a NetworkX Graph for the bonding network.

    The ``_species_labels`` list holds the original species labels (which
    may include R-group tags like 'H!' or 'C!'). The pymatgen Structure
    only stores clean element symbols.

    Attributes:
        name: Identifier for this crystal.
        structure: pymatgen Structure (lattice, clean species, fractional coords).
        bonds: NetworkX Graph with integer node IDs matching atom indices.
            Node attributes: 'species' (str — original label with tags).
            Edge attributes: 'distance' (float, Angstroms), 'cross_boundary' (bool).
        provenance: Optional provenance metadata tracking modifications.
    """

    def __init__(
        self,
        name: str,
        structure: Structure,
        bonds: nx.Graph | None = None,
        provenance: Provenance | None = None,
        species_labels: list[str] | None = None,
    ):
        self.name = name
        self.structure = structure
        self.bonds = bonds if bonds is not None else nx.Graph()
        self.provenance = provenance

        # If species_labels provided, use them; otherwise derive from structure
        if species_labels is not None:
            self._species_labels = list(species_labels)
        else:
            self._species_labels = [str(s) for s in self.structure.species]

        # Ensure bond graph has nodes for all atoms with species attributes
        if self.bonds.number_of_nodes() == 0 and self.n_atoms > 0:
            for i in range(self.n_atoms):
                self.bonds.add_node(i, species=self._species_labels[i])

    # -------------------------------------------------------------------------
    # Properties
    # -------------------------------------------------------------------------

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
        return self.structure.frac_coords.copy()

    @property
    def cart_coords(self) -> np.ndarray:
        """Cartesian coordinates, shape (N, 3) in Angstroms."""
        return self.structure.cart_coords.copy()

    @property
    def lattice(self) -> Lattice:
        """The pymatgen Lattice of this crystal."""
        return self.structure.lattice

    @property
    def n_bonds(self) -> int:
        """Number of bonds (edges in the bond graph)."""
        return self.bonds.number_of_edges()

    # -------------------------------------------------------------------------
    # Constructors
    # -------------------------------------------------------------------------

    @classmethod
    def from_cif(cls, filepath: str | Path, name: str | None = None) -> Crystal:
        """Load a Crystal from a CIF file.

        Args:
            filepath: Path to the CIF file.
            name: Optional name; defaults to the filename stem.

        Returns:
            Crystal with no bonds (use infer_bonds to add them).
        """
        filepath = Path(filepath)
        if name is None:
            name = filepath.stem
        structure = read_cif(filepath)
        return cls(name=name, structure=structure)

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

        Args:
            species: List of species strings (may include '!' tags).
            cart_coords: Cartesian coordinates, shape (N, 3).
            name: Name for the crystal.
            lattice: Optional lattice; defaults to a 100 A cube.
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

    # -------------------------------------------------------------------------
    # Indexing and slicing
    # -------------------------------------------------------------------------

    def __getitem__(self, indices: list[int] | np.ndarray) -> Crystal:
        """Extract a sub-crystal containing only the specified atom indices.

        Bonds between selected atoms are preserved. Node IDs are
        renumbered to 0..len(indices)-1.

        Args:
            indices: List of atom indices to extract.

        Returns:
            New Crystal containing only the specified atoms.
        """
        indices = list(indices)
        if not indices:
            return Crystal.empty(name=f"subset_{self.name}")
        if len(indices) != len(set(indices)):
            raise ValueError(
                "Duplicate atom indices are not supported in __getitem__. "
                f"Got {len(indices)} indices but only {len(set(indices))} are unique."
            )

        # Build new structure with clean species
        new_labels = [self._species_labels[i] for i in indices]
        new_clean = [_clean_species(s) for s in new_labels]
        new_frac_coords = self.frac_coords[indices]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            new_structure = Structure(
                self.lattice,
                new_clean,
                new_frac_coords,
            )

        # Build new bond graph with renumbered nodes
        old_to_new = {old: new for new, old in enumerate(indices)}
        new_bonds = nx.Graph()
        for new_idx, old_idx in enumerate(indices):
            new_bonds.add_node(new_idx, species=self._species_labels[old_idx])
        for u, v, data in self.bonds.edges(data=True):
            if u in old_to_new and v in old_to_new:
                new_bonds.add_edge(old_to_new[u], old_to_new[v], **data)

        return Crystal(
            name=f"subset_{self.name}",
            structure=new_structure,
            bonds=new_bonds,
            species_labels=new_labels,
        )

    def __contains__(self, query: Crystal) -> bool:
        """Check whether *query* is a substructure of this crystal.

        Args:
            query: The fragment to search for.

        Returns:
            True if at least one match is found.
        """
        from mofforge.search.search import find_pattern

        result = find_pattern(query, self)
        return result.nb_locations() > 0

    def find(self, query: Crystal) -> MatchResult:
        """Search for *query* as a substructure and return the full result.

        Args:
            query: The fragment to search for.

        Returns:
            MatchResult object with isomorphism results.
        """
        from mofforge.search.search import find_pattern

        return find_pattern(query, self)

    # -------------------------------------------------------------------------
    # Combination
    # -------------------------------------------------------------------------

    def __add__(self, other: Crystal) -> Crystal:
        """Combine two crystals into one (add atoms from other into self).

        Uses self's lattice. Other's atoms are converted to fractional
        coordinates in self's lattice. Bond graphs are merged with
        other's nodes offset by self.n_atoms.

        Args:
            other: Crystal whose atoms to add.

        Returns:
            New Crystal with combined atoms and bonds.
        """
        if other.n_atoms == 0:
            return self.copy()
        if self.n_atoms == 0:
            return other.copy()

        # Combine species labels
        combined_labels = self._species_labels + other._species_labels
        combined_clean = [_clean_species(s) for s in combined_labels]

        # Convert other's frac coords to self's lattice
        other_cart = other.cart_coords
        other_frac_in_self = self.lattice.get_fractional_coords(other_cart)
        combined_frac = np.vstack([self.frac_coords, other_frac_in_self])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            new_structure = Structure(
                self.lattice,
                combined_clean,
                combined_frac,
            )

        # Merge bond graphs
        offset = self.n_atoms
        new_bonds = self.bonds.copy()
        for node, data in other.bonds.nodes(data=True):
            new_bonds.add_node(node + offset, **data)
        for u, v, data in other.bonds.edges(data=True):
            new_bonds.add_edge(u + offset, v + offset, **data)

        # Warn if provenance is being lost from the other crystal
        if other.provenance is not None and self.provenance is not None:
            logger.debug(
                "Crystal.__add__: both operands have provenance; "
                "only the left operand's provenance is kept."
            )

        combined_provenance = self.provenance if self.provenance is not None else other.provenance

        return Crystal(
            name=self.name,
            structure=new_structure,
            bonds=new_bonds,
            provenance=combined_provenance,
            species_labels=combined_labels,
        )

    # -------------------------------------------------------------------------
    # Coordinate operations
    # -------------------------------------------------------------------------

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
        """Update fractional coordinates in-place.

        Args:
            new_coords: New fractional coordinates, shape (N, 3).
        """
        if new_coords.shape != (self.n_atoms, 3):
            raise ValueError(
                f"Expected coordinates shape ({self.n_atoms}, 3), got {new_coords.shape}."
            )
        # Rebuild the structure in one shot instead of O(N) item assignment
        clean = [_clean_species(s) for s in self._species_labels]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.structure = Structure(
                self.lattice,
                clean,
                new_coords,
            )

    def set_cart_coords(self, new_coords: np.ndarray) -> None:
        """Update Cartesian coordinates in-place (converted to fractional).

        Args:
            new_coords: New Cartesian coordinates, shape (N, 3).
        """
        frac = self.to_frac(new_coords)
        self.set_frac_coords(frac)

    # -------------------------------------------------------------------------
    # I/O
    # -------------------------------------------------------------------------

    def write_cif(self, filepath: str | Path) -> None:
        """Write this crystal to a CIF file."""
        write_cif(self.structure, filepath)

    def write_xyz(self, filepath: str | Path, comment: str = "") -> None:
        """Write this crystal's atoms to an XYZ file (Cartesian coordinates)."""
        write_xyz(self._species_labels, self.cart_coords, filepath, comment)

    # -------------------------------------------------------------------------
    # Copy and display
    # -------------------------------------------------------------------------

    def copy(self) -> Crystal:
        """Return a deep copy of this Crystal."""
        return Crystal(
            name=self.name,
            structure=self.structure.copy(),
            bonds=copy.deepcopy(self.bonds),
            provenance=copy.deepcopy(self.provenance),
            species_labels=list(self._species_labels),
        )

    def __repr__(self) -> str:
        return f"Crystal('{self.name}', n_atoms={self.n_atoms}, n_bonds={self.n_bonds})"

    def __len__(self) -> int:
        return self.n_atoms
