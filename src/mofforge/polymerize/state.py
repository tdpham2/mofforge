"""Engine-independent periodic topology, identity, and verified box persistence."""

from __future__ import annotations

import copy
import json
import os
import tempfile
from collections import Counter, deque
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import networkx as nx
import numpy as np
from pymatgen.core import Lattice, Structure

from mofforge.core.crystal import Crystal, PeriodicBond
from mofforge.polymerize.base import box_lengths as validate_lengths
from mofforge.polymerize.base import integer, positive
from mofforge.provenance import content_hash, file_hash

SCHEMA_VERSION = 1
AVOGADRO = 6.02214076e23


@dataclass(frozen=True)
class Atom:
    id: str
    template_id: str
    instance_id: str
    label: str
    species: str
    mass: float
    formal_charge: int | None = None


@dataclass(frozen=True)
class Bond:
    i: str
    j: str
    order: float
    image: tuple[int, int, int] = (0, 0, 0)
    length_range: tuple[float, float] | None = None

    def __post_init__(self):
        if len(self.image) != 3 or any(type(x) is not int for x in self.image):
            raise ValueError("Bond image must contain three integers.")
        image = tuple(self.image)
        if self.i > self.j or (self.i == self.j and image < (0, 0, 0)):
            i, j = self.j, self.i
            image = tuple(-x for x in image)
            object.__setattr__(self, "i", i)
            object.__setattr__(self, "j", j)
        if self.i == self.j and image == (0, 0, 0):
            raise ValueError("A bond cannot join an atom to itself in the same image.")
        if isinstance(self.order, bool) or self.order not in (1, 1.5, 2, 3):
            raise ValueError("Bond order must be 1, 1.5, 2, or 3.")
        object.__setattr__(self, "image", image)
        if self.length_range is not None:
            low, high = self.length_range
            positive(low, "minimum bond length")
            if positive(high, "maximum bond length") < low:
                raise ValueError("Invalid bond length range.")
            object.__setattr__(self, "length_range", (low, high))

    @property
    def key(self):
        return self.i, self.j, self.image


@dataclass
class Site:
    id: str
    instance_id: str
    label: str
    atom: str
    anchor: str
    role: str
    replaceable: list[str] = field(default_factory=list)
    consumed: bool = False


def adjacency(atoms, bonds):
    result = {a.id: [] for a in atoms}
    for bond in bonds:
        result[bond.i].append((bond.j, np.array(bond.image, dtype=int)))
        result[bond.j].append((bond.i, -np.array(bond.image, dtype=int)))
    return result


def components(atoms, bonds):
    """Return connected components, image potentials, and cycle-translation rank."""
    graph = adjacency(atoms, bonds)
    remaining = set(graph)
    result = []
    while remaining:
        root = min(remaining)
        offsets = {root: np.zeros(3, dtype=int)}
        queue = deque([root])
        cycles = []
        while queue:
            i = queue.popleft()
            for j, image in graph[i]:
                candidate = offsets[i] + image
                if j not in offsets:
                    offsets[j] = candidate
                    queue.append(j)
                else:
                    delta = candidate - offsets[j]
                    if np.any(delta):
                        cycles.append(delta)
        remaining.difference_update(offsets)
        rank = int(np.linalg.matrix_rank(np.array(cycles, dtype=float))) if cycles else 0
        result.append((set(offsets), offsets, rank))
    return result


def finite_cycle_too_small(atoms, bonds, i, j, image, minimum):
    """Bounded BFS in the lifted periodic graph; quotient cycles are not rings."""
    graph = adjacency(atoms, bonds)
    target = (j, tuple(image))
    queue = deque([(i, (0, 0, 0), 0)])
    seen = {(i, (0, 0, 0))}
    while queue:
        node, offset, depth = queue.popleft()
        if (node, offset) == target:
            return depth + 1 < minimum
        if depth >= minimum - 2:
            continue
        for neighbor, step in graph[node]:
            nxt = (neighbor, tuple(int(x) for x in np.array(offset) + step))
            if nxt not in seen:
                seen.add(nxt)
                queue.append((*nxt, depth + 1))
    return False


@dataclass
class ConstructionState:
    atoms: list[Atom]
    coordinates: np.ndarray
    box_lengths: tuple[float, float, float]
    bonds: list[Bond]
    connectors: list[Site] = field(default_factory=list)
    initial_connector_count: int = 0
    metadata: dict = field(default_factory=dict)
    events: list[dict] = field(default_factory=list)
    provenance: list[dict] = field(default_factory=list)
    periodicity: tuple[bool, bool, bool] = (True, True, True)

    def __post_init__(self):
        self.coordinates = np.array(self.coordinates, dtype=float, copy=True)
        self.box_lengths = validate_lengths(self.box_lengths)
        self.periodicity = tuple(self.periodicity)
        self.check()

    def check(self):
        if len(self.periodicity) != 3 or any(p is not True for p in self.periodicity):
            raise ValueError("Native boxes must be fully periodic and orthorhombic.")
        if not self.atoms or self.coordinates.shape != (len(self.atoms), 3):
            raise ValueError("A box needs atoms and exactly one coordinate per atom.")
        if not np.isfinite(self.coordinates).all():
            raise ValueError("Coordinates must be finite.")
        validate_lengths(self.box_lengths)
        ids = {a.id for a in self.atoms}
        if len(ids) != len(self.atoms):
            raise ValueError("Duplicate stable atom IDs.")
        if len({(a.instance_id, a.label) for a in self.atoms}) != len(self.atoms):
            raise ValueError("Duplicate template labels in a molecular instance.")
        for atom in self.atoms:
            if any(
                not isinstance(x, str) or not x
                for x in (atom.id, atom.template_id, atom.instance_id, atom.label, atom.species)
            ):
                raise ValueError("Atom identifiers and species must be nonempty strings.")
            positive(atom.mass, "atom mass")
            if atom.formal_charge is not None and type(atom.formal_charge) is not int:
                raise ValueError("Formal charges must be integer or unknown (null).")
        if len({b.key for b in self.bonds}) != len(self.bonds):
            raise ValueError("Duplicate periodic bonds.")
        if any(b.i not in ids or b.j not in ids for b in self.bonds):
            raise ValueError("Bond references a missing atom.")
        integer(self.initial_connector_count, "initial_connector_count", 0)
        if self.initial_connector_count != len(self.connectors):
            raise ValueError("Original connector count must match the connector ledger.")
        if len({s.id for s in self.connectors}) != len(self.connectors):
            raise ValueError("Duplicate connector IDs.")
        for site in self.connectors:
            if type(site.consumed) is not bool:
                raise ValueError("Connector consumption must be boolean.")
            if not site.consumed:
                referenced = {site.atom, site.anchor, *site.replaceable}
                if not referenced <= ids or site.atom == site.anchor:
                    raise ValueError("An active connector references missing or invalid atoms.")
                if any(a.instance_id != site.instance_id for a in self.atoms if a.id in referenced):
                    raise ValueError("Connector atoms must belong to its instance.")

    @property
    def index(self):
        return {a.id: i for i, a in enumerate(self.atoms)}

    @property
    def conversion(self):
        return (
            sum(s.consumed for s in self.connectors) / self.initial_connector_count
            if self.initial_connector_count
            else 0.0
        )

    @property
    def component_info(self):
        return components(self.atoms, self.bonds)

    @property
    def component_records(self):
        return [
            {"id": min(members), "atom_ids": sorted(members), "periodic_connectivity_rank": rank}
            for members, _, rank in self.component_info
        ]

    @property
    def statistics(self):
        mass = sum(a.mass for a in self.atoms)
        volume = float(np.prod(self.box_lengths))
        comps = self.component_info
        return {
            "composition": dict(sorted(Counter(a.species for a in self.atoms).items())),
            "mass_g_per_mol": mass,
            "volume_angstrom3": volume,
            "current_density": mass / AVOGADRO * 1e24 / volume,
            "conversion": self.conversion,
            "new_bonds": sum(e.get("committed", False) for e in self.events),
            "component_sizes": sorted((len(c[0]) for c in comps), reverse=True),
            "periodic_connectivity_ranks": [c[2] for c in comps],
        }

    @property
    def crystal(self) -> Crystal:
        index = self.index
        graph = nx.Graph()
        graph.add_nodes_from((i, {"species": a.species}) for i, a in enumerate(self.atoms))
        records = []
        for bond in self.bonds:
            i, j = index[bond.i], index[bond.j]
            image = bond.image
            if i > j:
                i, j = j, i
                image = tuple(-x for x in image)
            distance = float(
                np.linalg.norm(
                    self.coordinates[j] + np.array(image) * self.box_lengths - self.coordinates[i]
                )
            )
            records.append(PeriodicBond(i, j, image, distance, bond.order))
            if i != j and (not graph.has_edge(i, j) or distance < graph[i][j]["distance"]):
                graph.add_edge(
                    i,
                    j,
                    image=image,
                    order=bond.order,
                    distance=distance,
                    cross_boundary=any(image),
                )
        component_ids = {}
        for record in self.component_records:
            component_ids.update(dict.fromkeys(record["atom_ids"], record["id"]))
        structure = Structure(
            Lattice.orthorhombic(*self.box_lengths),
            [a.species for a in self.atoms],
            self.coordinates,
            coords_are_cartesian=True,
            site_properties={
                "atom_id": [a.id for a in self.atoms],
                "template_id": [a.template_id for a in self.atoms],
                "template_label": [a.label for a in self.atoms],
                "instance_id": [a.instance_id for a in self.atoms],
                "component_id": [component_ids[a.id] for a in self.atoms],
                "formal_charge": [a.formal_charge for a in self.atoms],
                "atom_mass": [a.mass for a in self.atoms],
            },
        )
        return Crystal("polymer_box", structure, bonds=graph, periodic_bonds=records)

    def copy(self):
        return copy.deepcopy(self)

    def wrap(self):
        result = self.copy()
        shifts = np.floor(result.coordinates / result.box_lengths).astype(int)
        result.coordinates -= shifts * result.box_lengths
        idx = result.index
        result.bonds = [
            replace(
                b,
                image=tuple(
                    int(x) for x in (np.array(b.image) + shifts[idx[b.j]] - shifts[idx[b.i]])
                ),
            )
            for b in result.bonds
        ]
        return result

    def to_dict(self):
        self.check()
        return {
            "schema_version": SCHEMA_VERSION,
            "atoms": [asdict(a) for a in self.atoms],
            "coordinates": self.coordinates.tolist(),
            "box_lengths": list(self.box_lengths),
            "periodicity": list(self.periodicity),
            "bonds": [asdict(b) for b in self.bonds],
            "connectors": [asdict(s) for s in self.connectors],
            "initial_connector_count": self.initial_connector_count,
            "metadata": copy.deepcopy(self.metadata),
            "events": copy.deepcopy(self.events),
            "provenance": copy.deepcopy(self.provenance),
            "statistics": self.statistics,
            "components": self.component_records,
        }

    @property
    def state_hash(self):
        return content_hash({"hash_version": 1, "state": self.to_dict()})

    @classmethod
    def from_dict(cls, data):
        data = copy.deepcopy(data)
        version = data.pop("schema_version", None)
        if type(version) is not int or version != SCHEMA_VERSION:
            raise ValueError("Unsupported native state schema_version.")
        expected_statistics = data.pop("statistics", None)
        expected_components = data.pop("components", None)
        data["atoms"] = [Atom(**a) for a in data["atoms"]]
        data["bonds"] = [Bond(**b) for b in data["bonds"]]
        data["connectors"] = [Site(**s) for s in data["connectors"]]
        state = cls(**data)
        if expected_statistics != state.statistics:
            raise ValueError("State statistics do not match its atoms and topology.")
        if expected_components != state.component_records:
            raise ValueError("State component identities do not match its topology.")
        return state

    def save(self, output_dir: str | Path) -> Path:
        """Publish a new bundle, with its verification manifest written last."""
        self.check()
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        bundle = Path(tempfile.mkdtemp(prefix="box-", dir=output_dir))
        native = bundle / "state.json"
        _atomic_json(native, self.to_dict())
        crystal = self.crystal
        crystal.write_cif(bundle / "box.cif")
        crystal.write_xyz(bundle / "box.xyz", comment="Interchange only; topology in state.json")
        files = {name: file_hash(bundle / name) for name in ("state.json", "box.cif", "box.xyz")}
        # Read the native state before publishing completion; catches serialization loss.
        roundtrip = self.from_dict(_read_json(native))
        if roundtrip.state_hash != self.state_hash:
            raise ValueError("Native state roundtrip changed its identity.")
        _atomic_json(
            bundle / "manifest.json",
            {
                "schema_version": SCHEMA_VERSION,
                "hash_version": 1,
                "bundle_complete": True,
                "state_hash": self.state_hash,
                "files": files,
            },
        )
        return native

    @classmethod
    def load(cls, path: str | Path):
        path = Path(path)
        if not path.is_dir() and path.name != "state.json":
            raise ValueError("Load a native bundle directory or its state.json file.")
        bundle = path if path.is_dir() else path.parent
        manifest = _read_json(bundle / "manifest.json")
        if (
            type(manifest.get("schema_version")) is not int
            or manifest.get("schema_version") != SCHEMA_VERSION
            or type(manifest.get("hash_version")) is not int
            or manifest.get("hash_version") != 1
            or manifest.get("bundle_complete") is not True
        ):
            raise ValueError("Unsupported or incomplete native box bundle.")
        if set(manifest.get("files", {})) != {"state.json", "box.cif", "box.xyz"}:
            raise ValueError("Incomplete native box bundle file list.")
        for name, digest in manifest["files"].items():
            if file_hash(bundle / name) != digest:
                raise ValueError(f"Box bundle file hash mismatch: {name}")
        data = _read_json(bundle / "state.json")
        if data.get("schema_version") != SCHEMA_VERSION:
            raise ValueError("Unsupported native state schema_version.")
        state = cls.from_dict(data)
        if state.state_hash != manifest["state_hash"]:
            raise ValueError("Native state hash mismatch.")
        return state

    def update_geometry(
        self, *, parent_hash: str, atoms: list[dict], box_lengths=None, unwrapped: bool = False
    ):
        """Import coordinates in the parent's image gauge, with explicit atom IDs.

        Each record has id, species, coordinates and (for wrapped data) an
        integer wrapping_offset. The latter reconstructs unwrapped coordinates
        relative to the parent's per-atom images, not a nearest-image guess.
        """
        if parent_hash != self.state_hash:
            raise ValueError("Geometry parent_hash does not identify this state.")
        if type(unwrapped) is not bool:
            raise ValueError("unwrapped must be boolean.")
        records = {a["id"]: a for a in atoms}
        if len(records) != len(atoms) or set(records) != set(self.index):
            raise ValueError("Geometry must map every atom ID exactly once.")
        result = self.copy()
        if box_lengths is not None:
            result.box_lengths = validate_lengths(box_lengths)
        positions = []
        for atom in self.atoms:
            rec = records[atom.id]
            if set(rec) - {"id", "species", "coordinates", "wrapping_offset"}:
                raise ValueError("Geometry updates cannot include chemical/topology edits.")
            if rec["species"] != atom.species:
                raise ValueError("Geometry update changed species.")
            pos = np.asarray(rec["coordinates"], dtype=float)
            if pos.shape != (3,) or not np.isfinite(pos).all():
                raise ValueError("Each atom needs three finite coordinates.")
            if unwrapped:
                if "wrapping_offset" in rec:
                    raise ValueError("Choose unwrapped positions or explicit wrapping offsets.")
            else:
                shift = rec.get("wrapping_offset")
                if shift is None or len(shift) != 3 or any(type(x) is not int for x in shift):
                    raise ValueError("Wrapped geometry requires explicit integer wrapping offsets.")
                pos = pos + np.array(shift) * result.box_lengths
            positions.append(pos)
        result.coordinates = np.array(positions)
        result.metadata.pop("validation", None)
        result.metadata.pop("derived_results", None)
        result.metadata["geometry_requires_validation"] = True
        result.provenance.append({"operation": "update_geometry", "parent_hash": parent_hash})
        result = result.wrap()
        result.check()
        return result


def _read_json(path):
    def invalid(value):
        raise ValueError(f"Nonfinite JSON value: {value}")

    return json.loads(Path(path).read_text(), parse_constant=invalid)


def _atomic_json(path, data):
    path = Path(path)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent, prefix=".pending-", delete=False
    ) as handle:
        json.dump(data, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(handle.name, path)
