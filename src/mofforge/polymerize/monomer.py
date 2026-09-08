"""Monomer preparation for simulated polymerization.

Turns a monomer SMILES (or XYZ/CIF path) into a relaxed 3-D structure with
tagged reactive sites, reusing the RDKit embed + UFF recipe from
:mod:`mofforge.build.smiles_to_bb`.  The output an MDL ``.mol`` file that pysimm
reads, plus the metadata the engine needs (reactive sites, molar mass) to size
the packing box.

This module deliberately does **not** import pysimm: monomer prep is pure
RDKit + numpy so it is fully unit-testable without the external engine.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

from mofforge.build.smiles_to_bb import _ensure_rdkit
from mofforge.io.xyz import read_xyz
from mofforge.polymerize.base import Monomer, ReactiveSite
from mofforge.polymerize.reactions import detect_reactive_sites
from mofforge.provenance import effective_seed

logger = logging.getLogger("mofforge")


@dataclass
class PreparedMonomer:
    """A monomer with 3-D geometry and resolved reactive sites.

    Attributes
    ----------
    name:
        Monomer name.
    mol_path:
        Path to the written MDL ``.mol`` file (what pysimm reads).
    sites:
        Resolved reactive sites (0-based heavy-atom indices into the mol).
    molar_mass:
        Molar mass in g/mol (for density-based box sizing).
    n_atoms:
        Atom count including explicit hydrogens.
    """

    name: str
    mol_path: Path
    sites: list[ReactiveSite]
    molar_mass: float
    n_atoms: int = 0
    species: list[str] = field(default_factory=list)


def prepare_monomer(
    monomer: Monomer,
    output_dir: str | Path,
    random_seed: int | None = None,
) -> PreparedMonomer:
    """Generate 3-D geometry + reactive sites for one monomer.

    SMILES monomers are embedded and UFF-optimized with RDKit; file-based
    monomers (XYZ/CIF) are read as-is.  Reactive sites are taken from
    ``monomer.sites`` when supplied, otherwise detected from the SMILES.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if monomer.is_smiles:
        return _prepare_from_smiles(monomer, output_dir, random_seed)
    return _prepare_from_file(monomer, output_dir)


def _prepare_from_smiles(
    monomer: Monomer,
    output_dir: Path,
    random_seed: int | None,
) -> PreparedMonomer:
    """Embed a SMILES monomer and resolve its reactive sites."""
    _ensure_rdkit()
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors

    smiles = str(monomer.source)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")

    # Resolve reactive sites on the implicit-H molecule (indices are stable
    # before AddHs appends hydrogens), mirroring smiles_to_bb.py.
    sites = _resolve_sites(monomer, smiles)

    molar_mass = float(Descriptors.MolWt(mol))

    mol_h = Chem.AddHs(mol)
    seed = effective_seed(random_seed)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol_h, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol_h, params) != 0:
            raise ValueError(f"RDKit could not generate 3-D coordinates for: {smiles!r}")
    try:
        AllChem.UFFOptimizeMolecule(mol_h, maxIters=2000)
    except Exception as exc:  # pragma: no cover - RDKit rarely raises here
        raise ValueError(f"UFF geometry optimization failed for {smiles!r}.") from exc

    mol_path = output_dir / f"{monomer.name}.mol"
    Chem.MolToMolFile(mol_h, str(mol_path))

    species = [atom.GetSymbol() for atom in mol_h.GetAtoms()]
    return PreparedMonomer(
        name=monomer.name,
        mol_path=mol_path,
        sites=sites,
        molar_mass=molar_mass,
        n_atoms=mol_h.GetNumAtoms(),
        species=species,
    )


def _prepare_from_file(monomer: Monomer, output_dir: Path) -> PreparedMonomer:
    """Prepare a monomer from an existing XYZ file.

    Reactive sites must be supplied explicitly (via ``monomer.sites``) or tagged
    with the ``!`` anchor convention in the XYZ; there is no SMILES to detect
    them from.
    """
    src = Path(str(monomer.source))
    if src.suffix.lower() != ".xyz":
        raise ValueError(
            f"File-based monomers must be .xyz (got {src.suffix!r}); "
            "convert to XYZ or provide a SMILES string."
        )
    species, coords = read_xyz(src)

    sites = list(monomer.sites)
    if not sites:
        # Fall back to '!' anchor tags (same convention as core.moiety).
        sites = [
            ReactiveSite(atom_idx=i, site_type="generic")
            for i, s in enumerate(species)
            if s.endswith("!")
        ]
    if not sites:
        raise ValueError(
            f"No reactive sites for monomer {monomer.name!r}: pass sites= or tag "
            "attachment atoms with '!' in the XYZ."
        )

    molar_mass = _molar_mass_from_species(species)

    # Re-emit as a clean XYZ mofforge/pysimm can read (strip '!' tags).
    from mofforge.utils.config import clean_species

    clean = [clean_species(s) for s in species]
    mol_path = output_dir / f"{monomer.name}.xyz"
    from mofforge.io.xyz import write_xyz

    write_xyz(clean, coords, mol_path, comment=monomer.name)

    return PreparedMonomer(
        name=monomer.name,
        mol_path=mol_path,
        sites=sites,
        molar_mass=molar_mass,
        n_atoms=len(species),
        species=clean,
    )


def _resolve_sites(monomer: Monomer, smiles: str) -> list[ReactiveSite]:
    """Return reactive sites for a monomer: explicit if given, else detected."""
    if monomer.sites:
        return list(monomer.sites)
    detected = detect_reactive_sites(smiles)
    sites = [
        ReactiveSite(atom_idx=d.atom_idx, anchor_idx=d.anchor_idx, site_type=d.site_type)
        for d in detected
    ]
    if not sites:
        raise ValueError(
            f"No reactive sites detected on monomer {monomer.name!r} ({smiles!r}). "
            "Supply sites= explicitly or use a monomer with a known reactive group."
        )
    return sites


def _molar_mass_from_species(species: list[str]) -> float:
    """Sum atomic masses for a list of (possibly tagged) element labels."""
    from pymatgen.core import Element

    from mofforge.utils.config import clean_species

    total = 0.0
    for s in species:
        total += float(Element(clean_species(s)).atomic_mass)
    return total


def box_length_for_density(
    prepared: list[PreparedMonomer],
    counts: list[int],
    target_density: float,
) -> float:
    """Return the cubic box edge (Angstrom) for a target mass density.

    ``density = total_mass / volume``; converts g/mol and g/cm^3 to a box edge
    in Angstrom using Avogadro's number.
    """
    if target_density <= 0:
        raise ValueError("target_density must be positive.")
    avogadro = 6.022_140_76e23
    total_mass_g = sum(
        p.molar_mass * n for p, n in zip(prepared, counts, strict=True)
    ) / avogadro
    volume_cm3 = total_mass_g / target_density
    volume_ang3 = volume_cm3 * 1.0e24  # 1 cm^3 = 1e24 Angstrom^3
    return float(volume_ang3 ** (1.0 / 3.0))
