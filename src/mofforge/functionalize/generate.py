"""Generate anchor-tagged query/replacement fragment pairs with RDKit.

Given a linker SMILES, a set of site indices (from
:func:`mofforge.functionalize.sites.find_functionalizable_sites`) and a
functional-group name (from :mod:`mofforge.functionalize.groups`), this module
produces two XYZ files that the existing find/replace pipeline consumes:

* the **query** — the linker's aromatic ring system with the chosen site
  hydrogen(s) tagged ``H!`` (the anchor tag).  Ring carbons that connect to the
  rest of the framework (carboxylates, other rings) are left "bare" (no H), so
  the query matches the ring as it sits inside the MOF.
* the **replacement** — the same ring system with the functional group grafted
  in place of the tagged hydrogen(s).

All 3-D geometry is produced by RDKit embedding + UFF optimisation of a
chemically valid capped molecule (e.g. benzene / nitrobenzene), after which the
cap hydrogens on the connection carbons are removed to expose the bare
attachment points.  The agent never authors coordinates or SMILES.

This procedure reproduces the hand-curated reference moieties: benzene with the
two para cap-hydrogens removed is the 10-atom ``2-!-p-phenylene`` query, and
nitrobenzene with the same two removed is the 12-atom ``2-nitro-p-phenylene``
replacement.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from mofforge.build.smiles_to_bb import _ensure_rdkit
from mofforge.functionalize.groups import FunctionalGroup, get_group
from mofforge.functionalize.sites import find_functionalizable_sites
from mofforge.io.xyz import write_xyz
from mofforge.utils.config import config

logger = logging.getLogger("mofforge")

# RDKit atom properties used to carry state across RWMol edits.
_PROP_SITE = "_mofforge_site"  # marks a chosen functionalization carbon
_PROP_CONNECTION = "_mofforge_connection"  # marks a bare-carbon (cap) position
_PROP_KEEP = "_mofforge_keep"  # marks an atom to retain (ring system or grafted group)


def _aromatic_ring_system(mol, seed_atom: int) -> set[int]:
    """Return the set of atoms fused to *seed_atom* through aromatic bonds.

    Biphenyl-type linkers (rings joined by a single non-aromatic bond) yield the
    single ring containing the seed; fused systems (naphthalene) yield all
    fused rings.
    """
    system = {seed_atom}
    stack = [seed_atom]
    while stack:
        idx = stack.pop()
        atom = mol.GetAtomWithIdx(idx)
        for bond in atom.GetBonds():
            if not bond.GetIsAromatic():
                continue
            other = bond.GetOtherAtomIdx(idx)
            if other not in system:
                system.add(other)
                stack.append(other)
    return system


def _build_core(mol, site_atoms: list[int], group: FunctionalGroup | None):
    """Build an embedded RDKit mol of the ring core.

    Parameters
    ----------
    mol:
        Parsed linker molecule (no explicit hydrogens).
    site_atoms:
        RDKit atom indices of the chosen functionalization carbons.  All must
        lie in a single aromatic ring system.
    group:
        The functional group to graft at each site (replacement), or ``None`` to
        keep the site hydrogens (query).

    Returns
    -------
    (mol_h, site_atom_map, connection_atoms)
        ``mol_h`` is the H-added, embedded, UFF-optimised molecule.  Indices are
        those of ``mol_h`` (post-edit).  ``site_atom_map`` maps each original
        site index to its atom index in ``mol_h``; ``connection_atoms`` is the
        set of bare-carbon indices in ``mol_h``.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    system = _aromatic_ring_system(mol, site_atoms[0])
    for s in site_atoms[1:]:
        if s not in system:
            raise ValueError(
                "All requested sites must lie in the same aromatic ring system. "
                f"Site atom {s} is not fused with site atom {site_atoms[0]}."
            )

    rw = Chem.RWMol(mol)

    # Flag connection carbons (ring atoms bonded to a heavy atom outside the
    # ring system) and the chosen sites, and mark the whole system to keep.
    site_set = set(site_atoms)
    for idx in range(rw.GetNumAtoms()):
        atom = rw.GetAtomWithIdx(idx)
        if idx in system:
            atom.SetProp(_PROP_KEEP, "1")
            atom.SetProp(_PROP_SITE, "1" if idx in site_set else "0")
            is_connection = any(
                nbr.GetIdx() not in system and nbr.GetAtomicNum() > 1
                for nbr in atom.GetNeighbors()
            )
            atom.SetProp(_PROP_CONNECTION, "1" if is_connection else "0")
        else:
            atom.SetProp(_PROP_KEEP, "0")

    # Graft the functional group at each site (replacement only).
    if group is not None:
        grp = Chem.MolFromSmiles(group.smiles)
        if grp is None:
            raise ValueError(f"RDKit could not parse group SMILES: {group.smiles!r}")
        for site_idx in site_atoms:
            offset = rw.GetNumAtoms()
            rw.InsertMol(grp)
            # Newly inserted group atoms are kept and are not sites/connections.
            for j in range(offset, rw.GetNumAtoms()):
                ga = rw.GetAtomWithIdx(j)
                ga.SetProp(_PROP_KEEP, "1")
                ga.SetProp(_PROP_SITE, "0")
            attach_in_combined = offset + group.attach_idx
            rw.AddBond(site_idx, attach_in_combined, Chem.BondType.SINGLE)
            # The aromatic site carbon now carries a substituent instead of H.
            rw.GetAtomWithIdx(site_idx).SetNoImplicit(False)

    # Remove everything outside the ring system (and outside grafted groups).
    for idx in sorted(range(rw.GetNumAtoms()), reverse=True):
        if rw.GetAtomWithIdx(idx).GetProp(_PROP_KEEP) != "1":
            rw.RemoveAtom(idx)

    core = rw.GetMol()
    Chem.SanitizeMol(core)

    mol_h = Chem.AddHs(core)

    params = AllChem.ETKDGv3()
    if AllChem.EmbedMolecule(mol_h, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol_h, params) != 0:
            raise ValueError("RDKit could not generate 3-D coordinates for the core.")
    try:
        AllChem.UFFOptimizeMolecule(mol_h, maxIters=2000)
    except Exception:
        logger.warning("UFF optimisation failed for core; using embedded coords")

    site_atom_map: dict[int, int] = {}
    connection_atoms: set[int] = set()
    # Track original site atoms by order: sites kept their identity through edits
    # because we only appended/removed atoms; recover via the SITE property.
    site_by_flag = [
        a.GetIdx()
        for a in mol_h.GetAtoms()
        if a.HasProp(_PROP_SITE) and a.GetProp(_PROP_SITE) == "1"
    ]
    for orig, new in zip(sorted(site_atoms), sorted(site_by_flag), strict=False):
        site_atom_map[orig] = new
    for a in mol_h.GetAtoms():
        if a.HasProp(_PROP_CONNECTION) and a.GetProp(_PROP_CONNECTION) == "1":
            connection_atoms.add(a.GetIdx())

    return mol_h, site_atom_map, connection_atoms


def _one_hydrogen_of(mol_h, carbon_idx: int) -> int:
    """Return the index of one hydrogen bonded to *carbon_idx*."""
    atom = mol_h.GetAtomWithIdx(carbon_idx)
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() == 1:
            return nbr.GetIdx()
    raise ValueError(f"Carbon atom {carbon_idx} has no hydrogen to substitute/tag.")


def _fragment_to_xyz(
    mol_h,
    site_atoms_new: list[int],
    connection_atoms: set[int],
    *,
    tag_sites: bool,
    r_tag: str,
) -> tuple[list[str], np.ndarray]:
    """Convert an embedded core into (species, coords) for an XYZ fragment.

    Removes the cap hydrogens on connection carbons.  When ``tag_sites`` is
    True (query), the hydrogen on each site carbon is kept and its label tagged
    with ``r_tag``; otherwise (replacement) the site hydrogen was already
    replaced by the grafted group.
    """
    conf = mol_h.GetConformer()

    # Hydrogens to drop entirely: the cap H on every connection carbon.
    drop: set[int] = set()
    for c in connection_atoms:
        drop.add(_one_hydrogen_of(mol_h, c))

    # Hydrogens to tag: the site H (query only).
    tagged: set[int] = set()
    if tag_sites:
        for c in site_atoms_new:
            tagged.add(_one_hydrogen_of(mol_h, c))

    species: list[str] = []
    coords: list[list[float]] = []
    for atom in mol_h.GetAtoms():
        idx = atom.GetIdx()
        if idx in drop:
            continue
        pos = conf.GetAtomPosition(idx)
        label = atom.GetSymbol()
        if idx in tagged:
            label = f"{label}{r_tag}"
        species.append(label)
        coords.append([pos.x, pos.y, pos.z])

    return species, np.array(coords, dtype=np.float64)


def make_query_replacement(
    linker_smiles: str,
    sites: int | list[int],
    group: str,
    output_dir: str | Path | None = None,
    r_tag: str | None = None,
) -> tuple[Path, Path]:
    """Generate an anchor-tagged query/replacement XYZ pair for the pipeline.

    Parameters
    ----------
    linker_smiles:
        SMILES of the linker to functionalize.
    sites:
        One site index, or a list of site indices (all in the same aromatic ring
        system), as returned by
        :func:`~mofforge.functionalize.sites.find_functionalizable_sites`.
        Multiple indices produce a multi-substitution pattern.
    group:
        Functional-group name from
        :func:`~mofforge.functionalize.groups.available_groups`.
    output_dir:
        Directory to write the two XYZ files into (created if needed).  Defaults
        to the current directory.
    r_tag:
        Anchor tag to use; defaults to :data:`config.r_tag` (``"!"``).

    Returns
    -------
    (query_path, replacement_path)
        Paths to the generated query and replacement XYZ files.
    """
    _ensure_rdkit()
    from rdkit import Chem

    if r_tag is None:
        r_tag = config.r_tag
    if isinstance(sites, int):
        sites = [sites]
    if not sites:
        raise ValueError("At least one site index is required.")

    fg = get_group(group)

    mol = Chem.MolFromSmiles(linker_smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {linker_smiles!r}")

    # Map public site indices -> RDKit atom indices via the site finder (which
    # uses the same SMILES parse and ordering).
    all_sites = find_functionalizable_sites(linker_smiles)
    by_index = {s.index: s for s in all_sites}
    site_atoms: list[int] = []
    for i in sites:
        if i not in by_index:
            raise ValueError(
                f"Site index {i} not found. Valid indices: "
                f"{[s.index for s in all_sites]}"
            )
        site_atoms.append(by_index[i].atom_idx)

    out_dir = Path(output_dir) if output_dir is not None else Path.cwd()
    out_dir.mkdir(parents=True, exist_ok=True)

    site_tag = "-".join(str(i) for i in sites)

    # --- Query ---------------------------------------------------------------
    q_mol, q_site_map, q_conn = _build_core(mol, site_atoms, group=None)
    q_species, q_coords = _fragment_to_xyz(
        q_mol, list(q_site_map.values()), q_conn, tag_sites=True, r_tag=r_tag
    )
    query_path = out_dir / f"query_site{site_tag}.xyz"
    write_xyz(q_species, q_coords, query_path, comment=f"query sites={sites}")

    # --- Replacement ---------------------------------------------------------
    r_mol, r_site_map, r_conn = _build_core(mol, site_atoms, group=fg)
    r_species, r_coords = _fragment_to_xyz(
        r_mol, list(r_site_map.values()), r_conn, tag_sites=False, r_tag=r_tag
    )
    replacement_path = out_dir / f"replacement_{group}_site{site_tag}.xyz"
    write_xyz(
        r_species, r_coords, replacement_path, comment=f"{group} at sites={sites}"
    )

    logger.debug(
        "Generated fragments: query=%s (%d atoms), replacement=%s (%d atoms)",
        query_path,
        len(q_species),
        replacement_path,
        len(r_species),
    )
    return query_path, replacement_path
