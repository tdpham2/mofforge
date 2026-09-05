"""Scientific regressions for the 0.2 reliability release."""

import json
import random

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from mofforge import Crystal, infer_bonds, validate_structure
from mofforge.adsorbate.placement import _filter_by_intermolecular_dist, place_adsorbate
from mofforge.adsorbate.sites import AdsorptionSite, _cluster_void_points
from mofforge.provenance import Provenance, structure_hash
from mofforge.utils.periodic import min_image_distance, nearest_image


def test_skewed_minimum_image():
    lattice = Lattice.from_parameters(10, 10, 10, 90, 90, 30)
    delta = np.array([0.49, 0.49, 0])
    expected = np.linalg.norm(lattice.get_cartesian_coords(delta - [1, 0, 0]))
    assert min_image_distance(np.zeros(3), delta, lattice) == pytest.approx(expected)
    assert np.linalg.norm(
        lattice.get_cartesian_coords(nearest_image(delta, lattice))
    ) == pytest.approx(expected)


def test_primitive_and_supercell_coordination():
    primitive = Structure(Lattice.cubic(1.5), ["C"], [[0, 0, 0]])
    small = infer_bonds(Crystal.from_structure(primitive))
    expanded = infer_bonds(Crystal.from_structure(primitive * [2, 2, 2]))
    assert small.coordination_number(0) == 6
    assert [expanded.coordination_number(i) for i in range(8)] == [6] * 8
    assert small.bonds.number_of_edges() == 0
    assert len(small.periodic_bonds) == 3


def test_reassembly_uses_skewed_lattice():
    import networkx as nx

    from mofforge.replace.conglomerate import reassemble

    lattice = Lattice.from_parameters(10, 10, 10, 90, 90, 30)
    structure = Structure(lattice, ["C", "C"], [[0, 0, 0], [0.49, 0.49, 0]])
    distance = lattice.get_distance_and_image(*structure.frac_coords)[0]
    bonds = nx.Graph()
    bonds.add_nodes_from([(0, {"species": "C"}), (1, {"species": "C"})])
    bonds.add_edge(0, 1, distance=distance, cross_boundary=True)
    result = reassemble(Crystal("skew", structure, bonds=bonds))
    assert np.linalg.norm(result.cart_coords[1] - result.cart_coords[0]) == pytest.approx(distance)


def test_addition_keeps_right_hand_periodic_images():
    lattice = Lattice.cubic(1.5)
    right = infer_bonds(Crystal.from_structure(Structure(lattice, ["C"], [[0, 0, 0]])))
    left = Crystal.from_structure(Structure(lattice, ["He"], [[0.5, 0.5, 0.5]]))
    result = left + right
    assert result.coordination_number(1) == 6


def test_transformations_preserve_compositions_and_properties():
    source = Structure(
        Lattice.cubic(10),
        [{"Fe2+": 0.5, "Mn2+": 0.5}, "O2-"],
        [[1.1, 0, 0], [0.3, 0, 0]],
        site_properties={"magmom": [4, 0]},
    )
    crystal = Crystal.from_structure(source)
    crystal.provenance = Provenance(operation="input")
    for transformed in (crystal.copy(), crystal.wrap(), crystal[[0, 1]]):
        assert [s.species for s in transformed.structure] == [s.species for s in source]
        assert transformed.structure.site_properties == {"magmom": [4, 0]}
        assert transformed.provenance.operation == "input"
    empty = crystal[[]]
    assert empty.lattice == crystal.lattice
    combined = crystal[[0]] + crystal[[1]]
    assert [s.species for s in combined.structure] == [s.species for s in source]
    assert combined.structure.site_properties == {"magmom": [4, 0]}


def test_bond_images_survive_wrap_reindex_and_displacement():
    crystal = infer_bonds(
        Crystal.from_structure(
            Structure(
                Lattice.cubic(10),
                ["C", "C"],
                [[0.95, 0, 0], [1.1, 0, 0]],
            )
        )
    )
    wrapped = crystal.wrap()[[1, 0]]
    assert wrapped.bonds[0][1]["distance"] == pytest.approx(1.5)
    assert wrapped.coordination_number(0) == 1
    assert np.linalg.norm(wrapped.bond_vectors(0)[0]) == pytest.approx(1.5)
    coords = wrapped.frac_coords
    coords[0, 0] += 0.01
    wrapped.set_frac_coords(coords)
    assert wrapped.bonds[0][1]["distance"] == pytest.approx(1.6)


@pytest.mark.parametrize("bonded", [False, True])
def test_severe_overlap_is_fatal_even_when_bonded(bonded):
    crystal = Crystal.from_xyz(["C", "C"], np.array([[0, 0, 0], [0.5, 0, 0]]))
    if bonded:
        crystal = infer_bonds(crystal)
    report = validate_structure(crystal, check_clashes=False)
    assert not report.is_valid
    assert any("Severe overlap" in error for error in report.errors)


def test_exact_duplicate_atoms_are_invalid():
    crystal = Crystal.from_xyz(["C", "C"], np.zeros((2, 3)))
    assert not validate_structure(crystal).is_valid


def test_empty_and_unchecked_reports_do_not_pass():
    from mofforge.validation import ValidationReport

    assert not validate_structure(Crystal.empty()).is_valid
    assert not ValidationReport().is_valid


def test_unknown_charge_is_not_neutral():
    crystal = Crystal.from_xyz(["C"], np.zeros((1, 3)))
    report = validate_structure(crystal, check_charges=True)
    assert report.is_valid
    assert report.charge_balance is None
    assert "charges" in report.checks_skipped
    charged = Crystal.from_structure(Structure(Lattice.cubic(10), ["Na+"], [[0, 0, 0]]))
    assert validate_structure(charged, check_charges=True).charge_balance == 1


def test_nonfinite_coordinates_rejected():
    crystal = Crystal.from_xyz(["C"], np.zeros((1, 3)))
    with pytest.raises(ValueError, match="finite"):
        crystal.set_frac_coords(np.array([[float("nan"), 0, 0]]))


def test_periodic_site_spacing_and_clustering():
    lattice = Lattice.cubic(10)
    frac = np.array([[0.01, 0.5, 0.5], [0.99, 0.5, 0.5]])
    sites = [AdsorptionSite(f, lattice.get_cartesian_coords(f), "void", 5) for f in frac]
    assert len(_filter_by_intermolecular_dist(sites, 3, lattice)) == 1
    assert (
        len(
            _cluster_void_points(
                frac, lattice.get_cartesian_coords(frac), np.array([5, 4]), lattice, 1
            )
        )
        == 1
    )


def test_co2_internal_contacts_are_not_clashes():
    host = Crystal.from_structure(Structure(Lattice.cubic(30), ["C"], [[0, 0, 0]]))
    site = AdsorptionSite(np.array([0.5] * 3), np.array([15.0] * 3), "void", 20)
    result = place_adsorbate(host, "CO2", site=site, orient="fixed")
    assert result.clashes == 0
    assert result.validation.is_valid
    assert result.crystal.n_bonds == 2


def test_omitted_seed_is_recorded_and_can_be_replayed():
    host = Crystal.from_structure(Structure(Lattice.cubic(30), ["C"], [[0, 0, 0]]))
    site = AdsorptionSite(np.array([0.5] * 3), np.array([15.0] * 3), "void", 20)
    first = place_adsorbate(host, "CO2", site=site)
    seed = first.crystal.provenance.parameters["random_seed"]
    assert isinstance(seed, int)
    replay = place_adsorbate(host, "CO2", site=site, random_seed=seed)
    assert structure_hash(first.crystal) == structure_hash(replay.crystal)


def test_provenance_growth_is_linear_and_reads_legacy():
    current = Provenance(operation="start")
    sizes = []
    for i in range(40):
        current = current.chain(Provenance(operation=str(i)))
        sizes.append(len(json.dumps(current.to_dict())))
    assert sizes[-1] < sizes[19] * 2.1
    assert all("history" not in entry for entry in current.history)
    legacy = {
        "operation": "third",
        "history": [
            {"operation": "first", "history": []},
            {"operation": "second", "history": [{"operation": "first"}]},
        ],
    }
    loaded = Provenance.from_dict(legacy)
    assert [entry["operation"] for entry in loaded.history] == ["first", "second"]
    assert all("history" not in entry for entry in loaded.to_dict()["history"])


def test_seeded_replacement_preserves_name_history_and_global_rng(
    crystal_dir, moiety_dir, tmp_path
):
    from mofforge import find_pattern, fragment, replace_pattern

    parent = infer_bonds(Crystal.from_cif(crystal_dir / "IRMOF-1.cif"))
    query = fragment("2-!-p-phenylene.xyz", fragment_path=moiety_dir)
    replacement = fragment("2-nitro-p-phenylene.xyz", fragment_path=moiety_dir)
    match = find_pattern(query, parent)
    state = random.getstate()
    first = replace_pattern(match, replacement, nb_loc=2, random_seed=27, name="wanted")
    second = replace_pattern(match, replacement, nb_loc=2, random_seed=27, name="wanted")
    assert random.getstate() == state
    assert first.name == "wanted"
    assert structure_hash(first) == structure_hash(second)
    assert first.provenance.history[0]["operation"] == "load"
    output = tmp_path / "result.cif"
    first.write_cif(output)
    manifest = json.loads(output.with_suffix(".cif.json").read_text())
    assert manifest["structure_sha256"] == structure_hash(first)
    assert manifest["provenance"]["parameters"]["random_seed"] == 27


@pytest.mark.parametrize("name, atoms", [("IRMOF-1", 424), ("UiO-66", 912), ("MOF-74", 162)])
def test_reference_framework_geometry(name, atoms, crystal_dir):
    crystal = infer_bonds(Crystal.from_cif(crystal_dir / f"{name}.cif"))
    report = validate_structure(crystal)
    assert crystal.n_atoms == atoms
    assert report.is_valid
    assert not report.errors
    assert not report.steric_clashes
    assert not report.unusual_bonds
    # Coordination and vdW proximity are advisory chemistry heuristics.
    assert report.close_contacts
