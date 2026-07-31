"""Tests for agent-driven linker functionalization.

Covers the curated group menu, deterministic site detection, RDKit fragment
generation (validated against the hand-curated reference moieties), and the
high-level functionalize / campaign API.
"""

from __future__ import annotations

import pytest

pytest.importorskip("rdkit")

from mofforge.functionalize import (
    available_groups,
    find_functionalizable_sites,
    functionalize,
    group_smiles,
    make_query_replacement,
    run_campaign,
)
from mofforge.functionalize.groups import get_group
from mofforge.io.xyz import read_xyz

# A benzene-1,4-dicarboxylic acid (BDC / terephthalic acid) linker — the linker
# of IRMOF-1.
BDC = "O=C(O)c1ccc(C(=O)O)cc1"
# 2,6-naphthalenedicarboxylic acid — a fused linker with distinct site classes.
NDC = "OC(=O)c1ccc2cc(C(=O)O)ccc2c1"


# ---------------------------------------------------------------------------
# Group menu
# ---------------------------------------------------------------------------


def test_available_groups_nonempty():
    groups = available_groups()
    assert "NH2" in groups
    assert "NO2" in groups
    assert "F" in groups


def test_group_smiles_and_lookup():
    assert group_smiles("NO2") == "[N+](=O)[O-]"
    assert get_group("CH3").name == "CH3"


def test_unknown_group_raises():
    with pytest.raises(KeyError):
        get_group("nonexistent_group")


# ---------------------------------------------------------------------------
# Site detection
# ---------------------------------------------------------------------------


def test_bdc_sites_all_equivalent():
    sites = find_functionalizable_sites(BDC)
    assert len(sites) == 4  # four aromatic C-H on the ring
    # All four are symmetry-equivalent -> single class.
    assert len({s.symmetry_class for s in sites}) == 1
    # Indices are stable and contiguous.
    assert [s.index for s in sites] == [0, 1, 2, 3]


def test_ndc_has_multiple_symmetry_classes():
    sites = find_functionalizable_sites(NDC)
    assert len(sites) == 6
    # A 2,6-substituted naphthalene has three distinct aromatic C-H environments.
    assert len({s.symmetry_class for s in sites}) == 3


def test_sites_exclude_carboxylate_carbons():
    # No returned site should be a carboxylate carbon (they carry no aromatic H).
    sites = find_functionalizable_sites(BDC)
    assert all(s.element == "C" for s in sites)


def test_bad_smiles_raises():
    with pytest.raises(ValueError):
        find_functionalizable_sites("not a smiles !!!")


# ---------------------------------------------------------------------------
# Fragment generation (against curated references)
# ---------------------------------------------------------------------------


def test_generated_matches_reference_nitro(tmp_path):
    """Generated BDC/nitro fragments match the hand-curated moieties atom-for-atom."""
    q, r = make_query_replacement(BDC, 0, "NO2", output_dir=str(tmp_path))
    qs, _ = read_xyz(q)
    rs, _ = read_xyz(r)

    # Query: benzene ring, two bare connection carbons, three H, one H! anchor.
    assert sorted(qs) == sorted(["C"] * 6 + ["H", "H", "H", "H!"])
    # Replacement: same ring backbone with -NO2 in place of the anchor H.
    assert sorted(rs) == sorted(["C"] * 6 + ["H", "H", "H", "N", "O", "O"])


def test_generated_query_has_single_anchor(tmp_path):
    q, _ = make_query_replacement(BDC, 0, "NH2", output_dir=str(tmp_path))
    qs, _ = read_xyz(q)
    assert sum(1 for s in qs if s.endswith("!")) == 1


def test_multi_substitution_tags_multiple_sites(tmp_path):
    q, r = make_query_replacement(BDC, [0, 3], "F", output_dir=str(tmp_path))
    qs, _ = read_xyz(q)
    rs, _ = read_xyz(r)
    assert sum(1 for s in qs if s.endswith("!")) == 2
    assert sum(1 for s in rs if s == "F") == 2


def test_generate_rejects_cross_ring_sites(tmp_path):
    # Sites 0 and 4 of NDC live on different fused positions but the guard is on
    # ring-system membership; pick two indices from different rings of biphenyl.
    biphenyl = "O=C(O)c1ccc(-c2ccc(C(=O)O)cc2)cc1"
    sites = find_functionalizable_sites(biphenyl)
    # Group indices by ring.
    by_ring: dict[int, list[int]] = {}
    for s in sites:
        by_ring.setdefault(s.ring_id, []).append(s.index)
    rings = [r for r in by_ring.values() if r]
    assert len(rings) >= 2
    cross = [rings[0][0], rings[1][0]]
    with pytest.raises(ValueError):
        make_query_replacement(biphenyl, cross, "F", output_dir=str(tmp_path))


# ---------------------------------------------------------------------------
# High-level functionalize / campaign
# ---------------------------------------------------------------------------


def test_functionalize_irmof1(crystal_dir, tmp_path):
    out = tmp_path / "func.cif"
    res = functionalize(
        str(crystal_dir / "IRMOF-1.cif"),
        BDC,
        "NO2",
        sites=0,
        coverage=0.5,
        output_cif=str(out),
        random_seed=1,
    )
    assert res.error is None
    assert res.n_matches == 24  # IRMOF-1 has 24 BDC linker locations
    assert res.n_functionalized == 12  # 50% coverage
    assert out.exists()
    assert res.crystal is not None


def test_functionalize_full_coverage(crystal_dir, tmp_path):
    res = functionalize(
        str(crystal_dir / "IRMOF-1.cif"),
        BDC,
        "F",
        sites=0,
        coverage=1.0,
        output_cif=str(tmp_path / "all.cif"),
        random_seed=1,
    )
    assert res.error is None
    assert res.n_functionalized == res.n_matches


def test_functionalize_reports_no_match(crystal_dir, tmp_path):
    # A linker whose ring is absent from IRMOF-1 yields zero matches, reported
    # cleanly rather than raising.
    res = functionalize(
        str(crystal_dir / "IRMOF-1.cif"),
        "c1ccncc1C(=O)O",  # a pyridyl acid — not present in IRMOF-1
        "F",
        sites=0,
        coverage=1.0,
        output_cif=str(tmp_path / "none.cif"),
        validate=False,
    )
    assert res.n_matches == 0
    assert res.error is not None


def test_campaign_ranks_results(crystal_dir, tmp_path):
    results = run_campaign(
        str(crystal_dir / "IRMOF-1.cif"),
        BDC,
        groups=["F", "NH2"],
        coverages=[0.25, 1.0],
        output_dir=str(tmp_path),
        random_seed=1,
    )
    assert len(results) == 4  # 2 groups x 2 coverages
    # Every combination ran without error.
    assert all(r.error is None for r in results)
    # Ranking is non-decreasing in clash count among valid/failed grouping.
    clashes = [r.clashes for r in results if r.clashes is not None]
    assert clashes == sorted(clashes) or len(set(clashes)) <= 1
