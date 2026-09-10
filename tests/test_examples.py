"""Executable cookbook contracts, using actual scripts and temporary artifacts."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import re
import struct
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import pytest

from mofforge import Crystal

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
CATALOG = json.loads((EXAMPLES / "catalog.json").read_text())
INTEGRATION_FLAGS = {"build": "MOFFORGE_RUN_BUILD_TESTS", "render": "MOFFORGE_RUN_RENDER_TESTS"}
CASES = [
    pytest.param(
        entry, id=entry["script"], marks=pytest.mark.integration if entry.get("integration") else ()
    )
    for entry in CATALOG
]
CORE_RUNNER = """
import importlib.abc
import runpy
import sys

class NoOptionalPackages(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {
            'rdkit', 'pormake', 'tobacco3', 'mcp', 'playwright',
            'chemgraph', 'langchain_mcp_adapters',
        }:
            raise ModuleNotFoundError(f'Optional package disabled: {fullname}', name=fullname)

def no_network(event, args):
    if event in {'socket.connect', 'socket.getaddrinfo'}:
        raise RuntimeError('Core lessons must execute without network access.')

sys.meta_path.insert(0, NoOptionalPackages())
sys.addaudithook(no_network)
script = sys.argv.pop(1)
sys.argv[0] = script
runpy.run_path(script, run_name='__main__')
"""


def run_example(script, out, *arguments, cwd=None, no_site=False, core_only=False):
    command = [sys.executable]
    if no_site:
        command.append("-S")
    if core_only:
        command.extend(["-c", CORE_RUNNER])
    command.extend([str(EXAMPLES / script), "--output-dir", str(out), *map(str, arguments)])
    result = subprocess.run(
        command,
        cwd=cwd or out.parent,
        capture_output=True,
        text=True,
        timeout=240,
        env={**os.environ, "PYTHONUNBUFFERED": "1"},
    )
    return result


def successful_example(script, out, *arguments, cwd=None, core_only=False):
    result = run_example(script, out, *arguments, cwd=cwd, core_only=core_only)
    assert result.returncode == 0, result.stdout[-4000:] + result.stderr[-4000:]
    summary = json.loads((out / "summary.json").read_text())
    assert summary["environment"]["python"]
    return summary


def test_catalog_covers_every_lesson_and_guide():
    scripts = {str(p.relative_to(EXAMPLES)) for p in EXAMPLES.glob("*/*.py")}
    listed = [entry["script"] for entry in CATALOG]
    assert len(listed) == len(set(listed)) == 37
    assert scripts == set(listed)
    for entry in CATALOG:
        assert (EXAMPLES / entry["guide"]).is_file(), entry
        assert entry["title"] and isinstance(entry["requires"], list)


def test_cookbook_local_links_resolve():
    # These guides use local relative links and simple external HTTPS URLs.
    for document in EXAMPLES.rglob("*.md"):
        if "_outputs" in document.parts:
            continue
        for target in re.findall(r"\]\(([^)]+)\)", document.read_text()):
            target = unquote(target.split("#", 1)[0])
            if not target or "://" in target:
                continue
            assert (document.parent / target).exists(), f"{document}: {target}"


@pytest.mark.parametrize("entry", CATALOG, ids=lambda entry: entry["script"])
def test_help_without_installed_dependencies(entry, tmp_path):
    # -S suppresses site-packages: even optional lessons must explain themselves.
    result = run_example(entry["script"], tmp_path / "unused", "--help", no_site=True)
    assert result.returncode == 0, result.stderr
    assert "--output-dir" in result.stdout
    assert not (tmp_path / "unused").exists()


@pytest.mark.parametrize("entry", CASES)
def test_lesson(entry, tmp_path):
    for dependency in entry["requires"]:
        if importlib.util.find_spec(dependency) is None:
            pytest.skip(f"Optional example dependency not installed: {dependency}")
    integration = entry.get("integration")
    if integration and os.environ.get(INTEGRATION_FLAGS[integration]) != "1":
        pytest.skip(f"Set {INTEGRATION_FLAGS[integration]}=1 to execute this integration.")
    # Alternate roots to cover both documented execution styles.
    cwd = ROOT if CATALOG.index(entry) % 2 == 0 else tmp_path
    before = {
        p: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (EXAMPLES / "data").rglob("*")
        if p.is_file() and p.suffix != ".db"
    }
    args = []
    if entry["script"] == "build/tobacco_build.py" and os.environ.get("MOFFORGE_TEST_TOBACCO_DATA"):
        args = ["--tobacco-data", os.environ["MOFFORGE_TEST_TOBACCO_DATA"]]
    out = tmp_path / "outputs"
    data = successful_example(entry["script"], out, *args, cwd=cwd, core_only=not entry["requires"])
    assert all(hashlib.sha256(p.read_bytes()).hexdigest() == digest for p, digest in before.items())
    name = Path(entry["script"]).stem
    expectations = {
        "structure_io": {"atoms": 424, "bonds": 512, "reloaded_atoms": 424, "fragment_atoms": 10},
        "periodic_geometry": {"periodic_distance_angstrom": pytest.approx(0.32)},
        "pattern_matching": {"locations": 24, "isomorphisms": 96},
        "linker_functionalization": {"parent_atoms": 424, "child_atoms": 466},
        "defect_engineering": {"parent_atoms": 912, "child_atoms": 896},
        "structure_repair": {"parent_atoms": 328, "child_atoms": 424, "added_hydrogens": 96},
        "cleanup_and_repair": {"repaired_atoms": 136, "child_atoms": 104, "removed_guests": 8},
        "cleanup_pipeline": {"child_atoms": 104},
        "symmetry_analysis": {"parent_atoms": 27, "child_atoms": 30, "supercell_atoms": 60},
        "solvent_removal": {"loaded_atoms": 430, "cleaned_atoms": 424, "removed_atoms": 6},
        "placement": {"requested": 2, "placed": 2, "loaded_atoms": 430},
        "provenance_replay": {"same_geometry": True, "seed": 42},
        "coremof_screening": {"synthetic_data": True, "deliberately_empty_count": 0},
    }
    for key, expected in expectations.get(name, {}).items():
        assert data[key] == expected, (name, key)
    if name == "bonding_rules":
        assert data["restricted_bonds"] < data["default_bonds"]
        assert data["restricted_zinc_coordination"] == data["removed_bonds"] == 0
    elif name == "string_pattern_search":
        assert data["patterns"]["C1-C-C-C-C-C-1"]["locations"] == 24
        assert data["patterns"]["[Xe]-[Xe]"]["locations"] == 0
        assert "Unclosed ring" in data["expected_parse_error"]
    elif name == "selective_modification":
        assert {key: value["atoms"] for key, value in data["modes"].items()} == {
            "all_optimal": 472,
            "random_locations": 440,
            "specific_locations": 432,
            "specific_orientations": 432,
            "random_orientations": 472,
        }
    elif name in {"agent_functionalization", "functionalization_campaign"}:
        assert data["results"] and all(row["error"] is None for row in data["results"])
        for row in data["results"]:
            assert row["n_matches"] == 24
            assert row["n_functionalized"] == (0 if row["coverage"] == 0 else 12)
            assert Path(row["output_cif"]).is_file()
    elif name == "site_discovery":
        assert len(data["void_sites"]) == 5
        assert data["framework_open_metal_count"] == data["impossible_clearance_count"] == 0
        assert len(data["illustrative_open_sites"]) == 1
    elif name == "structure_validation":
        assert data["reports"]["IRMOF-1"]["is_valid"]
        assert not data["expected_invalid"]["is_valid"]
        assert not data["expected_unchecked"]["is_valid"]
    elif name == "artifact_verification":
        assert data["complete"]["is_verified"]
        assert data["integrity_only"]["file_integrity"]
        assert not data["integrity_only"]["is_verified"]
        assert not data["altered"]["file_integrity"]
    elif name == "batch_processing":
        assert len(data["results"]) == 2
        assert len({row["output_path"] for row in data["results"]}) == 2
    elif name == "csd_lookup":
        assert data["records"][0]["refcode"] == "DEMOZN"
        assert set(data["searches"]["formula"]) == {"DEMOZN", "DEMOMIS"}
        assert data["searches"]["doi"] == ["DEMOZN"]
    elif name == "coremof_screening":
        assert len(data["candidates"]) == 2
    elif name == "bridge_and_resolve":
        matches = {row["refcode"]: row["variants"] for row in data["matches"]}
        assert len(matches["DEMOZN"]) == 2
        assert matches["DEMONONE"] == []
        assert matches["DEMOMIS"][0]["cif"] is None
    elif name == "interactive_html":
        for path in data["html_files"]:
            html = Path(path).read_text()
            assert "3Dmol" in html and "addCylinder" in html
    elif name == "png_rendering":
        for path in data["images"]:
            content = Path(path).read_bytes()
            assert content[:8] == b"\x89PNG\r\n\x1a\n"
            assert struct.unpack(">II", content[16:24]) == (640, 480)
    elif name == "local_client":
        assert data["advertised"] == ["mofforge_replace", "mofforge_search", "mofforge_validate"]
        assert data["replacement"]["success"] and data["validation"]["success"]
    elif name == "chemgraph_integration":
        assert data["tools"] == ["mofforge_validate"] and data["validation"]["success"]
    elif name in {"pormake_build", "tobacco_build"}:
        assert data["success"] and data["atoms"] > 0 and data["outputs"]
    elif name == "pcu_zn_random_linker":
        assert data["required_connectivity"] == 6 and data["attempts"][-1]["success"]
    elif name == "mcp_build_test":
        assert data["result"]["success"] and data["result"]["output_paths"]
    elif name == "build_and_modify":
        assert data["child_atoms"] == data["parent_atoms"] + 2
        assert data["matched_locations"] > 0
    elif name == "screen_and_place":
        statuses = [row["status"] for row in data["results"]]
        assert sorted(statuses) == ["missing_structure", "placed"]
    # Reparse complete framework artifacts, excluding deliberate overlap/modified-byte controls.
    if name in {"placement", "structure_repair", "linker_functionalization", "build_and_modify"}:
        for path in out.glob("*.cif"):
            assert Crystal.from_cif(path).n_atoms > 0


@pytest.mark.parametrize(
    "script,artifact",
    [
        ("modify/linker_functionalization.py", "acetylamido_IRMOF-1.cif"),
        ("adsorbate/placement.py", "CO2_loaded.cif"),
    ],
)
def test_seeded_lesson_repeats_geometry(script, artifact, tmp_path):
    first = successful_example(script, tmp_path / "first", "--seed", "7")
    second = successful_example(script, tmp_path / "second", "--seed", "7")
    a = Crystal.from_cif(tmp_path / "first" / artifact)
    b = Crystal.from_cif(tmp_path / "second" / artifact)
    assert a.species == b.species
    np.testing.assert_allclose(a.frac_coords, b.frac_coords, atol=1e-8)
    if "locations" in first:
        assert first["locations"] == second["locations"]


def test_batch_serial_and_parallel_agree(tmp_path):
    result = successful_example(
        "pipeline/batch_processing.py", tmp_path / "batch", "--compare-parallel"
    )
    assert result["parallel_matches_serial"]


def test_optional_dependency_error_is_actionable(tmp_path):
    result = run_example("build/smiles_building_blocks.py", tmp_path / "smiles", no_site=True)
    assert result.returncode != 0
    assert "rdkit" in result.stderr and "mofforge[chem]" in result.stderr


def test_expected_positive_search_does_not_silently_accept_no_matches(tmp_path):
    query = tmp_path / "xenon.xyz"
    query.write_text("2\nAbsent tutorial motif\nXe 0 0 0\nXe 1.2 0 0\n")
    result = run_example(
        "search/pattern_matching.py",
        tmp_path / "search",
        "--query",
        query.name,
        "--fragment-path",
        tmp_path,
    )
    assert result.returncode != 0 and "No matches" in result.stderr


def test_real_database_mode_requires_complete_paths(tmp_path):
    result = run_example(
        "databases/bridge_and_resolve.py", tmp_path / "db", "--csd-data", tmp_path / "example.tab"
    )
    assert result.returncode != 0 and "requires" in result.stderr


def test_real_build_failure_is_reported(tmp_path):
    if importlib.util.find_spec("pormake") is None:
        pytest.skip("Pormake is not installed.")
    result = run_example(
        "build/pormake_build.py", tmp_path / "build", "--node", "THIS_BLOCK_DOES_NOT_EXIST.xyz"
    )
    assert result.returncode != 0
    report = json.loads((tmp_path / "build" / "summary.json").read_text())
    assert not report["success"] and report["errors"]


@pytest.mark.integration
def test_external_database_paths(tmp_path):
    if os.environ.get("MOFFORGE_RUN_DATABASE_TESTS") != "1":
        pytest.skip("Set MOFFORGE_RUN_DATABASE_TESTS=1 with local licensed/database paths.")
    required = [
        "MOFFORGE_CSD_DATA_PATH",
        "MOFFORGE_COREMOF_DATA_PATH",
        "MOFFORGE_COREMOF_STRUCTURES_PATH",
    ]
    assert all(os.environ.get(key) for key in required), required
    summary = successful_example(
        "databases/bridge_and_resolve.py",
        tmp_path / "real",
        "--csd-data",
        os.environ[required[0]],
        "--coremof-data",
        os.environ[required[1]],
        "--structures-dir",
        os.environ[required[2]],
        "--name",
        os.environ.get("MOFFORGE_EXAMPLE_MOF_NAME", "HKUST"),
    )
    assert not summary["synthetic_data"]
    assert summary["matches"], "Choose MOFFORGE_EXAMPLE_MOF_NAME present in your export."
