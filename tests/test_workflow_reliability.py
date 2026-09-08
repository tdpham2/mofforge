"""Artifact, interface, and execution-order guarantees."""

import json
import shutil

import numpy as np
import pytest
import yaml
from click.testing import CliRunner
from pymatgen.core import Lattice, Structure

from mofforge import Crystal, Pipeline
from mofforge.batch import run_batch
from mofforge.build.base import BuildResult
from mofforge.build.results import finalize_build
from mofforge.cli import main
from mofforge.mcp._impl import validate_impl
from mofforge.provenance import file_hash, structure_hash


def test_invalid_geometry_has_same_report_in_cli_and_mcp(tmp_path):
    path = tmp_path / "overlap.cif"
    Crystal.from_structure(
        Structure(Lattice.cubic(10), ["C", "C"], [[0, 0, 0], [0.05, 0, 0]])
    ).write_cif(path)
    response = CliRunner().invoke(main, ["validate", str(path), "--json"])
    assert response.exit_code == 1
    report = json.loads(response.output)
    mcp = validate_impl(str(path))
    assert mcp["success"]  # validation executed successfully
    assert not mcp["is_valid"]
    assert json.loads(json.dumps(mcp["validation"])) == report


def test_builder_cannot_succeed_without_parseable_output(tmp_path):
    missing = tmp_path / "missing.cif"
    broken = tmp_path / "broken.cif"
    broken.write_text("not a CIF")
    for paths in ([], [missing], [broken]):
        result = finalize_build(BuildResult(success=True, output_paths=paths), {})
        assert not result.success
        assert result.errors


def test_builder_reports_validation_separately(tmp_path):
    path = tmp_path / "overlap.cif"
    Crystal.from_structure(
        Structure(Lattice.cubic(10), ["C", "C"], [[0, 0, 0], [0.05, 0, 0]])
    ).write_cif(path)
    result = finalize_build(BuildResult(success=True, output_paths=[path]), {"topology": "test"})
    assert result.success
    assert not result.validation.is_valid
    manifest = json.loads(path.with_suffix(".cif.json").read_text())
    assert manifest["output_file_sha256"] == file_hash(path)
    assert not manifest["validation"]["is_valid"]


def test_pipeline_seed_is_repeatable(crystal_dir, moiety_dir):
    pipeline = Pipeline(crystal_dir / "IRMOF-1.cif", fragment_path=moiety_dir, random_seed=17)
    pipeline.replace(query="2-!-p-phenylene.xyz", replacement="2-nitro-p-phenylene.xyz", nb_loc=2)
    pipeline.validate()
    first = pipeline.build()
    second = pipeline.build()
    assert structure_hash(first) == structure_hash(second)
    assert first.provenance.history[0]["operation"] == "load"
    assert first.provenance.validation == pipeline.validation_reports[0].to_dict()


@pytest.mark.parametrize("interface", ["cli", "mcp", "pipeline", "batch"])
def test_oversized_replacement_fails_without_artifacts(
    interface, tmp_path, crystal_dir, moiety_dir
):
    if interface == "mcp":
        pytest.importorskip("mcp")

    parent = crystal_dir / "IRMOF-1.cif"
    query = moiety_dir / "2-!-p-phenylene.xyz"
    replacement = moiety_dir / "2-nitro-p-phenylene.xyz"
    output = tmp_path / "outputs" / "modified.cif"
    output.parent.mkdir()
    expected = "Requested 25 replacement locations, but only 24 valid locations are available."
    parent_hash = file_hash(parent)

    if interface == "cli":
        result = CliRunner().invoke(main, [
            "replace", "-p", str(parent), "-q", str(query), "-r", str(replacement),
            "-o", str(output), "--nb-loc", "25",
        ])
        assert result.exit_code == 1
        assert f"Error: {expected}" in result.output
        assert "Output written" not in result.output
    elif interface == "mcp":
        import asyncio

        from mcp.server.fastmcp.exceptions import ToolError

        from mofforge.mcp.server import build_server

        server = build_server({"mofforge_replace"})
        with pytest.raises(ToolError) as error:
            asyncio.run(server.call_tool("mofforge_replace", {
                "parent_cif": str(parent), "query_xyz": str(query),
                "replacement_xyz": str(replacement), "output_cif": str(output), "nb_loc": 25,
            }))
        assert expected in str(error.value)
    elif interface == "pipeline":
        pipeline = Pipeline(parent, fragment_path=moiety_dir)
        pipeline.replace(query=query.name, replacement=replacement.name, nb_loc=25)
        with pytest.raises(ValueError) as error:
            pipeline.build().write_cif(output)
        assert str(error.value) == expected
    else:
        config = tmp_path / "batch.yaml"
        config.write_text(yaml.safe_dump({
            "parents": [str(parent)], "moiety_path": str(moiety_dir),
            "operations": [{"type": "replace", "query": query.name,
                            "replacement": replacement.name, "mode": "nb_loc_25"}],
            "output": {"directory": str(output.parent)},
        }))
        results = run_batch(config)
        assert len(results) == 1
        assert not results[0].success
        assert results[0].error == expected
        assert results[0].output_path is None

    assert not list(tmp_path.rglob("*.cif"))
    assert not list(tmp_path.rglob("*.json"))
    assert file_hash(parent) == parent_hash


def test_batch_parallel_is_repeatable_and_names_do_not_collide(tmp_path, crystal_dir, moiety_dir):
    parents = []
    for folder in ("a", "b"):
        source = tmp_path / folder / "same.cif"
        source.parent.mkdir()
        shutil.copy(crystal_dir / "IRMOF-1.cif", source)
        parents.append(str(source))
    config = {
        "parents": parents,
        "operations": [
            {
                "type": "replace",
                "query": "2-!-p-phenylene.xyz",
                "replacement": "2-nitro-p-phenylene.xyz",
                "mode": "nb_loc_2",
            }
        ],
        "moiety_path": str(moiety_dir),
        "random_seed": 31,
    }
    batches = []
    for parallel in (0, 2):
        config["parallel"] = parallel
        config["output"] = {"directory": str(tmp_path / f"outputs{parallel}")}
        yaml_path = tmp_path / "batch.yaml"
        yaml_path.write_text(yaml.safe_dump(config))
        results = run_batch(yaml_path)
        assert len(results) == 2
        assert all(r.success for r in results), [r.error for r in results]
        assert len({r.output_path for r in results}) == 2
        batches.append([(r.run_id, file_hash(r.output_path)) for r in results])
    assert batches[0] == batches[1]


def test_failed_and_empty_batches_return_nonzero(tmp_path):
    config = tmp_path / "batch.yaml"
    config.write_text(yaml.safe_dump({"parents": [str(tmp_path / "absent.cif")]}))
    assert CliRunner().invoke(main, ["batch", "-c", str(config)]).exit_code == 1


def test_seeded_functionalization_and_zero_coverage(tmp_path, crystal_dir):
    pytest.importorskip("rdkit")
    from mofforge.functionalize import functionalize, run_campaign

    source = str(crystal_dir / "IRMOF-1.cif")
    smiles = "O=C(O)c1ccc(C(=O)O)cc1"
    first = functionalize(source, smiles, "OCH3", coverage=0.25, random_seed=29)
    second = functionalize(source, smiles, "OCH3", coverage=0.25, random_seed=29)
    assert not first.error and not second.error
    np.testing.assert_allclose(first.crystal.frac_coords, second.crystal.frac_coords, atol=1e-12)
    untouched = functionalize(source, smiles, "F", coverage=0, random_seed=29)
    assert untouched.n_functionalized == 0
    assert structure_hash(untouched.crystal) == structure_hash(Crystal.from_cif(source))
    campaign = run_campaign(
        source, smiles, ["F"], coverages=[0.251, 0.254], output_dir=str(tmp_path), random_seed=3
    )
    assert all(not r.error for r in campaign)
    assert len({r.output_cif for r in campaign}) == 2
