"""CLI build option routing without optional construction dependencies."""

import json
from types import SimpleNamespace
from unittest.mock import Mock, patch

import pytest
from click.testing import CliRunner
from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter

from mofforge.build.base import BuildResult, Topology
from mofforge.build.pormake_backend import PormakeBackend
from mofforge.cli import main
from mofforge.core.crystal import Crystal
from mofforge.provenance import file_hash


@pytest.fixture
def pormake_engine(tmp_path, monkeypatch):
    """Stub Pormake itself while exercising the real facade and backend."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("mofforge.build.config._find_toml", lambda: None)
    monkeypatch.delenv("MOFFORGE_PORMAKE_OUTPUT_DIR", raising=False)

    topology = SimpleNamespace(unique_node_types=[0], unique_edge_types=[(0, 0)])
    database = Mock(topo_dir=tmp_path / "topologies", bb_dir=tmp_path / "blocks")
    database.get_topo.return_value = topology
    database.get_bb.side_effect = lambda name: name

    structure = Structure(Lattice.cubic(10), ["C", "C"], [[0, 0, 0], [0.15, 0, 0]])
    framework = SimpleNamespace(write_cif=CifWriter(structure).write_file, info={})
    engine = Mock()
    engine.build_by_type.return_value = framework
    module = SimpleNamespace(Database=lambda: database, Builder=lambda: engine)
    monkeypatch.setattr("mofforge.build.pormake_backend._pm", module)
    return engine


@pytest.fixture
def pormake_args(tmp_path):
    return [
        "build", "--backend", "pormake", "--topology", "pcu",
        "--node", "N108", "--edge", "E1", "--output", str(tmp_path / "output"),
    ]


@pytest.mark.parametrize("verbose", [False, True])
def test_pormake_cli_build_writes_artifacts(pormake_engine, pormake_args, tmp_path, verbose):
    args = pormake_args + (["--verbose"] if verbose else [])
    result = CliRunner().invoke(main, args)

    assert result.exit_code == 0, result.output
    assert "Build succeeded" in result.output
    pormake_engine.build_by_type.assert_called_once()
    assert "verbose" not in pormake_engine.build_by_type.call_args.kwargs

    outputs = list((tmp_path / "output").glob("*.cif"))
    assert len(outputs) == 1
    output = outputs[0]
    assert Crystal.from_cif(output).n_atoms == 2
    manifest = json.loads(output.with_suffix(".cif.json").read_text())
    assert manifest["provenance"]["operation"] == "build"
    assert manifest["provenance"]["parameters"]["backend"] == "pormake"
    assert manifest["output_file_sha256"] == file_hash(output)
    assert manifest["validation"]["is_valid"]


@pytest.mark.parametrize("verbose", [False, True])
def test_pormake_cli_reports_engine_failure(pormake_engine, pormake_args, tmp_path, verbose):
    pormake_engine.build_by_type.side_effect = RuntimeError("framework optimization failed")
    args = pormake_args + (["--verbose"] if verbose else [])
    result = CliRunner().invoke(main, args)

    assert result.exit_code == 1
    assert "Pormake build failed: framework optimization failed" in result.output
    pormake_engine.build_by_type.assert_called_once()
    assert not list((tmp_path / "output").glob("*.cif"))


@pytest.mark.parametrize("verbose", [False, True])
def test_tobacco_cli_forwards_verbosity(tmp_path, verbose):
    args = [
        "build", "--backend", "tobacco", "--topology", "pcu",
        "--output", str(tmp_path),
    ] + (["--verbose"] if verbose else [])

    with patch("mofforge.build.MOFBuilder", autospec=True) as builder_class:
        builder = builder_class.return_value
        builder.build.return_value = BuildResult(success=True, backend="tobacco")
        result = CliRunner().invoke(main, args)

    assert result.exit_code == 0, result.output
    builder_class.assert_called_once_with(backend="tobacco")
    builder.build.assert_called_once_with(
        topology="pcu", output_dir=str(tmp_path), verbose=verbose
    )


def test_pormake_still_rejects_unknown_options(pormake_engine, tmp_path):
    result = PormakeBackend(output_dir=tmp_path).build(
        Topology("pcu"), nodes=[], edges=[], unsupported_option=True
    )

    assert not result.success
    assert result.errors == ["Unknown build options: ['unsupported_option']"]
    pormake_engine.build_by_type.assert_not_called()
