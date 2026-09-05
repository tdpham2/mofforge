"""Input identity, byte integrity, and safe reuse regressions."""

import copy
import json
import shutil
from functools import wraps
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner
from pymatgen.core import Lattice, Structure

from mofforge import Crystal
from mofforge.artifacts import verify_artifact
from mofforge.batch import BatchConfig, _describe_batch_inputs, run_batch
from mofforge.cli import main
from mofforge.inputs import describe_inputs, input_identity, scientific_settings
from mofforge.provenance import Provenance, file_hash, record_operation
from mofforge.utils.config import config


@pytest.fixture
def artifact(tmp_path):
    parent = tmp_path / "parent.cif"
    crystal = Crystal.from_structure(Structure(Lattice.cubic(10), ["C"], [[0, 0, 0]]))
    crystal.write_cif(parent)
    descriptor = describe_inputs("example", {"random_seed": 17}, {"parent": parent})
    record_operation(crystal, crystal, "example", {}, input_descriptor=descriptor)
    output = tmp_path / "result.cif"
    crystal.write_cif(output)
    return output, descriptor


def _edit_manifest(output, edit):
    path = output.with_suffix(output.suffix + ".json")
    data = json.loads(path.read_text())
    updated = edit(data)
    path.write_text(json.dumps(data if updated is None else updated))


def test_verified_artifact_and_cli_are_read_only(artifact):
    output, descriptor = artifact
    paths = [output, output.with_suffix(".cif.json"), Path(descriptor["files"]["parent"]["path"])]
    before = [(p.read_bytes(), p.stat().st_mtime_ns) for p in paths]
    report = verify_artifact(output, expected_inputs=descriptor)
    assert report.is_verified
    result = CliRunner().invoke(main, ["verify", str(output), "--json"])
    assert result.exit_code == 0, result.output
    assert json.loads(result.output) == report.to_dict()
    assert [(p.read_bytes(), p.stat().st_mtime_ns) for p in paths] == before


@pytest.mark.parametrize("contents", ["", "partial output", "changed output\n"])
def test_changed_or_partial_artifact_never_passes(artifact, contents):
    output, _ = artifact
    output.write_text(contents)
    report = verify_artifact(output)
    assert not report.is_verified
    assert not report.file_integrity
    assert report.errors
    assert CliRunner().invoke(main, ["verify", str(output)]).exit_code == 1


@pytest.mark.parametrize("missing", ["artifact", "manifest", "input"])
def test_missing_files_are_explicit(artifact, missing):
    output, descriptor = artifact
    path = {
        "artifact": output,
        "manifest": output.with_suffix(".cif.json"),
        "input": Path(descriptor["files"]["parent"]["path"]),
    }[missing]
    path.unlink()
    report = verify_artifact(output)
    assert not report.is_verified
    assert report.errors


@pytest.mark.parametrize(
    "edit",
    [
        lambda data: data.update(schema_version=99),
        lambda data: data.update(schema_version=True),
        lambda data: data.update(output_file_sha256="missing"),
        lambda data: data.update(input_descriptor=[]),
        lambda data: data["input_descriptor"].update(schema_version=99),
        lambda data: data["input_descriptor"].update(identity_sha256="0" * 64),
        lambda data: data["input_descriptor"].update(files={"parent": {}}),
        lambda data: [],
    ],
)
def test_malformed_manifests_never_pass(artifact, edit):
    output, _ = artifact
    _edit_manifest(output, edit)
    report = verify_artifact(output)
    assert not report.is_verified
    assert report.errors


@pytest.mark.parametrize("raw", ['{"schema_version":', '{"schema_version":NaN}', "\xff"])
def test_unreadable_json_is_reported(artifact, raw):
    output, _ = artifact
    output.with_suffix(".cif.json").write_bytes(raw.encode("latin1"))
    result = CliRunner().invoke(main, ["verify", str(output), "--json"])
    assert result.exit_code == 1
    assert json.loads(result.output)["errors"]


def test_legacy_manifest_has_integrity_only(artifact):
    output, _ = artifact

    def remove_descriptor(data):
        del data["input_descriptor"]

    _edit_manifest(output, remove_descriptor)
    report = verify_artifact(output)
    assert report.file_integrity
    assert not report.inputs_verified
    assert not report.is_verified
    assert "Legacy" in report.warnings[0]
    assert CliRunner().invoke(main, ["verify", str(output)]).exit_code == 1


def test_provenance_alone_is_not_a_manifest(artifact):
    output, _ = artifact
    Provenance(operation="example").to_json(output.with_suffix(".cif.json"))
    report = verify_artifact(output)
    assert not report.is_verified
    assert "provenance alone" in report.errors[0]


def test_stale_inputs_and_requested_parameters_are_detected(artifact):
    output, descriptor = artifact
    expected = copy.deepcopy(descriptor)
    expected["parameters"]["random_seed"] = 18
    assert "requested workflow" in verify_artifact(output, expected_inputs=expected).errors[0]
    parent = Path(descriptor["files"]["parent"]["path"])
    parent.write_text(parent.read_text() + "\n# changed input\n")
    report = verify_artifact(output)
    assert report.file_integrity
    assert not report.inputs_verified
    assert "Input file changed" in report.errors[0]


def test_locations_do_not_define_content_identity(artifact, tmp_path):
    output, descriptor = artifact
    moved_input = tmp_path / "moved-input.cif"
    shutil.copyfile(descriptor["files"]["parent"]["path"], moved_input)
    relocated = describe_inputs("example", {"random_seed": 17}, {"parent": moved_input})
    assert relocated["identity_sha256"] == descriptor["identity_sha256"]
    moved_output = tmp_path / "moved-output.cif"
    shutil.copyfile(output, moved_output)
    report = verify_artifact(moved_output, manifest_path=output.with_suffix(".cif.json"))
    assert report.is_verified  # output_path in the manifest is informational


def test_scientific_settings_and_environment_are_checked(artifact, monkeypatch):
    output, _ = artifact
    monkeypatch.setattr(config, "bond_pad", config.bond_pad + 0.05)
    report = verify_artifact(output)
    assert report.file_integrity and not report.inputs_verified
    assert any("scientific settings" in message for message in report.errors)

    def edit(data):
        descriptor = data["input_descriptor"]
        descriptor["environment"]["implementation_sha256"] = "0" * 64
        descriptor["identity_sha256"] = input_identity(descriptor)

    _edit_manifest(output, edit)
    assert any("environment" in message for message in verify_artifact(output).errors)


def test_fragment_content_and_effective_defaults_change_identity(tmp_path, crystal_dir):
    query, replacement = tmp_path / "query.xyz", tmp_path / "replacement.xyz"
    query.write_text("query contents")
    replacement.write_text("replacement contents")
    batch = BatchConfig(
        operations=[{"type": "replace", "query": query.name, "replacement": replacement.name}],
        random_seed=17,
        moiety_path=str(tmp_path),
    )
    parent = crystal_dir / "IRMOF-1.cif"
    first = _describe_batch_inputs(parent, batch, scientific_settings())
    for path in (query, replacement):
        path.write_text(path.read_text() + " changed")
        changed = _describe_batch_inputs(parent, batch, scientific_settings())
        assert changed["identity_sha256"] != first["identity_sha256"]
        first = changed
    settings = scientific_settings()
    settings["validation"]["bond_tolerance"] += 0.1
    changed = _describe_batch_inputs(parent, batch, settings)
    assert changed["identity_sha256"] != first["identity_sha256"]
    assert changed["parameters"]["final_validation"]["check_coordination"] is True
    assert changed["files"]["operations.0.query"]["sha256"] == file_hash(query)


def test_batch_snapshots_global_settings_for_parallel_workers(
    tmp_path, crystal_dir, moiety_dir, monkeypatch
):
    monkeypatch.setattr(config, "bond_pad", 0.35)
    monkeypatch.setattr(config, "moiety_path", moiety_dir)
    yaml_path = tmp_path / "batch.yaml"
    data = {
        "parents": [str(crystal_dir / "IRMOF-1.cif")],
        "operations": [
            {
                "type": "replace",
                "query": "2-!-p-phenylene.xyz",
                "replacement": "2-nitro-p-phenylene.xyz",
                "mode": "nb_loc_1",
            }
        ],
        "random_seed": 17,
    }
    runs = []
    for workers in (0, 2):
        data.update(parallel=workers, output={"directory": str(tmp_path / str(workers))})
        yaml_path.write_text(yaml.safe_dump(data))
        (result,) = run_batch(yaml_path)
        assert result.success, result.error
        output = Path(result.output_path)
        report = verify_artifact(output)
        assert report.is_verified, report.errors
        runs.append((result.run_id, file_hash(output)))
    assert runs[0] == runs[1]
    assert config.bond_pad == 0.35


def test_replacement_edit_changes_batch_output_and_invalidates_previous_artifact(
    tmp_path, crystal_dir, moiety_dir
):
    for name in ("2-!-p-phenylene.xyz", "2-nitro-p-phenylene.xyz"):
        shutil.copyfile(moiety_dir / name, tmp_path / name)
    yaml_path = tmp_path / "batch.yaml"
    yaml_path.write_text(
        yaml.safe_dump(
            {
                "parents": [str(crystal_dir / "IRMOF-1.cif")],
                "moiety_path": str(tmp_path),
                "random_seed": 17,
                "operations": [
                    {
                        "type": "replace",
                        "query": "2-!-p-phenylene.xyz",
                        "replacement": "2-nitro-p-phenylene.xyz",
                        "mode": "nb_loc_1",
                    }
                ],
                "output": {"directory": str(tmp_path / "outputs")},
            }
        )
    )
    (first,) = run_batch(yaml_path)
    assert first.success, first.error
    replacement = tmp_path / "2-nitro-p-phenylene.xyz"
    lines = replacement.read_text().splitlines()
    lines[1] += " updated source comment"
    replacement.write_text("\n".join(lines) + "\n")
    (second,) = run_batch(yaml_path)
    assert second.success, second.error
    assert first.run_id != second.run_id
    assert first.output_path != second.output_path
    assert file_hash(first.output_path) == file_hash(second.output_path)
    assert not verify_artifact(first.output_path).is_verified
    assert verify_artifact(second.output_path).is_verified


def test_inputs_changed_during_execution_are_not_published(tmp_path, artifact, monkeypatch):
    import mofforge.batch as batch_module

    _, descriptor = artifact
    parent = Path(descriptor["files"]["parent"]["path"])
    original = batch_module.validate_structure

    @wraps(original)
    def change_input(*args, **kwargs):
        parent.write_text(parent.read_text() + "\n# concurrent edit\n")
        return original(*args, **kwargs)

    monkeypatch.setattr(batch_module, "validate_structure", change_input)
    yaml_path = tmp_path / "batch.yaml"
    output_dir = tmp_path / "outputs"
    yaml_path.write_text(
        yaml.safe_dump({"parents": [str(parent)], "output": {"directory": str(output_dir)}})
    )
    (result,) = run_batch(yaml_path)
    assert not result.success
    assert "Input file changed" in result.error
    assert not output_dir.exists()
