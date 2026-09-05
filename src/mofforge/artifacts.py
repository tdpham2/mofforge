"""Read-only verification of artifact integrity and recorded workflow inputs."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

from mofforge.inputs import (
    changed_inputs,
    execution_environment,
    input_identity,
    is_sha256,
    scientific_settings,
)
from mofforge.provenance import file_hash


@dataclass
class VerificationReport:
    artifact_path: str
    manifest_path: str
    schema_version: int | None = None
    file_integrity: bool = False
    inputs_verified: bool = False
    identity_sha256: str | None = None
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def is_verified(self) -> bool:
        return self.file_integrity and self.inputs_verified and not self.errors

    def to_dict(self) -> dict:
        return {**asdict(self), "is_verified": self.is_verified}


def _invalid_constant(value):
    raise ValueError(f"Nonfinite JSON value: {value}")


def verify_artifact(
    artifact_path: str | Path,
    *,
    manifest_path: str | Path | None = None,
    expected_inputs: dict | None = None,
) -> VerificationReport:
    """Verify bytes and inputs, without asserting geometric or physical validity.

    Legacy schema-2 artifacts without input descriptors can pass file integrity
    only. Pass expected_inputs when checking reuse for a particular requested
    workflow. The artifact path is supplied by the caller, so moved output
    bundles are not redirected to the manifest's original output location.
    """
    artifact = Path(artifact_path)
    manifest = (
        Path(manifest_path)
        if manifest_path is not None
        else artifact.with_suffix(artifact.suffix + ".json")
    )
    report = VerificationReport(str(artifact), str(manifest))
    try:
        data = json.loads(manifest.read_text(), parse_constant=_invalid_constant)
        if not isinstance(data, dict):
            raise ValueError("Manifest must be a JSON object.")
        schema = data.get("schema_version")
        if type(schema) is not int or schema != 2:
            raise ValueError("Unsupported artifact schema_version; expected 2.")
        report.schema_version = schema
        if not is_sha256(data.get("output_file_sha256")):
            raise ValueError(
                "Manifest requires output_file_sha256; "
                "provenance alone is not an artifact manifest."
            )
        if artifact.stat().st_size == 0:
            raise ValueError("Artifact is empty.")
        report.file_integrity = file_hash(artifact) == data["output_file_sha256"]
        if not report.file_integrity:
            report.errors.append("Artifact SHA-256 does not match the manifest.")

        descriptor = data.get("input_descriptor")
        if descriptor is None:
            report.warnings.append(
                "Legacy manifest has no input descriptor: file integrity only; "
                "workflow inputs and reuse cannot be verified."
            )
            return report
        if not isinstance(descriptor, dict):
            raise ValueError("input_descriptor must be a JSON object.")
        identity = input_identity(descriptor)
        if descriptor.get("identity_sha256") != identity:
            raise ValueError("Input descriptor identity_sha256 does not match its contents.")
        report.identity_sha256 = identity
        input_errors = changed_inputs(descriptor)
        if descriptor["scientific_settings"] != scientific_settings():
            input_errors.append("Effective scientific settings differ from the recorded inputs.")
        if descriptor["environment"] != execution_environment():
            input_errors.append("Execution environment differs from the recorded inputs.")
        if expected_inputs is not None and input_identity(expected_inputs) != identity:
            input_errors.append("Recorded inputs differ from the requested workflow.")
        report.inputs_verified = not input_errors
        report.errors.extend(input_errors)
    except (OSError, ValueError, TypeError, UnicodeError) as exc:
        report.errors.append(str(exc))
    return report
