"""Shared artifact verification for construction backends."""

from mofforge.build.base import BuildResult
from mofforge.core.bonding import infer_bonds
from mofforge.core.crystal import Crystal
from mofforge.provenance import Provenance, software_versions, structure_hash, write_manifest
from mofforge.validation import validate_structure


def finalize_build(result: BuildResult, parameters: dict) -> BuildResult:
    """Success means artifacts exist and parse; validation remains separate."""
    result.crystal = None
    reports = {}
    if not result.output_paths:
        result.errors.append("Builder produced no output structures.")
    for path in result.output_paths:
        try:
            crystal = Crystal.from_cif(path)
            if crystal.n_atoms == 0:
                raise ValueError("Output structure has no atoms.")
            report = validate_structure(crystal)
            if crystal.structure.is_ordered and "geometry" in report.checks_performed:
                crystal = infer_bonds(crystal)
            crystal.provenance = Provenance(
                operation="build",
                parameters={**parameters, "backend": result.backend},
                output_hash=structure_hash(crystal),
                software_versions=dict(software_versions()),
                validation=report.to_dict(),
            )
            write_manifest(crystal, path, report)
            reports[str(path)] = report.to_dict()
            if result.crystal is None:
                result.crystal = crystal
                result.validation = report
        except Exception as exc:
            result.errors.append(f"Unusable output {path}: {exc}")
    result.metadata["validation_by_output"] = reports
    result.success = result.crystal is not None and not result.errors
    return result
