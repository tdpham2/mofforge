"""Versioned input identities, independent of artifact locations."""

from __future__ import annotations

import copy
import inspect
import platform
from contextlib import contextmanager
from pathlib import Path

from mofforge.provenance import content_hash, file_hash, software_versions

INPUT_SCHEMA_VERSION = 1


def parameter_defaults(function, *excluded: str) -> dict:
    """Capture defaults from the function actually used for execution."""
    return {
        name: copy.deepcopy(parameter.default)
        for name, parameter in inspect.signature(function).parameters.items()
        if name not in excluded and parameter.default is not inspect.Parameter.empty
    }


def scientific_settings() -> dict:
    from mofforge.utils.config import COVALENT_RADII, VDW_RADII, config
    from mofforge.validation import EXPECTED_COORDINATION, validate_structure

    return {
        "bond_pad": config.bond_pad,
        "r_tag": config.r_tag,
        "covalent_radii": dict(COVALENT_RADII),
        "vdw_radii": dict(VDW_RADII),
        "expected_coordination": {k: list(v) for k, v in EXPECTED_COORDINATION.items()},
        "validation": parameter_defaults(validate_structure),
    }


def execution_environment() -> dict:
    """Include source changes in editable installs as well as package versions."""
    root = Path(__file__).parent
    return {
        "versions": dict(software_versions()),
        "system": platform.system(),
        "machine": platform.machine(),
        "python_implementation": platform.python_implementation(),
        "implementation_sha256": content_hash(
            {str(p.relative_to(root)): file_hash(p) for p in sorted(root.rglob("*.py"))}
        ),
    }


def input_identity(descriptor: dict) -> str:
    """Hash the full recipe and file contents, excluding their locations."""
    if type(descriptor.get("schema_version")) is not int or descriptor["schema_version"] != 1:
        raise ValueError("Unsupported input descriptor schema_version.")
    if not isinstance(descriptor.get("operation"), str) or not descriptor["operation"]:
        raise ValueError("Input descriptor requires an operation.")
    for key in ("files", "parameters", "scientific_settings", "environment"):
        if not isinstance(descriptor.get(key), dict):
            raise ValueError(f"Input descriptor requires a {key} object.")
    for role, entry in descriptor["files"].items():
        if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
            raise ValueError(f"Input file {role} requires a path.")
        if not is_sha256(entry.get("sha256")):
            raise ValueError(f"Input file {role} requires a SHA-256 hash.")
    payload = {k: v for k, v in descriptor.items() if k != "identity_sha256"}
    payload["files"] = {
        role: {k: v for k, v in entry.items() if k != "path"}
        for role, entry in descriptor["files"].items()
    }
    return content_hash(payload)


def is_sha256(value) -> bool:
    return (
        isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)
    )


def describe_inputs(
    operation: str,
    parameters: dict,
    files: dict[str, str | Path],
    *,
    settings: dict | None = None,
) -> dict:
    descriptor = {
        "schema_version": INPUT_SCHEMA_VERSION,
        "operation": operation,
        "parameters": copy.deepcopy(parameters),
        "files": {
            role: {"path": str(Path(path).resolve()), "sha256": file_hash(path)}
            for role, path in files.items()
        },
        "scientific_settings": copy.deepcopy(settings)
        if settings is not None
        else scientific_settings(),
        "environment": execution_environment(),
    }
    descriptor["identity_sha256"] = input_identity(descriptor)
    return descriptor


def changed_inputs(descriptor: dict) -> list[str]:
    """Check recorded source locations without modifying files."""
    errors = []
    for role, entry in descriptor["files"].items():
        path = Path(entry["path"])
        try:
            if file_hash(path) != entry["sha256"]:
                errors.append(f"Input file changed: {role} ({path}).")
        except OSError as exc:
            errors.append(f"Input file unavailable: {role} ({path}): {exc}")
    return errors


@contextmanager
def use_scientific_settings(settings: dict):
    """Apply a coordinator snapshot in a batch worker, restoring it on exit.

    Batch workers use separate processes. This is not a thread-local override
    for concurrent callers of the existing global configuration API.
    """
    from mofforge.core.bonding import default_bonding_rules
    from mofforge.utils.config import COVALENT_RADII, VDW_RADII, config
    from mofforge.validation import EXPECTED_COORDINATION

    previous = scientific_settings()

    def apply(snapshot):
        config.bond_pad, config.r_tag = snapshot["bond_pad"], snapshot["r_tag"]
        for target, key in (
            (COVALENT_RADII, "covalent_radii"),
            (VDW_RADII, "vdw_radii"),
            (EXPECTED_COORDINATION, "expected_coordination"),
        ):
            target.clear()
            target.update(
                {k: tuple(v) for k, v in snapshot[key].items()}
                if key == "expected_coordination"
                else copy.deepcopy(snapshot[key])
            )
        default_bonding_rules.cache_clear()

    try:
        apply(settings)
        yield
    finally:
        apply(previous)
