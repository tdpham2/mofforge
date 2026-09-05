"""Provenance tracking for crystal modifications."""

from __future__ import annotations

import copy
import hashlib
import json
import logging
import platform
import secrets
from dataclasses import dataclass, field
from datetime import datetime, timezone
from functools import lru_cache
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

logger = logging.getLogger("mofforge")


@dataclass
class Provenance:
    """Metadata tracking modifications made to a crystal structure."""

    parent: str | None = None
    query: str | None = None
    replacement: str | None = None
    operation: str | None = None
    parameters: dict = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    alignment_errors: list[float] = field(default_factory=list)
    history: list[dict] = field(default_factory=list)

    schema_version: int = 2
    input_hash: str | None = None
    output_hash: str | None = None
    software_versions: dict[str, str] = field(default_factory=dict)
    validation: dict | None = None

    def to_dict(self) -> dict:
        """Convert to a plain dictionary."""
        # Legacy records nested the whole earlier history inside each entry.
        # Their top-level list already contains those earlier operations.
        data = {
            key: copy.deepcopy(value) for key, value in self.__dict__.items() if key != "history"
        }
        data["history"] = [
            {key: copy.deepcopy(value) for key, value in entry.items() if key != "history"}
            for entry in self.history
        ]
        return data

    def to_json(self, filepath: str | Path) -> None:
        """Write provenance to a JSON file."""
        filepath = Path(filepath)
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, indent=2, default=_json_default, allow_nan=False)

    @classmethod
    def from_dict(cls, d: dict) -> Provenance:
        """Create a Provenance from a dictionary."""
        return cls(
            parent=d.get("parent"),
            query=d.get("query"),
            replacement=d.get("replacement"),
            operation=d.get("operation"),
            parameters=d.get("parameters", {}),
            timestamp=d.get("timestamp", datetime.now(timezone.utc).isoformat()),
            alignment_errors=d.get("alignment_errors", []),
            history=[
                {k: copy.deepcopy(v) for k, v in entry.items() if k != "history"}
                for entry in d.get("history", [])
            ],
            input_hash=d.get("input_hash"),
            output_hash=d.get("output_hash"),
            software_versions=d.get("software_versions", {}),
            validation=d.get("validation"),
        )

    @classmethod
    def from_json(cls, filepath: str | Path) -> Provenance:
        """Load provenance from a JSON file."""
        with open(filepath, encoding="utf-8") as f:
            d = json.load(f)
        return cls.from_dict(d)

    def chain(self, new_provenance: Provenance) -> Provenance:
        """Create a chained provenance record."""
        current = self.to_dict()
        history = current.pop("history")
        following = new_provenance.to_dict()["history"]
        new_provenance.history = [*history, current, *following]
        return new_provenance

    def summary(self) -> str:
        """Return a human-readable summary of the provenance chain."""
        lines = []
        for i, hist in enumerate(self.history):
            lines.append(
                f"  Step {i + 1}: {hist.get('operation', '?')} ({hist.get('timestamp', '?')})"
            )
        lines.append(f"  Current: {self.operation} ({self.timestamp})")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"Provenance(op='{self.operation}', "
            f"parent='{self.parent}', "
            f"query='{self.query}', "
            f"replacement='{self.replacement}')"
        )


def _json_default(value):
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "tolist"):
        return value.tolist()
    raise TypeError(f"Cannot serialize {type(value).__name__} in provenance.")


def content_hash(value) -> str:
    """Stable SHA-256 identity for JSON-serializable workflow inputs."""
    raw = json.dumps(
        value, sort_keys=True, separators=(",", ":"), default=_json_default, allow_nan=False
    )
    return hashlib.sha256(raw.encode()).hexdigest()


def file_hash(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def structure_hash(crystal) -> str:
    return content_hash(crystal.structure.as_dict())


@lru_cache(maxsize=1)
def software_versions() -> dict[str, str]:
    versions = {"python": platform.python_version()}
    for package in (
        "mofforge",
        "numpy",
        "scipy",
        "networkx",
        "pymatgen",
        "rdkit",
        "pormake",
        "tobacco3",
    ):
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            continue
    return versions


def effective_seed(random_seed: int | None) -> int:
    """Choose and retain a seed even when the caller leaves it unspecified."""
    if random_seed is None:
        return secrets.randbelow(2**31 - 1)
    if not isinstance(random_seed, int) or not 0 <= random_seed < 2**31 - 1:
        raise ValueError("random_seed must be an integer in [0, 2147483647).")
    return random_seed


def derive_seed(random_seed: int | None, identity) -> int | None:
    """Derive independent task seeds without touching process-global RNGs."""
    if random_seed is None:
        return None
    return int(content_hash([random_seed, identity])[:8], 16) % (2**31 - 1)


def record_operation(
    crystal, parent, operation: str, parameters: dict, *, validation=None
) -> Provenance:
    record = Provenance(
        parent=parent.name,
        operation=operation,
        parameters=copy.deepcopy(parameters),
        input_hash=structure_hash(parent),
        output_hash=structure_hash(crystal),
        software_versions=dict(software_versions()),
        validation=validation.to_dict() if validation is not None else None,
    )
    if parent.provenance is not None:
        record = parent.provenance.chain(record)
    crystal.provenance = record
    return record


def write_manifest(crystal, output_path: str | Path, validation=None) -> Path:
    """Write an artifact manifest beside a workflow-generated structure."""
    output = Path(output_path)
    manifest = output.with_suffix(output.suffix + ".json")
    data = {
        "schema_version": 2,
        "output_path": str(output),
        "output_file_sha256": file_hash(output),
        "structure_sha256": structure_hash(crystal),
        "software_versions": software_versions(),
        "provenance": crystal.provenance.to_dict() if crystal.provenance is not None else None,
        "validation": validation.to_dict()
        if validation is not None
        else (crystal.provenance.validation if crystal.provenance is not None else None),
    }
    manifest.write_text(json.dumps(data, indent=2, default=_json_default, allow_nan=False) + "\n")
    return manifest


def write_generation_manifest(output_path: str | Path, operation: str, parameters: dict) -> Path:
    """Record generated building blocks, which can contain nonphysical dummy atoms."""
    output = Path(output_path)
    provenance = Provenance(
        operation=operation,
        parameters=parameters,
        input_hash=content_hash(parameters),
        output_hash=file_hash(output),
        software_versions=dict(software_versions()),
    )
    manifest = output.with_suffix(output.suffix + ".json")
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 2,
                "output_path": str(output),
                "output_file_sha256": file_hash(output),
                "provenance": provenance.to_dict(),
                "software_versions": software_versions(),
                "validation": None,
            },
            indent=2,
            default=_json_default,
            allow_nan=False,
        )
        + "\n"
    )
    return manifest
