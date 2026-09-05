"""Batch processing for multiple crystal structures."""

from __future__ import annotations

import logging
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from glob import glob
from pathlib import Path

import yaml

from mofforge.core.bonding import infer_bonds
from mofforge.core.crystal import Crystal
from mofforge.core.moiety import fragment, resolve_fragment_path
from mofforge.inputs import (
    changed_inputs,
    describe_inputs,
    parameter_defaults,
    scientific_settings,
    use_scientific_settings,
)
from mofforge.provenance import (
    content_hash,
    derive_seed,
    effective_seed,
    file_hash,
    record_operation,
)
from mofforge.replace.replace import replace_pattern
from mofforge.search.search import find_pattern
from mofforge.validation import ValidationReport, validate_structure

logger = logging.getLogger("mofforge")


@dataclass
class BatchResult:
    """Result of processing a single structure in a batch."""

    parent_name: str
    output_path: str | None = None
    success: bool = True
    error: str | None = None
    validation: ValidationReport | None = None
    run_id: str = ""


@dataclass
class BatchConfig:
    """Configuration for a batch processing run."""

    parent_paths: list[str] = field(default_factory=list)
    operations: list[dict] = field(default_factory=list)
    output_dir: str = "results"
    output_format: str = "cif"
    naming: str = "{parent_name}_modified"
    parallel: int = 0
    moiety_path: str | None = None
    random_seed: int | None = None
    scientific_snapshot: dict | None = field(default=None, repr=False)
    resolved_moiety_path: str | None = field(default=None, repr=False)

    _VALID_FORMATS = ("cif", "xyz")

    def __post_init__(self) -> None:
        if self.output_format not in self._VALID_FORMATS:
            raise ValueError(
                f"Unsupported output format '{self.output_format}'. "
                f"Must be one of: {', '.join(self._VALID_FORMATS)}."
            )

    @classmethod
    def from_yaml(cls, filepath: str | Path) -> BatchConfig:
        """Load configuration from a YAML file.

        Expected YAML format::

            parents:
              - path: "structures/*.cif"
            operations:
              - type: replace
                query: BDC.xyz
                replacement: NH2-BDC.xyz
              - type: validate
            output:
              directory: results/
              format: cif
              naming: "{parent_name}_functionalized"
            parallel: 4
            moiety_path: ./data/moieties
        """
        with open(filepath, encoding="utf-8") as f:
            raw = yaml.safe_load(f) or {}

        parent_paths = []
        for p in raw.get("parents", []):
            if isinstance(p, str):
                parent_paths.append(p)
            elif isinstance(p, dict):
                parent_paths.append(p.get("path", ""))

        output = raw.get("output", {})

        # Support both 'moiety_path' (canonical) and 'fragment_path' (alias)
        moiety_path = raw.get("moiety_path")
        if moiety_path is None:
            fragment_path = raw.get("fragment_path")
            if fragment_path is not None:
                warnings.warn(
                    "The 'fragment_path' YAML key is deprecated. Use 'moiety_path' instead.",
                    DeprecationWarning,
                    stacklevel=2,
                )
                moiety_path = fragment_path

        return cls(
            parent_paths=parent_paths,
            operations=raw.get("operations", []),
            output_dir=output.get("directory", "results"),
            output_format=output.get("format", "cif"),
            naming=output.get("naming", "{parent_name}_modified"),
            parallel=raw.get("parallel", 0),
            moiety_path=moiety_path,
            random_seed=raw.get("random_seed"),
        )


def _resolve_parent_paths(patterns: list[str]) -> list[Path]:
    """Resolve glob patterns to actual file paths."""
    paths = []
    for pattern in patterns:
        matches = glob(pattern)
        if matches:
            paths.extend(Path(m) for m in sorted(matches))
        else:
            p = Path(pattern)
            if p.exists():
                paths.append(p)
            else:
                logger.warning("No files matched pattern: %s", pattern)
    return sorted(set(paths), key=lambda p: str(p.resolve()))


_VALID_OP_TYPES = ("replace", "remove", "validate", "desolvate")


def _describe_batch_inputs(parent_path: Path, config: BatchConfig, settings: dict) -> dict:
    """Capture consumed files and the effective options used by this batch."""
    from mofforge.solvent.removal import remove_solvent

    files = {"parent": parent_path}
    # Retain the 0.2 seed derivation so existing seeded geometries do not change.
    seed_identity = [
        str(parent_path.resolve()),
        file_hash(parent_path),
        config.operations,
        config.random_seed,
        config.moiety_path,
        config.output_format,
    ]
    operations = []
    for index, op in enumerate(config.operations):
        kind = op.get("type", "")
        if kind not in _VALID_OP_TYPES:
            raise ValueError(f"unknown operation type '{kind}'")
        entry = {"type": kind}
        if kind in {"replace", "remove"}:
            names = (
                {"query": op.get("query"), "replacement": op.get("replacement")}
                if kind == "replace"
                else {"guest": op.get("guest") or op.get("query")}
            )
            for role, name in names.items():
                key = f"operations.{index}.{role}"
                entry[role] = key if name is not None else None
                if name is not None:
                    files[key] = resolve_fragment_path(
                        name, config.resolved_moiety_path or config.moiety_path
                    )
            options = parameter_defaults(replace_pattern, "name", "verbose", "random_seed")
            mode = op.get("mode", "all_optimal") if kind == "replace" else "all_optimal"
            if mode == "random":
                options["random"] = True
            elif mode.startswith("nb_loc_"):
                options["nb_loc"] = int(mode.split("_")[-1])
            elif mode != "all_optimal":
                raise ValueError(f"Unknown replacement mode: {mode}")
            options["random_seed"] = effective_seed(
                op.get("random_seed", derive_seed(config.random_seed, [seed_identity, index]))
            )
            entry["options"] = options
        elif kind == "desolvate":
            entry["options"] = {
                **parameter_defaults(remove_solvent),
                **{k: v for k, v in op.items() if k != "type"},
            }
        else:
            entry["options"] = {
                **settings["validation"],
                **{k: v for k, v in op.items() if k != "type"},
            }
        operations.append(entry)
    return describe_inputs(
        "batch",
        {
            "operations": operations,
            "random_seed": config.random_seed,
            "seed_strategy": "batch-v1",
            "output_format": config.output_format,
            "fragment_options": {"presort": True, "periodic": False},
            "final_validation": settings["validation"],
        },
        files,
        settings=settings,
    )


def _process_single(
    parent_path: Path,
    config: BatchConfig,
) -> BatchResult:
    """Process a single parent structure through all operations."""
    settings = config.scientific_snapshot or scientific_settings()
    with use_scientific_settings(settings):
        return _execute_single(parent_path, config, settings)


def _execute_single(parent_path: Path, config: BatchConfig, settings: dict) -> BatchResult:
    parent_name = parent_path.stem
    result = BatchResult(parent_name=parent_name)

    try:
        descriptor = _describe_batch_inputs(parent_path, config, settings)
        # Keep source-instance filenames distinct when identical inputs occur
        # at two locations. The descriptor itself excludes file locations.
        result.run_id = content_hash([descriptor["identity_sha256"], str(parent_path.resolve())])[
            :12
        ]
        current = Crystal.from_cif(parent_path)
        current = infer_bonds(current, periodic=True)

        for step_index, op in enumerate(config.operations):
            op_type = op.get("type", "")
            effective_op = descriptor["parameters"]["operations"][step_index]
            options = effective_op["options"]

            if op_type not in _VALID_OP_TYPES:
                raise ValueError(f"unknown operation type '{op_type}'")

            if op_type == "replace":
                root = config.resolved_moiety_path or config.moiety_path
                q = fragment(op.get("query"), fragment_path=root)
                r = fragment(op.get("replacement"), fragment_path=root)

                match = find_pattern(q, current)
                current = replace_pattern(match, r, **options)

                # Re-infer bonds for subsequent steps
                if current.n_bonds == 0 and current.n_atoms > 0:
                    current = infer_bonds(current, periodic=True)

            elif op_type == "remove":
                guest_name = op.get("guest") or op.get("query")
                g = fragment(
                    guest_name, fragment_path=config.resolved_moiety_path or config.moiety_path
                )
                match = find_pattern(g, current, disconnected_component=True)
                current = replace_pattern(match, None, **options)

                if current.n_bonds == 0 and current.n_atoms > 0:
                    current = infer_bonds(current, periodic=True)

            elif op_type == "desolvate":
                from mofforge.solvent.removal import remove_solvent

                sol_result = remove_solvent(current, **options)
                current = sol_result.crystal

                if current.n_bonds == 0 and current.n_atoms > 0:
                    current = infer_bonds(current, periodic=True)

            elif op_type == "validate":
                report = validate_structure(current, **options)
                result.validation = report
                record_operation(current, current, "validate", op, validation=report)

        # Always describe the final artifact, even if validation also appeared
        # earlier in the pipeline.
        result.validation = validate_structure(current, **settings["validation"])
        input_errors = changed_inputs(descriptor)
        if input_errors:
            raise ValueError(" ".join(input_errors))
        record_operation(
            current,
            current,
            "batch_result",
            {"run_id": result.run_id, "random_seed": config.random_seed},
            validation=result.validation,
            input_descriptor=descriptor,
        )
        # Write output
        output_dir = Path(config.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        output_name = config.naming.format(parent_name=parent_name)
        output_path = output_dir / f"{output_name}_{result.run_id}.{config.output_format}"
        if config.output_format == "xyz":
            current.write_xyz(output_path)
        else:
            current.write_cif(output_path)
        result.output_path = str(output_path)

        logger.debug("Processed '%s' -> '%s'", parent_name, output_path)

    except Exception as e:
        result.success = False
        result.error = str(e)
        logger.error("Failed to process '%s': %s", parent_name, e)

    return result


def run_batch(config_path: str | Path) -> list[BatchResult]:
    """Run batch processing from a YAML configuration file."""
    config = BatchConfig.from_yaml(config_path)
    config.random_seed = effective_seed(config.random_seed)
    config.scientific_snapshot = scientific_settings()
    from mofforge.utils.config import config as global_config

    root = config.moiety_path if config.moiety_path is not None else global_config.moiety_path
    config.resolved_moiety_path = str(Path(root).resolve()) if root is not None else None
    parent_paths = _resolve_parent_paths(config.parent_paths)

    if not parent_paths:
        logger.warning("No parent structures found.")
        return []

    logger.debug("Batch processing %d structures", len(parent_paths))
    results: list[BatchResult] = []

    if config.parallel > 1:
        with ProcessPoolExecutor(max_workers=config.parallel) as executor:
            futures = {
                executor.submit(_process_single, path, config): path for path in parent_paths
            }
            for future in as_completed(futures):
                try:
                    results.append(future.result())
                except Exception as exc:
                    results.append(
                        BatchResult(parent_name=futures[future].stem, success=False, error=str(exc))
                    )
    else:
        for path in parent_paths:
            results.append(_process_single(path, config))

    # Summary
    successes = sum(1 for r in results if r.success)
    failures = sum(1 for r in results if not r.success)
    logger.debug("Batch complete: %d succeeded, %d failed", successes, failures)

    return sorted(results, key=lambda r: (r.parent_name, r.run_id))
