"""Batch processing for multiple crystal structures.

Processes multiple parent structures with configurable operations
from a YAML configuration file, with optional parallel execution.
"""

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
from mofforge.core.moiety import fragment
from mofforge.replace.replace import replace_pattern
from mofforge.search.search import find_pattern
from mofforge.validation import ValidationReport, validate_structure

logger = logging.getLogger("mofforge")


@dataclass
class BatchResult:
    """Result of processing a single structure in a batch.

    Attributes:
        parent_name: Name of the parent structure.
        output_path: Path to the output file (if written).
        success: Whether the operation succeeded.
        error: Error message if failed.
        validation: Validation report if validation was performed.
    """

    parent_name: str
    output_path: str | None = None
    success: bool = True
    error: str | None = None
    validation: ValidationReport | None = None


@dataclass
class BatchConfig:
    """Configuration for a batch processing run.

    Attributes:
        parent_paths: List of glob patterns or file paths for parent structures.
        operations: List of operation dicts.
        output_dir: Output directory for results.
        output_format: Output file format ('cif' or 'xyz').
        naming: Output naming template (e.g. '{parent_name}_modified').
        parallel: Number of parallel workers (0 = sequential).
        moiety_path: Path to fragment XYZ files.
    """

    parent_paths: list[str] = field(default_factory=list)
    operations: list[dict] = field(default_factory=list)
    output_dir: str = "results"
    output_format: str = "cif"
    naming: str = "{parent_name}_modified"
    parallel: int = 0
    moiety_path: str | None = None

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

        Args:
            filepath: Path to the YAML config file.

        Returns:
            BatchConfig object.
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
    return paths


def _process_single(
    parent_path: Path,
    config: BatchConfig,
) -> BatchResult:
    """Process a single parent structure through all operations."""
    parent_name = parent_path.stem
    result = BatchResult(parent_name=parent_name)

    try:
        current = Crystal.from_cif(parent_path)
        current = infer_bonds(current, periodic=True)

        for op in config.operations:
            op_type = op.get("type", "")

            _VALID_OP_TYPES = ("replace", "remove", "validate")
            if op_type not in _VALID_OP_TYPES:
                raise ValueError(
                    f"Unknown operation type '{op_type}' in batch config. "
                    f"Must be one of: {', '.join(_VALID_OP_TYPES)}."
                )

            if op_type == "replace":
                query_name = op.get("query")
                replacement_name = op.get("replacement")
                mode = op.get("mode", "all_optimal")

                q = fragment(query_name, fragment_path=config.moiety_path)
                r = fragment(replacement_name, fragment_path=config.moiety_path)

                match = find_pattern(q, current)
                kwargs = {}
                if mode == "random":
                    kwargs["random"] = True
                elif mode.startswith("nb_loc_"):
                    kwargs["nb_loc"] = int(mode.split("_")[-1])

                current = replace_pattern(match, r, **kwargs)

                # Re-infer bonds for subsequent steps
                if current.n_bonds == 0 and current.n_atoms > 0:
                    current = infer_bonds(current, periodic=True)

            elif op_type == "remove":
                guest_name = op.get("guest") or op.get("query")
                g = fragment(guest_name, fragment_path=config.moiety_path)
                match = find_pattern(g, current, disconnected_component=True)
                current = replace_pattern(match, None)

                if current.n_bonds == 0 and current.n_atoms > 0:
                    current = infer_bonds(current, periodic=True)

            elif op_type == "validate":
                report = validate_structure(current)
                result.validation = report

        # Write output
        output_dir = Path(config.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        output_name = config.naming.format(parent_name=parent_name)
        output_path = output_dir / f"{output_name}.{config.output_format}"
        if config.output_format == "xyz":
            current.write_xyz(output_path)
        else:
            current.write_cif(output_path)
        result.output_path = str(output_path)

        logger.info("Processed '%s' -> '%s'", parent_name, output_path)

    except Exception as e:
        result.success = False
        result.error = str(e)
        logger.error("Failed to process '%s': %s", parent_name, e)

    return result


def run_batch(config_path: str | Path) -> list[BatchResult]:
    """Run batch processing from a YAML configuration file.

    Args:
        config_path: Path to the YAML config file.

    Returns:
        List of BatchResult objects, one per parent structure.
    """
    config = BatchConfig.from_yaml(config_path)
    parent_paths = _resolve_parent_paths(config.parent_paths)

    if not parent_paths:
        logger.warning("No parent structures found.")
        return []

    logger.info("Batch processing %d structures", len(parent_paths))
    results: list[BatchResult] = []

    if config.parallel > 1:
        with ProcessPoolExecutor(max_workers=config.parallel) as executor:
            futures = {
                executor.submit(_process_single, path, config): path for path in parent_paths
            }
            for future in as_completed(futures):
                results.append(future.result())
    else:
        for path in parent_paths:
            results.append(_process_single(path, config))

    # Summary
    successes = sum(1 for r in results if r.success)
    failures = sum(1 for r in results if not r.success)
    logger.info("Batch complete: %d succeeded, %d failed", successes, failures)

    return results
