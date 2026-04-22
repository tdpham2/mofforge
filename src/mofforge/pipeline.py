"""Multi-step replacement pipeline API.

Provides a fluent interface for chaining multiple find-and-replace
operations on a crystal structure.  Can also start from scratch by
building a MOF with :meth:`Pipeline.build_mof`.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from mofforge.core.bonding import infer_bonds
from mofforge.core.crystal import Crystal
from mofforge.core.moiety import fragment
from mofforge.provenance import Provenance
from mofforge.replace.replace import replace_pattern
from mofforge.search.search import find_pattern
from mofforge.validation import ValidationReport, validate_structure

logger = logging.getLogger("mofforge")


@dataclass
class PipelineStep:
    """A single operation in a pipeline.

    Attributes:
        operation: Type of operation ('replace', 'remove', 'validate').
        kwargs: Keyword arguments for the operation.
    """

    operation: str
    kwargs: dict[str, Any] = field(default_factory=dict)


class Pipeline:
    """Fluent API for chaining crystal modification operations.

    Example::

        child = (Pipeline("parent.cif")
            .replace(query="BDC.xyz", replacement="NH2-BDC.xyz", nb_loc=6)
            .remove(guest="acetylene.xyz")
            .replace(query="OH.xyz", replacement="F.xyz")
            .build(name="functionalized_MOF"))
    """

    def __init__(
        self,
        parent: Crystal | str | Path,
        fragment_path: str | Path | None = None,
    ):
        """Initialize a pipeline with a parent crystal.

        Args:
            parent: A Crystal object or path to a CIF file.
            fragment_path: Optional path to directory containing fragment XYZ files.
        """
        if isinstance(parent, (str, Path)):
            parent = Crystal.from_cif(parent)
        self._parent = parent
        self._moiety_path = fragment_path
        self._steps: list[PipelineStep] = []
        self._intermediates: list[Crystal] = []
        self._reports: list[ValidationReport] = []

    def replace(
        self,
        query: str,
        replacement: str,
        **kwargs,
    ) -> Pipeline:
        """Queue a find-and-replace operation.

        Args:
            query: Query fragment XYZ filename.
            replacement: Replacement fragment XYZ filename.
            **kwargs: Additional arguments for replace_pattern.

        Returns:
            self (for chaining).
        """
        self._steps.append(
            PipelineStep(
                operation="replace",
                kwargs={"query": query, "replacement": replacement, **kwargs},
            )
        )
        return self

    def remove(self, guest: str, **kwargs) -> Pipeline:
        """Queue a guest removal operation (replace with nothing).

        Args:
            guest: Guest fragment XYZ filename.
            **kwargs: Additional arguments for replace_pattern.

        Returns:
            self (for chaining).
        """
        kwargs.setdefault("disconnected_component", True)
        self._steps.append(
            PipelineStep(
                operation="remove",
                kwargs={"guest": guest, **kwargs},
            )
        )
        return self

    def validate(self, **kwargs) -> Pipeline:
        """Queue a validation step.

        Args:
            **kwargs: Arguments for validate_structure.

        Returns:
            self (for chaining).
        """
        self._steps.append(
            PipelineStep(
                operation="validate",
                kwargs=kwargs,
            )
        )
        return self

    def build(self, name: str = "new_xtal") -> Crystal:
        """Execute all queued operations and return the final crystal.

        Args:
            name: Name for the final crystal.

        Returns:
            The modified Crystal.
        """
        # Reset state so repeated calls don't accumulate duplicates
        self._intermediates = []
        self._reports = []

        current = self._parent.copy()

        # Ensure parent has bonds
        if current.n_bonds == 0 and current.n_atoms > 0:
            current = infer_bonds(current, periodic=True)

        provenance_chain = Provenance(
            parent=current.name,
            operation="pipeline_start",
        )

        for i, step in enumerate(self._steps):
            logger.info("Pipeline step %d/%d: %s", i + 1, len(self._steps), step.operation)

            if step.operation == "replace":
                query_name = step.kwargs["query"]
                replacement_name = step.kwargs["replacement"]
                extra_kwargs = {
                    k: v for k, v in step.kwargs.items() if k not in ("query", "replacement")
                }
                q = fragment(query_name, fragment_path=self._moiety_path)
                r = fragment(replacement_name, fragment_path=self._moiety_path)

                match = find_pattern(q, current)
                step_name = f"step_{i + 1}_{name}"
                current = replace_pattern(match, r, name=step_name, **extra_kwargs)

                # Ensure bonds are inferred for the next step
                if current.n_bonds == 0 and current.n_atoms > 0:
                    current = infer_bonds(current, periodic=True)

                provenance_chain = provenance_chain.chain(
                    Provenance(
                        parent=current.name,
                        query=query_name,
                        replacement=replacement_name,
                        operation="replace",
                        parameters=step.kwargs,
                    )
                )

            elif step.operation == "remove":
                guest_name = step.kwargs["guest"]
                disconnected = step.kwargs.get("disconnected_component", True)
                extra_kwargs = {
                    k: v
                    for k, v in step.kwargs.items()
                    if k not in ("guest", "disconnected_component")
                }
                g = fragment(guest_name, fragment_path=self._moiety_path)

                match = find_pattern(g, current, disconnected_component=disconnected)
                step_name = f"step_{i + 1}_{name}"
                current = replace_pattern(match, None, name=step_name, **extra_kwargs)

                if current.n_bonds == 0 and current.n_atoms > 0:
                    current = infer_bonds(current, periodic=True)

                provenance_chain = provenance_chain.chain(
                    Provenance(
                        parent=current.name,
                        query=guest_name,
                        replacement=None,
                        operation="remove",
                    )
                )

            elif step.operation == "validate":
                report = validate_structure(current, **step.kwargs)
                self._reports.append(report)
                logger.info("Validation: %s", "PASSED" if report.is_valid else "ISSUES FOUND")

            self._intermediates.append(current.copy())

        current.name = name
        current.provenance = provenance_chain
        return current

    def build_all(self, name: str = "new_xtal") -> list[Crystal]:
        """Execute all operations and return all intermediate crystals.

        The last element is the final crystal (same object returned by
        ``build()``).

        Args:
            name: Name for the final crystal.

        Returns:
            List of Crystal objects (one per step).  The last element is
            the final named crystal.
        """
        final = self.build(name)
        # _intermediates already contains a copy after each step; replace the
        # last entry with the properly-named final crystal to avoid duplication.
        if self._intermediates:
            self._intermediates[-1] = final
            return list(self._intermediates)
        return [final]

    @property
    def validation_reports(self) -> list[ValidationReport]:
        """Access validation reports from completed pipeline runs."""
        return self._reports

    # ------------------------------------------------------------------ #
    # Factory: build a MOF from scratch, then modify
    # ------------------------------------------------------------------ #

    @classmethod
    def build_mof(
        cls,
        backend: str = "tobacco",
        topology: str = "pcu",
        nodes: list[str | Path] | None = None,
        edges: list[str | Path] | None = None,
        output_dir: str | Path = ".",
        fragment_path: str | Path | None = None,
        **backend_kwargs: Any,
    ) -> Pipeline:
        """Build a MOF from scratch and return a pipeline for further modifications.

        Combines construction (via TOBACCO or Pormake) with the
        existing find-and-replace pipeline in a single fluent chain.

        Example::

            child = (Pipeline.build_mof(
                        backend="tobacco", topology="pcu",
                        tobacco_path="/path/to/tobacco_3.0")
                .replace(query="BDC.xyz", replacement="NH2-BDC.xyz")
                .validate()
                .build(name="functionalized_pcu_MOF"))

        Args:
            backend: ``"tobacco"`` or ``"pormake"``.
            topology: Topology name (e.g. ``"pcu"`` or ``"pcu.cif"``).
            nodes: Paths to node building-block files to register
                before building.
            edges: Paths to edge building-block files to register
                before building.
            output_dir: Where to write generated CIF files.
            fragment_path: Optional path to directory containing
                fragment XYZ files (for downstream replace/remove steps).
            **backend_kwargs: Passed to :class:`~mofforge.build.MOFBuilder`
                constructor (e.g. ``tobacco_path="/path/to/tobacco_3.0"``).

        Returns:
            A :class:`Pipeline` initialised with the built crystal,
            ready for chaining ``.replace()`` / ``.remove()`` /
            ``.validate()`` / ``.build()`` calls.

        Raises:
            RuntimeError: If the MOF build fails.
        """
        from mofforge.build import MOFBuilder

        builder = MOFBuilder(backend=backend, **backend_kwargs)

        if nodes:
            for n in nodes:
                builder.add_node(n)
        if edges:
            for e in edges:
                builder.add_edge(e)

        result = builder.build(topology=topology, output_dir=output_dir)

        if not result.success or result.crystal is None:
            error_msg = "; ".join(result.errors) if result.errors else "unknown error"
            raise RuntimeError(f"MOF build failed: {error_msg}")

        logger.info(
            "MOF built with %s backend (topology=%s, %d output files)",
            backend,
            topology,
            len(result.output_paths),
        )

        return cls(parent=result.crystal, fragment_path=fragment_path)
