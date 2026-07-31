"""High-level functionalization: single runs and full campaigns.

These entry points wire the functionalization-specific pieces (group menu, site
detection, fragment generation) to mofforge's existing find/replace/validate
pipeline.  An AI agent drives everything through *external factors only*:

* which functional group (from :func:`mofforge.functionalize.groups.available_groups`),
* which site indices (from :func:`mofforge.functionalize.sites.find_functionalizable_sites`),
* the framework-wide coverage / concentration.

Coverage maps onto the existing ``nb_loc`` knob of
:func:`mofforge.replace.replace.replace_pattern`: a coverage of ``0.5`` on a
framework with 12 matched linkers functionalizes 6 of them (rounded).
"""

from __future__ import annotations

import logging
import tempfile
from dataclasses import dataclass
from itertools import product
from pathlib import Path

from mofforge.core.bonding import infer_bonds
from mofforge.core.crystal import Crystal
from mofforge.core.moiety import fragment as load_fragment
from mofforge.functionalize.generate import make_query_replacement
from mofforge.replace.replace import replace_pattern
from mofforge.search.search import find_pattern
from mofforge.validation import validate_structure

logger = logging.getLogger("mofforge")


@dataclass
class FunctionalizationResult:
    """Outcome of a single functionalization run.

    Attributes
    ----------
    group:
        Functional-group name applied.
    sites:
        Site indices functionalized on each modified linker.
    coverage:
        Requested fraction of matched linkers to functionalize (0-1).
    n_matches:
        Total matched linker locations found in the framework.
    n_functionalized:
        Number of linkers actually functionalized (``nb_loc``).
    output_cif:
        Path to the written structure, or ``None`` if not written.
    crystal:
        The resulting :class:`~mofforge.core.crystal.Crystal`.
    is_valid:
        Validation verdict (``None`` if validation was skipped).
    clashes:
        Number of steric clashes reported by validation (``None`` if skipped).
    validation_summary:
        Human-readable validation summary (``None`` if skipped).
    error:
        Error message if the run failed, else ``None``.
    """

    group: str
    sites: list[int]
    coverage: float
    n_matches: int = 0
    n_functionalized: int = 0
    output_cif: str | None = None
    crystal: Crystal | None = None
    is_valid: bool | None = None
    clashes: int | None = None
    validation_summary: str | None = None
    error: str | None = None


def _coverage_to_nb_loc(coverage: float, n_matches: int) -> int:
    """Translate a coverage fraction into a number of locations to replace."""
    if coverage >= 1.0:
        return 0  # 0 == "all locations" in replace_pattern
    if coverage <= 0.0:
        return 0
    return max(1, round(coverage * n_matches))


def functionalize(
    parent_cif: str,
    linker_smiles: str,
    group: str,
    sites: int | list[int] = 0,
    coverage: float = 1.0,
    output_cif: str | None = None,
    validate: bool = True,
    random_seed: int | None = None,
    name: str | None = None,
) -> FunctionalizationResult:
    """Functionalize a MOF linker with a chosen group at chosen sites.

    Parameters
    ----------
    parent_cif:
        Path to the MOF CIF to functionalize.
    linker_smiles:
        SMILES of the linker being modified (e.g. from MOFid).
    group:
        Functional-group name from the curated menu.
    sites:
        Site index or list of indices (same aromatic ring system) to
        functionalize on each linker.  Defaults to site ``0``.
    coverage:
        Fraction (0-1) of matched linkers to functionalize.  ``1.0`` = all.
    output_cif:
        Where to write the result.  If ``None``, no file is written.
    validate:
        Run structure validation on the result.
    random_seed:
        Seed for the random location selection when ``coverage < 1``.
    name:
        Name for the resulting crystal.

    Returns
    -------
    FunctionalizationResult
    """
    if isinstance(sites, int):
        sites = [sites]
    result = FunctionalizationResult(group=group, sites=list(sites), coverage=coverage)

    try:
        with tempfile.TemporaryDirectory() as tmp:
            query_path, replacement_path = make_query_replacement(
                linker_smiles, sites, group, output_dir=tmp
            )
            query = load_fragment(Path(query_path).name, fragment_path=tmp)
            replacement = load_fragment(Path(replacement_path).name, fragment_path=tmp)

            parent = Crystal.from_cif(parent_cif)
            parent = infer_bonds(parent, periodic=True)

            match = find_pattern(query, parent)
            result.n_matches = match.nb_locations()
            if result.n_matches == 0:
                result.error = (
                    "Query fragment not found in the framework. The linker SMILES "
                    "may not match this MOF, or its hydrogens may be absent in the CIF."
                )
                return result

            if random_seed is not None:
                import random as pyrandom

                pyrandom.seed(random_seed)

            nb_loc = _coverage_to_nb_loc(coverage, result.n_matches)
            child = replace_pattern(
                match,
                replacement,
                nb_loc=nb_loc,
                name=name or f"{group}_functionalized",
            )
            result.n_functionalized = nb_loc if nb_loc > 0 else result.n_matches
            result.crystal = child

            if output_cif is not None:
                child.write_cif(output_cif)
                result.output_cif = output_cif

            if validate:
                child_bonded = infer_bonds(child, periodic=True)
                report = validate_structure(child_bonded)
                result.is_valid = report.is_valid
                result.clashes = len(report.steric_clashes)
                result.validation_summary = report.summary()

    except Exception as exc:
        logger.warning("functionalize failed", exc_info=True)
        result.error = str(exc)

    return result


def run_campaign(
    parent_cif: str,
    linker_smiles: str,
    groups: list[str],
    coverages: list[float] | None = None,
    sites: int | list[int] = 0,
    output_dir: str | None = None,
    validate: bool = True,
    random_seed: int | None = None,
) -> list[FunctionalizationResult]:
    """Sweep groups x coverages, functionalize each, and rank the results.

    Ranking puts valid structures first, then orders by ascending steric-clash
    count so the agent sees the most promising functionalized structures on top.

    Parameters
    ----------
    parent_cif:
        Path to the MOF CIF to functionalize.
    linker_smiles:
        SMILES of the linker being modified.
    groups:
        Functional-group names to sweep.
    coverages:
        Coverage fractions to sweep.  Defaults to ``[0.25, 0.5, 1.0]``.
    sites:
        Site index or indices to functionalize (same for every combination).
    output_dir:
        Directory for the per-combination CIFs.  If ``None``, structures are
        computed but not written to disk.
    validate:
        Run validation on each result (needed for meaningful ranking).
    random_seed:
        Seed for reproducible location selection.

    Returns
    -------
    list[FunctionalizationResult]
        Ranked best-first.
    """
    if coverages is None:
        coverages = [0.25, 0.5, 1.0]

    out_dir = Path(output_dir) if output_dir is not None else None
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)

    site_list = [sites] if isinstance(sites, int) else list(sites)
    site_tag = "-".join(str(i) for i in site_list)

    results: list[FunctionalizationResult] = []
    for group, coverage in product(groups, coverages):
        cov_tag = f"{round(coverage * 100)}pct"
        output_cif = (
            str(out_dir / f"{group}_site{site_tag}_{cov_tag}.cif")
            if out_dir is not None
            else None
        )
        res = functionalize(
            parent_cif,
            linker_smiles,
            group,
            sites=site_list,
            coverage=coverage,
            output_cif=output_cif,
            validate=validate,
            random_seed=random_seed,
        )
        results.append(res)

    def rank_key(r: FunctionalizationResult) -> tuple:
        failed = r.error is not None
        invalid = not bool(r.is_valid) if r.is_valid is not None else True
        clashes = r.clashes if r.clashes is not None else float("inf")
        return (failed, invalid, clashes)

    results.sort(key=rank_key)
    return results
