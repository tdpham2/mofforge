"""Agent-driven MOF linker functionalization.

This subpackage lets an AI agent decorate a MOF's organic linker with a chosen
functional group without ever authoring 3-D geometry or SMILES:

* :mod:`~mofforge.functionalize.groups` — a curated menu of substituents
  (``NH2``, ``NO2``, ``F`` …), defined as SMILES attachments, not coordinates.
* :mod:`~mofforge.functionalize.sites` — deterministic detection of the
  functionalizable aromatic C-H positions on a linker, each tagged with a
  symmetry class.  This encodes the pre-defined "best position": the agent
  picks *indices*, never coordinates.
* :mod:`~mofforge.functionalize.generate` — RDKit-driven generation of an
  anchor-tagged query/replacement fragment pair for a chosen site + group.
* :mod:`~mofforge.functionalize.campaign` — the high-level
  :func:`functionalize` and :func:`run_campaign` entry points that drive the
  existing find/replace/validate pipeline.
"""

from __future__ import annotations

from mofforge.functionalize.campaign import (
    FunctionalizationResult,
    functionalize,
    run_campaign,
)
from mofforge.functionalize.generate import make_query_replacement
from mofforge.functionalize.groups import (
    available_groups,
    group_smiles,
)
from mofforge.functionalize.sites import (
    FunctionalizableSite,
    find_functionalizable_sites,
)

__all__ = [
    "FunctionalizableSite",
    "FunctionalizationResult",
    "available_groups",
    "find_functionalizable_sites",
    "functionalize",
    "group_smiles",
    "make_query_replacement",
    "run_campaign",
]
