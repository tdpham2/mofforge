"""mofforge: Find-and-replace tool for atomistic crystal structure models.

A Python package for searching and modifying crystal structures,
especially metal-organic frameworks (MOFs). Provides pattern
matching via VF2 graph isomorphism and fragment replacement via
SVD-based Procrustes alignment.

Features include VF2 algorithm, SMARTS-like pattern matching,
batch processing, multi-step pipelines, structure validation,
and provenance tracking.
"""

try:
    from importlib.metadata import version as _get_version

    __version__ = _get_version("mofforge")
except Exception:
    __version__ = "0.1.0"  # fallback for editable installs without metadata

# Core data structures
from mofforge.batch import run_batch
from mofforge.core.bonding import BondingRule, infer_bonds, remove_bonds
from mofforge.core.crystal import Crystal
from mofforge.core.moiety import (
    anchor_indices,
    fragment,
    subtract_anchor,
    untag_anchor,
)

# I/O
from mofforge.io.cif import read_cif, write_cif
from mofforge.io.xyz import read_xyz, write_xyz
from mofforge.pipeline import Pipeline
from mofforge.provenance import Provenance
from mofforge.replace.alignment import Alignment
from mofforge.replace.conglomerate import reassemble

# Replace
from mofforge.replace.replace import replace_pattern, swap

# Search
from mofforge.search.search import MatchResult, find_pattern
from mofforge.smarts import parse_smarts, smarts_search

# Configuration
from mofforge.utils.config import config, set_paths

# Validation
from mofforge.validation import ValidationReport, validate_structure

__all__ = [
    "Alignment",
    "BondingRule",
    "Crystal",
    "MatchResult",
    "Pipeline",
    "Provenance",
    "ValidationReport",
    "anchor_indices",
    "config",
    "find_pattern",
    "fragment",
    "infer_bonds",
    "parse_smarts",
    "read_cif",
    "read_xyz",
    "reassemble",
    "remove_bonds",
    "replace_pattern",
    "run_batch",
    "set_paths",
    "smarts_search",
    "subtract_anchor",
    "swap",
    "untag_anchor",
    "validate_structure",
    "write_cif",
    "write_xyz",
]
