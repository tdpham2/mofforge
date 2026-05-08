"""mofforge: Build and modify atomistic crystal structure models."""

try:
    from importlib.metadata import version as _get_version

    __version__ = _get_version("mofforge")
except Exception:
    __version__ = "0.1.0"  # fallback for editable installs without metadata

from mofforge.adsorbate import (
    AdsorbatePlacement,
    AdsorptionSite,
    available_molecules,
    find_adsorption_sites,
    get_molecule,
    place_adsorbate,
)
from mofforge.batch import run_batch
from mofforge.build import (
    BuildConfig,
    BuilderBackend,
    BuildingBlock,
    BuildResult,
    ConfigError,
    MOFBuilder,
    Topology,
)
from mofforge.core.bonding import BondingRule, infer_bonds, remove_bonds
from mofforge.core.crystal import Crystal
from mofforge.core.moiety import (
    anchor_indices,
    fragment,
    subtract_anchor,
    untag_anchor,
)
from mofforge.io.cif import read_cif, write_cif
from mofforge.io.xyz import read_xyz, write_xyz
from mofforge.pipeline import Pipeline
from mofforge.provenance import Provenance
from mofforge.replace.alignment import Alignment, apply_alignment, get_r2p_alignment
from mofforge.replace.conglomerate import reassemble
from mofforge.replace.replace import replace_pattern, swap
from mofforge.search.isomorphism import find_subgraph_isomorphisms
from mofforge.search.search import MatchResult, find_pattern
from mofforge.smarts import parse_smarts, smarts_search
from mofforge.solvent import RemovedMolecule, SolventRemovalResult, remove_solvent
from mofforge.utils.config import clean_species, config, set_paths
from mofforge.utils.periodic import min_image_distance, nearest_image, wrap_coords
from mofforge.validation import ValidationReport, validate_structure
from mofforge.vis import (
    DEFAULT_COLOR,
    JMOL_COLORS,
    METALS,
    async_render_file_to_png,
    async_render_to_png,
    build_html,
    generate_atom_labels,
    get_element_color,
    render_file_to_png,
    render_to_png,
)

__all__ = [
    "AdsorbatePlacement",
    "AdsorptionSite",
    "Alignment",
    "BondingRule",
    "BuildConfig",
    "BuildResult",
    "BuilderBackend",
    "BuildingBlock",
    "ConfigError",
    "Crystal",
    "DEFAULT_COLOR",
    "JMOL_COLORS",
    "METALS",
    "MOFBuilder",
    "MatchResult",
    "Pipeline",
    "Provenance",
    "Topology",
    "ValidationReport",
    "anchor_indices",
    "apply_alignment",
    "async_render_file_to_png",
    "async_render_to_png",
    "available_molecules",
    "build_html",
    "clean_species",
    "config",
    "find_adsorption_sites",
    "find_pattern",
    "find_subgraph_isomorphisms",
    "fragment",
    "generate_atom_labels",
    "get_element_color",
    "get_molecule",
    "get_r2p_alignment",
    "infer_bonds",
    "min_image_distance",
    "nearest_image",
    "parse_smarts",
    "place_adsorbate",
    "read_cif",
    "read_xyz",
    "reassemble",
    "RemovedMolecule",
    "remove_bonds",
    "remove_solvent",
    "render_file_to_png",
    "render_to_png",
    "replace_pattern",
    "run_batch",
    "set_paths",
    "smarts_search",
    "SolventRemovalResult",
    "subtract_anchor",
    "swap",
    "untag_anchor",
    "validate_structure",
    "wrap_coords",
    "write_cif",
    "write_xyz",
]
