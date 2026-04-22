"""Visualization subsystem for mofforge.

Renders crystal structures as PNG images using 3Dmol.js and Playwright.
Supports both periodic structures (MOFs with unit cell) and non-periodic
fragments.

Quick start::

    from mofforge import Crystal, render_to_png

    crystal = Crystal.from_cif("IRMOF-1.cif")
    render_to_png(crystal, "IRMOF-1.png", show_unit_cell=True)

    # Or from a file directly:
    from mofforge.vis import render_file_to_png
    render_file_to_png("IRMOF-1.cif", "IRMOF-1.png")
"""

from mofforge.vis.colors import (
    DEFAULT_COLOR,
    JMOL_COLORS,
    METALS,
    get_element_color,
)
from mofforge.vis.render import (
    async_render_file_to_png,
    async_render_to_png,
    build_html,
    generate_atom_labels,
    render_file_to_png,
    render_to_png,
)

__all__ = [
    "DEFAULT_COLOR",
    "JMOL_COLORS",
    "METALS",
    "async_render_file_to_png",
    "async_render_to_png",
    "build_html",
    "generate_atom_labels",
    "get_element_color",
    "render_file_to_png",
    "render_to_png",
]
