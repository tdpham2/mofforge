"""Render crystal structures as PNG images using 3Dmol.js and Playwright."""

from __future__ import annotations

import json
import logging
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from mofforge.utils.config import clean_species as _clean_symbol
from mofforge.vis.colors import METALS

if TYPE_CHECKING:
    from mofforge.core.crystal import Crystal

logger = logging.getLogger("mofforge")


def generate_atom_labels(
    species: list[str],
    mode: str = "sequential",
) -> list[str]:
    """Generate atom labels for a structure."""
    if mode == "none":
        return []

    clean = [_clean_symbol(s) for s in species]

    if mode == "sequential":
        return [f"{sym}{i + 1}" for i, sym in enumerate(clean)]
    elif mode == "per_element":
        counts: dict[str, int] = defaultdict(int)
        labels = []
        for sym in clean:
            counts[sym] += 1
            labels.append(f"{sym}{counts[sym]}")
        return labels
    else:
        raise ValueError(
            f"Unknown label mode: {mode!r}. Use 'sequential', 'per_element', or 'none'."
        )


def _crystal_to_xyz_string(crystal: Crystal) -> str:
    """Convert a Crystal's atoms to an XYZ-format string for 3Dmol.js."""
    species = crystal.species
    clean_species = [_clean_symbol(s) for s in species]
    coords = crystal.cart_coords

    lines = [str(len(species)), crystal.name]
    for sym, (x, y, z) in zip(clean_species, coords, strict=True):
        lines.append(f"{sym:8s} {x:14.5f} {y:14.5f} {z:14.5f}")
    return "\n".join(lines) + "\n"


def _crystal_to_cif_string(crystal: Crystal) -> str:
    """Convert a Crystal's structure to a CIF-format string for 3Dmol.js."""
    from pymatgen.io.cif import CifWriter

    writer = CifWriter(crystal.structure)
    return str(writer)


def _unit_cell_js(crystal: Crystal) -> str:
    """Generate JavaScript to draw unit cell edges in 3Dmol.js.

    Draws the 12 edges of the parallelepiped defined by the lattice
    vectors using ``viewer.addCylinder()`` calls.
    """
    lattice = crystal.lattice
    a, b, c = lattice.matrix  # 3x3, rows are lattice vectors

    # The 8 corners of the unit cell in Cartesian coordinates
    origin = np.zeros(3)
    corners = [
        origin,  # 0
        a,  # 1
        b,  # 2
        c,  # 3
        a + b,  # 4
        a + c,  # 5
        b + c,  # 6
        a + b + c,  # 7
    ]

    # The 12 edges (pairs of corner indices)
    edges = [
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 4),
        (1, 5),
        (2, 4),
        (2, 6),
        (3, 5),
        (3, 6),
        (4, 7),
        (5, 7),
        (6, 7),
    ]

    js_lines = []
    for i, j in edges:
        p1 = corners[i]
        p2 = corners[j]
        js_lines.append(
            f"viewer.addCylinder({{start: {{x:{p1[0]:.4f}, y:{p1[1]:.4f}, z:{p1[2]:.4f}}}, "
            f"end: {{x:{p2[0]:.4f}, y:{p2[1]:.4f}, z:{p2[2]:.4f}}}, "
            f"radius: 0.04, color: 'gray', fromCap: true, toCap: true}});"
        )
    return "\n        ".join(js_lines)


def build_html(
    crystal: Crystal,
    label_mode: str = "sequential",
    width: int = 800,
    height: int = 600,
    representation: str = "ball_stick",
    sphere_scale: float = 0.3,
    stick_scale: float = 0.25,
    metal_scale: float = 0.75,
    label_size: int = 14,
    show_formula: bool = True,
    show_unit_cell: bool = False,
    bg_color: str = "white",
    rotate: tuple[float, float, float] | None = None,
) -> str:
    """Build a self-contained HTML page that renders the structure with 3Dmol.js."""
    # Use XYZ format for 3Dmol.js (simpler, always works)
    xyz_string = _crystal_to_xyz_string(crystal)
    species = crystal.species
    clean_species = [_clean_symbol(s) for s in species]
    positions = crystal.cart_coords.tolist()

    # Identify metals
    metal_symbols = {sym for sym in clean_species if sym in METALS}

    # Find label anchor position (first metal, or centroid)
    metal_indices = [i for i, sym in enumerate(clean_species) if sym in METALS]
    if metal_indices:
        formula_pos = positions[metal_indices[0]]
    else:
        formula_pos = np.mean(positions, axis=0).tolist() if positions else [0, 0, 0]

    # Build atom labels
    atom_labels = generate_atom_labels(species, mode=label_mode)

    # Build the style commands as JS
    style_js_lines = []
    if representation == "ball_stick":
        style_js_lines.append(
            f"viewer.setStyle({{}}, {{stick: {{colorscheme: 'Jmol', radius: {stick_scale}}}}});"
        )
        style_js_lines.append(
            f"viewer.addStyle({{}}, {{sphere: {{colorscheme: 'Jmol', scale: {sphere_scale}}}}});"
        )
        # Scale metals larger
        for ms in metal_symbols:
            style_js_lines.append(
                f"viewer.setStyle({{elem: '{ms}'}}, "
                f"{{sphere: {{colorscheme: 'Jmol', scale: {metal_scale}}}, "
                f"stick: {{colorscheme: 'Jmol', radius: {stick_scale}}}}});"
            )
    elif representation == "stick":
        style_js_lines.append(
            f"viewer.setStyle({{}}, {{stick: {{colorscheme: 'Jmol', radius: {stick_scale}}}}});"
        )
    elif representation == "sphere":
        style_js_lines.append(
            f"viewer.setStyle({{}}, {{sphere: {{colorscheme: 'Jmol', scale: {sphere_scale}}}}});"
        )
    else:
        style_js_lines.append(
            f"viewer.setStyle({{}}, {{{representation}: {{colorscheme: 'Jmol'}}}});"
        )

    style_js = "\n        ".join(style_js_lines)

    # Build label JS
    label_js_lines = []

    # Atom index labels
    for i, lbl in enumerate(atom_labels):
        x, y, z = positions[i]
        label_js_lines.append(
            f'viewer.addLabel("{lbl}", {{'
            f"position: {{x: {x:.4f}, y: {y:.4f}, z: {z:.4f}}}, "
            f"backgroundColor: 'black', backgroundOpacity: 0.4, "
            f"fontOpacity: 1, fontSize: {label_size}, "
            f"fontColor: 'white', inFront: true"
            f"}});"
        )

    # Formula label
    if show_formula:
        formula = crystal.structure.composition.reduced_formula
        fx, fy, fz = formula_pos
        label_js_lines.append(
            f'viewer.addLabel("{formula}", {{'
            f"position: {{x: {fx:.4f}, y: {fy:.4f}, z: {fz:.4f}}}, "
            f"backgroundColor: 'black', backgroundOpacity: 0.5, "
            f"fontOpacity: 1, fontSize: {label_size + 2}, "
            f"fontColor: 'yellow', inFront: true"
            f"}});"
        )

    label_js = "\n        ".join(label_js_lines)

    # Escape the XYZ string for JS embedding
    xyz_js = json.dumps(xyz_string)

    # Unit cell JS (optional)
    unit_cell_js = ""
    if show_unit_cell:
        unit_cell_js = _unit_cell_js(crystal)

    # Build rotation JS
    rotate_js_lines = []
    if rotate is not None:
        rx, ry, rz = rotate
        if rx != 0:
            rotate_js_lines.append(f"viewer.rotate({rx}, 'x');")
        if ry != 0:
            rotate_js_lines.append(f"viewer.rotate({ry}, 'y');")
        if rz != 0:
            rotate_js_lines.append(f"viewer.rotate({rz}, 'z');")
    rotate_js = "\n        ".join(rotate_js_lines)

    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <style>
        * {{ margin: 0; padding: 0; }}
        #viewer {{ width: {width}px; height: {height}px; position: relative; }}
    </style>
</head>
<body>
    <div id="viewer"></div>
    <script>
        var viewer = $3Dmol.createViewer("viewer", {{
            backgroundColor: "{bg_color}"
        }});

        var xyzData = {xyz_js};
        viewer.addModel(xyzData, "xyz");

        {style_js}

        {label_js}

        {unit_cell_js}

        viewer.zoomTo();
        {rotate_js}
        viewer.render();

        // Signal that rendering is complete
        window.__3dmol_ready = true;
    </script>
</body>
</html>"""

    return html


def render_to_png(
    crystal: Crystal,
    output_file: str = "structure.png",
    label_mode: str = "sequential",
    width: int = 800,
    height: int = 600,
    representation: str = "ball_stick",
    sphere_scale: float = 0.3,
    stick_scale: float = 0.25,
    metal_scale: float = 0.75,
    label_size: int = 14,
    show_formula: bool = True,
    show_unit_cell: bool = False,
    bg_color: str = "white",
    timeout: int = 30000,
    rotate: tuple[float, float, float] | None = None,
) -> str:
    """Render a crystal structure to a PNG file."""
    try:
        from playwright.sync_api import sync_playwright
    except ImportError as exc:
        raise ImportError("playwright is required for PNG rendering") from exc

    html = build_html(
        crystal,
        label_mode=label_mode,
        width=width,
        height=height,
        representation=representation,
        sphere_scale=sphere_scale,
        stick_scale=stick_scale,
        metal_scale=metal_scale,
        label_size=label_size,
        show_formula=show_formula,
        show_unit_cell=show_unit_cell,
        bg_color=bg_color,
        rotate=rotate,
    )

    # Write HTML to a temp file
    tmp_dir = tempfile.mkdtemp(prefix="mofforge_render_")
    html_path = os.path.join(tmp_dir, "structure.html")
    with open(html_path, "w") as f:
        f.write(html)

    output_path = os.path.abspath(output_file)

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": width, "height": height})

        page.goto(f"file://{html_path}")
        page.wait_for_function("window.__3dmol_ready === true", timeout=timeout)

        # Small extra delay for WebGL to fully flush
        page.wait_for_timeout(500)

        # Screenshot the viewer div
        viewer_element = page.locator("#viewer")
        viewer_element.screenshot(path=output_path)

        browser.close()

    # Clean up temp file
    try:
        os.remove(html_path)
        os.rmdir(tmp_dir)
    except OSError:
        pass

    logger.debug("Rendered structure to %s", output_path)
    return output_path


async def async_render_to_png(
    crystal: Crystal,
    output_file: str = "structure.png",
    label_mode: str = "sequential",
    width: int = 800,
    height: int = 600,
    representation: str = "ball_stick",
    sphere_scale: float = 0.3,
    stick_scale: float = 0.25,
    metal_scale: float = 0.75,
    label_size: int = 14,
    show_formula: bool = True,
    show_unit_cell: bool = False,
    bg_color: str = "white",
    timeout: int = 30000,
    rotate: tuple[float, float, float] | None = None,
) -> str:
    """Async version of render_to_png."""
    try:
        from playwright.async_api import async_playwright
    except ImportError as exc:
        raise ImportError("playwright is required for PNG rendering") from exc

    html = build_html(
        crystal,
        label_mode=label_mode,
        width=width,
        height=height,
        representation=representation,
        sphere_scale=sphere_scale,
        stick_scale=stick_scale,
        metal_scale=metal_scale,
        label_size=label_size,
        show_formula=show_formula,
        show_unit_cell=show_unit_cell,
        bg_color=bg_color,
        rotate=rotate,
    )

    tmp_dir = tempfile.mkdtemp(prefix="mofforge_render_")
    html_path = os.path.join(tmp_dir, "structure.html")
    with open(html_path, "w") as f:
        f.write(html)

    output_path = os.path.abspath(output_file)

    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page(viewport={"width": width, "height": height})

        await page.goto(f"file://{html_path}")
        await page.wait_for_function("window.__3dmol_ready === true", timeout=timeout)
        await page.wait_for_timeout(500)

        viewer_element = page.locator("#viewer")
        await viewer_element.screenshot(path=output_path)

        await browser.close()

    try:
        os.remove(html_path)
        os.rmdir(tmp_dir)
    except OSError:
        pass

    logger.debug("Rendered structure to %s (async)", output_path)
    return output_path


def render_file_to_png(
    input_file: str,
    output_file: str = "structure.png",
    label_mode: str = "sequential",
    width: int = 800,
    height: int = 600,
    representation: str = "ball_stick",
    sphere_scale: float = 0.3,
    stick_scale: float = 0.25,
    metal_scale: float = 0.75,
    label_size: int = 14,
    show_formula: bool = True,
    show_unit_cell: bool | None = None,
    bg_color: str = "white",
    timeout: int = 30000,
    rotate: tuple[float, float, float] | None = None,
) -> str:
    """Read a structure file and render it to PNG.

    Accepts CIF and XYZ files.
    """
    crystal = _load_file(input_file)

    # Auto-enable unit cell for CIF files (periodic structures) unless
    # the caller explicitly passed show_unit_cell=False.
    ext = Path(input_file).suffix.lower()
    if ext == ".cif" and show_unit_cell is None:
        show_unit_cell = True
    elif show_unit_cell is None:
        show_unit_cell = False

    return render_to_png(
        crystal,
        output_file=output_file,
        label_mode=label_mode,
        width=width,
        height=height,
        representation=representation,
        sphere_scale=sphere_scale,
        stick_scale=stick_scale,
        metal_scale=metal_scale,
        label_size=label_size,
        show_formula=show_formula,
        show_unit_cell=show_unit_cell,
        bg_color=bg_color,
        timeout=timeout,
        rotate=rotate,
    )


async def async_render_file_to_png(
    input_file: str,
    output_file: str = "structure.png",
    label_mode: str = "sequential",
    width: int = 800,
    height: int = 600,
    representation: str = "ball_stick",
    sphere_scale: float = 0.3,
    stick_scale: float = 0.25,
    metal_scale: float = 0.75,
    label_size: int = 14,
    show_formula: bool = True,
    show_unit_cell: bool | None = None,
    bg_color: str = "white",
    timeout: int = 30000,
    rotate: tuple[float, float, float] | None = None,
) -> str:
    """Async version of render_file_to_png."""
    crystal = _load_file(input_file)

    # Auto-enable unit cell for CIF files (periodic structures) unless
    # the caller explicitly passed show_unit_cell=False.
    ext = Path(input_file).suffix.lower()
    if ext == ".cif" and show_unit_cell is None:
        show_unit_cell = True
    elif show_unit_cell is None:
        show_unit_cell = False

    return await async_render_to_png(
        crystal,
        output_file=output_file,
        label_mode=label_mode,
        width=width,
        height=height,
        representation=representation,
        sphere_scale=sphere_scale,
        stick_scale=stick_scale,
        metal_scale=metal_scale,
        label_size=label_size,
        show_formula=show_formula,
        show_unit_cell=show_unit_cell,
        bg_color=bg_color,
        timeout=timeout,
        rotate=rotate,
    )


def _load_file(filepath: str) -> Crystal:
    """Load a CIF or XYZ file into a Crystal object."""
    from mofforge.core.crystal import Crystal
    from mofforge.io.xyz import read_xyz

    ext = Path(filepath).suffix.lower()

    if ext == ".cif":
        return Crystal.from_cif(filepath)
    elif ext == ".xyz":
        species, coords = read_xyz(filepath)
        return Crystal.from_xyz(species, coords, name=Path(filepath).stem)
    else:
        raise ValueError(f"Unsupported file format: {ext!r}")
