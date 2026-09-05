# Visualize and export structures

Visualization is useful for inspecting replacements, unit-cell boundaries, and
guest placement. Rendering itself is not a geometry validation step.

## 1. Write interactive HTML

```bash
python examples/vis/interactive_html.py
```

[interactive_html.py](interactive_html.py) needs only core dependencies to create
three pages: `ball_stick.html`, `stick.html`, and `sphere.html`.

1. Load the IRMOF-1 CIF.
2. Call `build_html` for each representation.
3. Show unit-cell edges and rotate the camera by 15° and 30°.
4. Open a generated page in a browser to rotate and zoom interactively.

The pages load 3Dmol.js from its CDN, so **viewing requires network access** even
though generating the files is offline. They embed the structure coordinates;
they are not fully offline web bundles.

The defaults omit per-atom labels to keep a 424-atom cell readable. For a small
fragment, try `label_mode="sequential"` or `"per_element"`.
`show_formula`, dimensions, background, and sphere/stick scales are independent
display controls. Visual bonds are assigned by the renderer from geometry and
need not exactly match the library's inferred graph.

**Try:** load a replacement fragment and compare labeling styles. If a page is
blank, check CDN access and browser errors before changing the crystal data.

## 2. Render synchronous and asynchronous PNGs

Install the `vis` extra and its browser:

```bash
python -m playwright install chromium
python examples/vis/png_rendering.py
```

[png_rendering.py](png_rendering.py) writes `sync.png` and `async.png`, both
**640 × 480 pixels**. It calls `render_to_png` outside an event loop, then runs
`async_render_to_png` through `asyncio.run`. In an existing async application,
await the asynchronous function directly.

Both paths need a functioning Chromium installation and CDN access. A missing
browser, blocked launch, or network timeout is an execution failure rather than
an empty successful image. Automated live rendering is gated with
`MOFFORGE_RUN_RENDER_TESTS=1`; generating HTML is tested routinely.

**Try:** render the output of a modification lesson with labels disabled and
unit-cell edges enabled. PNGs are generated artifacts; retain source structures
and parameter summaries alongside exported figures.

Equivalent CLI, after creating an output directory with the lesson:

```bash
mofforge render --input examples/data/crystals/IRMOF-1.cif --output examples/_outputs/png_rendering/cli.png --label-mode none --show-unit-cell
```
