"""Cookbook helpers for paths, prerequisites, and report serialization."""

from __future__ import annotations

import argparse
import importlib.util
import json
import platform
import shutil
from contextlib import suppress
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

EXAMPLES = Path(__file__).resolve().parent
CRYSTALS = EXAMPLES / "data" / "crystals"
MOIETIES = EXAMPLES / "data" / "moieties"
BDC_SMILES = "O=C(O)c1ccc(C(=O)O)cc1"


def example_parser(name, description, *, seed=False):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--output-dir",
        "-o",
        type=Path,
        default=EXAMPLES / "_outputs" / name,
        help="Directory for generated artifacts and summary.json.",
    )
    if seed:
        parser.add_argument("--seed", "-s", type=int, default=42, help="Random seed (default: 42).")
    return parser


def output_dir(args):
    path = args.output_dir.expanduser().resolve()
    path.mkdir(parents=True, exist_ok=True)
    return path


def output_file(args, default):
    path = Path(args.output).expanduser().resolve() if args.output else output_dir(args) / default
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def require(module, extra):
    if importlib.util.find_spec(module) is None:
        raise SystemExit(
            f"This lesson requires {module}. Install with: pip install 'mofforge[{extra}]'"
        )


def json_default(value):
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "tolist"):
        return value.tolist()
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def write_report(directory, **payload):
    versions = {"python": platform.python_version()}
    for package in (
        "mofforge",
        "pymatgen",
        "numpy",
        "networkx",
        "scipy",
        "rdkit",
        "pormake",
        "tobacco3",
        "mcp",
        "playwright",
        "chemgraph",
        "langchain-mcp-adapters",
        "langchain-core",
    ):
        with suppress(PackageNotFoundError):
            versions[package] = version(package)
    text = json.dumps(
        {"environment": versions, **payload}, indent=2, default=json_default, allow_nan=False
    )
    (Path(directory) / "summary.json").write_text(text + "\n", encoding="utf-8")
    print(text)


def stage_demo_data(directory):
    """Copy fictional metadata and explicitly mapped geometry into outputs."""
    root = Path(directory) / "demo_data"
    root.mkdir(parents=True, exist_ok=True)
    source = EXAMPLES / "data" / "databases"
    for name in ("csd_demo.tab", "coremof_demo.csv"):
        shutil.copyfile(source / name, root / name)
    structures = root / "structures"
    structures.mkdir(exist_ok=True)
    mapping = json.loads((source / "structure_map.json").read_text())
    for identifier, filename in mapping.items():
        shutil.copyfile(CRYSTALS / filename, structures / f"{identifier}.cif")
    return root / "csd_demo.tab", root / "coremof_demo.csv", structures


def tool_payload(result):
    """Decode the SDK's structured or text JSON result without hiding tool errors."""
    if getattr(result, "isError", False):
        raise RuntimeError(f"MCP tool failed: {result.content}")
    payload = getattr(result, "structuredContent", None)
    if payload is None:
        payload = next(block.text for block in result.content if block.type == "text")
    if isinstance(payload, dict) and "result" in payload:
        payload = payload["result"]
    if isinstance(payload, str):
        payload = json.loads(payload)
    if not isinstance(payload, dict):
        raise TypeError("Expected a JSON object from the MCP tool.")
    return payload
