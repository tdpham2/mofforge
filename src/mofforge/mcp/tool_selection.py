"""Shared startup-time MCP tool selection.

The selector deliberately works before a server begins accepting requests.  It
therefore supports MCP Python SDK 1.16, which has public tool registration but
does not yet have the public ``FastMCP.remove_tool`` method.
"""

from __future__ import annotations

import importlib.util
from collections.abc import Callable, Collection, Mapping


class ToolSelectionError(ValueError):
    """Raised when an explicit MCP tool selection cannot be satisfied."""


def parse_tool_list(value: str | None) -> set[str] | None:
    """Parse a comma-separated tool allowlist.

    ``None`` means no explicit allowlist. An empty string intentionally selects
    no tools, which is valid in MCP.
    """
    if value is None:
        return None
    return {name.strip() for name in value.split(",") if name.strip()}


def capability_available(capability: str) -> bool:
    """Return whether an optional mofforge capability is importable."""
    modules = {
        "vis": ("playwright",),
        "chem": ("rdkit",),
        # Either construction backend makes the shared build tools useful.
        "build": ("pormake", "tobacco3"),
    }.get(capability)
    if modules is None:
        raise ValueError(f"Unknown MCP tool capability: {capability!r}")
    return any(importlib.util.find_spec(module) is not None for module in modules)


def select_tool_names(
    all_names: Collection[str],
    requirements: Mapping[str, str | None],
    *,
    requested: Collection[str] | None = None,
    available_only: bool = False,
    availability: Callable[[str], bool] | None = None,
) -> set[str]:
    """Resolve the names that should be registered for one server process.

    Explicit allowlists are strict: unknown names and names whose optional
    capability is unavailable are rejected. Without an explicit allowlist,
    dependency filtering occurs only when ``available_only`` is true, retaining
    the historical full catalog by default.
    """
    if availability is None:
        availability = capability_available

    known = set(all_names)
    requested_names = set(requested) if requested is not None else None

    if requested_names is not None:
        unknown = requested_names - known
        if unknown:
            raise ToolSelectionError(
                "Unknown MCP tool name(s): " + ", ".join(sorted(unknown))
            )

        unavailable = sorted(
            name
            for name in requested_names
            if (capability := requirements.get(name)) is not None
            and not availability(capability)
        )
        if unavailable:
            details = ", ".join(
                f"{name} ({requirements[name]})" for name in unavailable
            )
            raise ToolSelectionError(
                "Requested MCP tool(s) require unavailable optional dependencies: "
                + details
            )

    selected = known if requested_names is None else requested_names
    if available_only:
        selected = {
            name
            for name in selected
            if (capability := requirements.get(name)) is None
            or availability(capability)
        }
    return selected
