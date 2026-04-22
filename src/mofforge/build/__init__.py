"""MOF construction subsystem.

Provides a common interface for building MOFs from scratch using
different backends (TOBACCO, Pormake).

Quick start::

    from mofforge.build import MOFBuilder

    # Build with TOBACCO
    builder = MOFBuilder(backend="tobacco")
    result = builder.build(topology="pcu")

    # Build with Pormake
    builder = MOFBuilder(backend="pormake")
    builder.add_node("node.xyz")
    builder.add_edge("edge.xyz")
    result = builder.build(topology="pcu", output_dir="./output")
    crystal = result.crystal
"""

from mofforge.build.base import BuilderBackend, BuildingBlock, BuildResult, Topology
from mofforge.build.builder import MOFBuilder
from mofforge.build.config import BuildConfig, ConfigError

__all__ = [
    "BuildConfig",
    "BuilderBackend",
    "BuildingBlock",
    "BuildResult",
    "ConfigError",
    "MOFBuilder",
    "Topology",
]
