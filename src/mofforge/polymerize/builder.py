"""Unified facade for amorphous POP (porous organic polymer) generation.

Mirrors :class:`mofforge.build.builder.MOFBuilder`: collect monomers, then
``build()`` to run simulated polymerization.  The heavy lifting (packing +
bond formation + MD) is delegated to a backend that wraps :mod:`pysimm`; this
class only manages the monomer list and configuration.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from mofforge.polymerize.base import Monomer, PopConfig, POPResult, ReactiveSite

logger = logging.getLogger("mofforge")

_BACKENDS = ("pysimm",)


class PopBuilder:
    """Build amorphous porous organic polymers from reactive monomers."""

    def __init__(self, backend: str = "pysimm", **kwargs: Any) -> None:
        if backend not in _BACKENDS:
            raise ValueError(f"Unknown backend {backend!r}.  Choose from: {_BACKENDS}")
        self._backend_name = backend
        self._monomers: list[Monomer] = []
        self._backend_kwargs = kwargs

    @property
    def backend_name(self) -> str:
        """Short name of the active backend."""
        return self._backend_name

    def add_monomer(
        self,
        source: str | Path,
        name: str | None = None,
        functionality: int | None = None,
        sites: list[ReactiveSite] | None = None,
    ) -> Monomer:
        """Register a monomer (SMILES or XYZ/CIF path).

        Reactive sites and functionality are detected from the SMILES at build
        time when not supplied here.
        """
        if name is None:
            p = Path(str(source))
            if p.suffix in (".cif", ".xyz", ".mol2"):
                name = p.stem
            else:
                name = f"monomer{len(self._monomers) + 1}"
        monomer = Monomer(
            name=name,
            source=source,
            functionality=functionality if functionality is not None else 2,
            sites=list(sites) if sites is not None else [],
        )
        self._monomers.append(monomer)
        return monomer

    def list_monomers(self) -> list[str]:
        """Names of the registered monomers."""
        return [m.name for m in self._monomers]

    def build(self, output_dir: str | Path = ".", **options: Any) -> POPResult:
        """Run simulated polymerization and write a P1 CIF of the network.

        Any :class:`PopConfig` field may be passed as a keyword; the rest use
        defaults.  Requires ``pysimm`` plus the Packmol and LAMMPS binaries.
        """
        if not self._monomers:
            return POPResult(
                success=False,
                backend=self._backend_name,
                errors=["No monomers registered; call add_monomer() first."],
            )

        config = self._make_config(options)
        backend = self._load_backend()
        return backend.polymerize(
            monomers=list(self._monomers),
            config=config,
            output_dir=Path(output_dir),
        )

    def _load_backend(self):
        """Import the backend lazily so pysimm stays an optional dependency."""
        from mofforge.polymerize.engine import PysimmBackend

        return PysimmBackend(**self._backend_kwargs)

    @staticmethod
    def _make_config(options: dict[str, Any]) -> PopConfig:
        """Build a :class:`PopConfig` from build() keyword options."""
        fields = {
            "box_length",
            "target_density",
            "forcefield",
            "target_conversion",
            "n_monomers",
            "equilibrate",
            "random_seed",
            "md_settings",
        }
        kwargs = {k: v for k, v in options.items() if k in fields}
        unknown = set(options) - fields
        if unknown:
            raise TypeError(f"Unknown build option(s): {', '.join(sorted(unknown))}")
        return PopConfig(**kwargs)
