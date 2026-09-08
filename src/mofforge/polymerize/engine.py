"""pysimm-backed simulated-polymerization engine.

Orchestration only: mofforge prepares monomers (RDKit), sizes the packing box,
then drives :mod:`pysimm` — Packmol packing via ``apps.polymatic.pack`` and
Polymatic bond formation via ``apps.polymatic.polymatic`` (which itself runs
LAMMPS) — and finally converts the result back into a mofforge ``Crystal`` for
validation, provenance, and CIF output.  No packing, bond-formation, or MD logic
lives here; those belong to pysimm/Packmol/LAMMPS.

The external tools are resolved at run time (see
:mod:`mofforge.polymerize.config`); a missing one raises an actionable
:class:`ConfigError` rather than failing deep inside pysimm.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from mofforge.polymerize.base import Monomer, PopConfig, POPResult, Timer
from mofforge.polymerize.config import ConfigError, PopBuildConfig, validate_pysimm
from mofforge.polymerize.monomer import (
    PreparedMonomer,
    box_length_for_density,
    prepare_monomer,
)

logger = logging.getLogger("mofforge")


def _get_pysimm():
    """Import pysimm lazily with an actionable error message."""
    errors = validate_pysimm()
    if errors:
        raise ImportError(errors[0])
    import pysimm  # noqa: F401
    from pysimm import system
    from pysimm.apps import polymatic

    return system, polymatic


class PysimmBackend:
    """Backend that delegates POP construction to pysimm (Packmol + LAMMPS)."""

    name: str = "pysimm"

    def __init__(self, **kwargs: Any) -> None:
        # Binary locations may be passed through from PopBuilder kwargs.
        self._cfg = PopBuildConfig.load(**kwargs)

    def polymerize(
        self,
        monomers: list[Monomer],
        config: PopConfig,
        output_dir: Path,
        **options: Any,
    ) -> POPResult:
        """Run simulated polymerization and return a validated ``POPResult``."""
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        timer = Timer()

        try:
            with timer:
                crystal, metadata = self._run(monomers, config, output_dir)
        except (ConfigError, ImportError) as exc:
            # Missing tool / dependency: report cleanly, do not raise.
            return POPResult(
                success=False,
                backend=self.name,
                errors=[str(exc)],
                elapsed_seconds=round(timer.elapsed, 2),
            )
        except Exception as exc:
            logger.warning("polymerize failed", exc_info=True)
            return POPResult(
                success=False,
                backend=self.name,
                errors=[f"pysimm polymerization failed: {exc}"],
                elapsed_seconds=round(timer.elapsed, 2),
            )

        # Validate + write CIF + manifest (mofforge side).
        from mofforge.core.bonding import infer_bonds
        from mofforge.provenance import record_operation, write_manifest
        from mofforge.validation import validate_structure

        report = validate_structure(crystal)
        if crystal.n_bonds == 0:
            crystal = infer_bonds(crystal)

        cif_path = output_dir / f"{crystal.name}.cif"
        record_operation(
            crystal,
            parent=crystal,
            operation="polymerize",
            parameters={
                "backend": self.name,
                "monomers": [str(m.source) for m in monomers],
                "target_density": config.target_density,
                "forcefield": config.forcefield,
                "target_conversion": config.target_conversion,
                "equilibrate": config.equilibrate,
                "random_seed": config.random_seed,
            },
            validation=report,
        )
        crystal.write_cif(cif_path)
        write_manifest(crystal, cif_path, report)

        return POPResult(
            success=report.is_valid,
            output_paths=[cif_path],
            crystal=crystal,
            elapsed_seconds=round(timer.elapsed, 2),
            backend=self.name,
            metadata={**metadata, "forcefield": config.forcefield},
            validation=report,
        )

    def _run(
        self,
        monomers: list[Monomer],
        config: PopConfig,
        output_dir: Path,
    ):
        """Execute the pysimm pipeline; returns (Crystal, metadata dict)."""
        # Resolve binaries up front so failures are actionable, not buried.
        packmol_bin = self._cfg.resolve_packmol_binary()
        lammps_bin = self._cfg.resolve_lammps_binary()
        logger.info("Using Packmol=%s  LAMMPS=%s", packmol_bin, lammps_bin)

        system, polymatic = _get_pysimm()

        # 1. Prepare each monomer (RDKit geometry + reactive sites).
        prepared: list[PreparedMonomer] = [
            prepare_monomer(m, output_dir, random_seed=config.random_seed) for m in monomers
        ]

        # 2. Decide monomer counts and box size.
        counts = self._monomer_counts(len(prepared), config.n_monomers)
        boxl = config.box_length or box_length_for_density(
            prepared, counts, config.target_density
        )

        # 3. Load monomer reference systems and apply the force field.
        from pysimm import forcefield

        ff = _get_forcefield(forcefield, config.forcefield)
        ref_files: list[str] = []
        for prep in prepared:
            mono_sys = _read_monomer_system(system, prep)
            mono_sys.apply_forcefield(ff)
            ref_path = output_dir / f"{prep.name}_typed.lmps"
            mono_sys.write_lammps(str(ref_path))
            ref_files.append(str(ref_path))

        # 4. Pack the box (Packmol, via pysimm's polymatic.pack).
        packed_path = output_dir / "packed.lmps"
        pack_ok = polymatic.pack(
            script="pack",
            file_in=ref_files,
            nrep=counts,
            boxl=boxl,
            file_out=str(packed_path),
        )
        if pack_ok is False:
            raise RuntimeError("Packmol packing (polymatic.pack) failed.")

        # 5. Simulated polymerization (Polymatic bond loop + LAMMPS relaxation).
        #    polym.in and types.txt are expected in the working directory; the
        #    caller supplies them via config.md_settings or a template.
        final_path = output_dir / "polymerized.lmps"
        polym_ok = polymatic.polymatic(
            script="polym.pl",
            file_in=str(packed_path),
            file_out=str(final_path),
        )
        if polym_ok is False:
            raise RuntimeError("Simulated polymerization (polymatic.polymatic) failed.")

        # 6. Read the final system and convert to a Crystal.
        final_system = system.read_lammps(str(final_path))
        from mofforge.polymerize.convert import system_to_crystal

        name = "_".join(p.name for p in prepared) or "pop"
        crystal = system_to_crystal(final_system, name=name)

        metadata = {
            "n_monomers": sum(counts),
            "monomer_counts": dict(zip((p.name for p in prepared), counts, strict=True)),
            "box_length": round(boxl, 3),
            "target_density": config.target_density,
            "n_atoms": crystal.n_atoms,
            "n_bonds": crystal.n_bonds,
            "packmol_bin": str(packmol_bin),
            "lammps_bin": str(lammps_bin),
        }
        return crystal, metadata

    @staticmethod
    def _monomer_counts(n_types: int, n_monomers: int | None) -> list[int]:
        """Split a total monomer count evenly across monomer types."""
        total = n_monomers if n_monomers is not None else 20 * n_types
        base, extra = divmod(total, n_types)
        return [base + (1 if i < extra else 0) for i in range(n_types)]


def _get_forcefield(forcefield_mod, name: str):
    """Instantiate a pysimm force field by name (gaff2 | dreiding | pcff)."""
    classes = {
        "gaff2": "Gaff2",
        "dreiding": "Dreiding",
        "pcff": "Pcff",
    }
    cls_name = classes.get(name.lower())
    if cls_name is None:
        raise ValueError(
            f"Unknown force field {name!r}. Choose from: {', '.join(sorted(classes))}."
        )
    return getattr(forcefield_mod, cls_name)()


def _read_monomer_system(system_mod, prep: PreparedMonomer):
    """Read a prepared monomer file into a pysimm ``System``."""
    path = str(prep.mol_path)
    if prep.mol_path.suffix.lower() == ".mol":
        return system_mod.read_mol(path)
    return system_mod.read_xyz(path)
