"""End-to-end POP generation via pysimm + Packmol + LAMMPS.

Marked ``integration``: skipped unless pysimm is importable and both the Packmol
and LAMMPS binaries resolve.  Mirrors ``test_build_integration.py`` for the MOF
builders.  Run with:  pytest -m integration tests/test_polymerize_integration.py
"""

from __future__ import annotations

import pytest

from mofforge.polymerize import PopBuilder
from mofforge.polymerize.config import ConfigError, PopBuildConfig, validate_pysimm


def _binaries_available() -> bool:
    if validate_pysimm():
        return False
    cfg = PopBuildConfig.load()
    try:
        cfg.resolve_packmol_binary()
        cfg.resolve_lammps_binary()
    except ConfigError:
        return False
    return True


pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        not _binaries_available(),
        reason="requires pysimm plus Packmol and LAMMPS binaries",
    ),
]


def test_imine_network_end_to_end(tmp_path):
    """Diamine + dialdehyde -> connected imine network, clash-free CIF."""
    builder = PopBuilder()
    builder.add_monomer("NCCN", name="diamine", functionality=2)
    builder.add_monomer("O=Cc1ccc(C=O)cc1", name="dialdehyde", functionality=2)

    result = builder.build(
        output_dir=str(tmp_path),
        target_density=0.6,
        target_conversion=0.8,
        n_monomers=40,
        equilibrate=False,
        random_seed=42,
    )

    assert result.success, result.errors
    assert result.crystal is not None
    assert result.crystal.n_atoms > 0
    assert result.crystal.n_bonds > 0
    # A CIF and its manifest were written.
    cif = result.output_paths[0]
    assert cif.is_file()
    assert cif.with_suffix(cif.suffix + ".json").is_file()
    # No severe steric overlaps in the final structure.
    assert not result.validation.steric_clashes
    assert result.metadata["n_monomers"] == 40


def test_determinism_same_seed(tmp_path):
    """Two runs with the same seed give the same atom/bond counts."""
    counts = []
    for i in range(2):
        builder = PopBuilder()
        builder.add_monomer("NCCN", name="diamine")
        builder.add_monomer("O=Cc1ccc(C=O)cc1", name="dialdehyde")
        result = builder.build(
            output_dir=str(tmp_path / f"run{i}"),
            target_density=0.6,
            n_monomers=20,
            equilibrate=False,
            random_seed=7,
        )
        assert result.success, result.errors
        counts.append((result.crystal.n_atoms, result.crystal.n_bonds))
    assert counts[0] == counts[1]
