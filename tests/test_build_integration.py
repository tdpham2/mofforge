"""Real builds; enabled explicitly in the builder CI job."""

import json
import os
from pathlib import Path

import pytest

from mofforge.build import MOFBuilder

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        os.environ.get("MOFFORGE_RUN_BUILD_TESTS") != "1",
        reason="Set MOFFORGE_RUN_BUILD_TESTS=1 for real backend builds.",
    ),
]


def _assert_artifacts(result):
    assert result.success, result.errors
    assert result.crystal is not None
    assert result.crystal.n_atoms > 0
    assert result.validation is not None
    for output in result.output_paths:
        manifest = json.loads(output.with_suffix(".cif.json").read_text())
        assert manifest["provenance"]["operation"] == "build"
        assert manifest["validation"]["checks_performed"]


def test_real_pormake_build(tmp_path):
    pytest.importorskip("pormake")
    builder = MOFBuilder("pormake", output_dir=tmp_path)
    builder.add_node("N108")
    builder.add_edge("E1")
    result = builder.build("pcu", output_dir=tmp_path, accuracy=3)
    _assert_artifacts(result)


def test_real_tobacco_build(tmp_path):
    pytest.importorskip("tobacco3")
    # Tests can use a staged verified archive; CI exercises verified download.
    data = os.environ.get("MOFFORGE_TEST_TOBACCO_DATA")
    builder = MOFBuilder("tobacco", **({"tobacco_data_dir": data} if data else {}))
    root = Path(builder.backend.status()["data_dir"])
    fixture = root / "tests" / "fixtures"
    for node in ("4c_1Zn_Ch", "triazole"):
        builder.add_node(fixture / "nodes" / f"{node}.cif")
    for edge in ("ntn_edge", "oxalic_edge", "squOxa_ch"):
        builder.add_edge(fixture / "edges" / f"{edge}.cif")
    _assert_artifacts(builder.build("dmc", output_dir=tmp_path))
