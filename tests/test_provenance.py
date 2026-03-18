"""Tests for the Provenance class."""

import json
import tempfile
from pathlib import Path

from mofforge.provenance import Provenance


class TestProvenance:
    """Tests for Provenance serialization, chaining, and summary."""

    def test_to_dict(self):
        """to_dict should return all fields."""
        prov = Provenance(
            parent="IRMOF-1",
            query="BDC.xyz",
            replacement="NH2-BDC.xyz",
            operation="replace",
            parameters={"nb_loc": 3},
        )
        d = prov.to_dict()
        assert d["parent"] == "IRMOF-1"
        assert d["query"] == "BDC.xyz"
        assert d["replacement"] == "NH2-BDC.xyz"
        assert d["operation"] == "replace"
        assert d["parameters"]["nb_loc"] == 3
        assert "timestamp" in d
        assert isinstance(d["alignment_errors"], list)
        assert isinstance(d["history"], list)

    def test_from_dict_roundtrip(self):
        """from_dict(to_dict()) should produce an equivalent Provenance."""
        prov = Provenance(
            parent="MOF-74",
            query="linker.xyz",
            replacement="new_linker.xyz",
            operation="replace",
            parameters={"random": True},
            alignment_errors=[0.01, 0.02],
        )
        d = prov.to_dict()
        restored = Provenance.from_dict(d)
        assert restored.parent == prov.parent
        assert restored.query == prov.query
        assert restored.replacement == prov.replacement
        assert restored.operation == prov.operation
        assert restored.parameters == prov.parameters
        assert restored.alignment_errors == prov.alignment_errors
        assert restored.timestamp == prov.timestamp

    def test_to_json_from_json_roundtrip(self):
        """JSON file roundtrip should preserve all data."""
        prov = Provenance(
            parent="UiO-66",
            query="BDC.xyz",
            replacement=None,
            operation="remove",
            parameters={"disconnected": True},
        )
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            filepath = f.name

        try:
            prov.to_json(filepath)
            restored = Provenance.from_json(filepath)
            assert restored.parent == prov.parent
            assert restored.query == prov.query
            assert restored.replacement == prov.replacement
            assert restored.operation == prov.operation
            assert restored.parameters == prov.parameters
        finally:
            Path(filepath).unlink(missing_ok=True)

    def test_json_file_is_valid_json(self):
        """The written JSON file should be valid JSON."""
        prov = Provenance(parent="test", operation="search")
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            filepath = f.name

        try:
            prov.to_json(filepath)
            with open(filepath, encoding="utf-8") as f:
                data = json.load(f)
            assert isinstance(data, dict)
            assert data["parent"] == "test"
        finally:
            Path(filepath).unlink(missing_ok=True)

    def test_chain_appends_to_history(self):
        """chain() should append current provenance to new one's history."""
        prov1 = Provenance(parent="A", operation="step1")
        prov2 = Provenance(parent="B", operation="step2")

        result = prov1.chain(prov2)

        assert result is prov2
        assert len(result.history) == 1
        assert result.history[0]["parent"] == "A"
        assert result.history[0]["operation"] == "step1"

    def test_chain_preserves_existing_history(self):
        """chain() should merge histories, not overwrite."""
        prov1 = Provenance(parent="A", operation="step1")
        prov2 = Provenance(parent="B", operation="step2")
        prov3 = Provenance(
            parent="C",
            operation="step3",
            history=[{"parent": "Z", "operation": "pre-existing"}],
        )

        intermediate = prov1.chain(prov2)
        result = intermediate.chain(prov3)

        # prov3.history should now contain: prov1's history (empty) + prov1.to_dict()
        # + prov2's history (which now has prov1) + prov2.to_dict() + prov3's original history
        assert len(result.history) >= 2
        # The chain should include step1 and step2 at minimum
        operations = [h.get("operation") for h in result.history]
        assert "step1" in operations
        assert "step2" in operations

    def test_chain_three_steps(self):
        """Chaining three provenances should build up history correctly."""
        prov1 = Provenance(parent="A", operation="step1")
        prov2 = Provenance(parent="B", operation="step2")
        prov3 = Provenance(parent="C", operation="step3")

        p2 = prov1.chain(prov2)
        p3 = p2.chain(prov3)

        assert len(p3.history) == 2
        assert p3.history[0]["operation"] == "step1"
        assert p3.history[1]["operation"] == "step2"
        assert p3.operation == "step3"

    def test_summary_format(self):
        """summary() should produce a readable string."""
        prov = Provenance(
            parent="test",
            operation="replace",
            history=[
                {"operation": "load", "timestamp": "2026-01-01T00:00:00"},
                {"operation": "search", "timestamp": "2026-01-01T00:01:00"},
            ],
        )
        s = prov.summary()
        assert "Step 1" in s
        assert "Step 2" in s
        assert "Current: replace" in s
        assert "load" in s
        assert "search" in s

    def test_summary_empty_history(self):
        """summary() with no history should show only the current step."""
        prov = Provenance(parent="test", operation="validate")
        s = prov.summary()
        assert "Current: validate" in s

    def test_from_dict_defaults(self):
        """from_dict with minimal dict should use defaults for missing fields."""
        prov = Provenance.from_dict({"operation": "test"})
        assert prov.parent is None
        assert prov.query is None
        assert prov.replacement is None
        assert prov.parameters == {}
        assert prov.alignment_errors == []
        assert prov.history == []

    def test_repr(self):
        """repr should be informative."""
        prov = Provenance(parent="A", query="B", replacement="C", operation="replace")
        r = repr(prov)
        assert "replace" in r
        assert "A" in r
