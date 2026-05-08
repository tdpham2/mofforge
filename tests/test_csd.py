"""Tests for the CSD lookup module."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from mofforge.csd.database import (
    CSDDatabase,
    _extract_ccdc,
    _extract_doi,
    get_database,
)
from mofforge.csd.models import CSDRecord, CSDSearchResult

# ---------------------------------------------------------------------------
# Synthetic test data — mirrors the real MOF_subset.tab structure
# ---------------------------------------------------------------------------

# All 61 column headers from the real file
_HEADER = (
    "[REFCODE]\t[_publ_authors]\t[_journal_name_full]\t[_journal_volume]\t"
    "[_journal_pages]\t[_journal_year]\t[_chemical_formula_moiety]\t"
    "[_chemical_name_systematic]\t[_chemical_name_common]\t"
    "[_chemical_compound_source]\t[_chemical_melting_point]\t"
    "[_exptl_crystal_colour]\t[_chemical_QUAL]\t[_chemical_NOTE]\t"
    "[_symmetry_space_group_name]\t[_symmetry_Int_Tables_number]\t"
    "[_cell_length_a]\t[_cell_length_b]\t[_cell_length_c]\t"
    "[_cell_angle_alpha]\t[_cell_angle_beta]\t[_cell_angle_gamma]\t"
    "[_cell_volume]\t[_rcell_length_a]\t[_rcell_length_b]\t[_rcell_length_c]\t"
    "[_rcell_angle_alpha]\t[_rcell_angle_beta]\t[_rcell_angle_gamma]\t"
    "[_rcell_volume]\t[_molecular_volume]\t[_cell_RESIDUES]\t"
    "[_cell_formula_units_Z]\t[_cell_formula_units_Zprime]\t"
    "[_exptl_crystal_preparation]\t[_exptl_crystal_description]\t"
    "[POLYMORPH]\t[_exptl_crystal_NOTES]\t[_refine_ls_R_factor]\t"
    "[_cell_measurement_temperature]\t[_exptl_crystal_density_CCDC]\t"
    "[_exptl_crystal_density_diffrn]\t[_diffrn_measurement_device]\t"
    "[_exptl_DISORDER]\t[_exptl_QUAL]\t[_ccdc_REMARK]\t"
    "[PREF]\t[BATC]\t[ADAT]\t[MDAT]\t[NBSI]\t[CDRE]\t[MSDB]\t[NBIT]\t"
    "[MATF]\t[CATF]\t[CCOM]\t[_ccdc_NOTES]\t"
    "[_ccdc_RMARKS]\t[_ccdc_QUAL]\t[_ccdc_PROPS]"
)


def _row(refcode, name_sys, name_common, formula, year, space_group, remarks, **kw):
    """Build a tab-separated row with 61 fields."""
    fields = [""] * 61
    fields[0] = refcode
    fields[1] = kw.get("authors", "A.Test")
    fields[2] = kw.get("journal", "J.Test")
    fields[3] = kw.get("volume", "1")
    fields[4] = kw.get("pages", "100")
    fields[5] = year
    fields[6] = formula
    fields[7] = name_sys
    fields[8] = name_common
    fields[14] = space_group
    fields[38] = kw.get("r_factor", "0.05")
    fields[39] = kw.get("temperature", "293")
    fields[58] = remarks
    return "\t".join(fields)


_ROWS = [
    _row(
        "TESTMOF", "poly[test-systematic]", "TEST-MOF",
        "(C6 H12 Zn1 O4)n", "2020", "P-1",
        '"CCDC: CCDC 999999 "DOI: 10.1234/test.2020',
    ),
    _row(
        "TESTMOF01", "poly[test-systematic-v2]", "",
        "(C6 H12 Zn1 O4)n", "2021", "P-1",
        '"CCDC: CCDC 999998',
    ),
    _row(
        "IRMOF1", "catena-[Zinc 1,4-benzenedicarboxylate]", "IRMOF-1",
        "(C24 H12 O13 Zn4)n", "1999", "Fm-3m",
        '"CCDC: CCDC 100001 "DOI: 10.1038/46248',
        authors="O.Yaghi",
    ),
    _row(
        "HKUST1", "catena-[Copper trimesate]", "HKUST-1",
        "(C18 H6 Cu3 O12)n", "1999", "Fm-3m",
        '"CCDC: CCDC 100002 "DOI: 10.1126/science.283.5405.1148',
        authors="S.Chui",
    ),
    _row(
        "UIOMOF", "Zirconium 1,4-benzenedicarboxylate", "UiO-66",
        "(C48 H28 O32 Zr6)n", "2008", "Fm-3m",
        '"CCDC: CCDC 200003 "DOI: 10.1021/ja8057953',
    ),
]


@pytest.fixture
def sample_tsv(tmp_path):
    """Create a small synthetic CSD TSV file."""
    tsv = tmp_path / "test_subset.tab"
    content = _HEADER + "\n" + "\n".join(_ROWS) + "\n"
    tsv.write_text(content, encoding="utf-8")
    return tsv


@pytest.fixture
def db(sample_tsv):
    """Return a CSDDatabase backed by the synthetic fixture."""
    database = CSDDatabase(sample_tsv)
    yield database
    database.close()


# ---------------------------------------------------------------------------
# Extraction helpers
# ---------------------------------------------------------------------------


class TestExtraction:
    def test_extract_doi(self):
        assert _extract_doi('"DOI: 10.1234/test.2020') == "10.1234/test.2020"

    def test_extract_doi_with_trailing_quote(self):
        assert _extract_doi('"DOI: 10.1038/46248"') == "10.1038/46248"

    def test_extract_doi_none(self):
        assert _extract_doi("no doi here") is None

    def test_extract_ccdc(self):
        assert _extract_ccdc('"CCDC: CCDC 999999') == "999999"

    def test_extract_ccdc_none(self):
        assert _extract_ccdc("no ccdc") is None


# ---------------------------------------------------------------------------
# Auto-detection
# ---------------------------------------------------------------------------


class TestFieldDetection:
    def test_detect_refcode(self):
        assert CSDDatabase._detect_field("ABACUF") == "refcode"
        assert CSDDatabase._detect_field("TESTMOF01") == "refcode"

    def test_detect_doi(self):
        assert CSDDatabase._detect_field("10.1234/test") == "doi"

    def test_detect_ccdc(self):
        assert CSDDatabase._detect_field("999999") == "ccdc"

    def test_detect_name(self):
        assert CSDDatabase._detect_field("UiO-66") == "name"
        assert CSDDatabase._detect_field("copper trimesate") == "name"


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------


class TestModels:
    def test_record_summary(self):
        rec = CSDRecord(refcode="TEST", chemical_name_common="Test-MOF", year="2020")
        s = rec.summary()
        assert "TEST" in s
        assert "Test-MOF" in s
        assert "(2020)" in s

    def test_record_summary_with_doi(self):
        rec = CSDRecord(refcode="TEST", doi="10.1234/x")
        assert "DOI:10.1234/x" in rec.summary()

    def test_search_result(self):
        result = CSDSearchResult(
            query="test",
            field="name",
            records=[CSDRecord(refcode="A"), CSDRecord(refcode="B")],
        )
        assert result.n_matches == 2
        assert "2 match(es)" in result.summary()


# ---------------------------------------------------------------------------
# Database operations
# ---------------------------------------------------------------------------


class TestDatabase:
    def test_build_and_count(self, db):
        assert db.record_count() == 5

    def test_lookup_refcode(self, db):
        rec = db.lookup_refcode("TESTMOF")
        assert rec is not None
        assert rec.refcode == "TESTMOF"
        assert rec.doi == "10.1234/test.2020"
        assert rec.ccdc_number == "999999"
        assert rec.year == "2020"

    def test_lookup_refcode_case_insensitive(self, db):
        rec = db.lookup_refcode("testmof")
        assert rec is not None
        assert rec.refcode == "TESTMOF"

    def test_lookup_refcode_missing(self, db):
        assert db.lookup_refcode("NONEXIST") is None

    def test_search_name(self, db):
        results = db.search_name("IRMOF")
        assert len(results) == 1
        assert results[0].refcode == "IRMOF1"

    def test_search_name_systematic(self, db):
        results = db.search_name("trimesate")
        assert len(results) == 1
        assert results[0].refcode == "HKUST1"

    def test_search_doi(self, db):
        results = db.search_doi("10.1038/46248")
        assert len(results) == 1
        assert results[0].refcode == "IRMOF1"

    def test_search_formula(self, db):
        results = db.search_formula("C18 H6 Cu3")
        assert len(results) == 1
        assert results[0].refcode == "HKUST1"

    def test_search_ccdc(self, db):
        results = db.search_ccdc("200003")
        assert len(results) == 1
        assert results[0].refcode == "UIOMOF"

    def test_search_auto_refcode(self, db):
        result = db.search("HKUST1")
        assert result.field == "refcode"
        assert result.n_matches == 1

    def test_search_auto_doi(self, db):
        result = db.search("10.1234/test.2020")
        assert result.field == "doi"
        assert result.n_matches == 1

    def test_search_auto_name(self, db):
        result = db.search("UiO-66")
        assert result.field == "name"
        assert result.n_matches == 1

    def test_search_no_matches(self, db):
        result = db.search("nonexistent-mof-name")
        assert result.n_matches == 0

    def test_raw_field_populated(self, db):
        rec = db.lookup_refcode("IRMOF1")
        assert rec is not None
        assert rec.raw.get("REFCODE") == "IRMOF1"
        assert rec.raw.get("_publ_authors") == "O.Yaghi"

    def test_sqlite_cache_reused(self, sample_tsv):
        """Second instantiation should reuse the .db file without rebuilding."""
        db1 = CSDDatabase(sample_tsv)
        assert db1.record_count() == 5
        db1.close()

        db_path = sample_tsv.with_suffix(".db")
        assert db_path.exists()
        mtime = db_path.stat().st_mtime

        db2 = CSDDatabase(sample_tsv)
        assert db2.record_count() == 5
        db2.close()

        # The db file should not have been regenerated
        assert db_path.stat().st_mtime == mtime


# ---------------------------------------------------------------------------
# get_database() singleton
# ---------------------------------------------------------------------------


class TestGetDatabase:
    def test_explicit_path(self, sample_tsv):
        db = get_database(data_path=sample_tsv)
        assert db.record_count() == 5
        db.close()

    def test_missing_path_raises(self, monkeypatch):
        # Clear all config sources
        from mofforge.csd import database as db_module
        from mofforge.utils.config import config

        monkeypatch.setattr(config, "csd_data_path", None)
        monkeypatch.delenv("MOFFORGE_CSD_DATA_PATH", raising=False)
        monkeypatch.setattr(db_module, "_db", None)

        with pytest.raises(FileNotFoundError, match="CSD data path is not configured"):
            get_database()

    def test_env_var(self, sample_tsv, monkeypatch):
        from mofforge.csd import database as db_module
        from mofforge.utils.config import config

        monkeypatch.setattr(config, "csd_data_path", None)
        monkeypatch.setenv("MOFFORGE_CSD_DATA_PATH", str(sample_tsv))
        monkeypatch.setattr(db_module, "_db", None)

        db = get_database()
        assert db.record_count() == 5
        db.close()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


class TestCLI:
    def test_csd_command(self, sample_tsv):
        from click.testing import CliRunner

        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["csd", "TESTMOF", "--data-path", str(sample_tsv)])
        assert result.exit_code == 0
        assert "TESTMOF" in result.output
        assert "1 match" in result.output

    def test_csd_command_verbose(self, sample_tsv):
        from click.testing import CliRunner

        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(
            main, ["csd", "IRMOF1", "--data-path", str(sample_tsv), "-v"]
        )
        assert result.exit_code == 0
        assert "IRMOF1" in result.output
        assert "Fm-3m" in result.output

    def test_csd_command_no_match(self, sample_tsv):
        from click.testing import CliRunner

        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["csd", "ZZZNOPE", "--data-path", str(sample_tsv)])
        assert result.exit_code == 0
        assert "0 match" in result.output

    def test_csd_command_missing_config(self, monkeypatch):
        from click.testing import CliRunner

        from mofforge.csd import database as db_module
        from mofforge.cli import main
        from mofforge.utils.config import config

        monkeypatch.setattr(config, "csd_data_path", None)
        monkeypatch.delenv("MOFFORGE_CSD_DATA_PATH", raising=False)
        monkeypatch.setattr(db_module, "_db", None)

        runner = CliRunner()
        result = runner.invoke(main, ["csd", "TESTMOF"])
        assert result.exit_code != 0
