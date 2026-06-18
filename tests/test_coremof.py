"""Tests for the CoRE MOF database module."""

from __future__ import annotations

import csv
import io
import os
from pathlib import Path

import pytest

from mofforge.coremof.database import (
    CoreMOFDatabase,
    _derive_base_refcode,
    _safe_float,
    _safe_int,
    csd_to_coremof,
    get_database,
    search_csd_name,
)
from mofforge.coremof.models import BridgeResult, CoreMOFRecord, CoreMOFSearchResult

# ---------------------------------------------------------------------------
# Synthetic test data — mirrors the real CoRE MOF CSV structure
# ---------------------------------------------------------------------------

_HEADER = (
    "number,coreid,refcode,name,mofid-v1,mofid-v2,"
    "LCD (Å),PLD (Å),LFPD (Å),Density (g/cm3),"
    "ASA (A2),ASA (m2/cm3),ASA (m2/g),"
    "NASA (A2),NASA (m2/cm3),NASA (m2/g),"
    "PV (A3),VF,PV (cm3/g),"
    "NAV (A3),NAV_VF,NPV (cm3/g),"
    "structure_dimension,"
    "topology(SingleNodes),topology(SingleNodes)-zeo,"
    "topology(AllNodes),topology(AllNodes)-zeo,"
    "catenation,dimension_by_topo,hall,number_spacegroup,"
    "Metal Types,Has OMS,OMS Types,Charge,"
    "average_atomic_mass,"
    "Heat_capacity@300K (J/g/K),std @ 300 K (J/g/K),"
    "Heat_capacity@350K (J/g/K),std @ 350 K (J/g/K),"
    "Heat_capacity@400K (J/g/K),std @ 400 K (J/g/K),"
    "k_cp (J/g/K/K),cp0 (J/g/K),"
    "natoms,Source,DOI,Year,Time,Publication,"
    "Extension,unmodified,"
    "Thermal_stability (℃),Solvent_stability,Water_stability,KH_Classes"
)


def _row(
    number,
    coreid,
    refcode,
    name="-",
    lcd="5.0",
    pld="3.0",
    density="1.5",
    asa_m2g="100.0",
    pv_cm3g="0.3",
    vf="0.4",
    topo_single="pcu",
    topo_all="pcu",
    metals="Zn",
    has_oms="No",
    oms_types="N/A",
    thermal="350.0",
    solvent="0.8",
    water="0.6",
    kh="weak",
    doi="10.1234/test",
    year="2020",
    publication="ACS",
    extension="All Solvent Removed",
    natoms="200",
    spacegroup="225",
):
    """Build a CSV row with all 56 fields."""
    fields = [
        number, coreid, refcode, name,
        "smiles_v1", "smiles_v2",
        lcd, pld, lcd, density,
        "0", "0", asa_m2g,
        "0", "0", "0",
        "0", vf, pv_cm3g,
        "0", "0", "0",
        "3",
        topo_single, "unnamed",
        topo_all, "unnamed",
        "1", "3", "-P 4", spacegroup,
        metals, has_oms, oms_types, "PACMAN-DDEC6",
        "12.0",
        "0.8", "0.02",
        "0.9", "0.03",
        "1.0", "0.04",
        "0.001", "0.4",
        natoms, "CSD", doi, year, year, publication,
        extension, "FALSE",
        thermal, solvent, water, kh,
    ]
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(fields)
    return buf.getvalue().rstrip("\r\n")


_ROWS = [
    _row(
        "1", "2004[Co][rtl]3[ASR]1", "ABAVIJ_ASR_pacman",
        metals="Co", topo_single="rtl", topo_all="rtl",
        lcd="4.4", pld="2.5", density="1.52",
        thermal="333.7", water="0.46", kh="weak",
        doi="10.1039/b404485a", year="2004",
    ),
    _row(
        "2", "2016[Mn][tbo]3[ASR]1", "ABAYIO_ASR_pacman",
        metals="Mn", topo_single="tbo", topo_all="tbo",
        lcd="11.4", pld="4.3", density="0.95",
        thermal="428.7", water="0.60", kh="superstrong_high_loading",
        doi="10.1021/acs.cgd.6b00776", year="2016",
        has_oms="No",
    ),
    _row(
        "3", "2021[Zn][pts]3[ASR]1", "ABINIM_ASR_pacman",
        name="JOU-34",
        metals="Zn", topo_single="pts", topo_all="pts",
        lcd="7.3", pld="6.7", density="1.2",
        thermal="349.1", water="0.31", kh="weak",
        doi="10.1021/test.2021", year="2021",
        has_oms="Yes", oms_types="Zn",
    ),
    _row(
        "4", "2020[Cu][pcu]3[ASR]1", "HKUSTX_ASR_pacman",
        name="HKUST-variant",
        metals="Cu", topo_single="pcu", topo_all="pcu",
        lcd="9.0", pld="6.0", density="0.88",
        asa_m2g="1500.0", pv_cm3g="0.7", vf="0.65",
        thermal="400.0", water="0.85", kh="strong",
        doi="10.1126/test", year="2020",
        has_oms="Yes", oms_types="Cu",
    ),
    _row(
        "5", "2020[Cu][pcu]3[ION]1", "HKUSTX_ION_pacman",
        name="HKUST-variant",
        metals="Cu", topo_single="pcu", topo_all="pcu",
        lcd="8.5", pld="5.5", density="0.92",
        thermal="390.0", water="0.80", kh="strong",
        doi="10.1126/test", year="2020",
        extension="with ion",
        has_oms="Yes", oms_types="Cu",
    ),
    _row(
        "6", "2018[Mn,Cu][dia]3[ASR]1", "MULTIM_ASR_pacman",
        metals="Mn,Cu", topo_single="dia", topo_all="dia",
        lcd="6.0", pld="4.0", density="1.1",
        thermal="unknown", water="unknown", solvent="unknown",
        kh="unknown",
        doi="10.1007/test", year="2018",
    ),
]


@pytest.fixture
def sample_csv(tmp_path):
    """Create a small synthetic CoRE MOF CSV file."""
    csv_file = tmp_path / "test_coremof.csv"
    content = _HEADER + "\n" + "\n".join(_ROWS) + "\n"
    csv_file.write_text(content, encoding="utf-8")
    return csv_file


@pytest.fixture
def db(sample_csv):
    """Return a CoreMOFDatabase backed by the synthetic fixture."""
    database = CoreMOFDatabase(sample_csv)
    yield database
    database.close()


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------


class TestParsingHelpers:
    def test_safe_float(self):
        assert _safe_float("3.14") == pytest.approx(3.14)
        assert _safe_float("0") == 0.0

    def test_safe_float_none_cases(self):
        assert _safe_float(None) is None
        assert _safe_float("") is None
        assert _safe_float("unknown") is None
        assert _safe_float("-") is None
        assert _safe_float("N/A") is None

    def test_safe_int(self):
        assert _safe_int("42") == 42
        assert _safe_int("3.0") == 3

    def test_safe_int_none_cases(self):
        assert _safe_int(None) is None
        assert _safe_int("unknown") is None

    def test_derive_base_refcode(self):
        assert _derive_base_refcode("ABAVIJ_ASR_pacman") == "ABAVIJ"
        assert _derive_base_refcode("HKUST1") == "HKUST1"
        assert _derive_base_refcode("") == ""


# ---------------------------------------------------------------------------
# Auto-detection
# ---------------------------------------------------------------------------


class TestFieldDetection:
    def test_detect_coreid(self):
        assert CoreMOFDatabase._detect_field("2004[Co][rtl]3[ASR]1") == "coreid"

    def test_detect_doi(self):
        assert CoreMOFDatabase._detect_field("10.1234/test") == "doi"

    def test_detect_refcode(self):
        assert CoreMOFDatabase._detect_field("ABAVIJ") == "refcode"
        assert CoreMOFDatabase._detect_field("HKUST1") == "refcode"

    def test_detect_metal(self):
        assert CoreMOFDatabase._detect_field("Cu") == "metal"
        assert CoreMOFDatabase._detect_field("Zn") == "metal"

    def test_detect_topology(self):
        assert CoreMOFDatabase._detect_field("pcu") == "topology"
        assert CoreMOFDatabase._detect_field("dia") == "topology"

    def test_detect_name(self):
        assert CoreMOFDatabase._detect_field("HKUST-1") == "name"
        assert CoreMOFDatabase._detect_field("some mof name") == "name"


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------


class TestModels:
    def test_record_summary(self):
        rec = CoreMOFRecord(
            coreid="2020[Cu][pcu]3[ASR]1",
            metal_types="Cu",
            topology_single="pcu",
            year="2020",
        )
        s = rec.summary()
        assert "2020[Cu][pcu]3[ASR]1" in s
        assert "[Cu]" in s
        assert "pcu" in s

    def test_record_summary_with_name(self):
        rec = CoreMOFRecord(coreid="test", name="HKUST-1")
        assert "HKUST-1" in rec.summary()

    def test_record_summary_skips_dash_name(self):
        rec = CoreMOFRecord(coreid="test", name="-")
        assert "-" not in rec.summary().split("|")[0].strip()  # not in first part

    def test_properties_summary(self):
        rec = CoreMOFRecord(
            coreid="test", refcode="TEST_ASR", base_refcode="TEST",
            metal_types="Cu", lcd=9.0, pld=6.0, water_stability=0.85,
            extension="All Solvent Removed",
        )
        ps = rec.properties_summary()
        assert "LCD" in ps
        assert "9.000" in ps
        assert "H2O_stab" in ps

    def test_search_result(self):
        result = CoreMOFSearchResult(
            query="test",
            field="name",
            records=[CoreMOFRecord(coreid="A"), CoreMOFRecord(coreid="B")],
        )
        assert result.n_matches == 2
        assert "2 match(es)" in result.summary()


# ---------------------------------------------------------------------------
# Database operations
# ---------------------------------------------------------------------------


class TestDatabase:
    def test_build_and_count(self, db):
        assert db.record_count() == 6

    def test_lookup_coreid(self, db):
        rec = db.lookup_coreid("2004[Co][rtl]3[ASR]1")
        assert rec is not None
        assert rec.coreid == "2004[Co][rtl]3[ASR]1"
        assert rec.refcode == "ABAVIJ_ASR_pacman"
        assert rec.base_refcode == "ABAVIJ"
        assert rec.metal_types == "Co"
        assert rec.topology_single == "rtl"
        assert rec.lcd == pytest.approx(4.4)
        assert rec.year == "2004"

    def test_lookup_coreid_missing(self, db):
        assert db.lookup_coreid("NONEXIST") is None

    def test_lookup_refcode(self, db):
        results = db.lookup_refcode("ABAVIJ_ASR_pacman")
        assert len(results) == 1
        assert results[0].coreid == "2004[Co][rtl]3[ASR]1"

    def test_lookup_base_refcode(self, db):
        results = db.lookup_base_refcode("HKUSTX")
        assert len(results) == 2
        extensions = {r.extension for r in results}
        assert "All Solvent Removed" in extensions
        assert "with ion" in extensions

    def test_search_name(self, db):
        results = db.search_name("JOU-34")
        assert len(results) == 1
        assert results[0].name == "JOU-34"

    def test_search_name_hkust(self, db):
        results = db.search_name("HKUST")
        assert len(results) == 2

    def test_search_doi(self, db):
        results = db.search_doi("10.1039/b404485a")
        assert len(results) == 1
        assert results[0].metal_types == "Co"

    def test_search_metal_single(self, db):
        results = db.search_metal("Cu")
        assert len(results) == 3  # HKUSTX_ASR, HKUSTX_ION, MULTIM (multi-metal)
        coreids = {r.coreid for r in results}
        assert "2020[Cu][pcu]3[ASR]1" in coreids
        assert "2018[Mn,Cu][dia]3[ASR]1" in coreids

    def test_search_metal_no_false_positive(self, db):
        """Mn search should not return Cu-only MOFs."""
        results = db.search_metal("Mn")
        coreids = {r.coreid for r in results}
        assert "2016[Mn][tbo]3[ASR]1" in coreids
        assert "2018[Mn,Cu][dia]3[ASR]1" in coreids
        assert "2020[Cu][pcu]3[ASR]1" not in coreids

    def test_search_topology(self, db):
        results = db.search_topology("pcu")
        assert len(results) == 2

    def test_search_topology_all_nodes(self, db):
        results = db.search_topology("dia", nodes="all")
        assert len(results) == 1

    def test_search_kh_class(self, db):
        results = db.search_kh_class("strong")
        assert len(results) == 2

    def test_search_oms(self, db):
        results = db.search_oms(has_oms=True)
        assert len(results) == 3  # JOU-34, HKUSTX_ASR, HKUSTX_ION

    def test_search_auto_coreid(self, db):
        result = db.search("2004[Co][rtl]3[ASR]1")
        assert result.field == "coreid"
        assert result.n_matches == 1

    def test_search_auto_refcode(self, db):
        result = db.search("HKUSTX")
        assert result.field == "refcode"
        assert result.n_matches == 2

    def test_search_auto_doi(self, db):
        result = db.search("10.1039/b404485a")
        assert result.field == "doi"
        assert result.n_matches == 1

    def test_search_auto_name(self, db):
        result = db.search("JOU-34")
        assert result.field == "name"
        assert result.n_matches == 1

    def test_search_no_matches(self, db):
        result = db.search("nonexistent-mof")
        assert result.n_matches == 0

    def test_has_oms_field(self, db):
        rec = db.lookup_coreid("2021[Zn][pts]3[ASR]1")
        assert rec is not None
        assert rec.has_oms is True

        rec2 = db.lookup_coreid("2004[Co][rtl]3[ASR]1")
        assert rec2 is not None
        assert rec2.has_oms is False

    def test_unknown_values_are_none(self, db):
        rec = db.lookup_coreid("2018[Mn,Cu][dia]3[ASR]1")
        assert rec is not None
        assert rec.thermal_stability is None
        assert rec.water_stability is None

    def test_raw_field_populated(self, db):
        rec = db.lookup_coreid("2004[Co][rtl]3[ASR]1")
        assert rec is not None
        assert rec.raw.get("refcode") == "ABAVIJ_ASR_pacman"

    def test_sqlite_cache_reused(self, sample_csv):
        db1 = CoreMOFDatabase(sample_csv)
        assert db1.record_count() == 6
        db1.close()

        db_path = sample_csv.with_suffix(".db")
        assert db_path.exists()
        mtime = db_path.stat().st_mtime

        db2 = CoreMOFDatabase(sample_csv)
        assert db2.record_count() == 6
        db2.close()

        assert db_path.stat().st_mtime == mtime


# ---------------------------------------------------------------------------
# Property screening
# ---------------------------------------------------------------------------


class TestScreening:
    def test_screen_lcd_min(self, db):
        results = db.screen(lcd_min=8.0)
        assert all(r.lcd >= 8.0 for r in results)
        assert len(results) >= 2  # ABAYIO (11.4), HKUSTX_ASR (9.0), HKUSTX_ION (8.5)

    def test_screen_water_stability(self, db):
        results = db.screen(water_stability_min=0.7)
        assert all(r.water_stability is not None and r.water_stability >= 0.7 for r in results)
        assert len(results) >= 1

    def test_screen_metal(self, db):
        results = db.screen(metal="Co")
        assert len(results) == 1
        assert results[0].metal_types == "Co"

    def test_screen_topology(self, db):
        results = db.screen(topology="pcu")
        assert len(results) == 2

    def test_screen_has_oms(self, db):
        results = db.screen(has_oms=True)
        assert all(r.has_oms for r in results)

    def test_screen_combined(self, db):
        results = db.screen(
            lcd_min=7.0,
            has_oms=True,
            kh_class="strong",
        )
        assert len(results) >= 1
        for r in results:
            assert r.lcd >= 7.0
            assert r.has_oms
            assert r.kh_class == "strong"

    def test_screen_extension(self, db):
        results = db.screen(extension="with ion")
        assert len(results) == 1
        assert results[0].extension == "with ion"

    def test_screen_no_filters(self, db):
        results = db.screen(limit=100)
        assert len(results) == 6


# ---------------------------------------------------------------------------
# CSD bridge
# ---------------------------------------------------------------------------


class TestBridge:
    def test_csd_to_coremof(self, sample_csv):
        db = CoreMOFDatabase(sample_csv)
        results = csd_to_coremof("HKUSTX", db=db)
        assert len(results) == 2
        coreids = {r.coreid for r in results}
        assert "2020[Cu][pcu]3[ASR]1" in coreids
        assert "2020[Cu][pcu]3[ION]1" in coreids
        db.close()

    def test_csd_to_coremof_strips_suffixes(self, sample_csv):
        db = CoreMOFDatabase(sample_csv)
        # Pass a full CoreMOF refcode — should still find via base
        results = csd_to_coremof("HKUSTX_ASR_pacman", db=db)
        assert len(results) == 2
        db.close()

    def test_csd_to_coremof_no_match(self, sample_csv):
        db = CoreMOFDatabase(sample_csv)
        results = csd_to_coremof("NONEXIST", db=db)
        assert len(results) == 0
        db.close()


# ---------------------------------------------------------------------------
# get_database() singleton
# ---------------------------------------------------------------------------


class TestGetDatabase:
    def test_explicit_path(self, sample_csv):
        db = get_database(data_path=sample_csv)
        assert db.record_count() == 6
        db.close()

    def test_explicit_path_does_not_clobber_singleton(self, sample_csv, monkeypatch):
        """An explicit data_path must not overwrite the module singleton."""
        from mofforge.coremof import database as db_module

        sentinel = object()
        monkeypatch.setattr(db_module, "_db", sentinel)

        db = get_database(data_path=sample_csv)
        assert db is not sentinel
        # Singleton untouched by the explicit-path call.
        assert db_module._db is sentinel
        db.close()

    def test_missing_path_raises(self, monkeypatch):
        from mofforge.coremof import database as db_module
        from mofforge.utils.config import config

        monkeypatch.setattr(config, "coremof_data_path", None)
        monkeypatch.delenv("MOFFORGE_COREMOF_DATA_PATH", raising=False)
        monkeypatch.setattr(db_module, "_db", None)

        with pytest.raises(FileNotFoundError, match="CoRE MOF data path"):
            get_database()

    def test_env_var(self, sample_csv, monkeypatch):
        from mofforge.coremof import database as db_module
        from mofforge.utils.config import config

        monkeypatch.setattr(config, "coremof_data_path", None)
        monkeypatch.setenv("MOFFORGE_COREMOF_DATA_PATH", str(sample_csv))
        monkeypatch.setattr(db_module, "_db", None)

        db = get_database()
        assert db.record_count() == 6
        db.close()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


class TestCLI:
    def test_coremof_command(self, sample_csv):
        from click.testing import CliRunner

        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(
            main, ["coremof", "ABAVIJ", "--data-path", str(sample_csv)]
        )
        assert result.exit_code == 0
        assert "ABAVIJ" in result.output

    def test_coremof_bridge(self, sample_csv):
        from click.testing import CliRunner

        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(
            main,
            ["coremof", "HKUSTX", "--bridge", "--data-path", str(sample_csv)],
        )
        assert result.exit_code == 0
        assert "2 match" in result.output
        assert "2020[Cu][pcu]3[ASR]1" in result.output

    def test_coremof_no_match(self, sample_csv):
        from click.testing import CliRunner

        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(
            main, ["coremof", "ZZZNOPE", "--data-path", str(sample_csv)]
        )
        assert result.exit_code == 0
        assert "0 match" in result.output

    def test_coremof_missing_config(self, monkeypatch):
        from click.testing import CliRunner

        from mofforge.coremof import database as db_module
        from mofforge.cli import main
        from mofforge.utils.config import config

        monkeypatch.setattr(config, "coremof_data_path", None)
        monkeypatch.delenv("MOFFORGE_COREMOF_DATA_PATH", raising=False)
        monkeypatch.setattr(db_module, "_db", None)

        runner = CliRunner()
        result = runner.invoke(main, ["coremof", "TESTMOF"])
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# CSD fixture for bridge/lookup tests
# ---------------------------------------------------------------------------

# CSD TSV header (same 61 columns as test_csd.py)
_CSD_HEADER = (
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


def _csd_row(refcode, name_sys, name_common, formula, year, space_group, remarks):
    """Build a tab-separated CSD row with 61 fields."""
    fields = [""] * 61
    fields[0] = refcode
    fields[1] = "A.Test"
    fields[2] = "J.Test"
    fields[3] = "1"
    fields[4] = "100"
    fields[5] = year
    fields[6] = formula
    fields[7] = name_sys
    fields[8] = name_common
    fields[14] = space_group
    fields[38] = "0.05"
    fields[39] = "293"
    fields[58] = remarks
    return "\t".join(fields)


# CSD rows with refcodes matching the CoreMOF fixture (ABAVIJ, HKUSTX)
_CSD_ROWS = [
    _csd_row(
        "ABAVIJ", "poly[cobalt-test]", "Co-MOF-test",
        "(C8 H4 Co1 N2 O4)n", "2004", "P-1",
        '"DOI: 10.1039/b404485a',
    ),
    _csd_row(
        "HKUSTX", "catena-[Copper trimesate variant]", "HKUST-variant",
        "(C18 H6 Cu3 O12)n", "2020", "Fm-3m",
        '"DOI: 10.1126/test',
    ),
    _csd_row(
        "NOMATCH", "poly[no-match]", "NoMatch-MOF",
        "(C6 H12 Zr1 O4)n", "2022", "P21/c",
        '"DOI: 10.9999/nomatch',
    ),
]


@pytest.fixture
def sample_csd_tsv(tmp_path):
    """Create a small synthetic CSD TSV that shares refcodes with the CoreMOF fixture."""
    tsv = tmp_path / "test_csd.tab"
    content = _CSD_HEADER + "\n" + "\n".join(_CSD_ROWS) + "\n"
    tsv.write_text(content, encoding="utf-8")
    return tsv


@pytest.fixture
def csd_db(sample_csd_tsv):
    """Return a CSDDatabase backed by the synthetic CSD fixture."""
    from mofforge.csd.database import CSDDatabase
    database = CSDDatabase(sample_csd_tsv)
    yield database
    database.close()


# ---------------------------------------------------------------------------
# search_csd_name + BridgeResult
# ---------------------------------------------------------------------------


class TestSearchCsdName:
    def test_search_csd_name_with_matches(self, sample_csv, sample_csd_tsv):
        coremof_db = CoreMOFDatabase(sample_csv)
        from mofforge.csd.database import CSDDatabase
        csd_db = CSDDatabase(sample_csd_tsv)

        results = search_csd_name("HKUST", coremof_db=coremof_db, csd_db=csd_db)
        assert len(results) == 1  # one CSD match for "HKUST"
        assert results[0].csd_record.refcode == "HKUSTX"
        assert results[0].has_coremof
        assert len(results[0].coremof_records) == 2  # ASR + ION variants
        coreids = {r.coreid for r in results[0].coremof_records}
        assert "2020[Cu][pcu]3[ASR]1" in coreids
        assert "2020[Cu][pcu]3[ION]1" in coreids

        coremof_db.close()
        csd_db.close()

    def test_search_csd_name_partial_match(self, sample_csv, sample_csd_tsv):
        coremof_db = CoreMOFDatabase(sample_csv)
        from mofforge.csd.database import CSDDatabase
        csd_db = CSDDatabase(sample_csd_tsv)

        # "cobalt" matches ABAVIJ's systematic name
        results = search_csd_name("cobalt", coremof_db=coremof_db, csd_db=csd_db)
        assert len(results) == 1
        assert results[0].csd_record.refcode == "ABAVIJ"
        assert results[0].has_coremof
        assert len(results[0].coremof_records) == 1

        coremof_db.close()
        csd_db.close()

    def test_search_csd_name_no_coremof_match(self, sample_csv, sample_csd_tsv):
        coremof_db = CoreMOFDatabase(sample_csv)
        from mofforge.csd.database import CSDDatabase
        csd_db = CSDDatabase(sample_csd_tsv)

        results = search_csd_name("NoMatch", coremof_db=coremof_db, csd_db=csd_db)
        assert len(results) == 1
        assert results[0].csd_record.refcode == "NOMATCH"
        assert not results[0].has_coremof
        assert len(results[0].coremof_records) == 0

        coremof_db.close()
        csd_db.close()

    def test_search_csd_name_no_csd_match(self, sample_csv, sample_csd_tsv):
        coremof_db = CoreMOFDatabase(sample_csv)
        from mofforge.csd.database import CSDDatabase
        csd_db = CSDDatabase(sample_csd_tsv)

        results = search_csd_name("zzz-nonexistent", coremof_db=coremof_db, csd_db=csd_db)
        assert len(results) == 0

        coremof_db.close()
        csd_db.close()


class TestBridgeResultModel:
    def test_summary_with_matches(self):
        from mofforge.csd.models import CSDRecord

        csd_rec = CSDRecord(refcode="HKUST1", chemical_name_common="HKUST-1",
                            doi="10.1126/test", year="1999")
        coremof_recs = [
            CoreMOFRecord(coreid="2020[Cu]1", metal_types="Cu",
                          topology_single="pcu", extension="All Solvent Removed"),
        ]
        br = BridgeResult(csd_record=csd_rec, coremof_records=coremof_recs)
        s = br.summary()
        assert "HKUST1" in s
        assert "HKUST-1" in s
        assert "2020[Cu]1" in s
        assert "All Solvent Removed" in s

    def test_summary_no_matches(self):
        from mofforge.csd.models import CSDRecord

        csd_rec = CSDRecord(refcode="NOMATCH", year="2022")
        br = BridgeResult(csd_record=csd_rec, coremof_records=[])
        s = br.summary()
        assert "no CoreMOF match" in s

    def test_has_coremof(self):
        from mofforge.csd.models import CSDRecord

        br_yes = BridgeResult(
            csd_record=CSDRecord(refcode="X"),
            coremof_records=[CoreMOFRecord(coreid="Y")],
        )
        assert br_yes.has_coremof

        br_no = BridgeResult(csd_record=CSDRecord(refcode="X"))
        assert not br_no.has_coremof


# ---------------------------------------------------------------------------
# CLI: mofforge lookup
# ---------------------------------------------------------------------------


class TestLookupCLI:
    def test_lookup_command(self, sample_csv, sample_csd_tsv):
        from click.testing import CliRunner

        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "lookup", "HKUST",
                "--coremof-data-path", str(sample_csv),
                "--csd-data-path", str(sample_csd_tsv),
            ],
        )
        assert result.exit_code == 0
        assert "HKUSTX" in result.output
        assert "2020[Cu][pcu]3[ASR]1" in result.output
        assert "With CoreMOF matches" in result.output

    def test_lookup_no_coremof(self, sample_csv, sample_csd_tsv):
        from click.testing import CliRunner

        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "lookup", "NoMatch",
                "--coremof-data-path", str(sample_csv),
                "--csd-data-path", str(sample_csd_tsv),
            ],
        )
        assert result.exit_code == 0
        assert "Without CoreMOF matches" in result.output

    def test_lookup_no_csd_match(self, sample_csv, sample_csd_tsv):
        from click.testing import CliRunner

        from mofforge.cli import main

        runner = CliRunner()
        result = runner.invoke(
            main,
            [
                "lookup", "zzz-nonexistent",
                "--coremof-data-path", str(sample_csv),
                "--csd-data-path", str(sample_csd_tsv),
            ],
        )
        assert result.exit_code == 0
        assert "0 total CSD match" in result.output
