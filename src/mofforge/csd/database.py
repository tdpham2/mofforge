"""CSD dataset accessor with SQLite-backed search.

Parses a CSD MOF subset TSV file on first access and builds a local SQLite
cache for fast subsequent lookups.  The SQLite file is placed next to the
source TSV (e.g. ``MOF_subset.db`` alongside ``MOF_subset.tab``).
"""

from __future__ import annotations

import csv
import json
import logging
import os
import re
import sqlite3
from pathlib import Path

from mofforge.csd.models import CSDRecord, CSDSearchResult

logger = logging.getLogger("mofforge")

# ---------------------------------------------------------------------------
# Column name constants (brackets stripped during header parsing)
# ---------------------------------------------------------------------------
_COL_REFCODE = "REFCODE"
_COL_AUTHORS = "_publ_authors"
_COL_JOURNAL = "_journal_name_full"
_COL_VOLUME = "_journal_volume"
_COL_PAGES = "_journal_pages"
_COL_YEAR = "_journal_year"
_COL_FORMULA = "_chemical_formula_moiety"
_COL_NAME_SYS = "_chemical_name_systematic"
_COL_NAME_COMMON = "_chemical_name_common"
_COL_SPACE_GROUP = "_symmetry_space_group_name"
_COL_CELL_A = "_cell_length_a"
_COL_CELL_B = "_cell_length_b"
_COL_CELL_C = "_cell_length_c"
_COL_CELL_ALPHA = "_cell_angle_alpha"
_COL_CELL_BETA = "_cell_angle_beta"
_COL_CELL_GAMMA = "_cell_angle_gamma"
_COL_CELL_VOL = "_cell_volume"
_COL_TEMP = "_cell_measurement_temperature"
_COL_R_FACTOR = "_refine_ls_R_factor"
_COL_REMARKS = "_ccdc_RMARKS"

# Regex patterns for extracting DOI and CCDC number from _ccdc_RMARKS
_DOI_RE = re.compile(r"DOI:\s*(10\.\S+)")
_CCDC_RE = re.compile(r"CCDC[\s:]+(\d{5,8})")

# Regex for auto-detecting query type
_REFCODE_RE = re.compile(r"^[A-Z]{3,8}\d{0,2}$", re.IGNORECASE)

_SCHEMA = """\
CREATE TABLE IF NOT EXISTS records (
    refcode TEXT PRIMARY KEY,
    chemical_name_systematic TEXT,
    chemical_name_common TEXT,
    chemical_formula_moiety TEXT,
    doi TEXT,
    ccdc_number TEXT,
    authors TEXT,
    journal TEXT,
    volume TEXT,
    pages TEXT,
    year TEXT,
    space_group TEXT,
    cell_a TEXT,
    cell_b TEXT,
    cell_c TEXT,
    cell_alpha TEXT,
    cell_beta TEXT,
    cell_gamma TEXT,
    cell_volume TEXT,
    temperature TEXT,
    r_factor TEXT,
    raw_json TEXT
);

CREATE INDEX IF NOT EXISTS idx_name_sys ON records(chemical_name_systematic COLLATE NOCASE);
CREATE INDEX IF NOT EXISTS idx_name_common ON records(chemical_name_common COLLATE NOCASE);
CREATE INDEX IF NOT EXISTS idx_doi ON records(doi);
CREATE INDEX IF NOT EXISTS idx_ccdc ON records(ccdc_number);
CREATE INDEX IF NOT EXISTS idx_formula ON records(chemical_formula_moiety COLLATE NOCASE);
CREATE INDEX IF NOT EXISTS idx_year ON records(year);
"""

_INSERT = """\
INSERT OR REPLACE INTO records (
    refcode, chemical_name_systematic, chemical_name_common,
    chemical_formula_moiety, doi, ccdc_number, authors, journal,
    volume, pages, year, space_group,
    cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma,
    cell_volume, temperature, r_factor, raw_json
) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
"""


def _extract_doi(remarks: str | None) -> str | None:
    """Extract the first DOI from a _ccdc_RMARKS field."""
    if not remarks:
        return None
    m = _DOI_RE.search(remarks)
    return m.group(1).rstrip('"').rstrip() if m else None


def _extract_ccdc(remarks: str | None) -> str | None:
    """Extract the first CCDC deposition number from a _ccdc_RMARKS field."""
    if not remarks:
        return None
    m = _CCDC_RE.search(remarks)
    return m.group(1) if m else None


def _strip_brackets(name: str) -> str:
    """Strip ``[`` and ``]`` from a column header name."""
    return name.strip().removeprefix("[").removesuffix("]")


def _get(row: dict[str, str], key: str) -> str:
    """Get a value from a row dict, returning empty string if missing."""
    return row.get(key, "")


class CSDDatabase:
    """CSD dataset accessor with SQLite-backed search.

    Parameters
    ----------
    data_path : Path
        Path to the CSD MOF subset TSV file (e.g. ``MOF_subset.tab``).
    """

    def __init__(self, data_path: Path) -> None:
        self._data_path = Path(data_path)
        self._db_path = self._data_path.with_suffix(".db")
        self._conn: sqlite3.Connection | None = None

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _ensure_db(self) -> sqlite3.Connection:
        """Open or (re)build the SQLite cache and return a connection."""
        if self._conn is not None:
            return self._conn

        need_build = not self._db_path.exists()
        if not need_build:
            # Rebuild if the source TSV is newer than the cache.
            tsv_mtime = self._data_path.stat().st_mtime
            db_mtime = self._db_path.stat().st_mtime
            if tsv_mtime > db_mtime:
                logger.info("CSD TSV is newer than cache; rebuilding SQLite database.")
                need_build = True

        if need_build:
            self._build_db_from_tsv()

        self._conn = sqlite3.connect(str(self._db_path))
        self._conn.row_factory = sqlite3.Row

        # Guard against a stale empty DB from a previously failed build.
        try:
            count = self._conn.execute("SELECT COUNT(*) FROM records").fetchone()[0]
        except sqlite3.OperationalError:
            count = 0
        if count == 0:
            self._conn.close()
            self._conn = None
            self._build_db_from_tsv()
            self._conn = sqlite3.connect(str(self._db_path))
            self._conn.row_factory = sqlite3.Row

        return self._conn

    def _build_db_from_tsv(self) -> None:
        """Parse the TSV file and populate a fresh SQLite database."""
        if not self._data_path.exists():
            raise FileNotFoundError(
                f"CSD data file not found: {self._data_path}\n"
                "Configure the path via mofforge.toml, MOFFORGE_CSD_DATA_PATH "
                "env var, or set_paths(csd_data=...)."
            )

        logger.info("Building CSD SQLite cache from %s ...", self._data_path)

        # Remove stale db before rebuilding
        if self._db_path.exists():
            self._db_path.unlink()

        conn = sqlite3.connect(str(self._db_path))
        conn.executescript(_SCHEMA)

        # Try UTF-8 first, fall back to latin-1.
        for encoding in ("utf-8", "latin-1"):
            try:
                with open(self._data_path, newline="", encoding=encoding) as fh:
                    # Disable quoting — CSD fields use embedded quotes as
                    # sub-field delimiters, not CSV-style quoting.
                    reader = csv.DictReader(
                        fh, delimiter="\t", quoting=csv.QUOTE_NONE,
                    )
                    # Strip brackets from header names
                    if reader.fieldnames:
                        reader.fieldnames = [_strip_brackets(f) for f in reader.fieldnames]

                    batch: list[tuple] = []  # type: ignore[type-arg]
                    count = 0
                    for row in reader:
                        remarks = _get(row, _COL_REMARKS)
                        raw_json = json.dumps(row, ensure_ascii=False)
                        batch.append((
                            _get(row, _COL_REFCODE),
                            _get(row, _COL_NAME_SYS),
                            _get(row, _COL_NAME_COMMON),
                            _get(row, _COL_FORMULA),
                            _extract_doi(remarks),
                            _extract_ccdc(remarks),
                            _get(row, _COL_AUTHORS),
                            _get(row, _COL_JOURNAL),
                            _get(row, _COL_VOLUME),
                            _get(row, _COL_PAGES),
                            _get(row, _COL_YEAR),
                            _get(row, _COL_SPACE_GROUP),
                            _get(row, _COL_CELL_A),
                            _get(row, _COL_CELL_B),
                            _get(row, _COL_CELL_C),
                            _get(row, _COL_CELL_ALPHA),
                            _get(row, _COL_CELL_BETA),
                            _get(row, _COL_CELL_GAMMA),
                            _get(row, _COL_CELL_VOL),
                            _get(row, _COL_TEMP),
                            _get(row, _COL_R_FACTOR),
                            raw_json,
                        ))
                        count += 1
                        if len(batch) >= 5000:
                            conn.executemany(_INSERT, batch)
                            batch.clear()

                    if batch:
                        conn.executemany(_INSERT, batch)

                conn.commit()
                logger.info("CSD cache built: %d records indexed.", count)
                conn.close()
                return  # success
            except UnicodeDecodeError:
                conn.close()
                if self._db_path.exists():
                    self._db_path.unlink()
                if encoding == "latin-1":
                    raise
                logger.debug("UTF-8 decode failed, retrying with latin-1.")

    def _row_to_record(self, row: sqlite3.Row) -> CSDRecord:
        """Convert a sqlite3.Row to a CSDRecord."""
        raw = json.loads(row["raw_json"]) if row["raw_json"] else {}
        return CSDRecord(
            refcode=row["refcode"] or "",
            chemical_name_systematic=row["chemical_name_systematic"] or "",
            chemical_name_common=row["chemical_name_common"] or "",
            chemical_formula_moiety=row["chemical_formula_moiety"] or "",
            doi=row["doi"],
            ccdc_number=row["ccdc_number"],
            authors=row["authors"] or "",
            journal=row["journal"] or "",
            volume=row["volume"] or "",
            pages=row["pages"] or "",
            year=row["year"] or "",
            space_group=row["space_group"] or "",
            cell_a=row["cell_a"] or "",
            cell_b=row["cell_b"] or "",
            cell_c=row["cell_c"] or "",
            cell_alpha=row["cell_alpha"] or "",
            cell_beta=row["cell_beta"] or "",
            cell_gamma=row["cell_gamma"] or "",
            cell_volume=row["cell_volume"] or "",
            temperature=row["temperature"] or "",
            r_factor=row["r_factor"] or "",
            raw=raw,
        )

    def _query(self, sql: str, params: tuple = (), limit: int = 50) -> list[CSDRecord]:
        """Execute a query and return CSDRecord list."""
        conn = self._ensure_db()
        cursor = conn.execute(sql + " LIMIT ?", (*params, limit))
        return [self._row_to_record(r) for r in cursor.fetchall()]

    # ------------------------------------------------------------------
    # Public search methods
    # ------------------------------------------------------------------

    def lookup_refcode(self, refcode: str) -> CSDRecord | None:
        """Exact refcode lookup (case-insensitive)."""
        conn = self._ensure_db()
        cursor = conn.execute(
            "SELECT * FROM records WHERE refcode = ? COLLATE NOCASE",
            (refcode,),
        )
        row = cursor.fetchone()
        return self._row_to_record(row) if row else None

    def search_name(self, name: str, limit: int = 50) -> list[CSDRecord]:
        """Substring search across systematic and common chemical names."""
        pattern = f"%{name}%"
        return self._query(
            "SELECT * FROM records WHERE chemical_name_systematic LIKE ? COLLATE NOCASE "
            "OR chemical_name_common LIKE ? COLLATE NOCASE",
            (pattern, pattern),
            limit=limit,
        )

    def search_doi(self, doi: str, limit: int = 50) -> list[CSDRecord]:
        """Exact or partial DOI match."""
        pattern = f"%{doi}%"
        return self._query(
            "SELECT * FROM records WHERE doi LIKE ?",
            (pattern,),
            limit=limit,
        )

    def search_formula(self, formula: str, limit: int = 50) -> list[CSDRecord]:
        """Substring match on chemical formula."""
        pattern = f"%{formula}%"
        return self._query(
            "SELECT * FROM records WHERE chemical_formula_moiety LIKE ? COLLATE NOCASE",
            (pattern,),
            limit=limit,
        )

    def search_ccdc(self, ccdc_number: str, limit: int = 50) -> list[CSDRecord]:
        """Lookup by CCDC deposition number."""
        return self._query(
            "SELECT * FROM records WHERE ccdc_number = ?",
            (ccdc_number,),
            limit=limit,
        )

    def search(self, query: str, field: str = "auto", limit: int = 50) -> CSDSearchResult:
        """Unified search with auto-detection of query type.

        Parameters
        ----------
        query : str
            The search term.
        field : str
            One of ``"auto"``, ``"refcode"``, ``"name"``, ``"doi"``,
            ``"formula"``, ``"ccdc"``.
        limit : int
            Maximum number of results.
        """
        query = query.strip()

        if field == "auto":
            field = self._detect_field(query)

        if field == "refcode":
            rec = self.lookup_refcode(query)
            if rec:
                records = [rec]
            else:
                # Fallback: the query looks like a refcode but wasn't found;
                # try a name search instead (e.g. "HKUST", "MOF-5").
                field = "name"
                records = self.search_name(query, limit=limit)
        elif field == "name":
            records = self.search_name(query, limit=limit)
        elif field == "doi":
            records = self.search_doi(query, limit=limit)
        elif field == "formula":
            records = self.search_formula(query, limit=limit)
        elif field == "ccdc":
            records = self.search_ccdc(query, limit=limit)
        else:
            raise ValueError(f"Unknown search field: {field!r}")

        return CSDSearchResult(query=query, field=field, records=records)

    @staticmethod
    def _detect_field(query: str) -> str:
        """Guess which field to search based on the query string."""
        if query.startswith("10."):
            return "doi"
        if query.isdigit() and 5 <= len(query) <= 8:
            return "ccdc"
        if _REFCODE_RE.match(query):
            return "refcode"
        return "name"

    def record_count(self) -> int:
        """Return total number of records in the database."""
        conn = self._ensure_db()
        cursor = conn.execute("SELECT COUNT(*) FROM records")
        return cursor.fetchone()[0]

    def close(self) -> None:
        """Close the database connection."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None


# ---------------------------------------------------------------------------
# Module-level lazy accessor
# ---------------------------------------------------------------------------

_db: CSDDatabase | None = None


def get_database(data_path: str | Path | None = None) -> CSDDatabase:
    """Return the lazily-initialized CSD database singleton.

    The data path is resolved in order of priority:

    1. Explicit *data_path* argument.
    2. ``config.csd_data_path`` (set via :func:`set_paths`).
    3. ``MOFFORGE_CSD_DATA_PATH`` environment variable.
    4. ``[csd] data_path`` in ``mofforge.toml``.

    Raises
    ------
    FileNotFoundError
        If no data path is configured or the file does not exist.
    """
    global _db

    if data_path is not None:
        # Explicit path always creates a (new) database instance.
        _db = CSDDatabase(Path(data_path))
        return _db

    if _db is not None:
        return _db

    from mofforge.utils.config import config

    resolved: Path | None = config.csd_data_path

    # Try environment variable.
    if resolved is None:
        env = os.environ.get("MOFFORGE_CSD_DATA_PATH")
        if env:
            resolved = Path(env)

    # Try mofforge.toml.
    if resolved is None:
        resolved = _load_csd_path_from_toml()

    if resolved is None:
        raise FileNotFoundError(
            "CSD data path is not configured.\n"
            "Set it via one of:\n"
            "  1. set_paths(csd_data='/path/to/MOF_subset.tab')\n"
            "  2. MOFFORGE_CSD_DATA_PATH environment variable\n"
            "  3. [csd] data_path in mofforge.toml"
        )

    _db = CSDDatabase(resolved)
    return _db


def _load_csd_path_from_toml() -> Path | None:
    """Try to read ``[csd] data_path`` from ``mofforge.toml``."""
    from mofforge.build.config import _find_toml, _load_toml

    toml_path = _find_toml()
    if toml_path is None:
        return None
    data = _load_toml(toml_path)
    csd_section = data.get("csd", {})
    raw = csd_section.get("data_path")
    return Path(raw) if raw else None
