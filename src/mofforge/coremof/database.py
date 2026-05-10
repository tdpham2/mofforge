"""CoRE MOF dataset accessor with SQLite-backed search.

Parses a CoRE MOF CSV file on first access and builds a local SQLite
cache for fast subsequent lookups.  The SQLite file is placed next to the
source CSV (e.g. ``CR_data_CSD_modified_20250227.db``).
"""

from __future__ import annotations

import csv
import json
import logging
import os
import re
import sqlite3
from pathlib import Path

from mofforge.coremof.models import BridgeResult, CoreMOFRecord, CoreMOFSearchResult

logger = logging.getLogger("mofforge")

ZENODO_URL = "https://zenodo.org/records/14510695"

# ---------------------------------------------------------------------------
# CSV column name constants
# ---------------------------------------------------------------------------
_COL_COREID = "coreid"
_COL_REFCODE = "refcode"
_COL_NAME = "name"
_COL_MOFID_V1 = "mofid-v1"
_COL_MOFID_V2 = "mofid-v2"
_COL_LCD = "LCD (Å)"
_COL_PLD = "PLD (Å)"
_COL_LFPD = "LFPD (Å)"
_COL_DENSITY = "Density (g/cm3)"
_COL_ASA = "ASA (m2/g)"
_COL_PV = "PV (cm3/g)"
_COL_VF = "VF"
_COL_TOPO_SINGLE = "topology(SingleNodes)"
_COL_TOPO_ALL = "topology(AllNodes)"
_COL_CATENATION = "catenation"
_COL_STRUCT_DIM = "structure_dimension"
_COL_SPACEGROUP = "number_spacegroup"
_COL_HALL = "hall"
_COL_METALS = "Metal Types"
_COL_HAS_OMS = "Has OMS"
_COL_OMS_TYPES = "OMS Types"
_COL_CHARGE = "Charge"
_COL_THERMAL = "Thermal_stability (℃)"
_COL_SOLVENT = "Solvent_stability"
_COL_WATER = "Water_stability"
_COL_HC_300K = "Heat_capacity@300K (J/g/K)"
_COL_KH = "KH_Classes"
_COL_DOI = "DOI"
_COL_YEAR = "Year"
_COL_PUBLICATION = "Publication"
_COL_SOURCE = "Source"
_COL_EXTENSION = "Extension"
_COL_NATOMS = "natoms"

# Regex for auto-detecting query type
_COREID_RE = re.compile(r"^\d{4}\[")
_REFCODE_RE = re.compile(r"^[A-Z]{3,8}\d{0,2}$", re.IGNORECASE)
_ELEMENT_RE = re.compile(r"^[A-Z][a-z]?$")

# Known topology names for auto-detection (common ones)
_COMMON_TOPOLOGIES = frozenset({
    "pcu", "dia", "sra", "nbo", "rht", "tbo", "pts", "fcu", "bcu",
    "acs", "cds", "hxg", "kgd", "lvt", "ntt", "pyr", "qtz", "reo",
    "rtl", "she", "sod", "sql", "srs", "ths", "ubt", "unh",
})

_SCHEMA = """\
CREATE TABLE IF NOT EXISTS records (
    coreid TEXT PRIMARY KEY,
    refcode TEXT NOT NULL,
    base_refcode TEXT NOT NULL,
    name TEXT,
    mofid_v1 TEXT,
    mofid_v2 TEXT,

    lcd REAL,
    pld REAL,
    lfpd REAL,
    density REAL,
    asa REAL,
    pore_volume REAL,
    void_fraction REAL,
    topology_single TEXT,
    topology_all TEXT,
    catenation INTEGER,
    structure_dimension INTEGER,
    spacegroup_number INTEGER,
    hall TEXT,

    metal_types TEXT,
    has_oms INTEGER,
    oms_types TEXT,
    charge_method TEXT,

    thermal_stability REAL,
    solvent_stability REAL,
    water_stability REAL,
    heat_capacity_300k REAL,
    kh_class TEXT,

    doi TEXT,
    year TEXT,
    publication TEXT,
    source TEXT,
    extension TEXT,
    natoms INTEGER,

    raw_json TEXT
);

CREATE INDEX IF NOT EXISTS idx_refcode ON records(refcode);
CREATE INDEX IF NOT EXISTS idx_base_refcode ON records(base_refcode COLLATE NOCASE);
CREATE INDEX IF NOT EXISTS idx_name ON records(name COLLATE NOCASE);
CREATE INDEX IF NOT EXISTS idx_doi ON records(doi);
CREATE INDEX IF NOT EXISTS idx_topology_single ON records(topology_single);
CREATE INDEX IF NOT EXISTS idx_topology_all ON records(topology_all);
CREATE INDEX IF NOT EXISTS idx_lcd ON records(lcd);
CREATE INDEX IF NOT EXISTS idx_pld ON records(pld);
CREATE INDEX IF NOT EXISTS idx_water_stability ON records(water_stability);
CREATE INDEX IF NOT EXISTS idx_thermal_stability ON records(thermal_stability);
CREATE INDEX IF NOT EXISTS idx_has_oms ON records(has_oms);
CREATE INDEX IF NOT EXISTS idx_kh_class ON records(kh_class);
CREATE INDEX IF NOT EXISTS idx_metal_types ON records(metal_types);
CREATE INDEX IF NOT EXISTS idx_extension ON records(extension);

CREATE TABLE IF NOT EXISTS record_metals (
    coreid TEXT NOT NULL,
    metal TEXT NOT NULL,
    FOREIGN KEY (coreid) REFERENCES records(coreid)
);
CREATE INDEX IF NOT EXISTS idx_metal ON record_metals(metal);
CREATE INDEX IF NOT EXISTS idx_metal_coreid ON record_metals(coreid);
"""

_INSERT_RECORD = """\
INSERT OR REPLACE INTO records (
    coreid, refcode, base_refcode, name, mofid_v1, mofid_v2,
    lcd, pld, lfpd, density, asa, pore_volume, void_fraction,
    topology_single, topology_all, catenation, structure_dimension,
    spacegroup_number, hall,
    metal_types, has_oms, oms_types, charge_method,
    thermal_stability, solvent_stability, water_stability,
    heat_capacity_300k, kh_class,
    doi, year, publication, source, extension, natoms,
    raw_json
) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
"""

_INSERT_METAL = "INSERT INTO record_metals (coreid, metal) VALUES (?, ?)"


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------


def _safe_float(value: str | None) -> float | None:
    if not value or value.strip().lower() in ("", "unknown", "-", "n/a", "nan"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _safe_int(value: str | None) -> int | None:
    if not value or value.strip().lower() in ("", "unknown", "-", "n/a", "nan"):
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def _get(row: dict[str, str], key: str) -> str:
    return row.get(key, "")


def _derive_base_refcode(refcode: str) -> str:
    """Extract the base CSD refcode from a CoreMOF refcode.

    ``"ABAVIJ_ASR_pacman"`` -> ``"ABAVIJ"``
    """
    return refcode.split("_")[0] if refcode else ""


class CoreMOFDatabase:
    """CoRE MOF dataset accessor with SQLite-backed search.

    Parameters
    ----------
    data_path : Path
        Path to the CoRE MOF CSV file.
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
            csv_mtime = self._data_path.stat().st_mtime
            db_mtime = self._db_path.stat().st_mtime
            if csv_mtime > db_mtime:
                logger.info(
                    "CoreMOF CSV is newer than cache; rebuilding SQLite database."
                )
                need_build = True

        if need_build:
            self._build_db_from_csv()

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
            self._build_db_from_csv()
            self._conn = sqlite3.connect(str(self._db_path))
            self._conn.row_factory = sqlite3.Row

        return self._conn

    def _build_db_from_csv(self) -> None:
        """Parse the CSV file and populate a fresh SQLite database."""
        if not self._data_path.exists():
            raise FileNotFoundError(
                f"CoRE MOF data file not found: {self._data_path}\n"
                f"Download the dataset from: {ZENODO_URL}\n"
                "Then configure the path via mofforge.toml, "
                "MOFFORGE_COREMOF_DATA_PATH env var, or "
                "set_paths(coremof_data=...)."
            )

        logger.info("Building CoreMOF SQLite cache from %s ...", self._data_path)

        if self._db_path.exists():
            self._db_path.unlink()

        conn = sqlite3.connect(str(self._db_path))
        conn.executescript(_SCHEMA)

        with open(self._data_path, newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh)

            record_batch: list[tuple] = []
            metal_batch: list[tuple[str, str]] = []
            count = 0

            for row in reader:
                coreid = _get(row, _COL_COREID)
                if not coreid:
                    continue

                refcode = _get(row, _COL_REFCODE)
                base_refcode = _derive_base_refcode(refcode)
                metal_types = _get(row, _COL_METALS)
                has_oms = 1 if _get(row, _COL_HAS_OMS).strip().lower() == "yes" else 0
                raw_json = json.dumps(row, ensure_ascii=False)

                record_batch.append((
                    coreid,
                    refcode,
                    base_refcode,
                    _get(row, _COL_NAME),
                    _get(row, _COL_MOFID_V1),
                    _get(row, _COL_MOFID_V2),
                    _safe_float(_get(row, _COL_LCD)),
                    _safe_float(_get(row, _COL_PLD)),
                    _safe_float(_get(row, _COL_LFPD)),
                    _safe_float(_get(row, _COL_DENSITY)),
                    _safe_float(_get(row, _COL_ASA)),
                    _safe_float(_get(row, _COL_PV)),
                    _safe_float(_get(row, _COL_VF)),
                    _get(row, _COL_TOPO_SINGLE),
                    _get(row, _COL_TOPO_ALL),
                    _safe_int(_get(row, _COL_CATENATION)),
                    _safe_int(_get(row, _COL_STRUCT_DIM)),
                    _safe_int(_get(row, _COL_SPACEGROUP)),
                    _get(row, _COL_HALL),
                    metal_types,
                    has_oms,
                    _get(row, _COL_OMS_TYPES),
                    _get(row, _COL_CHARGE),
                    _safe_float(_get(row, _COL_THERMAL)),
                    _safe_float(_get(row, _COL_SOLVENT)),
                    _safe_float(_get(row, _COL_WATER)),
                    _safe_float(_get(row, _COL_HC_300K)),
                    _get(row, _COL_KH),
                    _get(row, _COL_DOI),
                    _get(row, _COL_YEAR),
                    _get(row, _COL_PUBLICATION),
                    _get(row, _COL_SOURCE),
                    _get(row, _COL_EXTENSION),
                    _safe_int(_get(row, _COL_NATOMS)),
                    raw_json,
                ))

                # Normalize metals into join table
                if metal_types:
                    for metal in metal_types.split(","):
                        metal = metal.strip()
                        if metal:
                            metal_batch.append((coreid, metal))

                count += 1
                if len(record_batch) >= 5000:
                    conn.executemany(_INSERT_RECORD, record_batch)
                    conn.executemany(_INSERT_METAL, metal_batch)
                    record_batch.clear()
                    metal_batch.clear()

            if record_batch:
                conn.executemany(_INSERT_RECORD, record_batch)
            if metal_batch:
                conn.executemany(_INSERT_METAL, metal_batch)

        conn.commit()
        logger.info("CoreMOF cache built: %d records indexed.", count)
        conn.close()

    def _row_to_record(self, row: sqlite3.Row) -> CoreMOFRecord:
        """Convert a sqlite3.Row to a CoreMOFRecord."""
        raw = json.loads(row["raw_json"]) if row["raw_json"] else {}
        return CoreMOFRecord(
            coreid=row["coreid"] or "",
            refcode=row["refcode"] or "",
            base_refcode=row["base_refcode"] or "",
            name=row["name"] or "",
            mofid_v1=row["mofid_v1"] or "",
            mofid_v2=row["mofid_v2"] or "",
            lcd=row["lcd"],
            pld=row["pld"],
            lfpd=row["lfpd"],
            density=row["density"],
            asa=row["asa"],
            pore_volume=row["pore_volume"],
            void_fraction=row["void_fraction"],
            topology_single=row["topology_single"] or "",
            topology_all=row["topology_all"] or "",
            catenation=row["catenation"],
            structure_dimension=row["structure_dimension"],
            spacegroup_number=row["spacegroup_number"],
            hall=row["hall"] or "",
            metal_types=row["metal_types"] or "",
            has_oms=bool(row["has_oms"]),
            oms_types=row["oms_types"] or "",
            charge_method=row["charge_method"] or "",
            thermal_stability=row["thermal_stability"],
            solvent_stability=row["solvent_stability"],
            water_stability=row["water_stability"],
            heat_capacity_300k=row["heat_capacity_300k"],
            kh_class=row["kh_class"] or "",
            doi=row["doi"] or "",
            year=row["year"] or "",
            publication=row["publication"] or "",
            source=row["source"] or "",
            extension=row["extension"] or "",
            natoms=row["natoms"],
            raw=raw,
        )

    def _query(
        self, sql: str, params: tuple = (), limit: int | None = None
    ) -> list[CoreMOFRecord]:
        """Execute a query and return CoreMOFRecord list."""
        conn = self._ensure_db()
        if limit is not None:
            cursor = conn.execute(sql + " LIMIT ?", (*params, limit))
        else:
            cursor = conn.execute(sql, params)
        return [self._row_to_record(r) for r in cursor.fetchall()]

    # ------------------------------------------------------------------
    # Public search methods
    # ------------------------------------------------------------------

    def lookup_coreid(self, coreid: str) -> CoreMOFRecord | None:
        """Exact coreid lookup."""
        conn = self._ensure_db()
        cursor = conn.execute(
            "SELECT * FROM records WHERE coreid = ?", (coreid,)
        )
        row = cursor.fetchone()
        return self._row_to_record(row) if row else None

    def lookup_refcode(self, refcode: str) -> list[CoreMOFRecord]:
        """Exact full CoreMOF refcode lookup (case-insensitive)."""
        return self._query(
            "SELECT * FROM records WHERE refcode = ? COLLATE NOCASE",
            (refcode,),
        )

    def lookup_base_refcode(self, base_refcode: str) -> list[CoreMOFRecord]:
        """Find all CoreMOF entries for a CSD-style base refcode.

        This is the primary CSD bridge method. Given ``"ABAVIJ"``, returns
        all processed variants (e.g. ``ABAVIJ_ASR_pacman``,
        ``ABAVIJ_ION_pacman``).
        """
        return self._query(
            "SELECT * FROM records WHERE base_refcode = ? COLLATE NOCASE",
            (base_refcode,),
        )

    def search_name(self, name: str, limit: int | None = None) -> list[CoreMOFRecord]:
        """Substring search on MOF name."""
        pattern = f"%{name}%"
        return self._query(
            "SELECT * FROM records WHERE name LIKE ? COLLATE NOCASE",
            (pattern,),
            limit=limit,
        )

    def search_doi(self, doi: str, limit: int | None = None) -> list[CoreMOFRecord]:
        """Exact or partial DOI match."""
        pattern = f"%{doi}%"
        return self._query(
            "SELECT * FROM records WHERE doi LIKE ?",
            (pattern,),
            limit=limit,
        )

    def search_metal(self, metal: str, limit: int | None = None) -> list[CoreMOFRecord]:
        """Find MOFs containing a specific metal element.

        Uses the normalized ``record_metals`` table for accurate matching
        (avoids substring false positives).
        """
        conn = self._ensure_db()
        sql = (
            "SELECT r.* FROM records r "
            "JOIN record_metals m ON r.coreid = m.coreid "
            "WHERE m.metal = ? COLLATE NOCASE"
        )
        if limit is not None:
            cursor = conn.execute(sql + " LIMIT ?", (metal.strip(), limit))
        else:
            cursor = conn.execute(sql, (metal.strip(),))
        return [self._row_to_record(r) for r in cursor.fetchall()]

    def search_topology(
        self, topology: str, nodes: str = "single", limit: int | None = None
    ) -> list[CoreMOFRecord]:
        """Find MOFs with a specific topology name.

        Parameters
        ----------
        topology : str
            Topology name (e.g. ``"pcu"``, ``"dia"``).
        nodes : str
            ``"single"`` for topology(SingleNodes) or ``"all"`` for
            topology(AllNodes).
        """
        col = "topology_single" if nodes == "single" else "topology_all"
        return self._query(
            f"SELECT * FROM records WHERE {col} = ? COLLATE NOCASE",
            (topology,),
            limit=limit,
        )

    def search_kh_class(self, kh_class: str, limit: int | None = None) -> list[CoreMOFRecord]:
        """Find MOFs by KH gas storage classification."""
        return self._query(
            "SELECT * FROM records WHERE kh_class = ? COLLATE NOCASE",
            (kh_class,),
            limit=limit,
        )

    def search_oms(
        self, has_oms: bool = True, limit: int | None = None
    ) -> list[CoreMOFRecord]:
        """Find MOFs with or without open metal sites."""
        return self._query(
            "SELECT * FROM records WHERE has_oms = ?",
            (1 if has_oms else 0,),
            limit=limit,
        )

    def screen(
        self,
        *,
        lcd_min: float | None = None,
        lcd_max: float | None = None,
        pld_min: float | None = None,
        pld_max: float | None = None,
        density_min: float | None = None,
        density_max: float | None = None,
        asa_min: float | None = None,
        asa_max: float | None = None,
        pore_volume_min: float | None = None,
        pore_volume_max: float | None = None,
        void_fraction_min: float | None = None,
        void_fraction_max: float | None = None,
        water_stability_min: float | None = None,
        solvent_stability_min: float | None = None,
        thermal_stability_min: float | None = None,
        metal: str | None = None,
        topology: str | None = None,
        has_oms: bool | None = None,
        kh_class: str | None = None,
        extension: str | None = None,
        limit: int | None = None,
    ) -> list[CoreMOFRecord]:
        """Screen MOFs by property ranges and categorical filters.

        All parameters are optional — only non-None filters are applied.
        Numeric filters exclude records with NULL values for that property.
        """
        conditions: list[str] = []
        params: list[float | int | str] = []

        range_filters = [
            ("lcd", lcd_min, lcd_max),
            ("pld", pld_min, pld_max),
            ("density", density_min, density_max),
            ("asa", asa_min, asa_max),
            ("pore_volume", pore_volume_min, pore_volume_max),
            ("void_fraction", void_fraction_min, void_fraction_max),
        ]
        for col, lo, hi in range_filters:
            if lo is not None:
                conditions.append(f"{col} >= ?")
                params.append(lo)
            if hi is not None:
                conditions.append(f"{col} <= ?")
                params.append(hi)

        min_only_filters = [
            ("water_stability", water_stability_min),
            ("solvent_stability", solvent_stability_min),
            ("thermal_stability", thermal_stability_min),
        ]
        for col, lo in min_only_filters:
            if lo is not None:
                conditions.append(f"{col} >= ?")
                params.append(lo)

        if has_oms is not None:
            conditions.append("has_oms = ?")
            params.append(1 if has_oms else 0)

        if kh_class is not None:
            conditions.append("kh_class = ? COLLATE NOCASE")
            params.append(kh_class)

        if extension is not None:
            conditions.append("extension = ? COLLATE NOCASE")
            params.append(extension)

        if topology is not None:
            conditions.append("topology_single = ? COLLATE NOCASE")
            params.append(topology)

        # Metal filter uses the join table
        use_metal_join = metal is not None

        if use_metal_join:
            base = (
                "SELECT r.* FROM records r "
                "JOIN record_metals rm ON r.coreid = rm.coreid "
                "WHERE rm.metal = ? COLLATE NOCASE"
            )
            all_params: list[float | int | str] = [metal]  # type: ignore[list-item]
            if conditions:
                base += " AND " + " AND ".join(conditions)
                all_params.extend(params)
        else:
            if conditions:
                base = "SELECT * FROM records WHERE " + " AND ".join(conditions)
            else:
                base = "SELECT * FROM records"
            all_params = params

        return self._query(base, tuple(all_params), limit=limit)

    def search(
        self, query: str, field: str = "auto", limit: int | None = None
    ) -> CoreMOFSearchResult:
        """Unified search with auto-detection of query type.

        Parameters
        ----------
        query : str
            The search term.
        field : str
            One of ``"auto"``, ``"coreid"``, ``"refcode"``, ``"name"``,
            ``"doi"``, ``"metal"``, ``"topology"``.
        limit : int
            Maximum number of results.
        """
        query = query.strip()

        if field == "auto":
            field = self._detect_field(query)

        if field == "coreid":
            rec = self.lookup_coreid(query)
            records = [rec] if rec else []
        elif field == "refcode":
            records = self.lookup_base_refcode(query)
            if not records:
                # Try full refcode match
                records = self.lookup_refcode(query)
            if not records:
                # Fallback to name search
                field = "name"
                records = self.search_name(query, limit=limit)
        elif field == "name":
            records = self.search_name(query, limit=limit)
        elif field == "doi":
            records = self.search_doi(query, limit=limit)
        elif field == "metal":
            records = self.search_metal(query, limit=limit)
        elif field == "topology":
            records = self.search_topology(query, limit=limit)
        else:
            raise ValueError(f"Unknown search field: {field!r}")

        return CoreMOFSearchResult(query=query, field=field, records=records)

    @staticmethod
    def _detect_field(query: str) -> str:
        """Guess which field to search based on the query string."""
        if _COREID_RE.match(query):
            return "coreid"
        if query.startswith("10."):
            return "doi"
        # Check topology before refcode — short lowercase strings like "pcu"
        # would otherwise match the refcode regex.
        if query.lower() in _COMMON_TOPOLOGIES:
            return "topology"
        if _ELEMENT_RE.match(query):
            return "metal"
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

_db: CoreMOFDatabase | None = None


def get_database(data_path: str | Path | None = None) -> CoreMOFDatabase:
    """Return the lazily-initialized CoRE MOF database singleton.

    The data path is resolved in order of priority:

    1. Explicit *data_path* argument.
    2. ``config.coremof_data_path`` (set via :func:`set_paths`).
    3. ``MOFFORGE_COREMOF_DATA_PATH`` environment variable.
    4. ``[coremof] data_path`` in ``mofforge.toml``.

    Raises
    ------
    FileNotFoundError
        If no data path is configured or the file does not exist.
    """
    global _db

    if data_path is not None:
        _db = CoreMOFDatabase(Path(data_path))
        return _db

    if _db is not None:
        return _db

    from mofforge.utils.config import config

    resolved: Path | None = config.coremof_data_path

    if resolved is None:
        env = os.environ.get("MOFFORGE_COREMOF_DATA_PATH")
        if env:
            resolved = Path(env)

    if resolved is None:
        resolved = _load_coremof_path_from_toml()

    if resolved is None:
        raise FileNotFoundError(
            "CoRE MOF data path is not configured.\n"
            f"Download the dataset from: {ZENODO_URL}\n"
            "Then set the path via one of:\n"
            "  1. set_paths(coremof_data='/path/to/CR_data.csv')\n"
            "  2. MOFFORGE_COREMOF_DATA_PATH environment variable\n"
            "  3. [coremof] data_path in mofforge.toml"
        )

    _db = CoreMOFDatabase(resolved)
    return _db


def csd_to_coremof(
    refcode: str, db: CoreMOFDatabase | None = None
) -> list[CoreMOFRecord]:
    """Find CoreMOF entries for a CSD refcode.

    Given a CSD refcode like ``"ABAVIJ"``, returns all CoreMOF entries
    matching that base refcode (e.g. ``ABAVIJ_ASR_pacman``,
    ``ABAVIJ_ION_pacman``).
    """
    if db is None:
        db = get_database()
    base = _derive_base_refcode(refcode)
    return db.lookup_base_refcode(base)


def search_csd_name(
    name: str,
    coremof_db: CoreMOFDatabase | None = None,
    csd_db: "CSDDatabase | None" = None,
    limit: int | None = None,
) -> list[BridgeResult]:
    """Search for a MOF name in CSD and return CoreMOF entries for each match.

    Chains CSD name search with CoreMOF refcode bridge: for each CSD
    record found, looks up all corresponding CoreMOF entries via the
    base refcode.

    Parameters
    ----------
    name : str
        MOF name to search in CSD (substring match).
    coremof_db : CoreMOFDatabase, optional
        CoreMOF database instance. Uses singleton if not provided.
    csd_db : CSDDatabase, optional
        CSD database instance. Uses singleton if not provided.
    limit : int
        Maximum number of CSD records to return.

    Returns
    -------
    list[BridgeResult]
        Each entry pairs a CSD record with its CoreMOF matches (may be empty).
    """
    from mofforge.csd import get_database as get_csd_database

    if csd_db is None:
        csd_db = get_csd_database()
    if coremof_db is None:
        coremof_db = get_database()

    csd_records = csd_db.search_name(name, limit=limit)
    results: list[BridgeResult] = []
    for csd_rec in csd_records:
        coremof_recs = coremof_db.lookup_base_refcode(csd_rec.refcode)
        results.append(BridgeResult(csd_record=csd_rec, coremof_records=coremof_recs))
    return results


def _load_coremof_path_from_toml() -> Path | None:
    """Try to read ``[coremof] data_path`` from ``mofforge.toml``."""
    from mofforge.build.config import _find_toml, _load_toml

    toml_path = _find_toml()
    if toml_path is None:
        return None
    data = _load_toml(toml_path)
    section = data.get("coremof", {})
    raw = section.get("data_path")
    return Path(raw) if raw else None
