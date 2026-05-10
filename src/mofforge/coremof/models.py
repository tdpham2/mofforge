"""Data models for CoRE MOF database lookup results."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from mofforge.csd.models import CSDRecord


@dataclass
class CoreMOFRecord:
    """A single CoRE MOF database entry."""

    # Identifiers
    coreid: str
    refcode: str = ""
    base_refcode: str = ""
    name: str = ""
    mofid_v1: str = ""
    mofid_v2: str = ""

    # Structural properties
    lcd: float | None = None
    pld: float | None = None
    lfpd: float | None = None
    density: float | None = None
    asa: float | None = None
    pore_volume: float | None = None
    void_fraction: float | None = None
    topology_single: str = ""
    topology_all: str = ""
    catenation: int | None = None
    structure_dimension: int | None = None
    spacegroup_number: int | None = None
    hall: str = ""

    # Chemistry
    metal_types: str = ""
    has_oms: bool = False
    oms_types: str = ""
    charge_method: str = ""

    # Properties
    thermal_stability: float | None = None
    solvent_stability: float | None = None
    water_stability: float | None = None
    heat_capacity_300k: float | None = None
    kh_class: str = ""

    # Metadata
    doi: str = ""
    year: str = ""
    publication: str = ""
    source: str = ""
    extension: str = ""
    natoms: int | None = None

    raw: dict[str, str] = field(default_factory=dict, repr=False)

    def summary(self) -> str:
        """Human-readable one-line summary."""
        parts = [self.coreid]
        if self.name and self.name != "-":
            parts.append(self.name[:40])
        if self.metal_types:
            parts.append(f"[{self.metal_types}]")
        if self.topology_single and self.topology_single not in ("nan", "unnamed"):
            parts.append(self.topology_single)
        if self.doi:
            parts.append(f"DOI:{self.doi}")
        if self.year:
            parts.append(f"({self.year})")
        return " | ".join(parts)

    def properties_summary(self) -> str:
        """Multi-line summary of key properties."""
        lines = [f"CoreMOF: {self.coreid}"]
        lines.append(f"  Refcode:    {self.refcode} (base: {self.base_refcode})")
        if self.name and self.name != "-":
            lines.append(f"  Name:       {self.name}")
        lines.append(f"  Metals:     {self.metal_types}")
        lines.append(f"  Extension:  {self.extension}")

        if self.lcd is not None:
            lines.append(f"  LCD:        {self.lcd:.3f} A")
        if self.pld is not None:
            lines.append(f"  PLD:        {self.pld:.3f} A")
        if self.density is not None:
            lines.append(f"  Density:    {self.density:.4f} g/cm3")
        if self.asa is not None:
            lines.append(f"  ASA:        {self.asa:.2f} m2/g")
        if self.pore_volume is not None:
            lines.append(f"  Pore Vol:   {self.pore_volume:.4f} cm3/g")
        if self.void_fraction is not None:
            lines.append(f"  Void Frac:  {self.void_fraction:.4f}")

        topo = self.topology_single
        if topo and topo not in ("nan", "unnamed"):
            lines.append(f"  Topology:   {topo}")
        lines.append(f"  Has OMS:    {'Yes' if self.has_oms else 'No'}")

        if self.thermal_stability is not None:
            lines.append(f"  T_stab:     {self.thermal_stability:.1f} C")
        if self.water_stability is not None:
            lines.append(f"  H2O_stab:   {self.water_stability:.3f}")
        if self.solvent_stability is not None:
            lines.append(f"  Solv_stab:  {self.solvent_stability:.3f}")
        if self.kh_class:
            lines.append(f"  KH Class:   {self.kh_class}")

        if self.doi:
            lines.append(f"  DOI:        {self.doi}")
        if self.year:
            lines.append(f"  Year:       {self.year}")
        return "\n".join(lines)


@dataclass
class BridgeResult:
    """A CSD record paired with its CoreMOF matches."""

    csd_record: CSDRecord
    coremof_records: list[CoreMOFRecord] = field(default_factory=list)

    @property
    def has_coremof(self) -> bool:
        return len(self.coremof_records) > 0

    def summary(self) -> str:
        """Multi-line summary showing CSD entry and CoreMOF matches."""
        csd = self.csd_record
        name = csd.chemical_name_common or csd.chemical_name_systematic
        parts = [f"CSD: {csd.refcode}"]
        if name:
            parts.append(name[:50])
        if csd.doi:
            parts.append(f"DOI:{csd.doi}")
        if csd.year:
            parts.append(f"({csd.year})")
        lines = [" | ".join(parts)]

        if self.coremof_records:
            for rec in self.coremof_records:
                topo = rec.topology_single
                topo_str = f"  {topo}" if topo and topo not in ("nan", "unnamed") else ""
                lines.append(
                    f"  -> {rec.coreid}  ({rec.extension})  "
                    f"{rec.metal_types}{topo_str}"
                )
        else:
            lines.append("  (no CoreMOF match)")

        return "\n".join(lines)


@dataclass
class CoreMOFSearchResult:
    """Result of a CoRE MOF database query."""

    query: str
    field: str
    records: list[CoreMOFRecord] = field(default_factory=list)

    @property
    def n_matches(self) -> int:
        return len(self.records)

    def summary(self) -> str:
        """Multi-line human-readable summary."""
        lines = [
            f"CoreMOF lookup: {self.n_matches} match(es) "
            f"for '{self.query}' (field: {self.field})"
        ]
        for rec in self.records:
            lines.append(f"  {rec.summary()}")
        return "\n".join(lines)
