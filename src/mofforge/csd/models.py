"""Data models for CSD (Cambridge Structural Database) lookup results."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class CSDRecord:
    """A single CSD database entry."""

    refcode: str
    chemical_name_systematic: str = ""
    chemical_name_common: str = ""
    chemical_formula_moiety: str = ""
    doi: str | None = None
    ccdc_number: str | None = None
    authors: str = ""
    journal: str = ""
    volume: str = ""
    pages: str = ""
    year: str = ""
    space_group: str = ""
    cell_a: str = ""
    cell_b: str = ""
    cell_c: str = ""
    cell_alpha: str = ""
    cell_beta: str = ""
    cell_gamma: str = ""
    cell_volume: str = ""
    temperature: str = ""
    r_factor: str = ""
    raw: dict[str, str] = field(default_factory=dict, repr=False)

    def summary(self) -> str:
        """Human-readable one-line summary."""
        parts = [self.refcode]
        name = self.chemical_name_common or self.chemical_name_systematic
        if name:
            parts.append(name[:60])
        if self.doi:
            parts.append(f"DOI:{self.doi}")
        if self.year:
            parts.append(f"({self.year})")
        return " | ".join(parts)


@dataclass
class CSDSearchResult:
    """Result of a CSD lookup query."""

    query: str
    field: str
    records: list[CSDRecord] = field(default_factory=list)

    @property
    def n_matches(self) -> int:
        return len(self.records)

    def summary(self) -> str:
        """Multi-line human-readable summary."""
        lines = [f"CSD lookup: {self.n_matches} match(es) for '{self.query}' (field: {self.field})"]
        for rec in self.records:
            lines.append(f"  {rec.summary()}")
        return "\n".join(lines)
