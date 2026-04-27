"""Provenance tracking for crystal modifications."""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path

logger = logging.getLogger("mofforge")


@dataclass
class Provenance:
    """Metadata tracking modifications made to a crystal structure."""

    parent: str | None = None
    query: str | None = None
    replacement: str | None = None
    operation: str | None = None
    parameters: dict = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())
    alignment_errors: list[float] = field(default_factory=list)
    history: list[dict] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Convert to a plain dictionary."""
        return asdict(self)

    def to_json(self, filepath: str | Path) -> None:
        """Write provenance to a JSON file."""
        filepath = Path(filepath)
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def from_dict(cls, d: dict) -> Provenance:
        """Create a Provenance from a dictionary."""
        return cls(
            parent=d.get("parent"),
            query=d.get("query"),
            replacement=d.get("replacement"),
            operation=d.get("operation"),
            parameters=d.get("parameters", {}),
            timestamp=d.get("timestamp", datetime.now().isoformat()),
            alignment_errors=d.get("alignment_errors", []),
            history=d.get("history", []),
        )

    @classmethod
    def from_json(cls, filepath: str | Path) -> Provenance:
        """Load provenance from a JSON file."""
        with open(filepath, encoding="utf-8") as f:
            d = json.load(f)
        return cls.from_dict(d)

    def chain(self, new_provenance: Provenance) -> Provenance:
        """Create a chained provenance record."""
        new_provenance.history = [*self.history, self.to_dict(), *new_provenance.history]
        return new_provenance

    def summary(self) -> str:
        """Return a human-readable summary of the provenance chain."""
        lines = []
        for i, hist in enumerate(self.history):
            lines.append(
                f"  Step {i + 1}: {hist.get('operation', '?')} ({hist.get('timestamp', '?')})"
            )
        lines.append(f"  Current: {self.operation} ({self.timestamp})")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"Provenance(op='{self.operation}', "
            f"parent='{self.parent}', "
            f"query='{self.query}', "
            f"replacement='{self.replacement}')"
        )
