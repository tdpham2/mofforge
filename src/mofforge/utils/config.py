"""Global configuration for mofforge: paths, bonding radii, R-group tag."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger("mofforge")

# Covalent radii (Angstroms) from Cordero et al., Dalton Trans. 2008, 2832-2838.
# With overrides from Alvarez, Dalton Trans. 2008, 2832.
COVALENT_RADII: dict[str, float] = {
    "H": 0.31,
    "He": 0.28,
    "Li": 1.28,
    "Be": 0.96,
    "B": 0.84,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "Ne": 0.58,
    "Na": 1.66,
    "Mg": 1.41,
    "Al": 1.21,
    "Si": 1.11,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Ar": 1.06,
    "K": 2.03,
    "Ca": 1.76,
    "Sc": 1.70,
    "Ti": 1.60,
    "V": 1.53,
    "Cr": 1.39,
    "Mn": 1.39,
    "Fe": 1.32,
    "Co": 1.26,
    "Ni": 1.24,
    "Cu": 1.32,
    "Zn": 1.22,
    "Ga": 1.22,
    "Ge": 1.20,
    "As": 1.19,
    "Se": 1.20,
    "Br": 1.20,
    "Kr": 1.16,
    "Rb": 2.20,
    "Sr": 1.95,
    "Y": 1.90,
    "Zr": 1.75,
    "Nb": 1.64,
    "Mo": 1.54,
    "Tc": 1.47,
    "Ru": 1.46,
    "Rh": 1.42,
    "Pd": 1.39,
    "Ag": 1.45,
    "Cd": 1.44,
    "In": 1.42,
    "Sn": 1.39,
    "Sb": 1.39,
    "Te": 1.38,
    "I": 1.39,
    "Xe": 1.40,
    "Cs": 2.44,
    "Ba": 2.15,
    "La": 2.07,
    "Ce": 2.04,
    "Pr": 2.03,
    "Nd": 2.01,
    "Pm": 1.99,
    "Sm": 1.98,
    "Eu": 1.98,
    "Gd": 1.96,
    "Tb": 1.94,
    "Dy": 1.92,
    "Ho": 1.92,
    "Er": 1.89,
    "Tm": 1.90,
    "Yb": 1.87,
    "Lu": 1.87,
    "Hf": 1.75,
    "Ta": 1.70,
    "W": 1.62,
    "Re": 1.51,
    "Os": 1.44,
    "Ir": 1.41,
    "Pt": 1.36,
    "Au": 1.36,
    "Hg": 1.32,
    "Tl": 1.45,
    "Pb": 1.46,
    "Bi": 1.48,
    "Po": 1.40,
    "At": 1.50,
    "Rn": 1.50,
    "Fr": 2.60,
    "Ra": 2.21,
    "Ac": 2.15,
    "Th": 2.06,
    "Pa": 2.00,
    "U": 1.96,
    "Np": 1.90,
    "Pu": 1.87,
    "Am": 1.80,
    "Cm": 1.69,
}

# Van der Waals radii (Angstroms) for steric clash detection
VDW_RADII: dict[str, float] = {
    "H": 1.20,
    "He": 1.40,
    "Li": 1.82,
    "Be": 1.53,
    "B": 1.92,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "Ne": 1.54,
    "Na": 2.27,
    "Mg": 1.73,
    "Al": 1.84,
    "Si": 2.10,
    "P": 1.80,
    "S": 1.80,
    "Cl": 1.75,
    "Ar": 1.88,
    "K": 2.75,
    "Ca": 2.31,
    "Sc": 2.11,
    "Ti": 1.87,
    "V": 1.79,
    "Cr": 1.89,
    "Mn": 1.97,
    "Fe": 1.94,
    "Co": 1.92,
    "Ni": 1.63,
    "Cu": 1.40,
    "Zn": 1.39,
    "Ga": 1.87,
    "Ge": 2.11,
    "As": 1.85,
    "Se": 1.90,
    "Br": 1.85,
    "Kr": 2.02,
    "Rb": 3.03,
    "Sr": 2.49,
    "Y": 2.19,
    "Zr": 1.86,
    "Nb": 2.07,
    "Mo": 2.09,
    "Pd": 1.63,
    "Ag": 1.72,
    "Cd": 1.58,
    "In": 1.93,
    "Sn": 2.17,
    "Sb": 2.06,
    "Te": 2.06,
    "I": 1.98,
    "Xe": 2.16,
    "Cs": 3.43,
    "Ba": 2.68,
    "Pt": 1.75,
    "Au": 1.66,
    "Hg": 1.55,
    "Tl": 1.96,
    "Pb": 2.02,
    "Bi": 2.07,
    "U": 1.86,
}


@dataclass
class MofforgeConfig:
    """Global configuration for mofforge."""

    # R-group tag character
    r_tag: str = "!"

    # Default data paths (None means no default; user must set explicitly)
    crystal_path: Path | None = field(default=None)
    moiety_path: Path | None = field(default=None)
    csd_data_path: Path | None = field(default=None)
    coremof_data_path: Path | None = field(default=None)

    @property
    def fragment_path(self) -> Path | None:
        """Alias for ``moiety_path`` (preferred name)."""
        return self.moiety_path

    @fragment_path.setter
    def fragment_path(self, value: Path | None) -> None:
        self.moiety_path = value

    # Bonding parameters
    bond_pad: float = 0.25  # padding added to sum of covalent radii (Angstroms)

    def get_covalent_radius(self, species: str) -> float:
        """Get covalent radius for a species, stripping R-group tag if present."""
        clean = species.removesuffix(self.r_tag)
        if clean in COVALENT_RADII:
            return COVALENT_RADII[clean]
        raise ValueError(f"No covalent radius for species '{clean}'.")

    def get_vdw_radius(self, species: str) -> float:
        """Get van der Waals radius for a species, stripping R-group tag if present."""
        clean = species.removesuffix(self.r_tag)
        if clean in VDW_RADII:
            return VDW_RADII[clean]
        # Fall back to covalent radius + 0.5 if vdW not available
        return self.get_covalent_radius(species) + 0.5

    def max_bond_distance(self, species_i: str, species_j: str) -> float:
        """Compute maximum bond distance from sum of covalent radii + padding."""
        return (
            self.get_covalent_radius(species_i)
            + self.get_covalent_radius(species_j)
            + self.bond_pad
        )


# Global config singleton
config = MofforgeConfig()


def set_paths(
    crystals: str | Path | None = None,
    moieties: str | Path | None = None,
    csd_data: str | Path | None = None,
    coremof_data: str | Path | None = None,
) -> None:
    """Update global data directory paths."""
    if crystals is not None:
        config.crystal_path = Path(crystals)
    if moieties is not None:
        config.moiety_path = Path(moieties)
    if csd_data is not None:
        config.csd_data_path = Path(csd_data)
    if coremof_data is not None:
        config.coremof_data_path = Path(coremof_data)


# Regex to extract bare element symbol (e.g. "Zn2+" -> "Zn", "H!" -> "H")
_ELEMENT_RE = re.compile(r"^([A-Z][a-z]?)")


def clean_species(species: str) -> str:
    """Extract the bare element symbol from a species label (e.g. ``'H!'`` -> ``'H'``)."""
    cleaned = species.removesuffix(config.r_tag)
    m = _ELEMENT_RE.match(cleaned)
    return m.group(1) if m else cleaned
