"""Built-in adsorbate molecule geometries for MOF adsorption studies."""

from __future__ import annotations

import numpy as np

_MOLECULES: dict[str, list[tuple[str, float, float, float]]] = {
    # --- Monatomic / noble gases ---
    "He": [("He", 0.0, 0.0, 0.0)],
    "Ne": [("Ne", 0.0, 0.0, 0.0)],
    "Ar": [("Ar", 0.0, 0.0, 0.0)],
    "Kr": [("Kr", 0.0, 0.0, 0.0)],
    "Xe": [("Xe", 0.0, 0.0, 0.0)],
    # --- Diatomic ---
    "H2": [
        ("H", 0.0, 0.0, -0.37),
        ("H", 0.0, 0.0, 0.37),
    ],
    "N2": [
        ("N", 0.0, 0.0, -0.5488),
        ("N", 0.0, 0.0, 0.5488),
    ],
    "O2": [
        ("O", 0.0, 0.0, -0.6037),
        ("O", 0.0, 0.0, 0.6037),
    ],
    "CO": [
        ("C", 0.0, 0.0, -0.5641),
        ("O", 0.0, 0.0, 0.5641),
    ],
    "NO": [
        ("N", 0.0, 0.0, -0.5756),
        ("O", 0.0, 0.0, 0.5756),
    ],
    # --- Triatomic linear ---
    "CO2": [
        ("O", 0.0, 0.0, -1.162),
        ("C", 0.0, 0.0, 0.0),
        ("O", 0.0, 0.0, 1.162),
    ],
    # --- Triatomic bent ---
    "H2O": [
        # H-O-H angle = 104.5 deg, O-H = 0.9572 A, centered
        ("O", 0.0, 0.0, 0.3907),
        ("H", 0.0, 0.7568, -0.1953),
        ("H", 0.0, -0.7568, -0.1953),
    ],
    "H2S": [
        # H-S-H angle = 92.1 deg, S-H = 1.336 A, centered
        ("S", 0.0, 0.0, 0.6181),
        ("H", 0.0, 0.9618, -0.3091),
        ("H", 0.0, -0.9618, -0.3091),
    ],
    "SO2": [
        # O-S-O angle = 119.3 deg, S-O = 1.4321 A, centered
        ("S", 0.0, 0.0, 0.4824),
        ("O", 0.0, 1.2358, -0.2412),
        ("O", 0.0, -1.2358, -0.2412),
    ],
    "NO2": [
        # O-N-O angle = 134.1 deg, N-O = 1.1934 A, centered
        ("N", 0.0, 0.0, 0.3102),
        ("O", 0.0, 1.0989, -0.1551),
        ("O", 0.0, -1.0989, -0.1551),
    ],
    # --- Tetrahedral ---
    "CH4": [
        # Td symmetry, C-H = 1.087 A
        ("C", 0.0, 0.0, 0.0),
        ("H", 0.6276, 0.6276, 0.6276),
        ("H", -0.6276, -0.6276, 0.6276),
        ("H", -0.6276, 0.6276, -0.6276),
        ("H", 0.6276, -0.6276, -0.6276),
    ],
    # --- Trigonal pyramidal ---
    "NH3": [
        # C3v symmetry, H-N-H angle = 107.8 deg, N-H = 1.012 A, centered
        ("N", 0.0, 0.0, 0.2732),
        ("H", 0.9442, 0.0, -0.0911),
        ("H", -0.4721, 0.8177, -0.0911),
        ("H", -0.4721, -0.8177, -0.0911),
    ],
    # --- Other common adsorbates ---
    "C2H2": [
        # Acetylene, linear, C-C = 1.203 A, C-H = 1.060 A
        ("H", 0.0, 0.0, -1.6615),
        ("C", 0.0, 0.0, -0.6015),
        ("C", 0.0, 0.0, 0.6015),
        ("H", 0.0, 0.0, 1.6615),
    ],
    "C2H4": [
        # Ethylene, planar, C=C = 1.339 A, C-H = 1.086 A
        ("C", 0.0, 0.0, -0.6695),
        ("C", 0.0, 0.0, 0.6695),
        ("H", 0.0, 0.9289, -1.2321),
        ("H", 0.0, -0.9289, -1.2321),
        ("H", 0.0, 0.9289, 1.2321),
        ("H", 0.0, -0.9289, 1.2321),
    ],
}

# Case-insensitive aliases for convenience
_ALIASES: dict[str, str] = {
    "water": "H2O",
    "carbon_dioxide": "CO2",
    "carbon dioxide": "CO2",
    "carbondioxide": "CO2",
    "methane": "CH4",
    "nitrogen": "N2",
    "oxygen": "O2",
    "hydrogen": "H2",
    "ammonia": "NH3",
    "acetylene": "C2H2",
    "ethylene": "C2H4",
    "carbon_monoxide": "CO",
    "carbon monoxide": "CO",
    "carbonmonoxide": "CO",
    "nitric_oxide": "NO",
    "nitric oxide": "NO",
    "nitrogen_dioxide": "NO2",
    "nitrogen dioxide": "NO2",
    "sulfur_dioxide": "SO2",
    "sulfur dioxide": "SO2",
    "hydrogen_sulfide": "H2S",
    "hydrogen sulfide": "H2S",
    "argon": "Ar",
    "helium": "He",
    "neon": "Ne",
    "krypton": "Kr",
    "xenon": "Xe",
}


def available_molecules() -> list[str]:
    """Return sorted list of available molecule names."""
    return sorted(_MOLECULES.keys())


def get_molecule(name: str) -> tuple[list[str], np.ndarray]:
    """Get the geometry of a predefined adsorbate molecule."""
    # Try direct lookup first (case-sensitive)
    if name in _MOLECULES:
        atoms = _MOLECULES[name]
        species = [a[0] for a in atoms]
        coords = np.array([[a[1], a[2], a[3]] for a in atoms])
        return species, coords

    # Try case-insensitive lookup
    name_lower = name.lower()
    for key in _MOLECULES:
        if key.lower() == name_lower:
            atoms = _MOLECULES[key]
            species = [a[0] for a in atoms]
            coords = np.array([[a[1], a[2], a[3]] for a in atoms])
            return species, coords

    # Try aliases
    if name_lower in _ALIASES:
        return get_molecule(_ALIASES[name_lower])

    raise ValueError(f"Unknown molecule '{name}'")
