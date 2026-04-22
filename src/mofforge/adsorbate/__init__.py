"""Adsorbate initialization for MOF structures.

Provides tools for identifying adsorption sites in metal-organic
frameworks and placing adsorbate molecules at those sites, prior to
running MLIP or other atomistic calculations.

Submodules:
    - :mod:`sites`: Adsorption site identification (void detection,
      open metal sites).
    - :mod:`placement`: Adsorbate molecule placement and combination
      with the host framework.
    - :mod:`molecules`: Built-in geometries for common adsorbate molecules.
"""

from mofforge.adsorbate.molecules import available_molecules, get_molecule
from mofforge.adsorbate.placement import AdsorbatePlacement, place_adsorbate
from mofforge.adsorbate.sites import AdsorptionSite, find_adsorption_sites

__all__ = [
    "AdsorbatePlacement",
    "AdsorptionSite",
    "available_molecules",
    "find_adsorption_sites",
    "get_molecule",
    "place_adsorbate",
]
