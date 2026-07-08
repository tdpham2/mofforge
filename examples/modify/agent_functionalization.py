#!/usr/bin/env python3
"""Agent-driven linker functionalization.

Demonstrates the workflow an AI agent uses to functionalize a MOF linker without
authoring any geometry or SMILES:

    1. Start from the linker SMILES (e.g. obtained from MOFid).
    2. Enumerate functionalizable aromatic C-H sites (each with a symmetry class).
    3. Pick a functional group from the curated menu.
    4. Functionalize at a chosen site with a chosen coverage.
    5. (Optionally) run a full campaign sweeping groups x coverages.

Input:
    data/crystals/IRMOF-1.cif   - Parent MOF (linker is BDC / terephthalate)

Usage:
    python modify/agent_functionalization.py
    python modify/agent_functionalization.py --group NH2 --coverage 0.5
    python modify/agent_functionalization.py --campaign
"""

import argparse
from pathlib import Path

from mofforge import (
    available_groups,
    find_functionalizable_sites,
    functionalize,
    run_campaign,
)

SCRIPT_DIR = Path(__file__).parent
CRYSTAL_DIR = SCRIPT_DIR.parent / "data" / "crystals"

# BDC (benzene-1,4-dicarboxylic acid), the linker of IRMOF-1.
BDC_SMILES = "O=C(O)c1ccc(C(=O)O)cc1"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", default="NO2", help="Functional group to add.")
    parser.add_argument("--site", type=int, default=0, help="Site index to functionalize.")
    parser.add_argument("--coverage", type=float, default=0.5, help="Fraction of linkers (0-1).")
    parser.add_argument("--output", default="functionalized.cif", help="Output CIF path.")
    parser.add_argument("--campaign", action="store_true", help="Run a full group x coverage sweep.")
    args = parser.parse_args()

    cif = str(CRYSTAL_DIR / "IRMOF-1.cif")

    # ------------------------------------------------------------------ #
    # Step 1-2: inspect the linker's functionalizable sites.
    # ------------------------------------------------------------------ #
    print(f"Linker SMILES: {BDC_SMILES}")
    sites = find_functionalizable_sites(BDC_SMILES)
    print(f"\nFunctionalizable sites ({len(sites)}):")
    for s in sites:
        print(f"  index={s.index}  symmetry_class={s.symmetry_class}  {s.description}")

    print(f"\nAvailable functional groups: {', '.join(available_groups())}")

    if args.campaign:
        # -------------------------------------------------------------- #
        # Autonomous campaign: sweep groups x coverages, ranked best-first.
        # -------------------------------------------------------------- #
        print("\nRunning campaign (NH2, F, NO2) x (0.25, 0.5, 1.0) ...")
        results = run_campaign(
            cif,
            BDC_SMILES,
            groups=["NH2", "F", "NO2"],
            coverages=[0.25, 0.5, 1.0],
            output_dir="campaign_output",
            random_seed=1,
        )
        print("\nRanked results (best first):")
        for r in results:
            print(
                f"  {r.group:10s} coverage={r.coverage:<4} "
                f"n_func={r.n_functionalized:<3} valid={r.is_valid} "
                f"clashes={r.clashes} -> {r.output_cif}"
            )
        return

    # ------------------------------------------------------------------ #
    # Step 3-4: functionalize one group at one site with a coverage.
    # ------------------------------------------------------------------ #
    print(
        f"\nFunctionalizing: group={args.group} site={args.site} "
        f"coverage={args.coverage}"
    )
    result = functionalize(
        cif,
        BDC_SMILES,
        args.group,
        sites=args.site,
        coverage=args.coverage,
        output_cif=args.output,
        random_seed=1,
    )
    if result.error:
        print(f"  ERROR: {result.error}")
        return
    print(f"  matched linkers:     {result.n_matches}")
    print(f"  functionalized:      {result.n_functionalized}")
    print(f"  valid:               {result.is_valid}  (clashes={result.clashes})")
    print(f"  wrote:               {result.output_cif}")


if __name__ == "__main__":
    main()
