#!/usr/bin/env python3
"""Multi-step Pipeline

Demonstrates the Pipeline API for chaining multiple find-and-replace
operations on a crystal structure in a single fluent call.

Usage:
    python pipeline/multi_step_pipeline.py
"""

from pathlib import Path

from mofforge import Crystal, infer_bonds
from mofforge.pipeline import Pipeline

SCRIPT_DIR = Path(__file__).parent
EXAMPLES_DIR = SCRIPT_DIR.parent
CRYSTAL_DIR = EXAMPLES_DIR / "data" / "crystals"
MOIETY_DIR = EXAMPLES_DIR / "data" / "moieties"


def run_pipeline_demo():
    """Demonstrate multi-step pipeline operations."""

    # =========================================================================
    # Example 1: Single-step pipeline (equivalent to direct replace)
    # =========================================================================
    print("PIPELINE EXAMPLE 1: Single-step functionalization")
    print("-" * 50)

    parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
    parent = infer_bonds(parent, periodic=True)

    child = (
        Pipeline(parent, fragment_path=str(MOIETY_DIR))
        .replace(
            query="2-!-p-phenylene.xyz",
            replacement="2-acetylamido-p-phenylene.xyz",
            nb_loc=4,
        )
        .build(name="pipeline_single_step")
    )

    print(f"  Parent: {parent.n_atoms} atoms")
    print(f"  Child:  {child.n_atoms} atoms")
    print(f"  Name:   {child.name}")
    print()

    # =========================================================================
    # Example 2: Multi-step pipeline with validation
    #
    # Step 1: Functionalize 3 linkers with acetylamido
    # Step 2: Validate the resulting structure
    # =========================================================================
    print("PIPELINE EXAMPLE 2: Functionalize + validate")
    print("-" * 50)

    pipeline = Pipeline(parent, fragment_path=str(MOIETY_DIR))
    pipeline.replace(
        query="2-!-p-phenylene.xyz",
        replacement="2-acetylamido-p-phenylene.xyz",
        nb_loc=3,
    )
    pipeline.validate()  # Check structure quality after modification

    child = pipeline.build(name="pipeline_validated")

    print(f"  Parent: {parent.n_atoms} atoms")
    print(f"  Child:  {child.n_atoms} atoms")

    # Access validation results
    reports = pipeline.validation_reports
    if reports:
        print(f"\n  Validation Report:")
        print(f"    Valid: {reports[0].is_valid}")
        print(f"    Steric clashes: {len(reports[0].steric_clashes)}")
        print(f"    Unusual bonds:  {len(reports[0].unusual_bonds)}")
    print()

    # =========================================================================
    # Example 3: Pipeline from CIF file path
    #
    # You can pass a file path directly instead of pre-loading the Crystal.
    # =========================================================================
    print("PIPELINE EXAMPLE 3: Pipeline from file path")
    print("-" * 50)

    child = (
        Pipeline(CRYSTAL_DIR / "IRMOF-1.cif", fragment_path=str(MOIETY_DIR))
        .replace(
            query="2-!-p-phenylene.xyz",
            replacement="2-nitro-p-phenylene.xyz",
            nb_loc=2,
        )
        .validate()
        .build(name="nitro_IRMOF-1")
    )

    print(f"  Child: {child.n_atoms} atoms, name='{child.name}'")
    print()

    # =========================================================================
    # Example 4: Access intermediates with build_all()
    # =========================================================================
    print("PIPELINE EXAMPLE 4: Access all intermediates")
    print("-" * 50)

    pipeline = Pipeline(parent, fragment_path=str(MOIETY_DIR))
    pipeline.replace(
        query="2-!-p-phenylene.xyz",
        replacement="2-acetylamido-p-phenylene.xyz",
        nb_loc=2,
    )
    pipeline.replace(
        query="2-!-p-phenylene.xyz",
        replacement="2-nitro-p-phenylene.xyz",
        nb_loc=2,
    )

    all_steps = pipeline.build_all(name="multi_step")

    print(f"  Number of intermediate structures: {len(all_steps)}")
    for i, step_xtal in enumerate(all_steps):
        print(f"    Step {i}: {step_xtal.n_atoms} atoms")

    # Check provenance
    final = all_steps[-1]
    if final.provenance:
        print(f"\n  Provenance chain:")
        print(f"    {final.provenance.summary()}")


if __name__ == "__main__":
    run_pipeline_demo()
