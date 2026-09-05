#!/usr/bin/env python3
"""Execute YAML batches and optionally compare serial and parallel results.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import (
    CRYSTALS,
    EXAMPLES,
    MOIETIES,
    example_parser,
    output_dir,
    write_report,
)


def main():
    parser = example_parser("batch_processing", __doc__, seed=True)
    parser.add_argument("--parallel", type=int, default=0, help="Workers; 0 runs serially.")
    parser.add_argument(
        "--compare-parallel",
        action="store_true",
        help="Also run two workers and compare with a serial reference.",
    )
    args = parser.parse_args()
    out = output_dir(args)
    import shutil
    from dataclasses import asdict

    import numpy as np
    import yaml

    from mofforge import Crystal, run_batch

    # Two same-named sources show why callers must consume returned output paths.
    inputs = []
    for name in ["first", "second"]:
        directory = out / "inputs" / name
        directory.mkdir(parents=True, exist_ok=True)
        source = directory / "parent.cif"
        shutil.copyfile(CRYSTALS / "IRMOF-1.cif", source)
        inputs.append(str(source))
    recipe = yaml.safe_load((EXAMPLES / "pipeline" / "batch.yaml").read_text())
    recipe.update(
        parents=inputs,
        moiety_path=str(MOIETIES),
        random_seed=args.seed,
        parallel=0 if args.compare_parallel else args.parallel,
    )
    recipe["operations"][0]["random_seed"] = args.seed
    recipe["output"]["directory"] = str(out / "results")
    config = out / "batch.yaml"
    config.write_text(yaml.safe_dump(recipe), encoding="utf-8")
    results = run_batch(config)
    if not results or any(not r.success for r in results):
        raise SystemExit(f"Batch failed: {[r.error for r in results]}")
    comparison = None
    if args.compare_parallel:
        recipe["parallel"] = 2
        recipe["output"]["directory"] = str(out / "parallel_results")
        parallel_config = out / "batch_parallel.yaml"
        parallel_config.write_text(yaml.safe_dump(recipe), encoding="utf-8")
        parallel = run_batch(parallel_config)
        if len(results) != len(parallel) or any(not r.success for r in parallel):
            raise SystemExit("Parallel batch did not complete.")
        for first, second in zip(results, parallel, strict=True):
            if first.run_id != second.run_id:
                raise RuntimeError("Run identities differ between serial and parallel batches.")
            a, b = Crystal.from_cif(first.output_path), Crystal.from_cif(second.output_path)
            if a.species != b.species:
                raise RuntimeError("Serial and parallel species differ.")
            np.testing.assert_allclose(a.frac_coords, b.frac_coords, atol=1e-8)
        comparison = True
    write_report(
        out,
        seed=args.seed,
        config=config,
        results=[asdict(r) for r in results],
        parallel_matches_serial=comparison,
    )


if __name__ == "__main__":
    main()
