#!/usr/bin/env python3
"""Chain replacements and validation while saving intermediate structures.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, MOIETIES, example_parser, output_dir, write_report


def main():
    args = example_parser("multi_step_pipeline", __doc__, seed=True).parse_args()
    out = output_dir(args)
    from mofforge import Pipeline

    # build_all executes the queued steps; it does not resume an earlier run.
    pipeline = Pipeline(CRYSTALS / "IRMOF-1.cif", fragment_path=MOIETIES, random_seed=args.seed)
    pipeline.replace("2-!-p-phenylene.xyz", "2-acetylamido-p-phenylene.xyz", nb_loc=2)
    pipeline.replace("2-!-p-phenylene.xyz", "2-nitro-p-phenylene.xyz", nb_loc=2)
    pipeline.validate()
    intermediates = pipeline.build_all(name="mixed_linker")
    for index, crystal in enumerate(intermediates):
        crystal.write_cif(out / f"step_{index + 1}.cif")
    final = intermediates[-1]
    final.provenance.to_json(out / "provenance.json")
    write_report(
        out,
        seed=args.seed,
        intermediate_atoms=[c.n_atoms for c in intermediates],
        validation=[r.to_dict() for r in pipeline.validation_reports],
        provenance=final.provenance.to_dict(),
    )


if __name__ == "__main__":
    main()
