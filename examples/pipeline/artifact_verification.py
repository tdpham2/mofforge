#!/usr/bin/env python3
"""Distinguish complete batch verification from integrity-only manifests.

Run this file directly; see the topic README for the scientific walkthrough.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Only cookbook helpers are added to the path; install mofforge separately.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import CRYSTALS, MOIETIES, example_parser, output_dir, write_report


def main():
    args = example_parser("artifact_verification", __doc__, seed=True).parse_args()
    out = output_dir(args)
    import shutil

    import yaml

    from mofforge import Crystal, run_batch, verify_artifact

    # Full verification requires batch input descriptors with current source files.
    recipe = {
        "parents": [str(CRYSTALS / "IRMOF-1.cif")],
        "operations": [
            {
                "type": "replace",
                "query": "2-!-p-phenylene.xyz",
                "replacement": "2-nitro-p-phenylene.xyz",
                "mode": "nb_loc_2",
                "random_seed": args.seed,
            }
        ],
        "moiety_path": str(MOIETIES),
        "random_seed": args.seed,
        "output": {"directory": str(out / "batch"), "format": "cif"},
    }
    config = out / "batch.yaml"
    config.write_text(yaml.safe_dump(recipe), encoding="utf-8")
    results = run_batch(config)
    if not results or not results[0].success:
        raise SystemExit("Batch output could not be generated.")
    artifact = Path(results[0].output_path)
    verified = verify_artifact(artifact)
    # A plain export has file integrity but no complete batch input descriptor.
    plain = out / "plain.cif"
    Crystal.from_cif(CRYSTALS / "IRMOF-1.cif").write_cif(plain)
    integrity_only = verify_artifact(plain)
    # Alter a COPY. Its original manifest must reject the changed bytes.
    altered = out / "altered_copy.cif"
    shutil.copyfile(artifact, altered)
    with altered.open("a", encoding="utf-8") as stream:
        stream.write("\n# Deliberate tutorial modification\n")
    mismatch = verify_artifact(altered, manifest_path=artifact.with_suffix(".cif.json"))
    write_report(
        out,
        artifact=artifact,
        complete=verified.to_dict(),
        integrity_only=integrity_only.to_dict(),
        altered=mismatch.to_dict(),
    )
    if not verified.is_verified or not integrity_only.file_integrity or mismatch.file_integrity:
        raise SystemExit("Verification results differ from the tutorial's expected controls.")


if __name__ == "__main__":
    main()
