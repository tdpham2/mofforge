#!/usr/bin/env python3
"""Prepare, pack, or connect a generic aromatic construction fixture.

The default mode prepares labeled templates without Packmol. --mode pack uses
real Packmol; --mode connect additionally creates a finite chain, saves/reloads
its partial state, imports mapped geometry, and resumes connection. The fixture
is a geometric construction example, not a polymer synthesis recipe.
"""

from __future__ import annotations

import json
import sys
from dataclasses import asdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import example_parser, output_dir, require, write_report


def main():
    parser = example_parser("amorphous_pop", __doc__, seed=True)
    parser.add_argument("--mode", choices=["prepare", "pack", "connect"], default="prepare")
    parser.add_argument("--initial-packing-density", type=float)
    args = parser.parse_args()
    require("rdkit", "pop")
    out = output_dir(args)

    from mofforge.polymerize import (
        ConnectionRule,
        Connector,
        ConstructionState,
        PopBuilder,
        connect,
    )

    connectors = [
        Connector("left", "a0", "h:a0:1", "aryl", ("h:a0:1",)),
        Connector("right", "a3", "h:a3:1", "aryl", ("h:a3:1",)),
    ]
    rule = ConnectionRule(
        "aryl_construction_fixture",
        ("aryl", "aryl"),
        1,
        (1.45, 1.60),
        {"a": ["$replaceable"], "b": ["$replaceable"]},
    )
    packing = {"random_seed": args.seed, "material_state": "construction_fixture"}
    if args.initial_packing_density is None:
        packing["box_lengths"] = [25, 25, 25]  # This example's volume; not an API default.
    else:
        packing["initial_packing_density"] = args.initial_packing_density
    config = {
        "components": [
            {"source": "c1ccccc1", "count": 3, "connectors": [asdict(c) for c in connectors]}
        ],
        "packing": packing,
        "connection": {
            "rules": [asdict(rule)],
            "target_conversion": 2 / 3,
            "candidate_attempt_budget": 30,
            "min_nonbonded_distance": 1.5,
        },
    }
    (out / "construction.json").write_text(json.dumps(config, indent=2) + "\n")
    (out / "packing.json").write_text(
        json.dumps({k: v for k, v in config.items() if k != "connection"}, indent=2) + "\n"
    )
    builder = PopBuilder()
    builder.add_monomer("c1ccccc1", count=3, connectors=connectors)
    templates = builder.prepare(random_seed=args.seed)
    results = []
    if args.mode == "pack":
        result = builder.pack(output_dir=out, **packing)
        results.append(result.to_dict())
        if not result.success:
            raise RuntimeError(result.errors)
    elif args.mode == "connect":
        first = builder.build(
            output_dir=out,
            rules=[rule],
            target_conversion=2 / 3,
            candidate_attempt_budget=1,
            min_nonbonded_distance=1.5,
            **packing,
        )
        results.append(first.to_dict())
        if first.state_path is None or first.status == "failed":
            raise RuntimeError(first.errors)
        state = ConstructionState.load(first.state_path)
        records = [
            {
                "id": a.id,
                "species": a.species,
                "coordinates": [float(xyz[0] + 0.1), float(xyz[1] + 0.2), float(xyz[2] + 0.3)],
            }
            for a, xyz in zip(state.atoms, state.coordinates, strict=True)
        ]
        # A mapped translation exercises the handoff; no external simulation is claimed.
        updated = state.update_geometry(parent_hash=state.state_hash, atoms=records, unwrapped=True)
        final = connect(
            updated,
            [rule],
            target_conversion=2 / 3,
            candidate_attempt_budget=30,
            min_nonbonded_distance=1.5,
            output_dir=out,
        )
        results.append(final.to_dict())
        if not final.success:
            raise RuntimeError(final.errors)
    write_report(
        out,
        seed=args.seed,
        mode=args.mode,
        templates=[template.to_dict() for template in templates],
        results=results,
        configuration=str(out / "construction.json"),
    )


if __name__ == "__main__":
    main()
