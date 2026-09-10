"""Direct Packmol execution and strict topology-preserving output validation."""

from __future__ import annotations

import json
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist

from mofforge.io.xyz import read_xyz, write_xyz
from mofforge.polymerize.base import POPResult, positive, reject_removed
from mofforge.polymerize.base import box_lengths as validate_lengths
from mofforge.polymerize.config import PopBuildConfig, probe_packmol
from mofforge.polymerize.geometry import validate_box
from mofforge.polymerize.provision import instantiate
from mofforge.polymerize.state import AVOGADRO
from mofforge.provenance import effective_seed, file_hash, software_versions


def persist_result(result, output_dir):
    if result.state is not None:
        result.crystal = result.state.crystal
        result.metadata = {**result.state.metadata, **result.state.statistics, **result.metadata}
        try:
            result.state_path = result.state.save(output_dir)
        except (OSError, ValueError) as exc:
            result.success = False
            result.status = "failed"
            result.errors.append(f"Could not save native box bundle in {output_dir}: {exc}")
            result.output_paths.append(Path(output_dir))
            return result
        bundle = result.state_path.parent
        result.output_paths.extend(
            bundle / name for name in ("state.json", "box.cif", "box.xyz", "manifest.json")
        )
    return result


def pack(
    monomers,
    prepared,
    *,
    output_dir=".",
    packmol_bin=None,
    box_lengths=None,
    initial_packing_density=None,
    minimum_distance=2.0,
    random_seed=None,
    timeout=300.0,
    material_state="unspecified",
    experimental_density=None,
    **options,
):
    reject_removed(options)
    if options:
        raise TypeError(f"Unknown packing options: {sorted(options)}")
    if not monomers:
        raise ValueError("No monomers registered; call add_monomer with an explicit count.")
    if (box_lengths is None) == (initial_packing_density is None):
        raise ValueError("Supply exactly one of box_lengths or initial_packing_density.")
    minimum_distance = positive(minimum_distance, "minimum_distance")
    timeout = positive(timeout, "timeout")
    if type(random_seed) is bool:
        raise ValueError("random_seed must be an integer.")
    seed = effective_seed(random_seed)
    if not isinstance(material_state, str) or not material_state:
        raise ValueError("material_state must be a nonempty declaration.")
    if experimental_density is not None:
        if set(experimental_density) != {"value", "unit", "measurement"}:
            raise ValueError("Experimental density requires value, unit, measurement.")
        positive(experimental_density["value"], "experimental density")
        if not experimental_density["unit"] or not experimental_density["measurement"]:
            raise ValueError("Experimental density requires units and measurement definition.")
    mass = sum(p.molar_mass * m.count for m, p in zip(monomers, prepared, strict=True))
    if initial_packing_density is not None:
        density = positive(initial_packing_density, "initial_packing_density")
        length = (mass / AVOGADRO * 1e24 / density) ** (1 / 3)
        lengths = (length, length, length)
    else:
        lengths = validate_lengths(box_lengths)
        density = mass / AVOGADRO * 1e24 / float(np.prod(lengths))
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix="pack-", dir=out)).resolve()
    started = time.monotonic()
    metadata = {
        "initial_packing_density": density,
        "packing_separation": minimum_distance,
        "material_state": material_state,
        "random_seed": seed,
        "experimental_density": experimental_density,
        "requested_counts": [m.count for m in monomers],
        "component_names": [m.name for m in monomers],
    }
    result = POPResult(
        False, operation="pack", backend="packmol", output_paths=[run], metadata=metadata.copy()
    )
    try:
        lines = [
            "tolerance " + str(minimum_distance),
            "filetype xyz",
            "output packed.xyz",
            "pbc " + " ".join(str(v) for v in lengths),
            f"seed {seed}",
        ]
        state = instantiate(monomers, prepared, lengths)
        for monomer, template in zip(monomers, prepared, strict=True):
            filename = template.template_id + ".xyz"
            write_xyz(template.species, template.coordinates, run / filename, "Mofforge template")
            (run / (template.template_id + ".json")).write_text(
                json.dumps(template.to_dict(), indent=2, allow_nan=False) + "\n"
            )
            lines.extend([f"structure {filename}", f"  number {monomer.count}", "end structure"])
        text = "\n".join(lines) + "\n"
        (run / "packmol.inp").write_text(text)
        binary = PopBuildConfig.load(packmol_bin=packmol_bin).resolve_packmol_binary()
        version, transcript = probe_packmol(binary)
        (run / "packmol.version.txt").write_text(transcript)
        try:
            # Packmol rewinds stdin; a pipe fails with Illegal seek on some builds.
            with (run / "packmol.inp").open() as input_file:
                proc = subprocess.run(
                    [str(binary)],
                    stdin=input_file,
                    cwd=run,
                    capture_output=True,
                    text=True,
                    timeout=timeout,
                )
            stdout, stderr = proc.stdout, proc.stderr
        except subprocess.TimeoutExpired as exc:
            stdout = exc.stdout or b""
            stderr = exc.stderr or b""
            (run / "stdout.log").write_text(
                stdout.decode(errors="replace") if isinstance(stdout, bytes) else stdout
            )
            (run / "stderr.log").write_text(
                stderr.decode(errors="replace") if isinstance(stderr, bytes) else stderr
            )
            raise RuntimeError(
                f"Packmol timed out after {timeout} seconds; diagnostics in {run}"
            ) from exc
        (run / "stdout.log").write_text(stdout)
        (run / "stderr.log").write_text(stderr)
        if proc.returncode != 0:
            raise RuntimeError(f"Packmol exited {proc.returncode}; diagnostics in {run}")
        packed = run / "packed.xyz"
        species, coords = read_xyz(packed)
        if (
            species != [a.species for a in state.atoms]
            or len(packed.read_text().splitlines()) != len(state.atoms) + 2
        ):
            raise ValueError(
                "Packmol output atom count, element sequence, or record boundaries changed."
            )
        if coords.shape != (len(state.atoms), 3) or not np.isfinite(coords).all():
            raise ValueError("Malformed or nonfinite Packmol coordinates.")
        offset = 0
        for monomer, template in zip(monomers, prepared, strict=True):
            reference = pdist(template.coordinates)
            for _ in range(monomer.count):
                block = coords[offset : offset + template.n_atoms]
                if not np.allclose(pdist(block), reference, atol=0.002, rtol=0):
                    raise ValueError(
                        "Packmol altered intramolecular geometry or instance boundaries."
                    )
                offset += template.n_atoms
        state.coordinates = coords
        state.metadata = {**metadata, "actual_counts": [m.count for m in monomers]}
        state.provenance = [
            {
                "operation": "pack",
                "seed": seed,
                "packmol_version": version,
                "packmol_sha256": file_hash(binary),
                "software_versions": software_versions(),
                "templates": [p.provenance for p in prepared],
            }
        ]
        state = state.wrap()
        report, checks = validate_box(state, separation=minimum_distance, packing=True)
        result.validation = report
        result.metadata.update(checks)
        if not report.is_valid:
            raise ValueError("Invalid periodic packing: " + "; ".join(report.errors[:10]))
        state.metadata["validation"] = report.to_dict()
        result.state = state
        result.status = "completed"
        result.success = True
        persist_result(result, run)
    except (ValueError, OSError, RuntimeError, subprocess.SubprocessError, ImportError) as exc:
        result.success = False
        result.status = "failed"
        result.errors.append(str(exc))
    result.elapsed_seconds = time.monotonic() - started
    (run / "result.json").write_text(json.dumps(result.to_dict(), indent=2, allow_nan=False) + "\n")
    return result
