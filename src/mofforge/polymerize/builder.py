"""Public facade for explicit molecular packing and native connections."""

from __future__ import annotations

from pathlib import Path

from mofforge.core.crystal import Crystal
from mofforge.polymerize.base import Connector, Monomer, reject_removed
from mofforge.polymerize.provision import instantiate
from mofforge.provenance import derive_seed, effective_seed


class PopBuilder:
    def __init__(self, backend="native", *, packmol_bin=None, **options):
        reject_removed(options)
        if options:
            raise TypeError(f"Unknown builder options: {sorted(options)}")
        if backend != "native":
            raise ValueError(f"Unknown backend {backend!r}; use native (direct Packmol).")
        self.packmol_bin = packmol_bin
        self._monomers = []
        self._seed = None

    @property
    def backend_name(self):
        return "native"

    def add_monomer(
        self,
        source: str | Path | Crystal,
        *,
        count: int,
        name=None,
        connectors=None,
        functionality=None,
        graph=None,
        **options,
    ):
        reject_removed(options)
        if options:
            raise TypeError(f"Unknown component options: {sorted(options)}")
        monomer = Monomer(
            name or f"component{len(self._monomers) + 1}",
            source,
            count,
            connectors or [],
            functionality,
            graph,
        )
        self._monomers.append(monomer)
        return monomer

    def set_connectors(self, monomer, connectors):
        """Attach explicitly mapped sites after inspecting prepare() template labels."""
        if not any(m is monomer for m in self._monomers):
            raise ValueError("Monomer is not registered in this builder.")
        monomer.connectors = [Connector(**c) if isinstance(c, dict) else c for c in connectors]
        monomer.functionality = len(monomer.connectors)
        monomer.__post_init__()

    def list_monomers(self):
        return [m.name for m in self._monomers]

    def prepare(self, *, random_seed=None):
        """Return prepared templates and addressable labels, without running Packmol."""
        from mofforge.polymerize.monomer import prepare_monomer

        if not self._monomers:
            raise ValueError("No monomers registered; call add_monomer() first.")
        if type(random_seed) is bool:
            raise ValueError("random_seed must be an integer.")
        self._seed = effective_seed(self._seed if random_seed is None else random_seed)
        result = []
        for i, monomer in enumerate(self._monomers):
            monomer.__post_init__()
            identity = f"t{i:04d}"
            result.append(
                prepare_monomer(
                    monomer, template_id=identity, random_seed=derive_seed(self._seed, identity)
                )
            )
        return result

    def pack(self, output_dir=".", **options):
        from mofforge.polymerize.engine import pack

        reject_removed(options)
        prepared = self.prepare(random_seed=options.pop("random_seed", None))
        return pack(
            self._monomers,
            prepared,
            output_dir=output_dir,
            packmol_bin=self.packmol_bin,
            random_seed=self._seed,
            **options,
        )

    def build(
        self,
        output_dir=".",
        *,
        rules=None,
        target_conversion=None,
        candidate_attempt_budget=None,
        min_nonbonded_distance=None,
        allow_cycles=False,
        **packing_options,
    ):
        from mofforge.polymerize.connect import (
            check_options,
            check_rules,
            connect,
            rules_from_input,
        )
        from mofforge.polymerize.engine import pack

        reject_removed(packing_options)
        rules = rules_from_input(rules)
        if any(
            v is None for v in (target_conversion, candidate_attempt_budget, min_nonbonded_distance)
        ):
            raise ValueError(
                "build requires target_conversion, candidate_attempt_budget, "
                "and min_nonbonded_distance."
            )
        check_options(
            rules, target_conversion, candidate_attempt_budget, min_nonbonded_distance, allow_cycles
        )
        prepared = self.prepare(random_seed=packing_options.pop("random_seed", None))
        # Resolve chemistry edits against one copy of each template before invoking Packmol.
        preview = instantiate(self._monomers, prepared, (100.0, 100.0, 100.0), one_copy=True)
        check_rules(preview, rules)
        packed = pack(
            self._monomers,
            prepared,
            output_dir=output_dir,
            packmol_bin=self.packmol_bin,
            random_seed=self._seed,
            **packing_options,
        )
        if not packed.success:
            return packed
        result = connect(
            packed.state,
            rules,
            target_conversion=target_conversion,
            candidate_attempt_budget=candidate_attempt_budget,
            min_nonbonded_distance=min_nonbonded_distance,
            allow_cycles=allow_cycles,
            output_dir=output_dir,
        )
        result.output_paths = packed.output_paths + result.output_paths
        result.elapsed_seconds += packed.elapsed_seconds
        return result


def run_config(config: dict, *, operation="pack", output_dir="."):
    """One JSON request shape shared by Python, CLI, and MCP.

    Keys: components, packing, connection, packmol_bin. A saved state_path may
    replace components/packing for connection-only requests.
    """
    from mofforge.polymerize.connect import connect
    from mofforge.polymerize.state import ConstructionState

    if not isinstance(config, dict):
        raise ValueError("Box configuration must be a JSON object.")
    reject_removed(config)
    unknown = set(config) - {"components", "packing", "connection", "packmol_bin", "state_path"}
    if unknown:
        raise ValueError(f"Unknown box configuration keys: {sorted(unknown)}")
    if operation not in ("pack", "connect"):
        raise ValueError("operation must be pack or connect.")
    if "state_path" in config:
        if operation != "connect" or set(config) & {"components", "packing", "packmol_bin"}:
            raise ValueError("state_path is for connection-only requests without packing inputs.")
        return connect(
            ConstructionState.load(config["state_path"]),
            output_dir=output_dir,
            **config.get("connection", {}),
        )
    if operation == "pack" and "connection" in config:
        raise ValueError("Packing requests cannot include connection options; use polymerize.")
    builder = PopBuilder(packmol_bin=config.get("packmol_bin"))
    for component in config.get("components", []):
        builder.add_monomer(**component)
    if operation == "pack":
        return builder.pack(output_dir=output_dir, **config.get("packing", {}))
    return builder.build(
        output_dir=output_dir, **config.get("packing", {}), **config.get("connection", {})
    )
