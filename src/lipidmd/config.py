"""Load and sanity-check the YAML files in config/.

Everything is kept as plain dicts (what PyYAML returns) so you can print any
piece of it and see exactly what is in the file. `Config` adds a few helpers
on top: the list of systems, the list of runs, and `unresolved()`, which finds
every value still set to `null` (an open decision or a missing fact).
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]

CONFIG_FILES = ("lipids", "membranes", "systems", "md", "hpc", "analysis", "validation", "paths")

# nulls that legitimately mean "not applicable / use the default" and are
# therefore not reported as unresolved.
OPTIONAL_NULLS = {
    "hpc.account",
    "hpc.gpu_type",
    "hpc.mem",
}


class ConfigError(ValueError):
    """A config file is inconsistent or a required value is missing."""


@dataclass(frozen=True)
class SystemSpec:
    """One cell of the 3x3 matrix, e.g. MC3 at 50% protonation."""

    lipid: str
    fraction: float

    @property
    def name(self) -> str:
        return f"{self.lipid}_p{round(self.fraction * 100):03d}"


@dataclass(frozen=True)
class RunSpec:
    """One replica of one system -> one directory, one SLURM job chain."""

    system: SystemSpec
    replica: int  # 1-based
    velocity_seed: int
    protonation_seed: int

    @property
    def name(self) -> str:
        return f"{self.system.name}/rep{self.replica}"


# ---------------------------------------------------------------------------
# Values used ONLY by --smoke (1000-step pipeline test) when the real value is
# still null. They let you test the plumbing before the science is decided.
# They are printed loudly whenever used and never apply to real runs.
# ---------------------------------------------------------------------------
SMOKE_FALLBACKS: dict[str, Any] = {
    "md.salt_basis": "water",
    "md.thermostat": "V-rescale",
    "md.tau_t_ps": 1.0,
    "md.coupling_groups": ["MEMB_AGG", "SOLV"],
    "md.barostat_equilibration": "C-rescale",
    "md.barostat_production": "C-rescale",
    "md.tau_p_ps": 5.0,
    "md.compressibility_per_bar": 4.5e-5,
    "md.output.xtc_interval_ps": 0.2,
    "md.output.energy_interval_ps": 0.2,
    "md.output.log_interval_ps": 0.2,
    "md.minimization.restraints.lipid_pos": 1000.0,
    "md.minimization.restraints.lipid_dih": 1000.0,
    "md.minimization.restraints.aggregate_pos": 1000.0,
    "md.equilibration": [
        {"name": "eq1", "ensemble": "nvt", "length_ps": 2.0,
         "restraints": {"lipid_pos": 1000.0, "lipid_dih": 1000.0, "aggregate_pos": 1000.0}},
        {"name": "eq2", "ensemble": "npt", "length_ps": 2.0,
         "restraints": {"lipid_pos": 100.0, "lipid_dih": 100.0, "aggregate_pos": 100.0}},
    ],
    "md.production.length_ns": 0.002,
    "systems.replica_variation": "velocities",
    "systems.aggregate.gap_nm": 1.0,
    "systems.aggregate.water_padding_nm": 2.0,
    "lipids.ECO.protonation_order": ["amine_primary", "amine_tertiary"],
}


def _get(d: Any, dotted: str) -> Any:
    for part in dotted.split("."):
        if not isinstance(d, dict) or part not in d:
            return None
        d = d[part]
    return d


def _set(d: dict, dotted: str, value: Any) -> None:
    parts = dotted.split(".")
    for part in parts[:-1]:
        d = d.setdefault(part, {})
    d[parts[-1]] = value


def _find_nulls(node: Any, prefix: str) -> list[str]:
    """Dotted paths of every None leaf below `node`."""
    out: list[str] = []
    if node is None:
        return [prefix]
    if isinstance(node, dict):
        for k, v in node.items():
            out += _find_nulls(v, f"{prefix}.{k}" if prefix else str(k))
    elif isinstance(node, list):
        for i, v in enumerate(node):
            out += _find_nulls(v, f"{prefix}[{i}]")
    return out


class Config:
    def __init__(self, config_dir: Path | str | None = None, root: Path | str | None = None):
        self.root = Path(root) if root else REPO_ROOT
        self.config_dir = Path(config_dir) if config_dir else self.root / "config"
        self.data: dict[str, dict] = {}
        for name in CONFIG_FILES:
            path = self.config_dir / f"{name}.yaml"
            if not path.exists():
                raise ConfigError(f"missing config file: {path}")
            with open(path) as fh:
                self.data[name] = yaml.safe_load(fh) or {}
        self.smoke_fallbacks_used: list[str] = []
        self.check()

    # convenient attribute access: cfg.md, cfg.lipids, ...
    def __getattr__(self, item: str) -> dict:
        data = self.__dict__.get("data", {})
        if item in data:
            return data[item]
        raise AttributeError(item)

    def path(self, rel: str | Path) -> Path:
        """Resolve a repo-relative path from config."""
        p = Path(rel)
        return p if p.is_absolute() else self.root / p

    def get(self, dotted: str) -> Any:
        """cfg.get("md.output.xtc_interval_ps")"""
        head, _, rest = dotted.partition(".")
        return _get(self.data[head], rest) if rest else self.data[head]

    def require(self, dotted: str) -> Any:
        """Like get(), but a null is an error that names the file to edit."""
        value = self.get(dotted)
        if value is None:
            head = dotted.split(".")[0]
            raise ConfigError(
                f"'{dotted}' is not set (null) in config/{head}.yaml. This is an open "
                f"decision -- see the comment next to it and docs/DECISIONS.md.")
        return value

    # ---------------------------------------------------------------- smoke
    def with_smoke_fallbacks(self) -> "Config":
        """Copy of this config with SMOKE_FALLBACKS filling any nulls."""
        new = copy.copy(self)
        new.data = copy.deepcopy(self.data)
        new.smoke_fallbacks_used = []
        for dotted, value in SMOKE_FALLBACKS.items():
            head, _, rest = dotted.partition(".")
            if head == "lipids" and rest.split(".")[0] not in new.data["lipids"]:
                continue
            if _get(new.data[head], rest) is None:
                _set(new.data[head], rest, copy.deepcopy(value))
                new.smoke_fallbacks_used.append(dotted)
        return new

    # --------------------------------------------------------------- matrix
    def systems_list(self) -> list[SystemSpec]:
        s = self.systems
        return [SystemSpec(lip, float(f)) for lip in s["lipids"] for f in s["protonation_levels"]]

    def find_system(self, name: str) -> SystemSpec:
        for spec in self.systems_list():
            if spec.name == name:
                return spec
        raise ConfigError(f"unknown system {name!r}; known: {[s.name for s in self.systems_list()]}")

    def runs(self, systems: list[str] | None = None, replicas: list[int] | None = None) -> list[RunSpec]:
        s = self.systems
        base = int(s["seeds"]["protonation"])
        variation = s.get("replica_variation")
        out = []
        specs = [self.find_system(n) for n in systems] if systems else self.systems_list()
        for spec in specs:
            for r in range(1, int(s["replicas"]) + 1):
                if replicas and r not in replicas:
                    continue
                # Same protonation pattern for every replica unless the config
                # says replicas should also differ in which molecules are charged.
                pseed = base + (r - 1) if variation == "velocities+protonation" else base
                out.append(RunSpec(spec, r, int(s["seeds"]["velocities"][r - 1]), pseed))
        return out

    # ------------------------------------------------------------ membranes
    def target_membrane(self) -> str:
        return self.membranes["target"]

    def membrane(self, name: str | None = None) -> dict:
        name = name or self.target_membrane()
        try:
            return self.membranes["membranes"][name]
        except KeyError:
            raise ConfigError(f"membrane {name!r} not defined in membranes.yaml") from None

    def membrane_counts(self, name: str | None = None) -> dict[str, int]:
        """Lipids of each type per leaflet. Refuses to round a ratio for you."""
        m = self.membrane(name)
        if m.get("counts_per_leaflet"):
            counts = {k: int(v) for k, v in m["counts_per_leaflet"].items()}
        else:
            total = sum(m["ratio"].values())
            n = m["lipids_per_leaflet"]
            exact = {k: n * v / total for k, v in m["ratio"].items()}
            if any(abs(x - round(x)) > 1e-9 for x in exact.values()):
                raise ConfigError(
                    f"membrane {name or self.target_membrane()!r}: ratio {m['ratio']} of {n} lipids "
                    f"gives non-integer counts {exact}. Set counts_per_leaflet explicitly.")
            counts = {k: int(round(x)) for k, x in exact.items()}
        if sum(counts.values()) != m["lipids_per_leaflet"]:
            raise ConfigError(f"counts_per_leaflet {counts} do not sum to {m['lipids_per_leaflet']}")
        return counts

    # --------------------------------------------------------------- lipids
    def lipid(self, key: str) -> dict:
        try:
            return self.lipids[key]
        except KeyError:
            raise ConfigError(f"lipid {key!r} not in lipids.yaml") from None

    def n_sites(self, key: str) -> int:
        return len(self.lipid(key)["protonatable_sites"])

    # ------------------------------------------------------------ checking
    def check(self) -> None:
        """Structural consistency. Nulls are allowed here (see unresolved())."""
        s = self.systems
        for lip in s["lipids"]:
            info = self.lipid(lip)
            sites = info.get("protonatable_sites") or []
            if not sites:
                raise ConfigError(f"lipid {lip}: protonatable_sites is empty")
            ids = [site["id"] for site in sites]
            if len(set(ids)) != len(ids):
                raise ConfigError(f"lipid {lip}: duplicate site ids {ids}")
            states = {int(k) for k in info["states"]}
            if states != set(range(len(sites) + 1)):
                raise ConfigError(
                    f"lipid {lip}: has {len(sites)} protonatable sites so needs states "
                    f"{list(range(len(sites) + 1))}, found {sorted(states)}")
            order = info.get("protonation_order")
            if order is not None and sorted(order) != sorted(ids):
                raise ConfigError(f"lipid {lip}: protonation_order {order} must list each of {ids} once")
        for f in s["protonation_levels"]:
            if not 0.0 <= float(f) <= 1.0:
                raise ConfigError(f"protonation level {f} is not a fraction between 0 and 1")
        if int(s["replicas"]) < 1:
            raise ConfigError("replicas must be >= 1")
        if len(s["seeds"]["velocities"]) != int(s["replicas"]):
            raise ConfigError(
                f"systems.yaml: seeds.velocities has {len(s['seeds']['velocities'])} entries "
                f"but replicas = {s['replicas']}")
        if s.get("replica_variation") not in (None, "velocities", "velocities+protonation"):
            raise ConfigError(f"replica_variation must be 'velocities' or 'velocities+protonation'")
        if self.md.get("salt_basis") not in (None, "water", "box"):
            raise ConfigError("md.yaml: salt_basis must be 'water' or 'box'")
        self.membrane()  # target exists
        names = [spec.name for spec in self.systems_list()]
        if len(set(names)) != len(names):
            raise ConfigError(f"duplicate system names: {names}")

    def unresolved(self) -> list[str]:
        """Every null value, as 'file.dotted.path'. Empty list = nothing open."""
        out = []
        for name in CONFIG_FILES:
            for p in _find_nulls(self.data[name], ""):
                full = f"{name}.{p}"
                if full in OPTIONAL_NULLS:
                    continue
                # Literature pKa fields are informational and don't block anything.
                if name == "lipids" and ".pka_apparent." in full:
                    continue
                out.append(full)
        return out
