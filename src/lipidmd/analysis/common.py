"""Shared analysis plumbing: loading a run, selections, membrane centre.

Trajectory conventions (set up by the job script's final trjconv step):
    prod.tpr            topology incl. masses, charges and bonds
    prod_centered.xtc   molecules made whole, membrane centred in the box

All distances in the output CSVs are in nm and times in ps (GROMACS units);
MDAnalysis works in Angstrom internally, so we convert at the edges.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import Config

A_PER_NM = 10.0


@dataclass
class AnalysisContext:
    cfg: Config
    run_dir: Path
    out_dir: Path
    system: str
    replica: int
    lipid: str | None = None          # None for a membrane-only (validation) run
    membrane: str | None = None
    universe: object = None
    manifest: dict = field(default_factory=dict)

    # --- selections -------------------------------------------------------
    @property
    def membrane_resnames(self) -> list[str]:
        return list(self.cfg.membrane(self.membrane)["ratio"])

    def membrane_atoms(self):
        return self.universe.select_atoms("resname " + " ".join(self.membrane_resnames))

    def lipid_resnames(self) -> list[str]:
        if not self.lipid:
            return []
        return [s["resname"] for s in self.cfg.lipid(self.lipid)["states"].values()]

    def aggregate_atoms(self):
        names = self.lipid_resnames()
        return self.universe.select_atoms("resname " + " ".join(names)) if names else None

    def amine_atoms(self):
        """The protonatable nitrogens of every ionizable lipid molecule."""
        sites = self.cfg.lipid(self.lipid)["protonatable_sites"]
        names = [s["atom"] for s in sites]
        if any(n is None for n in names):
            raise ValueError(f"lipids.yaml: protonatable atom names for {self.lipid} are not set")
        return self.universe.select_atoms(
            f"resname {' '.join(self.lipid_resnames())} and name {' '.join(names)}")

    def charge_of_resname(self) -> dict[str, int]:
        if not self.lipid:
            return {}
        return {s["resname"]: int(q) for q, s in self.cfg.lipid(self.lipid)["states"].items()}

    # --- helpers ------------------------------------------------------------
    def frames(self):
        stride = int(self.cfg.analysis.get("stride", 1))
        return self.universe.trajectory[::stride]

    def write(self, tables: dict[str, pd.DataFrame]) -> list[Path]:
        self.out_dir.mkdir(parents=True, exist_ok=True)
        paths = []
        for name, df in tables.items():
            df.insert(0, "replica", self.replica)
            df.insert(0, "system", self.system)
            p = self.out_dir / f"{name}.csv"
            df.to_csv(p, index=False)
            paths.append(p)
        return paths


def load_context(cfg: Config, run_dir: Path, out_dir: Path, system: str, replica: int,
                 lipid: str | None, membrane: str | None = None,
                 trajectory: str = "prod_centered.xtc", topology: str = "prod.tpr") -> AnalysisContext:
    import MDAnalysis as mda  # imported here so non-analysis code doesn't need it

    tpr, xtc = run_dir / topology, run_dir / trajectory
    for p in (tpr, xtc):
        if not p.exists():
            raise FileNotFoundError(f"{p} not found (copy it from the cluster; see README step 6)")
    u = mda.Universe(str(tpr), str(xtc))
    manifest = {}
    if (run_dir / "build_manifest.json").exists():
        manifest = json.loads((run_dir / "build_manifest.json").read_text())
    membrane = membrane or manifest.get("membrane") or cfg.target_membrane()
    return AnalysisContext(cfg, run_dir, out_dir, system, replica, lipid, membrane, u, manifest)


def membrane_center_z(memb_atoms) -> float:
    """z of the membrane's centre of mass (Angstrom). Valid because trjconv
    centred the membrane, so it is not split across the periodic boundary."""
    return float(memb_atoms.center_of_mass()[2])


def min_image_dz(dz: np.ndarray, lz: float) -> np.ndarray:
    """Wrap z differences into [-Lz/2, Lz/2): a molecule that drifts out of
    the top of the box is really approaching the membrane from below."""
    return (dz + lz / 2) % lz - lz / 2


def headgroup_atoms(ctx: AnalysisContext):
    """One headgroup atom per membrane lipid (P for phospholipids, O3 for
    cholesterol). Their z relative to the membrane centre assigns leaflets."""
    hg = ctx.cfg.analysis["membrane"]["headgroup_atoms"]
    sel = " or ".join(f"(resname {r} and name {a})" for r, a in hg.items() if r in ctx.membrane_resnames)
    heads = ctx.universe.select_atoms(sel)
    if heads.n_atoms != ctx.membrane_atoms().n_residues:
        raise ValueError(f"found {heads.n_atoms} headgroup atoms for {ctx.membrane_atoms().n_residues} "
                         f"membrane lipids -- check analysis.yaml membrane.headgroup_atoms")
    return heads
