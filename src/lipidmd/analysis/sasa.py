"""Hydrophobic solvent-accessible surface area (SASA) of the ionizable-lipid aggregate.

Why: the working hypothesis is that protonation switches the lipids'
amphiphilicity. If charged headgroups pull the aggregate apart or reshape it,
more of its oily tail surface becomes exposed to water -- this is the direct
readout.

Method: Shrake-Rupley. Each atom is inflated to (vdW radius + probe radius);
points are spread evenly on that sphere; a point is "accessible" if it lies
inside no other inflated atom. Atom SASA = sphere area x accessible fraction.
The hydrophobic SASA is the sum over atoms matching
analysis.yaml sasa.hydrophobic_selection. Implemented here (numpy + scipy)
rather than adding another dependency.

Output: sasa.csv (time_ps, observable, value) with total and hydrophobic SASA
in nm^2.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .common import A_PER_NM, AnalysisContext

# Bondi (1964) van der Waals radii, nm, by element.
BONDI_NM = {"H": 0.120, "C": 0.170, "N": 0.155, "O": 0.152, "S": 0.180, "P": 0.180,
            "NA": 0.227, "CL": 0.175}


def sphere_points(n: int) -> np.ndarray:
    """n roughly evenly spaced unit vectors (golden-spiral construction)."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5**0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])


def shrake_rupley(xyz: np.ndarray, radii: np.ndarray, probe: float, n_points: int = 240,
                  box: np.ndarray | None = None, occluder_xyz: np.ndarray | None = None,
                  occluder_radii: np.ndarray | None = None, chunk: int = 200_000) -> np.ndarray:
    """Per-atom SASA (units of xyz squared).

    xyz, radii        atoms whose SASA we want
    occluder_*        extra atoms that can bury surface but aren't scored
    box               orthorhombic box lengths for periodic systems (or None)
    """
    all_xyz, all_r = xyz, radii
    if occluder_xyz is not None and len(occluder_xyz):
        all_xyz = np.vstack([xyz, occluder_xyz])
        all_r = np.concatenate([radii, occluder_radii])
    R = all_r + probe
    if box is not None:
        all_xyz = np.mod(all_xyz, box)
    tree = cKDTree(all_xyz, boxsize=box)
    unit = sphere_points(n_points)
    n = len(xyz)
    accessible = np.zeros(n)
    rmax = R.max()
    # build points for scored atoms in chunks to bound memory
    atoms_per_chunk = max(1, chunk // n_points)
    for a0 in range(0, n, atoms_per_chunk):
        a1 = min(n, a0 + atoms_per_chunk)
        owner = np.repeat(np.arange(a0, a1), n_points)
        pts = (all_xyz[a0:a1, None, :] + R[a0:a1, None, None] * unit[None]).reshape(-1, 3)
        if box is not None:
            pts = np.mod(pts, box)
        ptree = cKDTree(pts, boxsize=box)
        sdm = ptree.sparse_distance_matrix(tree, rmax, output_type="coo_matrix")
        buried_pair = (sdm.data < R[sdm.col] - 1e-9) & (sdm.col != owner[sdm.row])
        buried = np.zeros(len(pts), bool)
        buried[sdm.row[buried_pair]] = True
        accessible[a0:a1] = (~buried).reshape(a1 - a0, n_points).mean(axis=1)
    return 4 * np.pi * R[:n] ** 2 * accessible


def radii_for(atoms) -> np.ndarray:
    """Radius per atom from its element (tpr masses -> element via name)."""
    out = []
    for a in atoms:
        el = getattr(a, "element", "") or ""
        el = el.upper() if el else ("H" if a.mass < 1.5 else a.name[:1].upper())
        if el not in BONDI_NM:
            el = el[:1]
        if el not in BONDI_NM:
            raise ValueError(f"no radius for atom {a.name} ({el})")
        out.append(BONDI_NM[el])
    return np.array(out)


def compute(ctx: AnalysisContext) -> dict[str, pd.DataFrame]:
    s = ctx.cfg.analysis["sasa"]
    for k in ("probe_radius_nm", "hydrophobic_selection", "membrane_occludes"):
        if s.get(k) is None:
            raise ValueError(f"analysis.yaml sasa.{k} is null (open decision)")
    agg = ctx.aggregate_atoms()
    hyd = agg.select_atoms(s["hydrophobic_selection"])
    hyd_mask = np.isin(agg.indices, hyd.indices)
    r = radii_for(agg)
    occ = ctx.membrane_atoms() if s["membrane_occludes"] else None
    occ_r = radii_for(occ) if occ is not None else None
    probe = float(s["probe_radius_nm"])

    rows = []
    for ts in ctx.frames():
        box = ts.dimensions[:3] / A_PER_NM
        per_atom = shrake_rupley(agg.positions / A_PER_NM, r, probe, int(s["n_sphere_points"]), box,
                                 occ.positions / A_PER_NM if occ is not None else None, occ_r)
        rows += [(ts.time, "sasa_total_nm2", per_atom.sum()),
                 (ts.time, "sasa_hydrophobic_nm2", per_atom[hyd_mask].sum())]
    return {"sasa": pd.DataFrame(rows, columns=["time_ps", "observable", "value"])}
