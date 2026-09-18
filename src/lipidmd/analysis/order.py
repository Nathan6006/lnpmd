"""Deuterium order parameters S_CD per tail carbon.

    S_CD = < (3 cos^2 theta - 1) / 2 >

theta = angle between a C-H bond and the membrane normal (z). The average is
over every C-H bond on that carbon, every lipid, and (in the report) time.
S_CD is what 2H-NMR measures, so it is the standard check that lipid tails
have the right degree of order. Fully ordered, all-trans tails give
S_CD = -0.5; isotropic (liquid-like) tails give 0. Values are usually
reported as -S_CD (positive); we keep the sign.

Hydrogens are found from the bonds in the .tpr, so this works regardless of
hydrogen naming. Output: order.csv (time_ps, resname, chain, carbon, index, scd).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .common import AnalysisContext


def _ch_pairs(ag_lipids, tails: dict[str, list[str]]):
    """(carbon index, hydrogen index, chain, carbon name, position) for every C-H bond."""
    out = []
    for res in ag_lipids.residues:
        for chain, carbons in tails.items():
            for pos, cname in enumerate(carbons, start=1):
                c = res.atoms.select_atoms(f"name {cname}")
                if c.n_atoms != 1:
                    raise ValueError(f"{res.resname} {res.resid}: carbon {cname} not found")
                for h in c[0].bonded_atoms:
                    if h.mass < 1.5:  # hydrogen
                        out.append((c[0].index, h.index, chain, cname, pos))
    return out


def compute(ctx: AnalysisContext) -> dict[str, pd.DataFrame]:
    u = ctx.universe
    tails_cfg = ctx.cfg.analysis["order_parameters"]["tails"]
    frames = []
    for resname, tails in tails_cfg.items():
        if resname not in ctx.membrane_resnames:
            continue
        pairs = _ch_pairs(u.select_atoms(f"resname {resname}"), tails)
        if not pairs:
            continue
        ci = np.array([p[0] for p in pairs])
        hi = np.array([p[1] for p in pairs])
        keys = pd.DataFrame([(p[2], p[3], p[4]) for p in pairs], columns=["chain", "carbon", "index"])
        group = keys.groupby(["chain", "carbon", "index"], sort=False).ngroup().to_numpy()
        labels = keys.drop_duplicates().reset_index(drop=True)
        n_groups = group.max() + 1
        for ts in ctx.frames():
            v = u.atoms.positions[hi] - u.atoms.positions[ci]
            cos2 = v[:, 2] ** 2 / np.einsum("ij,ij->i", v, v)
            s = 0.5 * (3 * cos2 - 1)
            mean = np.bincount(group, weights=s, minlength=n_groups) / np.bincount(group, minlength=n_groups)
            df = labels.copy()
            df.insert(0, "resname", resname)
            df.insert(0, "time_ps", ts.time)
            df["scd"] = mean
            frames.append(df)
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
        columns=["time_ps", "resname", "chain", "carbon", "index", "scd"])
    return {"order": out}
