"""Area per lipid (APL) and bilayer thickness of the target membrane.

APL       = box_x * box_y / (lipids per leaflet)
            The membrane spans the periodic box, so the box's xy area is the
            membrane area. Whether cholesterol counts in the denominator is
            set in analysis.yaml (apl_counts_cholesterol).
thickness = mean z of upper-leaflet P atoms - mean z of lower-leaflet P atoms
            ("D_PP", phosphate-to-phosphate distance).

Both are tracked frame by frame so their convergence can be inspected.
Output: membrane.csv with columns time_ps, observable, value.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .common import A_PER_NM, AnalysisContext, headgroup_atoms, membrane_center_z


def compute(ctx: AnalysisContext) -> dict[str, pd.DataFrame]:
    cfg_m = ctx.cfg.analysis["membrane"]
    count_chol = cfg_m.get("apl_counts_cholesterol")
    has_chol = "CHL1" in ctx.membrane_resnames
    if has_chol and count_chol is None:
        raise ValueError("analysis.yaml membrane.apl_counts_cholesterol is null (open decision)")
    memb = ctx.membrane_atoms()
    heads = headgroup_atoms(ctx)
    is_chol = heads.resnames == "CHL1"
    phos = heads.names == cfg_m["thickness_atom"]

    rows = []
    for ts in ctx.frames():
        zc = membrane_center_z(memb)
        z = heads.positions[:, 2]
        upper = z > zc
        counted = np.ones_like(upper) if (count_chol or not has_chol) else ~is_chol
        n_leaf = (counted & upper).sum(), (counted & ~upper).sum()
        area = ts.dimensions[0] * ts.dimensions[1] / A_PER_NM**2
        apl = area / np.mean(n_leaf)
        thick = (z[phos & upper].mean() - z[phos & ~upper].mean()) / A_PER_NM
        rows += [(ts.time, "area_per_lipid_nm2", apl), (ts.time, "thickness_dpp_nm", thick),
                 (ts.time, "box_area_nm2", area)]
    return {"membrane": pd.DataFrame(rows, columns=["time_ps", "observable", "value"])}
