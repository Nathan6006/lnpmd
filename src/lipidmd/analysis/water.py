"""Water density profile across the membrane, and full permeation events.

Density profile: number density of water molecules (oxygen positions) in
slabs along z, measured from the membrane centre, averaged over frames.
Bulk water is ~33.4 molecules/nm^3; the dip toward zero in the middle is the
hydrophobic core. A membrane being disrupted shows water leaking into it.

Permeation: z of each water relative to the membrane centre is classified as
    +1 upper bulk  (dz >  bulk_boundary)
    -1 lower bulk  (dz < -bulk_boundary)
     0 core        (|dz| < core_half_width)
and "in between" otherwise. A full permeation is: leave one bulk side, visit
the core, arrive in the OPPOSITE bulk side, without first returning. A water
that wraps through the periodic boundary (top of box -> bottom) jumps between
bulks without visiting the core and is correctly not counted.

Outputs: water_density.csv (z_nm, density_per_nm3; averaged over ALL frames),
water_core.csv (time_ps, observable=n_water_in_core, value) -- a per-frame
number the report can average after the equilibration cut, and
water_permeation.csv (one row per event: time_ps, water_residue_index,
direction +1 = lower->upper).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .common import A_PER_NM, AnalysisContext, membrane_center_z, min_image_dz


class PermeationCounter:
    """Vectorised state machine over all waters. Unit-testable without MDAnalysis."""

    def __init__(self, n_waters: int, core_half_width: float, bulk_boundary: float):
        if not 0 < core_half_width < bulk_boundary:
            raise ValueError("need 0 < core_half_width < bulk_boundary")
        self.core, self.bulk = core_half_width, bulk_boundary
        self.origin = np.zeros(n_waters, dtype=np.int8)       # last bulk side visited (0 = none yet)
        self.visited_core = np.zeros(n_waters, dtype=bool)
        self.events: list[tuple[float, int, int]] = []

    def update(self, time: float, dz: np.ndarray) -> None:
        in_core = np.abs(dz) < self.core
        side = np.where(dz > self.bulk, 1, np.where(dz < -self.bulk, -1, 0)).astype(np.int8)
        self.visited_core |= in_core & (self.origin != 0)
        arrived = side != 0
        crossed = arrived & (self.origin == -side) & self.visited_core
        for i in np.flatnonzero(crossed):
            self.events.append((time, int(i), int(side[i])))
        # entering any bulk region resets the journey
        self.origin[arrived] = side[arrived]
        self.visited_core[arrived] = False


def compute(ctx: AnalysisContext) -> dict[str, pd.DataFrame]:
    w = ctx.cfg.analysis["water"]
    for k in ("core_half_width_nm", "bulk_boundary_nm"):
        if w.get(k) is None:
            raise ValueError(f"analysis.yaml water.{k} is null (open decision)")
    u = ctx.universe
    memb = ctx.membrane_atoms()
    # water oxygens: first atom of each TIP3 residue
    ow = u.select_atoms("resname TIP3 and name OH2")
    bw = float(w["bin_width_nm"])
    half = u.dimensions[2] / A_PER_NM / 2
    edges = np.arange(-half, half + bw, bw)
    hist = np.zeros(len(edges) - 1)
    counter = PermeationCounter(ow.n_atoms, float(w["core_half_width_nm"]), float(w["bulk_boundary_nm"]))
    n_frames = 0
    core_rows = []
    for ts in ctx.frames():
        zc = membrane_center_z(memb)
        dz = min_image_dz(ow.positions[:, 2] - zc, ts.dimensions[2]) / A_PER_NM
        slab_vol = ts.dimensions[0] * ts.dimensions[1] / A_PER_NM**2 * bw
        hist += np.histogram(dz, edges)[0] / slab_vol
        counter.update(ts.time, dz)
        core_rows.append((ts.time, "n_water_in_core", int((np.abs(dz) < counter.core).sum())))
        n_frames += 1
    centers = (edges[:-1] + edges[1:]) / 2
    density = pd.DataFrame({"z_nm": centers, "density_per_nm3": hist / max(n_frames, 1)})
    res_ix = ow.resindices
    perm = pd.DataFrame([(t, int(res_ix[i]), d) for t, i, d in counter.events],
                        columns=["time_ps", "water_residue_index", "direction"])
    core = pd.DataFrame(core_rows, columns=["time_ps", "observable", "value"])
    return {"water_density": density, "water_core": core, "water_permeation": perm}
