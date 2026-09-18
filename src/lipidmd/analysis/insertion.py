"""Insertion depth of each ionizable lipid relative to the membrane centre.

For every ionizable-lipid molecule and frame:
    com_dz_nm    z(molecule centre of mass) - z(membrane centre of mass)
    amine_dz_nm  same for the mean position of its protonatable N atoms
Both are wrapped into [-Lz/2, Lz/2). The aggregate starts on the +z side, so
values shrinking from ~+3 nm toward ~+2 nm (the headgroup region) or below
mean it is entering the upper leaflet.

Each row also carries the molecule's charge, so at 50% protonation the
charged and neutral MC3/SM-102 molecules can be compared directly.
Output: insertion.csv (time_ps, molecule, charge, com_dz_nm, amine_dz_nm).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .common import A_PER_NM, AnalysisContext, membrane_center_z, min_image_dz


def compute(ctx: AnalysisContext) -> dict[str, pd.DataFrame]:
    memb = ctx.membrane_atoms()
    agg = ctx.aggregate_atoms()
    amines = ctx.amine_atoms()
    residues = agg.residues
    charge = np.array([ctx.charge_of_resname()[r] for r in residues.resnames])
    # map each amine atom to its molecule's position in `residues`
    pos_of_res = {rix: i for i, rix in enumerate(residues.ix)}
    amine_mol = np.array([pos_of_res[a.residue.ix] for a in amines])
    counts = np.bincount(amine_mol, minlength=len(residues))

    rows = []
    for ts in ctx.frames():
        zc = membrane_center_z(memb)
        lz = ts.dimensions[2]
        com_z = agg.center_of_mass(compound="residues")[:, 2]
        am_z = np.bincount(amine_mol, weights=amines.positions[:, 2], minlength=len(residues)) / counts
        com_dz = min_image_dz(com_z - zc, lz) / A_PER_NM
        am_dz = min_image_dz(am_z - zc, lz) / A_PER_NM
        rows.append(pd.DataFrame({"time_ps": ts.time, "molecule": np.arange(len(residues)),
                                  "charge": charge, "com_dz_nm": com_dz, "amine_dz_nm": am_dz}))
    return {"insertion": pd.concat(rows, ignore_index=True)}
