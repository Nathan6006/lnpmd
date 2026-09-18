"""Every analysis module run end-to-end on a tiny synthetic system.

Skipped automatically if MDAnalysis isn't installed. Checks the plumbing
(selections, units, CSV shapes) and a few values that are known exactly.
"""
import copy

import numpy as np
import pandas as pd
import pytest

mda = pytest.importorskip("MDAnalysis")

from lipidmd.analysis import contacts, convergence, insertion, membrane, order, sasa, water  # noqa: E402
from lipidmd.analysis.common import AnalysisContext  # noqa: E402
from lipidmd.config import Config  # noqa: E402

BOX = [40.0, 40.0, 80.0, 90, 90, 90]  # Angstrom
ZC = 40.0


def _system():
    """8 POPC (4 per leaflet), 2 POPS (upper), 2 MC3 + 2 MC3H above, 3 waters."""
    res = []   # (resname, [(name, mass, xyz)], bonds as local index pairs)
    for i in range(8):
        up = i < 4
        zp = ZC + 20 if up else ZC - 20
        x, y = 5 + 8 * (i % 4), 10
        sgn = -1 if up else 1   # tails point toward the centre
        res.append(("POPC", [("P", 30.97, (x, y, zp)), ("C22", 12.0, (x, y, zp + sgn * 5)),
                             ("H2R", 1.008, (x, y, zp + sgn * 6)),     # C-H parallel to z -> S = 1
                             ("H2S", 1.008, (x + 1, y, zp + sgn * 5))],  # C-H perpendicular -> S = -0.5
                    [(1, 2), (1, 3)]))
    for i in range(2):
        res.append(("POPS", [("P", 30.97, (10 + 20 * i, 30, ZC + 20)), ("O13A", 16.0, (10 + 20 * i, 30, ZC + 21))], []))
    for i, rn in enumerate(["MC3", "MC3H", "MC3", "MC3H"]):
        x = 10 + 20 * (i % 2)
        # molecule 1 (MC3H) sits with its N 3 A from a POPS carboxylate O -> contact
        zN = ZC + 24 if i == 1 else ZC + 35
        res.append((rn, [("N1", 14.0, (x, 30, zN)), ("C1", 12.0, (x, 30, zN + 1.5))], [(0, 1)]))
    for i in range(3):
        res.append(("TIP3", [("OH2", 16.0, (5 + 10 * i, 35, 70)), ("H1", 1.008, (6 + 10 * i, 35, 70)),
                             ("H2", 1.008, (5 + 10 * i, 36, 70))], []))

    atoms, resindex, bonds, resnames = [], [], [], []
    for ri, (rn, ats, bs) in enumerate(res):
        base = len(atoms)
        atoms += ats
        resindex += [ri] * len(ats)
        bonds += [(base + a, base + b) for a, b in bs]
        resnames.append(rn)
    u = mda.Universe.empty(len(atoms), n_residues=len(res), atom_resindex=resindex, trajectory=True)
    u.add_TopologyAttr("names", [a[0] for a in atoms])
    u.add_TopologyAttr("masses", [a[1] for a in atoms])
    u.add_TopologyAttr("resnames", resnames)
    u.add_TopologyAttr("resids", list(range(1, len(res) + 1)))
    u.add_TopologyAttr("bonds", bonds)
    xyz = np.array([a[2] for a in atoms], float)
    from MDAnalysis.coordinates.memory import MemoryReader
    u.load_new(np.stack([xyz, xyz]), format=MemoryReader, dimensions=BOX, dt=100.0)
    return u


@pytest.fixture
def ctx(tmp_path):
    cfg = copy.copy(Config())
    cfg.data = copy.deepcopy(cfg.data)
    a = cfg.data["analysis"]
    a["membrane"]["apl_counts_cholesterol"] = True
    a["order_parameters"]["tails"] = {"POPC": {"sn2": ["C22"]}}
    a["contacts"].update(cutoff_nm=0.4, gap_tolerance_frames=0)
    a["sasa"].update(probe_radius_nm=0.14, hydrophobic_selection="name C*", membrane_occludes=False,
                     n_sphere_points=100)
    a["water"].update(core_half_width_nm=1.0, bulk_boundary_nm=2.5)
    cfg.data["lipids"]["MC3"]["protonatable_sites"][0]["atom"] = "N1"
    c = AnalysisContext(cfg, tmp_path, tmp_path, "MC3_p050", 1, lipid="MC3", membrane="popc_pops_chol",
                        universe=_system())
    return c


def test_membrane(ctx):
    df = membrane.compute(ctx)["membrane"]
    apl = df[df.observable == "area_per_lipid_nm2"].value
    assert np.allclose(apl, 16.0 / 5)               # 16 nm^2 box / mean(6 upper, 4 lower)
    assert np.allclose(df[df.observable == "thickness_dpp_nm"].value, 4.0)


def test_order(ctx):
    df = order.compute(ctx)["order"]
    assert set(df.carbon) == {"C22"}
    assert np.allclose(df.scd, (1.0 + -0.5) / 2)


def test_insertion(ctx):
    df = insertion.compute(ctx)["insertion"]
    assert len(df) == 2 * 4
    assert set(df.charge) == {0, 1}
    zc = ctx.membrane_atoms().center_of_mass()[2]   # not exactly 40: POPS only in the upper leaflet
    first = df[(df.time_ps == 0) & (df.molecule == 1)].iloc[0]
    assert first.amine_dz_nm == pytest.approx((ZC + 24 - zc) / 10)


def test_contacts(ctx):
    out = contacts.compute(ctx)
    ts = out["contacts"]
    charged = ts[(ts.charge == 1)]
    assert charged.n_amines_in_contact.tolist() == [1, 1]
    assert ts[ts.charge == 0].n_amines_in_contact.sum() == 0
    ev = out["contact_events"]
    assert len(ev) == 1 and bool(ev.censored.iloc[0])


def test_sasa(ctx):
    df = sasa.compute(ctx)["sasa"]
    tot = df[df.observable == "sasa_total_nm2"].value.iloc[0]
    hyd = df[df.observable == "sasa_hydrophobic_nm2"].value.iloc[0]
    assert 0 < hyd < tot


def test_water_and_convergence(ctx):
    w = water.compute(ctx)
    dens = w["water_density"]
    assert dens.density_per_nm3.sum() * 0.1 * 16 == pytest.approx(3.0, rel=1e-6)  # 3 waters
    assert len(w["water_permeation"]) == 0
    ctx.write({**membrane.compute(ctx), **w})
    conv = convergence.compute(ctx)
    assert "area_per_lipid_nm2" in set(conv["convergence_running"].observable)
