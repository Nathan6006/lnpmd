"""build.build_run end to end, with a fake `gmx` and toy topologies.

Checks the bookkeeping that is easy to get wrong: protonated residues get
their hydrogen and new name, the topology's [ molecules ] matches the
coordinates exactly, the box ends up neutral, water is removed from the
bilayer core, and index groups cover every atom.
"""
import copy
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from lipidmd import topology as top
from lipidmd.build import build_run
from lipidmd.config import Config

FIX = Path(__file__).parent / "fixtures"


def itp(name, atoms, bonds=(), extra=""):
    lines = ["[ moleculetype ]", f"{name} 3", "[ atoms ]"]
    for i, (an, at, q, m) in enumerate(atoms, 1):
        lines.append(f"{i} {at} 1 {name} {an} {i} {q} {m}")
    if bonds:
        lines.append("[ bonds ]")
        lines += [f"{a} {b} 1" for a, b in bonds]
    return "\n".join(lines) + "\n" + extra


RESTR = "#ifdef POSRES\n; POSRES_FC_LIPID\n#endif\n#ifdef DIHRES\n; DIHRES_FC\n#endif\n"


@pytest.fixture
def world(tmp_path, monkeypatch):
    # --- fake GROMACS data dir with a small water box
    monkeypatch.setenv("FAKE_GMX_PREFIX", str(tmp_path / "gmxdata"))
    wat_dir = tmp_path / "gmxdata" / "share" / "gromacs" / "top"
    wat_dir.mkdir(parents=True)
    g = np.arange(3) * 0.62 + 0.1
    waters = [top.Residue("SOL", ["OW", "HW1", "HW2"], np.array([[x, y, z], [x + .09, y, z], [x, y + .09, z]]))
              for x in g for y in g for z in g]
    top.write_gro(wat_dir / "spc216.gro", waters, np.array([1.86, 1.86, 1.86]))

    # --- fake CHARMM-GUI output: 2x2 lipids per leaflet (POPC 2, POPS 1, CHL1 1)
    cg = tmp_path / "charmm-gui" / "gromacs"
    (cg / "toppar").mkdir(parents=True)
    (cg / "toppar" / "forcefield.itp").write_text(
        "[ defaults ]\n1 2 yes 1.0 1.0\n[ atomtypes ]\n" +
        "\n".join(f"{t} 6 12.0 0.0 A 0.3 0.3" for t in ["PL", "OS", "CH", "OW", "HW", "SOD", "CLA", "NG301", "NG3P1", "CG331", "HGP2"]) + "\n")
    (cg / "toppar" / "POPC.itp").write_text(itp("POPC", [("P", "PL", 0.0, 31.0), ("C2", "CH", 0.0, 12.0)], [(1, 2)], RESTR))
    (cg / "toppar" / "POPS.itp").write_text(itp("POPS", [("P", "PL", -1.0, 31.0), ("C2", "CH", 0.0, 12.0)], [(1, 2)], RESTR))
    (cg / "toppar" / "CHL1.itp").write_text(itp("CHL1", [("O3", "OS", 0.0, 16.0), ("C2", "CH", 0.0, 12.0)], [(1, 2)], RESTR))
    (cg / "toppar" / "TIP3.itp").write_text(itp("TIP3", [("OH2", "OW", -0.834, 16.0), ("H1", "HW", 0.417, 1.0), ("H2", "HW", 0.417, 1.0)]))
    (cg / "toppar" / "SOD.itp").write_text(itp("SOD", [("SOD", "SOD", 1.0, 23.0)]))
    (cg / "toppar" / "CLA.itp").write_text(itp("CLA", [("CLA", "CLA", -1.0, 35.0)]))
    (cg / "topol.top").write_text("\n".join(f'#include "toppar/{n}.itp"' for n in
                                            ["forcefield", "POPC", "POPS", "CHL1", "SOD", "CLA", "TIP3"]) + "\n")
    memb = []
    for zP, sgn in ((6.0, -1), (2.0, 1)):
        for k, rn in enumerate(["POPC", "POPC", "POPS", "CHL1"]):
            x, y = 0.8 + 1.6 * (k % 2), 0.8 + 1.6 * (k // 2)
            head = "O3" if rn == "CHL1" else "P"
            memb.append(top.Residue(rn, [head, "C2"], np.array([[x, y, zP], [x, y, zP + sgn * 1.5]])))
    memb.append(top.Residue("TIP3", ["OH2", "H1", "H2"], np.array([[1, 1, 7.5], [1.1, 1, 7.5], [1, 1.1, 7.5]])))
    top.write_gro(cg / "step5_input.gro", memb, np.array([3.2, 3.2, 8.0]))

    # --- MC3 neutral / protonated topologies (N1 bonded to three methyls)
    par = tmp_path / "params" / "MC3"
    par.mkdir(parents=True)
    neutral = [("N1", "NG301", -0.3, 14.0), ("C1", "CG331", 0.1, 12.0), ("C2", "CG331", 0.1, 12.0), ("C3", "CG331", 0.1, 12.0)]
    prot = [("N1", "NG3P1", -0.2, 14.0), ("HN1", "HGP2", 0.3, 1.0), ("C1", "CG331", 0.3, 12.0), ("C2", "CG331", 0.3, 12.0), ("C3", "CG331", 0.3, 12.0)]
    (par / "MC3.itp").write_text(itp("MC3", neutral, [(1, 2), (1, 3), (1, 4)]))
    (par / "MC3H.itp").write_text(itp("MC3H", prot, [(1, 2), (1, 3), (1, 4), (1, 5)]))
    (par / "MC3.prm").write_text("[ bondtypes ]\nNG301 CG331 1 0.147 200000\n")
    (par / "MC3H.prm").write_text("[ bondtypes ]\nNG3P1 CG331 1 0.149 200000\nCG331 NG301 1 0.147 200000\n")

    # --- aggregate: 4 neutral MC3 in a 2x2 grid, atoms deliberately out of order
    lines = []
    n = 0
    for m in range(4):
        cx, cy = 1.0 + 1.2 * (m % 2), 1.0 + 1.2 * (m // 2)
        for name, d in (("C2", (-0.7, 1.2, -0.5)), ("N1", (0, 0, 0)), ("C1", (1.4, 0, -0.5)), ("C3", (-0.7, -1.2, -0.5))):
            n += 1
            x, y, z = cx * 10 + d[0], cy * 10 + d[1], 5 + d[2]
            lines.append(f"ATOM  {n:5d} {name:<4s} MC3  {m + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00")
    agg = tmp_path / "agg.pdb"
    agg.write_text("\n".join(lines) + "\nEND\n")

    cfg = Config().with_smoke_fallbacks()
    cfg.root = tmp_path
    cfg.data["paths"]["params_dir"] = str(tmp_path / "params")
    cfg.data["membranes"]["membranes"]["popc_pops_chol"].update(
        charmm_gui_dir=str(cg), counts_per_leaflet={"POPC": 2, "POPS": 1, "CHL1": 1}, lipids_per_leaflet=4)
    cfg.data["systems"]["n_molecules"] = 4
    cfg.data["systems"]["aggregate"]["file"] = str(agg)
    cfg.data["lipids"]["MC3"]["protonatable_sites"][0].update(atom="N1", added_h="HN1")
    gmx = f"{sys.executable} {FIX / 'fake_gmx.py'}"
    wrapper = tmp_path / "gmx"
    wrapper.write_text(f"#!/bin/sh\nexec {gmx} \"$@\"\n")
    wrapper.chmod(0o755)
    return cfg, str(wrapper), tmp_path


def test_build_mc3_half_protonated(world):
    cfg, gmx, tmp = world
    run = cfg.runs(systems=["MC3_p050"], replicas=[1])[0]
    out = tmp / "run"
    m = build_run(cfg, run, out, gmx=gmx)

    assert m["protonation"]["charge_counts"] == {0: 2, 1: 2}
    res, box = top.read_gro(out / "start.gro")
    mc3h = [r for r in res if r.resname == "MC3H"]
    assert len(mc3h) == 2 and all(r.names == ["N1", "HN1", "C1", "C2", "C3"] for r in mc3h)
    assert all(r.names == ["N1", "C1", "C2", "C3"] for r in res if r.resname == "MC3")

    # [ molecules ] must reproduce the coordinate order exactly
    text = (out / "topol.top").read_text()
    mol_lines = text.split("[ molecules ]")[1].strip().splitlines()[1:]
    expanded = [name for line in mol_lines for name in [line.split()[0]] * int(line.split()[1])]
    assert expanded == [r.resname for r in res]

    # includes: forcefield, then merged prm, before any moleculetype
    inc = top.parse_top_includes(out / "topol.top")
    assert inc[0].endswith("forcefield.itp") and inc[1] == "ionizable/ionizable.prm"
    assert "ionizable/MC3H_posre.itp" in inc

    # neutral overall: 2 POPS (-2) + 2 MC3H (+2) + ions
    charge = {"POPC": 0, "POPS": -1, "CHL1": 0, "MC3": 0, "MC3H": 1, "TIP3": 0, "SOD": 1, "CLA": -1}
    assert sum(charge[r.resname] for r in res) == 0
    assert m["ions"]["solute_charge"] == 0

    # no water inside the bilayer (between the phosphate planes)
    zP = sorted(r.xyz[0, 2] for r in res if r.resname in ("POPC", "POPS"))
    lo, hi = np.mean(zP[:3]), np.mean(zP[3:])
    assert not any(lo < r.xyz[0, 2] < hi for r in res if r.resname == "TIP3")

    # aggregate sits gap_nm above the upper phosphate plane
    agg_z = np.vstack([r.xyz for r in res if r.resname.startswith("MC3")])[:, 2]
    assert agg_z.min() - hi == pytest.approx(cfg.systems["aggregate"]["gap_nm"], abs=2e-3)

    # index groups partition the system
    ndx = (out / "index.ndx").read_text()
    for g in ("[ MEMB ]", "[ AGG ]", "[ MEMB_AGG ]", "[ SOLV ]"):
        assert g in ndx
    manifest = json.loads((out / "build_manifest.json").read_text())
    assert manifest["n_atoms"] == sum(len(r.names) for r in res)


def test_build_membrane_only(world):
    cfg, gmx, tmp = world
    m = build_run(cfg, None, tmp / "val", gmx=gmx, membrane="popc_pops_chol", with_aggregate=False)
    res, _ = top.read_gro(tmp / "val" / "start.gro")
    assert not any(r.resname.startswith("MC3") for r in res)
    assert m["ions"]["SOD"] - m["ions"]["CLA"] == 2    # counterions for 2 POPS
