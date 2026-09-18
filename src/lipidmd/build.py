"""Assemble one simulation box: membrane + ionizable-lipid aggregate + water + ions.

Inputs (made by the manual steps in docs/MANUAL_STEPS.md):
    * CHARMM-GUI Membrane Builder output (gromacs/ folder): membrane
      coordinates (step5_input.gro), topology (topol.top) and the CHARMM36m
      parameter files (toppar/).
    * For each ionizable lipid and charge state: GROMACS .itp/.prm converted
      from the CGenFF .str.
    * A pre-equilibrated aggregate of NEUTRAL ionizable lipids (.pdb).

Steps (each numbered in build_run below):
    1  take the lipids (only) out of the CHARMM-GUI box, centre them in x/y
    2  read the aggregate, check its atom names against the neutral topology
    3  decide which molecules are protonated (protonation.py) and add those
       hydrogens geometrically; rename residues to the charged state
    4  place the aggregate `gap_nm` above the upper leaflet, size the box
    5  write topology; solvate with CHARMM TIP3P; remove water that landed
       inside the hydrophobic core of the bilayer
    6  add Na+/Cl-: 150 mM plus counterions to make the box neutral
    7  write index groups and a manifest of every choice made

Why not reuse CHARMM-GUI's water and ions? Adding the aggregate changes the
box height and the net charge, so both must be redone anyway.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import asdict
from pathlib import Path

import numpy as np

from . import topology as top
from .config import Config, ConfigError, RunSpec
from .mdp import ION_MDP
from .params import cgenff_files, find_param_file
from .protonation import MoleculeState, assign_protonation, summarize

N_H_BOND_NM = 0.101          # starting N-H length; minimization refines it
WATER_MOLARITY = 55.5        # mol/L of pure water, for the "water" salt basis


class BuildError(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# Geometry (pure functions -- unit tested)
# ---------------------------------------------------------------------------
def add_proton(res: top.Residue, n_name: str, neighbor_names: list[str], h_name: str) -> top.Residue:
    """Return a copy of `res` with a hydrogen added to nitrogen `n_name`.

    An amine N with three bonded neighbours has its lone pair pointing away
    from them; the new H goes along that direction: -(sum of unit vectors from
    N to each neighbour). Works for tertiary (3 C) and primary (1 C + 2 H)
    amines alike. Energy minimization fixes the exact geometry afterwards.
    """
    if len(neighbor_names) != 3:
        raise BuildError(f"{res.resname}:{n_name} has {len(neighbor_names)} bonded neighbours; "
                         f"expected 3 for a neutral amine N")
    n_xyz = res.xyz[res.index(n_name)]
    s = np.zeros(3)
    for nb in neighbor_names:
        v = res.xyz[res.index(nb)] - n_xyz
        s += v / np.linalg.norm(v)
    if np.linalg.norm(s) < 1e-3:
        raise BuildError(f"{res.resname}:{n_name} is planar; cannot place the proton")
    h_xyz = n_xyz - N_H_BOND_NM * s / np.linalg.norm(s)
    out = res.copy()
    out.names.append(h_name)
    out.xyz = np.vstack([out.xyz, h_xyz])
    return out


def reorder_to(res: top.Residue, atom_names: list[str], resname: str) -> top.Residue:
    """Put atoms in the exact order of the topology (GROMACS matches by order)."""
    if sorted(res.names) != sorted(atom_names):
        extra = set(res.names) - set(atom_names)
        missing = set(atom_names) - set(res.names)
        raise BuildError(f"atom names of {res.resname} don't match topology {resname}: "
                         f"in coordinates only {sorted(extra)}, in topology only {sorted(missing)}")
    idx = [res.names.index(n) for n in atom_names]
    return top.Residue(resname, list(atom_names), res.xyz[idx])


def headgroup_planes(residues: list[top.Residue], headgroup_atoms: dict[str, str]) -> dict:
    """Mean z of the upper and lower leaflet headgroup atoms, and the centre."""
    z = np.array([r.xyz[r.index(headgroup_atoms[r.resname])][2] for r in residues])
    center = z.mean()
    upper, lower = z[z > center], z[z <= center]
    return {"center": float(center), "upper": float(upper.mean()), "lower": float(lower.mean()),
            "n_upper": int(len(upper)), "n_lower": int(len(lower)), "z": z}


def leaflet_counts(residues: list[top.Residue], headgroup_atoms: dict[str, str]) -> dict[str, dict[str, int]]:
    planes = headgroup_planes(residues, headgroup_atoms)
    out: dict[str, dict[str, int]] = {"upper": {}, "lower": {}}
    for r, z in zip(residues, planes["z"]):
        side = "upper" if z > planes["center"] else "lower"
        out[side][r.resname] = out[side].get(r.resname, 0) + 1
    return out


def place_aggregate(agg: list[top.Residue], box_xy: np.ndarray, z_upper: float, gap: float) -> list[top.Residue]:
    """Centre the aggregate over the box in x/y; lowest atom `gap` above z_upper."""
    allxyz = np.vstack([r.xyz for r in agg])
    shift = np.array([box_xy[0] / 2, box_xy[1] / 2, 0.0]) - np.array([*allxyz[:, :2].mean(0), 0.0])
    shift[2] = z_upper + gap - allxyz[:, 2].min()
    out = []
    for r in agg:
        c = r.copy()
        c.xyz = c.xyz + shift
        out.append(c)
    return out


# ---------------------------------------------------------------------------
# GROMACS calls
# ---------------------------------------------------------------------------
def _gmx(gmx: str, args: list[str], cwd: Path, log: Path, stdin: str | None = None) -> None:
    cmd = [gmx, *args]
    proc = subprocess.run(cmd, cwd=cwd, input=stdin, capture_output=True, text=True)
    with open(log, "a") as fh:
        fh.write(f"\n$ {' '.join(cmd)}\n{proc.stdout}\n{proc.stderr}\n")
    if proc.returncode != 0:
        raise BuildError(f"`{' '.join(cmd)}` failed (see {log}):\n{proc.stderr[-3000:]}")


def gmx_data_dir(gmx: str) -> Path:
    """Where GROMACS keeps its data files (for spc216.gro)."""
    proc = subprocess.run([gmx, "-version"], capture_output=True, text=True)
    for line in (proc.stdout + proc.stderr).splitlines():
        if line.strip().startswith("Data prefix:"):
            return Path(line.split(":", 1)[1].strip()) / "share" / "gromacs" / "top"
    raise BuildError(f"could not find the GROMACS data directory from `{gmx} -version` "
                     f"-- is GROMACS installed and on PATH?")


def make_tip3_box(gmx: str, water: top.MoleculeType, out: Path) -> None:
    """GROMACS ships an equilibrated 216-water box (spc216.gro). Any 3-site
    water box works as starting coordinates for TIP3P; we only rename atoms
    and residue so they match CHARMM-GUI's TIP3 topology."""
    res, box = top.read_gro(gmx_data_dir(gmx) / "spc216.gro")
    if len(water.atom_names) != 3:
        raise BuildError(f"expected a 3-site water model, got {water.atom_names}")
    renamed = [top.Residue(water.name, list(water.atom_names), r.xyz) for r in res]
    top.write_gro(out, renamed, box, "TIP3 box from spc216")


# ---------------------------------------------------------------------------
# Parameter bookkeeping
# ---------------------------------------------------------------------------
def load_charmm_gui(cfg: Config, membrane: str) -> dict:
    d = cfg.path(cfg.membrane(membrane)["charmm_gui_dir"])
    gro, topfile, toppar = d / "step5_input.gro", d / "topol.top", d / "toppar"
    for p in (gro, topfile, toppar):
        if not p.exists():
            raise BuildError(f"{p} not found. Download and unzip the CHARMM-GUI output as described "
                             f"in docs/MANUAL_STEPS.md step 1.")
    includes = top.parse_top_includes(topfile)
    moltypes: dict[str, top.MoleculeType] = {}
    for inc in includes:
        moltypes.update(top.parse_itp(d / inc))
    return {"dir": d, "gro": gro, "includes": includes, "moltypes": moltypes}


def resname_map(moltypes: dict[str, top.MoleculeType]) -> dict[str, str]:
    """residue name -> moleculetype name (they're usually identical)."""
    return {m.resnames[0]: name for name, m in moltypes.items() if m.resnames}


def _find_by_resname(moltypes: dict[str, top.MoleculeType], resname: str, what: str) -> top.MoleculeType:
    for m in moltypes.values():
        if m.resnames and m.resnames[0] == resname:
            return m
    raise BuildError(f"no {what} topology with residue name {resname} in the CHARMM-GUI toppar. "
                     f"When building the membrane, keep the water and choose NaCl ions "
                     f"(MANUAL_STEPS.md) so TIP3/SOD/CLA topologies are included.")


def check_restraint_macros(cfg: Config, cg: dict) -> None:
    rm = cfg.md["restraint_macros"]
    text = "".join((cg["dir"] / inc).read_text() for inc in cg["includes"])
    missing = [rm[k] for k in ("lipid_pos_switch", "lipid_pos_fc", "lipid_dih_switch", "lipid_dih_fc")
               if rm[k] not in text]
    if missing:
        raise BuildError(f"restraint macros {missing} (md.yaml restraint_macros) do not appear in the "
                         f"CHARMM-GUI lipid topologies; open toppar/POPC.itp and copy the names used "
                         f"in its #ifdef blocks into md.yaml.")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
def build_run(cfg: Config, run: RunSpec | None, out_dir: Path, gmx: str = "gmx",
              membrane: str | None = None, with_aggregate: bool = True) -> dict:
    """Build start.gro, topol.top and index.ndx in out_dir. Returns the manifest.

    run=None / with_aggregate=False builds a membrane-only box (validation).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    log = out_dir / "build.log"
    log.write_text("lipidmd build log\n")
    membrane = membrane or cfg.target_membrane()
    hg = cfg.analysis["membrane"]["headgroup_atoms"]
    manifest: dict = {"membrane": membrane}

    # 1 ---- membrane lipids from CHARMM-GUI --------------------------------
    cg = load_charmm_gui(cfg, membrane)
    check_restraint_macros(cfg, cg)
    counts_cfg = cfg.membrane(membrane)["counts_per_leaflet"]
    residues, box = top.read_gro(cg["gro"])
    lipid_names = set(cfg.membrane(membrane)["ratio"])
    memb = [r for r in residues if r.resname in lipid_names]
    unknown = {r.resname for r in residues} - lipid_names - {"TIP3", "SOD", "CLA", "POT"}
    if unknown:
        raise BuildError(f"CHARMM-GUI box contains residues {unknown} that are not in the membrane "
                         f"definition {sorted(lipid_names)}")
    shift = np.array([box[0] / 2, box[1] / 2, 0]) - np.array([*np.vstack([r.xyz for r in memb])[:, :2].mean(0), 0])
    for r in memb:
        r.xyz = r.xyz + shift
    found = leaflet_counts(memb, hg)
    manifest["leaflet_counts"] = found
    if counts_cfg:
        for side in ("upper", "lower"):
            if found[side] != {k: int(v) for k, v in counts_cfg.items()}:
                raise BuildError(f"{side} leaflet of the CHARMM-GUI membrane has {found[side]}, "
                                 f"membranes.yaml says {counts_cfg}")
    planes = headgroup_planes(memb, hg)

    # 2-3 ---- aggregate + protonation ---------------------------------------
    agg: list[top.Residue] = []
    lipid_itps: list[Path] = []
    lipid_moltypes: dict[str, top.MoleculeType] = {}
    if with_aggregate:
        lip = run.system.lipid
        info = cfg.lipid(lip)
        sites = {s["id"]: s for s in info["protonatable_sites"]}
        for sid, s in sites.items():
            if not s.get("atom") or not s.get("added_h"):
                raise ConfigError(f"lipids.yaml {lip}: site {sid} needs `atom` and `added_h` names")
        states = assign_protonation(int(cfg.systems["n_molecules"]), info.get("protonation_order"),
                                    run.system.fraction, run.protonation_seed, n_sites=len(sites))
        used_charges = sorted({s.charge for s in states} | {0})
        prm_files = []
        for q in used_charges:
            f = cgenff_files(cfg, lip, q)
            itp, prm = find_param_file(f["itp"]), find_param_file(f["prm"])
            mt = next(iter(top.parse_itp(itp).values()))
            lipid_moltypes[f["resname"]] = mt
            if q in {s.charge for s in states}:
                lipid_itps.append(itp)
            prm_files.append(prm)
            if abs(mt.charge - q) > 1e-3:
                raise BuildError(f"{itp}: net charge {mt.charge:+.3f} but it is meant to be the {q:+d} state")

        neutral = lipid_moltypes[info["states"][0]["resname"]]
        agg_file = cfg.path(cfg.require("systems.aggregate.file").format(lipid=lip))
        if not agg_file.exists():
            raise BuildError(f"aggregate structure {agg_file} not found (see MANUAL_STEPS.md)")
        raw = top.read_pdb(agg_file)
        if len(raw) != int(cfg.systems["n_molecules"]):
            raise BuildError(f"{agg_file} has {len(raw)} molecules, systems.yaml n_molecules = "
                             f"{cfg.systems['n_molecules']}")
        for st, res in zip(states, raw):
            res = reorder_to(res, neutral.atom_names, info["states"][0]["resname"])
            for sid in st.protonated_sites:
                s = sites[sid]
                res = add_proton(res, s["atom"], neutral.neighbors(s["atom"]), s["added_h"])
            target = info["states"][st.charge]["resname"]
            agg.append(reorder_to(res, lipid_moltypes[target].atom_names, target))
        agg = place_aggregate(agg, box[:2], planes["upper"], float(cfg.require("systems.aggregate.gap_nm")))
        manifest["protonation"] = {"lipid": lip, "fraction": run.system.fraction,
                                   "seed": run.protonation_seed, "charge_counts": summarize(states),
                                   "molecules": [asdict(s) for s in states]}

    # 4 ---- box height -----------------------------------------------------
    pad = float(cfg.require("systems.aggregate.water_padding_nm")) if with_aggregate else \
        float(box[2] - (planes["upper"] - planes["lower"])) / 2
    top_z = np.vstack([r.xyz for r in agg])[:, 2].max() if agg else planes["upper"]
    z0 = planes["lower"] - pad
    new_box = np.array([box[0], box[1], top_z + pad - z0])
    solute = memb + agg
    for r in solute:
        r.xyz[:, 2] -= z0
    zP = {"lower": planes["lower"] - z0, "upper": planes["upper"] - z0}
    manifest["box_nm"] = new_box.round(4).tolist()

    # 5 ---- topology + water ------------------------------------------------
    toppar_dst = out_dir / "toppar"
    if toppar_dst.exists():
        shutil.rmtree(toppar_dst)
    shutil.copytree(cg["dir"] / "toppar", toppar_dst)
    lipdir = out_dir / "ionizable"
    lipdir.mkdir(exist_ok=True)
    ff = [i for i in cg["includes"] if "forcefield" in Path(i).name]
    if len(ff) != 1:
        raise BuildError(f"expected exactly one forcefield .itp include in the CHARMM-GUI topol.top, found {ff}")
    rest = [i for i in cg["includes"] if i != ff[0]]
    includes = list(ff)                                # force field first (defaults + atom types)
    extra = cfg.path(cfg.paths["params_dir"]) / "extra_atomtypes.itp"
    if with_aggregate and extra.exists():
        # CGenFF atom types missing from CHARMM-GUI's forcefield.itp, copied by
        # hand from charmm36-jul2022.ff (MANUAL_STEPS.md troubleshooting).
        shutil.copy(extra, lipdir / extra.name)
        includes.append(f"ionizable/{extra.name}")
    if with_aggregate:
        (lipdir / "ionizable.prm").write_text(top.merge_prm(prm_files))
        includes.append("ionizable/ionizable.prm")       # bonded types must precede moleculetypes
    includes += rest
    rm = cfg.md["restraint_macros"]
    for itp in lipid_itps:
        shutil.copy(itp, lipdir / itp.name)
        mt = next(iter(top.parse_itp(itp).values()))
        top.write_posre_itp(lipdir / f"{mt.name}_posre.itp", mt, rm["aggregate_pos_switch"], rm["aggregate_pos_fc"])
        includes += [f"ionizable/{itp.name}", f"ionizable/{mt.name}_posre.itp"]
    all_types = cg["moltypes"] | {m.name: m for m in lipid_moltypes.values()}
    _check_atomtypes(out_dir, includes, lipid_moltypes)

    water = _find_by_resname(cg["moltypes"], "TIP3", "water")
    sod = _find_by_resname(cg["moltypes"], "SOD", "sodium")
    cla = _find_by_resname(cg["moltypes"], "CLA", "chloride")
    rmap = resname_map(all_types)

    top.write_gro(out_dir / "solute.gro", solute, new_box, "membrane + aggregate")
    make_tip3_box(gmx, water, out_dir / "tip3_216.gro")
    _gmx(gmx, ["solvate", "-cp", "solute.gro", "-cs", "tip3_216.gro", "-o", "solvated_raw.gro"], out_dir, log)
    solvated, sbox = top.read_gro(out_dir / "solvated_raw.gro")
    # gmx solvate fills every gap, including inside the bilayer's oily core.
    # Remove any water whose oxygen lies between the two phosphate planes.
    o = water.atom_names[0]
    keep = [r for r in solvated if r.resname != water.name or not (zP["lower"] < r.xyz[r.index(o)][2] < zP["upper"])]
    n_removed = len(solvated) - len(keep)
    top.write_gro(out_dir / "solvated.gro", keep, sbox, "solvated")
    n_water = sum(r.resname == water.name for r in keep)
    manifest["waters"] = {"added": n_water, "removed_from_core": n_removed}

    # 6 ---- ions -------------------------------------------------------------
    q = sum(all_types[rmap[r.resname]].charge for r in keep)
    q_int = int(round(q))
    if abs(q - q_int) > 1e-3:
        raise BuildError(f"total charge {q:.4f} is not an integer -- a topology is inconsistent")
    conc = float(cfg.md["salt_molar"])
    basis = cfg.require("md.salt_basis")
    if basis == "water":
        n_pairs = int(round(conc / WATER_MOLARITY * n_water))
    else:  # whole box volume, what `gmx genion -conc` would do
        n_pairs = int(round(conc * 6.02214076e23 * float(np.prod(sbox)) * 1e-24))
    n_na, n_cl = n_pairs + max(0, -q_int), n_pairs + max(0, q_int)
    manifest["ions"] = {"solute_charge": q_int, "salt_basis": basis, "salt_pairs": n_pairs,
                        "SOD": n_na, "CLA": n_cl}
    top.write_top(out_dir / "topol.top", includes, top.run_length(keep, rmap))
    (out_dir / "ions.mdp").write_text(ION_MDP)
    wat_idx = _atom_indices(keep, {water.name})
    top.write_ndx(out_dir / "ions.ndx", {"WATER": wat_idx})
    _gmx(gmx, ["grompp", "-f", "ions.mdp", "-c", "solvated.gro", "-p", "topol.top", "-o", "ions.tpr",
               "-maxwarn", "0"], out_dir, log)
    _gmx(gmx, ["genion", "-s", "ions.tpr", "-n", "ions.ndx", "-o", "ionized.gro", "-pname", sod.name,
               "-nname", cla.name, "-np", str(n_na), "-nn", str(n_cl),
               "-seed", str(cfg.systems["seeds"]["ions"])], out_dir, log, stdin="WATER\n")
    final, fbox = top.read_gro(out_dir / "ionized.gro")
    top.write_gro(out_dir / "start.gro", final, fbox, "lipidmd start structure")
    top.write_top(out_dir / "topol.top", includes, top.run_length(final, rmap))
    q_final = sum(all_types[rmap[r.resname]].charge for r in final)
    if abs(q_final) > 1e-3:
        raise BuildError(f"final system charge is {q_final:+.4f}, expected 0")

    # 7 ---- index groups -------------------------------------------------------
    agg_names = set(lipid_moltypes)
    groups = {
        "System": list(range(sum(len(r.names) for r in final))),
        "MEMB": _atom_indices(final, lipid_names),
        "AGG": _atom_indices(final, agg_names),
        "MEMB_AGG": _atom_indices(final, lipid_names | agg_names),
        "SOLV": _atom_indices(final, {water.name, sod.name, cla.name}),
    }
    if not groups["AGG"]:
        del groups["AGG"]
        groups["MEMB_AGG"] = groups["MEMB"]
    top.write_ndx(out_dir / "index.ndx", groups)
    manifest["n_atoms"] = len(groups["System"])
    (out_dir / "build_manifest.json").write_text(json.dumps(manifest, indent=2, default=_jsonable))
    for tmp in ("solvated_raw.gro", "tip3_216.gro", "ions.tpr", "mdout.mdp"):
        (out_dir / tmp).unlink(missing_ok=True)
    return manifest


def _atom_indices(residues: list[top.Residue], resnames: set[str]) -> list[int]:
    out, i = [], 0
    for r in residues:
        n = len(r.names)
        if r.resname in resnames:
            out.extend(range(i, i + n))
        i += n
    return out


def _check_atomtypes(out_dir: Path, includes: list[str], lipid_moltypes: dict[str, top.MoleculeType]) -> None:
    """Every atom type used by the ionizable lipids must be defined somewhere,
    otherwise grompp fails with a much less helpful message."""
    defined: set[str] = set()
    for inc in includes:
        p = out_dir / inc
        if p.exists():
            defined |= top.parse_atomtype_names(p)
    for p in (out_dir / "toppar").glob("*.itp"):
        defined |= top.parse_atomtype_names(p)
    for name, mt in lipid_moltypes.items():
        missing = sorted(set(mt.atom_types) - defined)
        if missing:
            raise BuildError(
                f"{name} uses atom types {missing} that the CHARMM-GUI force field files do not define. "
                f"See MANUAL_STEPS.md 'Troubleshooting: missing atom types'.")


def _jsonable(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, MoleculeState):
        return asdict(o)
    return str(o)
