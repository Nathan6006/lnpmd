#!/usr/bin/env python3
"""Check and convert the CGenFF output for one or all ionizable lipids.

    python scripts/prepare_params.py --lipid MC3
    python scripts/prepare_params.py --all

For every charge state of the lipid:
  1. penalty report (same as check_cgenff.py)
  2. formula in the .mol2 == expected formula (neutral + added protons)
  3. convert .str -> GROMACS .itp/.prm with cgenff_charmm2gmx (if configured
     in paths.yaml; otherwise tells you the command to run)
  4. the .itp's atom names == neutral names + the added hydrogens, and its
     net charge == the state's charge
"""
import _setup  # noqa: F401
import argparse
import sys

from lipidmd import topology as top
from lipidmd.config import Config
from lipidmd.params import (cgenff_files, check_penalties, convert_to_gromacs, expected_formula,
                            find_param_file, format_penalty_report, mol2_formula, parse_str)


def prepare(cfg: Config, lip: str, convert: bool) -> bool:
    info = cfg.lipid(lip)
    sites = info["protonatable_sites"]
    amines = [s["atom"] for s in sites if s.get("atom")]
    order = info.get("protonation_order") or [s["id"] for s in sites]
    by_id = {s["id"]: s for s in sites}
    ok = True
    neutral_names = None
    for q in sorted(int(k) for k in info["states"]):
        f = cgenff_files(cfg, lip, q)
        print(f"==== {lip} charge {q:+d} (resname {f['resname']}) ====")
        if not f["str"].exists():
            print(f"  MISSING {f['str']} -- see MANUAL_STEPS.md step 2\n")
            ok = False
            continue
        mol = parse_str(f["str"])
        print(format_penalty_report(mol, check_penalties(mol, amines), amines))
        if mol.resname != f["resname"]:
            print(f"  ERROR: RESI in .str is {mol.resname}, lipids.yaml expects {f['resname']}")
            ok = False
        if abs(mol.net_charge - q) > 1e-3:
            print(f"  ERROR: .str net charge {mol.net_charge:+.3f} != {q:+d}")
            ok = False
        if f["mol2"].exists() and info.get("formula_neutral"):
            got, want = mol2_formula(f["mol2"]), expected_formula(info["formula_neutral"], q)
            print(f"  formula: mol2 {got}, expected {want}  {'OK' if got == want else 'MISMATCH'}")
            ok &= got == want
        if convert:
            try:
                convert_to_gromacs(cfg, lip, q)
                print(f"  converted -> {f['itp']}")
            except (RuntimeError, FileNotFoundError) as e:
                print(f"  conversion not done: {e}")
                ok = False
        try:
            itp = find_param_file(f["itp"])
        except FileNotFoundError as e:
            print(f"  {e}\n")
            ok = False
            continue
        mt = next(iter(top.parse_itp(itp).values()))
        if abs(mt.charge - q) > 1e-3:
            print(f"  ERROR: {itp.name} net charge {mt.charge:+.3f} != {q:+d}")
            ok = False
        if q == 0:
            neutral_names = set(mt.atom_names)
        elif neutral_names is not None and all(by_id[s].get("added_h") for s in order[:q]):
            want = neutral_names | {by_id[s]["added_h"] for s in order[:q]}
            if set(mt.atom_names) != want:
                print(f"  ERROR: atom names differ from neutral + added H: "
                      f"extra {sorted(set(mt.atom_names) - want)}, missing {sorted(want - set(mt.atom_names))}")
                ok = False
            else:
                print("  atom names consistent with the neutral form + added protons: OK")
        print()
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lipid")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--no-convert", action="store_true", help="only check, don't run cgenff_charmm2gmx")
    a = ap.parse_args()
    cfg = Config()
    lipids = cfg.systems["lipids"] if a.all else [a.lipid]
    if not lipids or lipids == [None]:
        ap.error("give --lipid or --all")
    ok = all([prepare(cfg, lip, not a.no_convert) for lip in lipids])
    print("ALL OK" if ok else "Problems found (see above).")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
