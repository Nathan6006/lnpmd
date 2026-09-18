#!/usr/bin/env python3
"""CGenFF penalty checker.

    python scripts/check_cgenff.py inputs/cgenff/MC3/MC3H.str --lipid MC3
    python scripts/check_cgenff.py some.str --amine N1 --threshold 10
    python scripts/check_cgenff.py --all          # every .str configured in lipids.yaml

Lists every parameter and partial charge whose CGenFF penalty is above the
threshold (default 10), flags > 50 as ATTENTION, and marks with '*' anything
within --near-bonds bonds of an ionizable nitrogen. Exit code 1 if any
ATTENTION item touches an ionizable N.
"""
import _setup  # noqa: F401
import argparse
import sys
from pathlib import Path

from lipidmd.params import (ATTENTION_THRESHOLD, REVIEW_THRESHOLD, cgenff_files, check_penalties,
                            format_penalty_report, parse_str)


def report(path: Path, amines: list[str], threshold: float, near: int) -> bool:
    mol = parse_str(path)
    items = check_penalties(mol, amines, threshold, near)
    print(format_penalty_report(mol, items, amines, threshold, near))
    print()
    return any(it.near_amine and it.penalty > ATTENTION_THRESHOLD for it in items)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("str_file", nargs="?", type=Path)
    ap.add_argument("--lipid", help="lipid key in lipids.yaml (to look up amine atom names)")
    ap.add_argument("--amine", action="append", default=[], help="amine N atom name (repeatable)")
    ap.add_argument("--threshold", type=float, default=REVIEW_THRESHOLD)
    ap.add_argument("--near-bonds", type=int, default=3)
    ap.add_argument("--all", action="store_true", help="check every configured lipid/charge state")
    a = ap.parse_args()

    bad = False
    if a.all:
        from lipidmd.config import Config
        cfg = Config()
        for lip in cfg.systems["lipids"]:
            amines = [s["atom"] for s in cfg.lipid(lip)["protonatable_sites"] if s.get("atom")]
            for q in cfg.lipid(lip)["states"]:
                f = cgenff_files(cfg, lip, int(q))["str"]
                if not f.exists():
                    print(f"[{lip} {int(q):+d}] {f} not found -- skipped\n")
                    continue
                bad |= report(f, amines, a.threshold, a.near_bonds)
        return int(bad)
    if not a.str_file:
        ap.error("give a .str file or --all")
    amines = list(a.amine)
    if a.lipid:
        from lipidmd.config import Config
        amines += [s["atom"] for s in Config().lipid(a.lipid)["protonatable_sites"] if s.get("atom")]
    return int(report(a.str_file, amines, a.threshold, a.near_bonds))


if __name__ == "__main__":
    sys.exit(main())
