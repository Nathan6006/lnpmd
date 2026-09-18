#!/usr/bin/env python3
"""Validation harness: pure POPC vs published CHARMM36. Must pass first.

    python scripts/validate.py prepare [--dry-run] [--smoke] [--gmx gmx]
    python scripts/validate.py submit  [--smoke]
    python scripts/validate.py status  [--smoke]
    python scripts/validate.py analyze [--smoke]
    python scripts/validate.py compare

`compare` writes validate/report.txt and, only on PASS, validate/PASSED.
Production submission (scripts/run.py) checks for validate/PASSED.
"""
import _setup  # noqa: F401
import argparse
import sys

from lipidmd import validation as V
from lipidmd.analysis import membrane, order
from lipidmd.analysis.common import load_context
from lipidmd.build import BuildError, build_run
from lipidmd.config import Config, ConfigError
from lipidmd.run import prepare_run_dir, status_for_dir, submit


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["prepare", "submit", "status", "analyze", "compare"])
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--gmx", default="gmx")
    a = ap.parse_args()

    cfg = Config()
    if a.smoke:
        cfg = cfg.with_smoke_fallbacks()
    d = V.run_dir(cfg, a.smoke)
    mem = cfg.validation["membrane"]
    try:
        if a.step == "prepare":
            vcfg = cfg if a.smoke else V.validation_config(cfg)
            T = cfg.validation.get("temperature_K")
            if T is None and not a.smoke:
                raise ConfigError("validation.yaml temperature_K is null (open decision)")
            m = build_run(vcfg, None, d, gmx=a.gmx, membrane=mem, with_aggregate=False)
            prepare_run_dir(vcfg, d, f"{'smk' if a.smoke else 'md'}_validation", cfg.systems["seeds"]["velocities"][0],
                            a.smoke, placeholders=a.dry_run, temperature_K=T)
            print(f"built {d}: {m['n_atoms']} atoms. Next: python scripts/validate.py submit")
        elif a.step == "submit":
            print(f"submitted job {submit(d)}")
        elif a.step == "status":
            print(status_for_dir(d, "validation").row())
        elif a.step == "analyze":
            out = V.results_dir(cfg)
            ctx = load_context(cfg, d, out, "validation", 1, lipid=None, membrane=mem)
            for mod in (membrane, order):
                ctx.write(mod.compute(ctx))
            print(f"wrote {out}. Next: python scripts/validate.py compare")
        elif a.step == "compare":
            ok, text = V.compare(cfg)
            print(text)
            return 0 if ok else 1
    except (BuildError, ConfigError, FileNotFoundError, RuntimeError, ValueError) as e:
        print(f"{a.step} failed: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
