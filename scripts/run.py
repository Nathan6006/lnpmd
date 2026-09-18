#!/usr/bin/env python3
"""Build systems, write every input file, and submit to SLURM.

    python scripts/run.py --dry-run                   # all 27 runs: build + inputs, no submission
    python scripts/run.py --smoke --systems MC3_p050 --replicas 1
    python scripts/run.py --systems MC3_p050 ECO_p100 # real submission (needs validation passed)

Flags
  --dry-run     generate everything, submit nothing (unset hpc.yaml values
                become TODO_ placeholders in job.sh so you can read it)
  --smoke       1000 steps per stage, written to runs_smoke/. Open decisions
                are filled with clearly-labelled placeholders (config.py
                SMOKE_FALLBACKS) so the plumbing can be tested early.
  --skip-build  only (re)write mdp files + job.sh in already-built directories
  --rebuild     rebuild even if start.gro already exists
  --gmx         GROMACS executable used for BUILDING (default: gmx)
  --skip-validation-check   submit production even without validate/PASSED
"""
import _setup  # noqa: F401
import argparse
import shutil
import sys

from lipidmd.build import BuildError, build_run
from lipidmd.config import Config, ConfigError
from lipidmd.run import prepare_run_dir, run_dir_for, submit
from lipidmd.validation import require_passed


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--systems", nargs="*", help="e.g. MC3_p050 (default: all)")
    ap.add_argument("--replicas", nargs="*", type=int, help="e.g. 1 2 (default: all)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--skip-build", action="store_true")
    ap.add_argument("--rebuild", action="store_true")
    ap.add_argument("--gmx", default="gmx")
    ap.add_argument("--skip-validation-check", action="store_true")
    a = ap.parse_args()

    cfg = Config()
    if a.smoke:
        cfg = cfg.with_smoke_fallbacks()
        if cfg.smoke_fallbacks_used:
            print("*" * 78)
            print("SMOKE TEST: these open decisions are filled with PLACEHOLDER values that are")
            print("not scientific choices and must never be used for real runs:")
            for k in cfg.smoke_fallbacks_used:
                print(f"   {k}")
            print("*" * 78)
    if not a.dry_run and not a.smoke and not a.skip_validation_check:
        require_passed(cfg)
    if not a.skip_build and not shutil.which(a.gmx):
        print(f"'{a.gmx}' not found: building needs GROMACS (see README). Use --skip-build to only "
              f"write mdp/job files for already-built runs.")
        return 1

    runs = cfg.runs(a.systems, a.replicas)
    failures = skipped = 0
    for run in runs:
        d = run_dir_for(cfg, run, a.smoke)
        label = f"{run.name:20s}"
        try:
            if not a.skip_build and (a.rebuild or not (d / "start.gro").exists()):
                m = build_run(cfg, run, d, gmx=a.gmx)
                print(f"{label} built: {m['n_atoms']} atoms, charges {m['protonation']['charge_counts']}, "
                      f"ions {m['ions']['SOD']} Na+ / {m['ions']['CLA']} Cl-")
            elif not (d / "start.gro").exists():
                print(f"{label} not built yet -- skipping (drop --skip-build)")
                skipped += 1
                continue
            prepare_run_dir(cfg, d, f"{'smk' if a.smoke else 'md'}_{run.system.name}_r{run.replica}",
                            run.velocity_seed, a.smoke, placeholders=a.dry_run)
            if a.dry_run:
                print(f"{label} inputs written to {d} (dry run: not submitted)")
            else:
                print(f"{label} submitted as job {submit(d)}")
        except (BuildError, ConfigError, FileNotFoundError, RuntimeError, ValueError) as e:
            failures += 1
            print(f"{label} FAILED: {e}")
    done = len(runs) - failures - skipped
    print(f"\n{done}/{len(runs)} runs prepared" + (" (dry run)" if a.dry_run else "")
          + (f", {skipped} skipped" if skipped else "") + (f", {failures} failed" if failures else ""))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
