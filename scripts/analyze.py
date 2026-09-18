#!/usr/bin/env python3
"""Run the analysis modules on finished runs (locally, with MDAnalysis).

    python scripts/analyze.py --systems MC3_p050 --replicas 1
    python scripts/analyze.py                         # everything that has a trajectory
    python scripts/analyze.py --observables membrane order

Needs, per run directory: prod.tpr, prod_centered.xtc, build_manifest.json
(see README step 6 for the rsync command). Writes results/<system>/rep<N>/*.csv.
Convergence diagnostics always run last, on whatever the other modules wrote.
"""
import _setup  # noqa: F401
import argparse
import importlib
import sys

from lipidmd.analysis import OBSERVABLES
from lipidmd.analysis.common import load_context
from lipidmd.config import Config
from lipidmd.run import run_dir_for


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--systems", nargs="*")
    ap.add_argument("--replicas", nargs="*", type=int)
    ap.add_argument("--observables", nargs="*", choices=OBSERVABLES, default=OBSERVABLES)
    ap.add_argument("--smoke", action="store_true", help="analyze runs_smoke/ instead of runs/")
    a = ap.parse_args()
    cfg = Config()
    if a.smoke:
        cfg = cfg.with_smoke_fallbacks()
    results = cfg.path(cfg.systems["results_dir"]) / ("smoke" if a.smoke else "")
    failures = 0
    for run in cfg.runs(a.systems, a.replicas):
        d = run_dir_for(cfg, run, a.smoke)
        if not (d / "prod_centered.xtc").exists():
            continue
        out = results / run.system.name / f"rep{run.replica}"
        ctx = load_context(cfg, d, out, run.system.name, run.replica, run.system.lipid)
        print(f"{run.name}: {ctx.universe.trajectory.n_frames} frames")
        for obs in a.observables + ["convergence"]:
            try:
                mod = importlib.import_module(f"lipidmd.analysis.{obs}")
                paths = ctx.write(mod.compute(ctx))
                print(f"   {obs:12s} -> {', '.join(p.name for p in paths)}")
            except ValueError as e:
                failures += 1
                print(f"   {obs:12s} SKIPPED: {e}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
