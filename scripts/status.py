#!/usr/bin/env python3
"""Per run: state (not-built / built / QUEUED / RUNNING / DONE / FAILED /
stalled), current stage, and simulated production time.

    python scripts/status.py
    python scripts/status.py --smoke

"stalled" = the job is gone from the queue but never recorded finishing or
failing (usually killed by SLURM before mdrun could checkpoint). Look at the
newest slurm-*.out in the run directory, then resubmit with
    sbatch runs/<system>/rep<N>/job.sh
It will continue from the last checkpoint.
"""
import _setup  # noqa: F401
import argparse
import sys

from lipidmd.config import Config
from lipidmd.run import run_dir_for, status_for_dir
from lipidmd.validation import run_dir as validation_dir


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    cfg = Config()
    print(f"{'run':22s} {'state':9s} {'stage':6s} {'production time':>18s}  job")
    vd = validation_dir(cfg, a.smoke)
    if vd.exists():
        print(status_for_dir(vd, "validation").row())
    counts = {}
    for run in cfg.runs():
        st = status_for_dir(run_dir_for(cfg, run, a.smoke), run.name)
        counts[st.state] = counts.get(st.state, 0) + 1
        print(st.row())
    print("\n" + ", ".join(f"{k}: {v}" for k, v in sorted(counts.items())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
