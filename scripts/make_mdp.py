#!/usr/bin/env python3
"""Write the reference .mdp files to mdp/ from config/md.yaml.

    python scripts/make_mdp.py            # needs all md.yaml decisions filled in
    python scripts/make_mdp.py --smoke    # 1000-step versions, placeholders for open decisions

Each run directory gets its own copy (with its replica's velocity seed) when
scripts/run.py prepares it; the files in mdp/ are for reading and review.
"""
import _setup  # noqa: F401
import argparse
import sys

from lipidmd.config import Config, ConfigError
from lipidmd.mdp import write_all


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    cfg = Config()
    if a.smoke:
        cfg = cfg.with_smoke_fallbacks()
    out = cfg.path("mdp") / ("smoke" if a.smoke else "")
    try:
        stages = write_all(cfg, out, gen_seed=cfg.systems["seeds"]["velocities"][0], smoke=a.smoke)
    except ConfigError as e:
        print(f"Cannot generate mdp files: {e}")
        return 1
    for st in stages:
        print(f"wrote {out / st.mdp_name}   ({st.kind}, {st.nsteps} steps)")
    if cfg.smoke_fallbacks_used:
        print("\nSMOKE placeholders used for:", ", ".join(cfg.smoke_fallbacks_used))
    return 0


if __name__ == "__main__":
    sys.exit(main())
